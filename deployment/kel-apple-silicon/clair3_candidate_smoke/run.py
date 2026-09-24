"""Inject three known synthetic SNVs and smoke-test native Clair3 execution."""
import datetime
import gzip
import hashlib
from itertools import zip_longest
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time

import pysam

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'setup/clair3_candidate_smoke'
PREFIX = ROOT / 'setup/envs/clair3-arm64'
SOURCE = ROOT / 'setup/smoke/run/results/mapping/genes/SYNTHETIC__TEST/SYNTHETIC__TEST.variant.bam'
SOURCE_REF = ROOT / 'setup/smoke/run/results/reference/SYNTHETIC__TEST/reference.fasta'
env = os.environ.copy()
env.update({
    'CONDA_PREFIX': str(PREFIX),
    'PATH': str(PREFIX / 'bin') + ':/opt/homebrew/bin:/usr/bin:/bin',
    'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'VECLIB_MAXIMUM_THREADS': '1',
    'PYTHONDONTWRITEBYTECODE': '1', 'TMPDIR': str(OUT / 'tmp'),
})
(OUT / 'tmp').mkdir(exist_ok=True)

def sha256(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()

summary_path = OUT / 'summary.json'
report = {
    'purpose': 'Synthetic candidate/inference smoke only; not biological or clinical validation.',
    'started_at': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'source_bam': str(SOURCE), 'source_bam_sha256_before': sha256(SOURCE),
    'source_reference': str(SOURCE_REF), 'source_reference_sha256': sha256(SOURCE_REF),
    'status': 'preparing',
}

def save():
    summary_path.write_text(json.dumps(report, indent=2) + '\n')

save()
ref = OUT / 'reference.fasta'
shutil.copy2(SOURCE_REF, ref)
subprocess.run([str(PREFIX / 'bin/samtools'), 'faidx', str(ref)], check=True, env=env)
bed = OUT / 'targets.bed'
bed.write_text('SYNTHETIC_REF\t10999\t12000\n')
positions = [11200, 11500, 11800]  # 1-based reference coordinates.
with pysam.AlignmentFile(str(SOURCE), 'rb') as source:
    covering_ids = sorted({read.query_name for read in source.fetch('SYNTHETIC_REF', 10999, 12000)
                           if read.reference_start <= positions[0]-1 and read.reference_end >= positions[-1]})
alternate_haplotype_ids = set(covering_ids[::2])
(OUT / 'alternate_haplotype_read_ids.txt').write_text('\n'.join(sorted(alternate_haplotype_ids)) + '\n')
truth = []
with pysam.FastaFile(str(ref)) as fasta:
    for index, pos in enumerate(positions):
        base = fasta.fetch('SYNTHETIC_REF', pos-1, pos).upper()
        alt = {'A': 'G', 'G': 'A', 'C': 'T', 'T': 'C'}[base]
        truth.append({'chrom': 'SYNTHETIC_REF', 'pos': pos, 'ref': base, 'alt': alt,
                      'expected_genotype': '1/1' if index == 0 else '0/1',
                      'intended_alt_read_ids': len(covering_ids) if index == 0 else len(alternate_haplotype_ids),
                      'assigned_alt_reads': 0, 'edited_reads': 0, 'baseline_allele_counts': {}})
before_calmd = OUT / 'candidates.before_calmd.bam'
edited_read_count = 0
with pysam.AlignmentFile(str(SOURCE), 'rb') as source, pysam.AlignmentFile(str(before_calmd), 'wb', template=source) as dest:
    for read in source:
        sequence = list(read.query_sequence)
        qualities = read.query_qualities
        aligned_positions = {rpos: qpos for qpos, rpos in read.get_aligned_pairs(matches_only=True)}
        changed = False
        for index, variant in enumerate(truth):
            query_position = aligned_positions.get(variant['pos']-1)
            if query_position is not None:
                existing = sequence[query_position].upper()
                variant['baseline_allele_counts'][existing] = variant['baseline_allele_counts'].get(existing,0)+1
                alt_read = index == 0 or read.query_name in alternate_haplotype_ids
                intended_base = variant['alt'] if alt_read else variant['ref']
                variant['assigned_alt_reads'] += int(alt_read)
                if existing != intended_base:
                    sequence[query_position] = intended_base
                    variant['edited_reads'] += 1
                    changed = True
        if changed:
            read.query_sequence = ''.join(sequence)
            read.query_qualities = qualities
            for tag in ('NM', 'MD'):
                if read.has_tag(tag):
                    read.set_tag(tag, None)
            edited_read_count += 1
        dest.write(read)
candidate_bam = OUT / 'candidates.bam'
with candidate_bam.open('wb') as bam_out, (OUT / 'calmd.log').open('w') as log:
    subprocess.run([str(PREFIX / 'bin/samtools'), 'calmd', '-b', str(before_calmd), str(ref)],
                   env=env, stdout=bam_out, stderr=log, check=True)
subprocess.run([str(PREFIX / 'bin/samtools'), 'index', str(candidate_bam)], check=True, env=env)
with pysam.AlignmentFile(str(candidate_bam), 'rb') as bam:
    for variant in truth:
        counts = {}
        for read in bam.fetch('SYNTHETIC_REF', variant['pos']-1, variant['pos']):
            query_pos = dict((rpos,qpos) for qpos,rpos in read.get_aligned_pairs(matches_only=True)).get(variant['pos']-1)
            if query_pos is not None:
                observed = read.query_sequence[query_pos]
                counts[observed] = counts.get(observed,0)+1
        variant['bam_allele_counts'] = counts
checked_reads = 0
with pysam.AlignmentFile(str(SOURCE), 'rb') as original, pysam.AlignmentFile(str(candidate_bam), 'rb') as modified:
    for before, after in zip_longest(original, modified):
        assert before is not None and after is not None
        assert before.query_name == after.query_name
        assert before.flag == after.flag
        assert before.reference_start == after.reference_start
        assert before.cigarstring == after.cigarstring
        assert before.mapping_quality == after.mapping_quality
        assert before.query_qualities == after.query_qualities
        assert before.query_length == after.query_length
        checked_reads += 1
assert all(v['assigned_alt_reads'] == v['intended_alt_read_ids'] for v in truth)
(OUT / 'synthetic_truth.json').write_text(json.dumps(truth, indent=2) + '\n')
report.update({'status': 'calling', 'covering_read_ids': len(covering_ids),
               'alternate_haplotype_read_ids': len(alternate_haplotype_ids),
               'edited_read_count': edited_read_count, 'truth': truth,
               'preserved_cigar_qualities_mapping_verified_reads': checked_reads,
               'candidate_bam_sha256': sha256(candidate_bam),
               'nm_md_recomputed': 'samtools calmd -b',
               'reference_identical_to_source': sha256(ref) == sha256(SOURCE_REF)})
cmd = [str(ROOT / 'setup/run_clair3.sh'), '--bam_fn=' + str(candidate_bam),
       '--ref_fn=' + str(ref), '--bed_fn=' + str(bed), '--threads=2', '--platform=ont',
       '--model_path=' + str(PREFIX / 'bin/models/r1041_e82_400bps_hac_v520'),
       '--output=' + str(OUT / 'output'), '--sample_name=SYNTHETIC__TEST',
       '--var_pct_full=1', '--ref_pct_full=1', '--var_pct_phasing=1',
       '--enable_variant_calling_at_sequence_head_and_tail']
report['command'] = cmd
save()
start = time.monotonic()
with (OUT / 'run.log').open('w') as log:
    result = subprocess.run(cmd, cwd=OUT, env=env, stdout=log, stderr=subprocess.STDOUT)
report['exit_code'] = result.returncode
report['elapsed_calling_seconds'] = round(time.monotonic()-start,3)
report['source_bam_sha256_after'] = sha256(SOURCE)
report['original_bam_unchanged'] = report['source_bam_sha256_before'] == report['source_bam_sha256_after']
run_text = (OUT / 'run.log').read_text()
report['stage_messages'] = [line for line in run_text.splitlines() if re.search(r'\[INFO\] [1-7]/7',line)]
report['processed_position_messages'] = [line for line in run_text.splitlines() if 'Total processed positions' in line]
report['full_alignment_stage_reached'] = '[INFO] 6/7' in run_text and 'full-alignment model' in run_text
vcf = OUT / 'output/merge_output.vcf.gz'
report['merged_vcf_exists'] = vcf.exists()
report['merged_vcf_tbi_exists'] = Path(str(vcf)+'.tbi').exists()
if vcf.exists():
    with gzip.open(vcf,'rt') as handle:
        lines = handle.readlines()
    report['gzip_valid'] = True
    report['vcf_chrom_header'] = next((line.strip() for line in lines if line.startswith('#CHROM')),None)
    records = [line.strip().split('\t') for line in lines if not line.startswith('#')]
    report['variant_record_count'] = len(records)
    report['variants'] = [{'chrom': f[0], 'pos': int(f[1]), 'ref': f[3], 'alt': f[4],
                            'qual': f[5], 'filter': f[6], 'format': f[8], 'sample': f[9]} for f in records]
    def genotype(sample):
        return sorted(sample.split(':')[0].replace('|','/').split('/'))
    comparisons = []
    for expected in truth:
        matches = [v for v in report['variants'] if all(v[k] == expected[k] for k in ('chrom','pos','ref','alt'))]
        comparisons.append({'pos': expected['pos'], 'expected_genotype': expected['expected_genotype'],
                            'exact_allele_found': len(matches) == 1,
                            'genotype_matches': len(matches) == 1 and genotype(matches[0]['sample']) == genotype(expected['expected_genotype'])})
    report['synthetic_truth_comparison'] = comparisons
    queried = subprocess.run([str(PREFIX/'bin/bcftools'),'view','-H','-r','SYNTHETIC_REF:11000-12000',str(vcf)],env=env,text=True,capture_output=True)
    report['indexed_query_exit_code'] = queried.returncode
    report['indexed_query_records'] = len(queried.stdout.splitlines())
report['status'] = 'passed' if (result.returncode == 0 and report.get('variant_record_count',0)>0
    and report.get('gzip_valid') and report.get('merged_vcf_tbi_exists')
    and report.get('indexed_query_exit_code')==0 and report.get('full_alignment_stage_reached')) else 'failed_or_partial'
report['completed_at'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
save()
print(json.dumps(report,indent=2))
