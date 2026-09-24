"""Small execution-only Clair3 smoke on synthetic data."""
import datetime
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'setup/clair3_smoke'
PREFIX = ROOT / 'setup/envs/clair3-arm64'
env = os.environ.copy()
env.update({
    'CONDA_PREFIX': str(PREFIX),
    'PATH': str(PREFIX / 'bin') + ':/opt/homebrew/bin:/usr/bin:/bin',
    'OMP_NUM_THREADS': '1',
    'OPENBLAS_NUM_THREADS': '1',
    'VECLIB_MAXIMUM_THREADS': '1',
    'PYTHONDONTWRITEBYTECODE': '1',
    'TMPDIR': str(OUT / 'tmp'),
})
(OUT / 'tmp').mkdir(exist_ok=True)
ref = OUT / 'reference.fasta'
shutil.copy2(ROOT / 'setup/smoke/run/results/reference/SYNTHETIC__TEST/reference.fasta', ref)
subprocess.run([str(PREFIX / 'bin/samtools'), 'faidx', str(ref)], env=env, check=True)
bed = OUT / 'targets.bed'
bed.write_text('SYNTHETIC_REF\t10999\t12000\n')
bam = ROOT / 'setup/smoke/run/results/mapping/genes/SYNTHETIC__TEST/SYNTHETIC__TEST.variant.bam'
cmd = [
    str(ROOT / 'setup/run_clair3.sh'),
    '--bam_fn=' + str(bam),
    '--ref_fn=' + str(ref),
    '--bed_fn=' + str(bed),
    '--threads=2', '--platform=ont',
    '--model_path=' + str(PREFIX / 'bin/models/r1041_e82_400bps_hac_v520'),
    '--output=' + str(OUT / 'output'),
    '--sample_name=SYNTHETIC__TEST',
    '--var_pct_full=1', '--ref_pct_full=1', '--var_pct_phasing=1',
    '--enable_variant_calling_at_sequence_head_and_tail',
]
report = {
    'started_at': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'purpose': 'Synthetic execution smoke only, not variant-truth or clinical validation.',
    'command': cmd,
    'cwd': str(OUT),
    'target_bed': bed.read_text().strip(),
    'status': 'running',
}
summary = OUT / 'summary.json'
summary.write_text(json.dumps(report, indent=2) + '\n')
start = time.monotonic()
with (OUT / 'run.log').open('w') as log:
    result = subprocess.run(cmd, cwd=OUT, env=env, stdout=log, stderr=subprocess.STDOUT)
report['exit_code'] = result.returncode
report['elapsed_seconds'] = round(time.monotonic() - start, 3)
report['completed_at'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
vcf = OUT / 'output/merge_output.vcf.gz'
report['merged_vcf_exists'] = vcf.exists()
report['vcf_tbi_exists'] = Path(str(vcf) + '.tbi').exists()
run_text = (OUT / 'run.log').read_text()
report['pileup_processed_zero_candidate_positions'] = 'Total processed positions in SYNTHETIC_REF (chunk 1/1) : 0' in run_text
report['full_alignment_stage_reached'] = 'Call variants using full-alignment model' in run_text
report['execution_limitations'] = [
    'Unmodified synthetic reads match the reference in this target; zero candidate positions may cause early exit after pileup.',
    'A valid empty VCF establishes operational empty-call handling, not successful variant detection or full-alignment inference.',
    'Separate setup/clair3_preflight.json records CPU inference of both models on synthetic tensors.',
]
if vcf.exists():
    try:
        with gzip.open(vcf, 'rt') as handle:
            lines = handle.readlines()
        report['gzip_valid'] = True
        report['vcf_chrom_header'] = next((line.strip() for line in lines if line.startswith('#CHROM')), None)
        records = [line.strip().split('\t') for line in lines if not line.startswith('#')]
        report['variant_record_count'] = len(records)
        report['variants'] = [{
            'chrom': fields[0], 'pos': int(fields[1]), 'ref': fields[3], 'alt': fields[4],
            'qual': fields[5], 'filter': fields[6], 'format': fields[8], 'sample': fields[9],
        } for fields in records]
        queried = subprocess.run([str(PREFIX / 'bin/bcftools'), 'view', '-H', '-r', 'SYNTHETIC_REF:11000-12000', str(vcf)], env=env, text=True, capture_output=True)
        report['indexed_region_query_exit'] = queried.returncode
        report['indexed_region_query_records'] = len(queried.stdout.splitlines())
        report['indexed_region_query_stderr'] = queried.stderr
    except Exception as error:
        report['output_validation_error'] = repr(error)
report['status'] = 'passed' if (
    result.returncode == 0 and report.get('gzip_valid') and report.get('vcf_chrom_header')
    and report.get('vcf_tbi_exists') and report.get('indexed_region_query_exit') == 0
) else 'failed'
summary.write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
