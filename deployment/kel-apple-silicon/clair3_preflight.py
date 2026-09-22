"""Inspect installed Clair3 and exercise both bundled models on CPU only."""
import datetime
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import shlex
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
PREFIX = ROOT / 'setup/envs/clair3-arm64'
BIN = PREFIX / 'bin'
sys.path.insert(0, str(BIN))
report = {
    'checked_at': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'prefix': str(PREFIX),
    'python': sys.version,
    'machine': platform.machine(),
    'conda_prefix': os.environ.get('CONDA_PREFIX'),
    'scope': 'Read-only installed-package/model checks; no sequencing data or variant calling.',
}

def save():
    (ROOT / 'setup/clair3_preflight.json').write_text(json.dumps(report, indent=2) + '\n')

start = time.monotonic()
report['status'] = 'importing_torch'
save()
import torch
torch.set_num_threads(1)
torch.set_num_interop_threads(1)
report['torch'] = {
    'version': torch.__version__,
    'import_seconds': round(time.monotonic() - start, 3),
    'path': torch.__file__,
    'threads': torch.get_num_threads(),
    'interop_threads': torch.get_num_interop_threads(),
    'cuda_available': torch.cuda.is_available(),
}
report['status'] = 'checking_models'
save()

from clair3.model import Clair3_P, Clair3_F
from clair3.CallVariants import _load_torch_checkpoint, _select_device
import shared.param_p as param_p
import shared.param_f as param_f

model_dir = BIN / 'models/r1041_e82_400bps_hac_v520'
report['models'] = []
for basename, model_class, params in (
    ('pileup', Clair3_P, param_p),
    ('full_alignment', Clair3_F, param_f),
):
    checkpoint_path = model_dir / (basename + '.pt')
    checkpoint = torch.load(str(checkpoint_path), map_location='cpu', weights_only=True)
    state_dict = checkpoint['state_dict'] if 'state_dict' in checkpoint else checkpoint
    indel_length = any(key.startswith('Y_indel_length') for key in state_dict)
    shape = list(params.ont_input_shape)
    model = model_class(add_indel_length=indel_length, predict=True, input_channels=shape[-1])
    model.to('cpu').eval()
    _load_torch_checkpoint(model, str(checkpoint_path), torch.device('cpu'))
    with torch.inference_mode():
        output = model(torch.zeros([1] + shape, dtype=torch.float32))
    report['models'].append({
        'name': basename,
        'path': str(checkpoint_path),
        'sha256': hashlib.file_digest(checkpoint_path.open('rb'), 'sha256').hexdigest(),
        'bytes': checkpoint_path.stat().st_size,
        'state_dict_keys': len(state_dict),
        'parameters': sum(p.numel() for p in model.parameters()),
        'add_indel_length': indel_length,
        'checkpoint_loaded_with_installed_clair3_loader': True,
        'cpu_input_shape': [1] + shape,
        'cpu_output_shape': list(output.shape),
        'cpu_output_all_finite': bool(torch.isfinite(output).all()),
    })
    save()

spec = importlib.util.spec_from_file_location('run_clair3_preflight', BIN / 'run_clair3.py')
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)
extra = '--var_pct_full=1 --ref_pct_full=1 --var_pct_phasing=1 --enable_variant_calling_at_sequence_head_and_tail'
original_argv = sys.argv
sys.argv = ['run_clair3.py', '--bam_fn=/not-executed/input.bam', '--ref_fn=/not-executed/ref.fa',
            '--model_path=' + str(model_dir), '--threads=1', '--platform=ont',
            '--output=/not-executed/output', '--sample_name=KEL_PREFLIGHT',
            '--bed_fn=/not-executed/targets.bed'] + shlex.split(extra)
args = runner.parse_args()
sys.argv = original_argv
report['clair3_version'] = runner.VERSION
report['pipeline_cli_compatible'] = True
report['parsed_pipeline_args'] = {k: getattr(args, k) for k in (
    'var_pct_full', 'ref_pct_full', 'var_pct_phasing',
    'enable_variant_calling_at_sequence_head_and_tail',
    'bed_fn', 'sample_name', 'use_gpu', 'enable_dwell_time')}
report['selected_device'] = str(_select_device(args.use_gpu))
report['architecture_checks'] = subprocess.run(
    ['/usr/bin/file', str(BIN/'python'), str(BIN/'libclair3.so'), str(BIN/'pypy3.11/bin/pypy3')],
    text=True, capture_output=True, check=True).stdout.splitlines()
report['caveats'] = [
    'Successful checkpoint loading and one synthetic zero-tensor inference do not validate calling accuracy.',
    'Installed CallVariants._select_device selects CUDA when requested and available, otherwise CPU; no MPS selection here.',
    'hac_v520 model suitability still requires verified sequencing chemistry/basecaller metadata.',
    'No move-table model is selected; FASTQ-remapped BAMs lack original Dorado mv tags.',
    'Full workflow shell execution and read preprocessing are outside this preflight.',
]
report['elapsed_seconds'] = round(time.monotonic() - start, 3)
report['status'] = 'passed'
save()
print(json.dumps(report, indent=2))
