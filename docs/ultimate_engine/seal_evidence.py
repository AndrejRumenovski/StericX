#!/usr/bin/env python3
"""Bind final non-timing receipts without copying large or failed raw artifacts away."""
import argparse, hashlib, json, subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / '.stericx/ultimate_engine'

def digest(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    assert not a.output.exists(), 'Refuse to overwrite a seal'
    paths = '''baseline/manifest.json
final_scientific_build/manifest.json
final_profile/manifest.json
final_profile/complete_build.json
final_exact_comparison/summary.json
final_science_review/review.json
final_science_review/differences.json
final_science_review/inferred/comparison.json
final_science_review/topology/comparison.json
final_science/volume_reference/complete.json
final_science/all_bins/complete.json
final_science/fixed_union/complete.json
final_science/references_inferred/result.json
final_science/references_topology/result.json
final_science/thermodynamics/summary.json
final_science/models/manifest.json
final_science/model_supplement/manifest.json
final_science/cli_comparison_2/manifest.json
final_science/cli_comparison_4/manifest.json
final_science/cli_comparison_6/manifest.json
final_engineering_resume/receipt.json
final_engineering_interruption.json
final_python_tests/receipt.json
final_combined_python_tests.json
final_combined_python_tests.log
final_python_reference_tests.log
final_dbstep_tests.log
final_advanced/receipt.json
final_advanced/volume/summary.json
final_advanced/sterimol/summary.json
final_advanced/density/summary.json
final_advanced/profiles/summary.json
final_advanced/ensemble_reference/summary.json
final_point_domains/summary.json
final_numerical_replay/summary.json
final_allocations/summary.json
final_simd/summary.json
arm_validation/summary.json
competitive_gate/summary.json
competitive_build/receipt.json
competitive_build_pinned/receipt.json
competitive_preflight/summary.json
competitive_warm/receipt.json
competitive_warm_preflight/summary.json
competitive_memory/summary.json
paired_python_preflight/summary.json
paired_python_preflight_v2/summary.json
final_cli_ab_preflight/summary.json'''.splitlines()
    files = [BASE / name for name in paths]
    files += sorted((BASE / 'final_engineering_resume').glob('*.log'))
    files += sorted((BASE / 'final_valgrind').glob('*.log'))
    files += [ROOT / 'docs/ultimate_engine' / name for name in [
        'msrv_validation.json', 'arm_source_provenance.json', 'final_profile_evidence.json',
        'final_flamegraph.json', 'final_simd.json', 'n1_accuracy_measurements.json',
        'n1_external_comparisons.json', 'n2_accuracy_measurements.json',
        'n2_molecular_comparisons.json', 'n2_allocations.json']]
    records = []
    for f in files:
        records.append({'path': str(f.relative_to(ROOT)), 'bytes': f.stat().st_size, 'sha256': digest(f)})
    changed = subprocess.check_output(['git', 'diff', 'f9c952e2', '--', 'src', 'tests', 'examples',
        'Cargo.toml', 'Cargo.lock', 'pyproject.toml', 'python'], cwd=ROOT)
    assert not changed, 'Production source changed after final replay'
    review = json.loads((BASE / 'final_science_review/review.json').read_text())
    assert review['legacy_reference_preservation_pass'] and review['changed_legacy_reference_values'] == 0
    result = {'scope': 'Final non-timing evidence; failed attempts and historical negative statuses retained.',
        'production_commit': 'f9c952e2', 'production_diff_empty': True,
        'controller_sha256': digest(Path(__file__)), 'records': records,
        'legacy_reference_preservation_pass': True, 'final_performance_complete': False}
    a.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'files_bound': len(records), 'output': str(a.output)}))

if __name__ == '__main__':
    main()
