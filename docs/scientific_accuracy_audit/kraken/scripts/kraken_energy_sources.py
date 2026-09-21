#!/usr/bin/env python3
"""Document primary-source checks for recoverable per-conformer thermochemistry."""
import concurrent.futures,json
from pathlib import Path
from kraken_acquire import ROOT,API,acquire
out=ROOT/'energy_source_search';out.mkdir(exist_ok=True)
repo=json.loads((ROOT/'raw_sources/official_tree_master.body').read_text())['sha']
sigman=json.loads((ROOT/'raw_sources/sigman_tree.body').read_text())['sha']
framework=json.loads((ROOT/'raw_sources/framework_tree.body').read_text())['sha']
jobs={
 'official_readme_lower':f'https://raw.githubusercontent.com/the-matter-lab/kraken/{repo}/readme.md',
 'sigman_data_readme':f'https://raw.githubusercontent.com/SigmanGroup/kraken/{sigman}/data/README.md',
 'sigman_confdata_sample':f'https://raw.githubusercontent.com/SigmanGroup/kraken/{sigman}/validation/new_dft_yamls/00000401_confdata.yml',
 'official_readme':f'https://raw.githubusercontent.com/the-matter-lab/kraken/{repo}/README.md',
 'official_DFT_readme':f'https://raw.githubusercontent.com/the-matter-lab/kraken/{repo}/conf_selection_and_DFT/readme.md',
 'official_releases':'https://api.github.com/repos/the-matter-lab/kraken/releases',
 'framework_model':f'https://raw.githubusercontent.com/Descriptor-Libraries/descriptor-libraries-framework/{framework}/backend/app/app/models/conformer.py',
 'framework_conformer_v2':f'https://raw.githubusercontent.com/Descriptor-Libraries/descriptor-libraries-framework/{framework}/backend/app/app/api/v2/endpoints/conformer.py',
 'framework_import_readme':f'https://raw.githubusercontent.com/Descriptor-Libraries/descriptor-libraries-framework/{framework}/data_import/README.md',
 'framework_schema':f'https://raw.githubusercontent.com/Descriptor-Libraries/descriptor-libraries-framework/{framework}/data_import/schemas/schema_dump.sql',
 'framework_production_config':f'https://raw.githubusercontent.com/Descriptor-Libraries/descriptor-libraries-framework/{framework}/data_import/configs/production.template.yaml',
}
for cid in [31676,38989,54241,54318,63296]:
 for version in ['', '/v2']:jobs[f'api{version.replace("/", "_")}_data_{cid}']=f'{API}{version}/conformers/data/{cid}'
for cid in [31676,54318]:
 for fmt in ['xyz','json']:jobs[f'export_{fmt}_{cid}']=f'{API}/conformers/export/{fmt}/{cid}'
def run(item):
 name,url=item;return name,acquire(url,out/name)
with concurrent.futures.ThreadPoolExecutor(max_workers=5) as pool:results=dict(pool.map(run,jobs.items()))
(out/'request_summary.json').write_text(json.dumps(results,indent=2,sort_keys=True)+'\n')
print(json.dumps({k:(v.get('status'),v.get('bytes'),v.get('error')) for k,v in results.items()},indent=2))
