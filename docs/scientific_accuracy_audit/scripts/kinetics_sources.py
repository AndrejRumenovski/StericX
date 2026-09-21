"""Acquire and hash primary reference documents, including failed responses."""
from pathlib import Path
from datetime import datetime, UTC
import hashlib, json
import requests

ROOT = Path(__file__).resolve().parents[1] / 'kinetics' / 'sources'
SOURCES = {
    'iupac_transition_state.html': 'https://goldbook.iupac.org/terms/view/T06470',
    'iupac_transition_state.pdf': 'https://goldbook.iupac.org/terms/view/T06470/pdf',
    'eyring_1935.pdf': 'https://pubs.aip.org/aip/jcp/article-pdf/3/2/107/18788362/107_1_online.pdf',
    'nist_2022_constants.txt': 'https://physics.nist.gov/cuu/Constants/Table/allascii.txt',
    'iupac_greenbook.pdf': 'https://iupac.org/wp-content/uploads/2019/05/IUPAC-GB3-2012-2ndPrinting-PDFsearchable.pdf',
    'crest_keywords.html': 'https://crest-lab.github.io/crest-docs/page/documentation/keywords.html',
    'crest_v2.12_cregen.f90': 'https://raw.githubusercontent.com/crest-lab/crest/v2.12/src/cregen.f90',
}
ROOT.mkdir(parents=True, exist_ok=True)
manifest = []
for name, url in SOURCES.items():
    p = ROOT / name
    if p.exists():
        raise SystemExit(f'Refusing to overwrite source {p}')
    row = {'url': url, 'retrieved_utc': datetime.now(UTC).isoformat(), 'path': name}
    try:
        response = requests.get(url, timeout=90)
        p.write_bytes(response.content)
        row.update(status=response.status_code, final_url=response.url, content_type=response.headers.get('Content-Type'), sha256=hashlib.sha256(response.content).hexdigest(), bytes=len(response.content))
    except Exception as exc:
        row['error'] = str(exc)
    manifest.append(row)
    (ROOT / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(name, row.get('status', row.get('error')), flush=True)
