#!/usr/bin/env python3
"""Read original SI DFT worksheet without rewriting the source workbook."""
import csv, datetime, hashlib, json, zipfile, xml.etree.ElementTree as E
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 src=ROOT/'raw_sources/descriptors.xlsx';ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
 with zipfile.ZipFile(src) as z:
  strings=[''.join(e.itertext()) for e in E.fromstring(z.read('xl/sharedStrings.xml'))]
  sheet=E.fromstring(z.read('xl/worksheets/sheet1.xml'));rows=[]
  for row in sheet.find('s:sheetData',ns):
   result=['']*192
   for cell in row:
    col=''.join(c for c in cell.attrib['r'] if c.isalpha());i=0
    for c in col:i=i*26+ord(c)-64
    v=cell.find('s:v',ns)
    if v is not None:result[i-1]=strings[int(v.text)] if cell.attrib.get('t')=='s' else v.text
   rows.append(result)
 out=ROOT/'primary_si';out.mkdir(exist_ok=False)
 with (out/'DFT_data.csv').open('w',newline='') as f:csv.writer(f).writerows(rows)
 a=pd.read_csv(out/'DFT_data.csv').set_index('ID');b=pd.read_csv(ROOT/'raw_sources/ni_hda_raw.body').set_index('Unnamed: 0')
 common=a.index.intersection(b.index);cols=[c for c in a.columns.intersection(b.columns) if pd.api.types.is_numeric_dtype(a[c]) and pd.api.types.is_numeric_dtype(b[c])];diff=[]
 for c in cols:
  for i in common:
   av,bv=a.loc[i,c],b.loc[i,c]
   if pd.notna(av) and pd.notna(bv) and av!=bv:diff.append({'id':int(i),'property':c,'SI':float(av),'Ni_hDA':float(bv),'difference':float(bv-av)})
 (out/'reference_differences.json').write_text(json.dumps(diff,indent=2,allow_nan=False)+'\n')
 manifest={'created_utc':datetime.datetime.now(datetime.UTC).isoformat(),'source_workbook_sha256':sha(src),'script_sha256':sha(Path(__file__)),'sheet':'DFT_data (sheet1.xml)','rows':len(a),'columns':len(a.columns)+1,'current_reference_rows':len(b),'overlap_IDs':len(common),'only_SI_IDs':list(map(int,a.index.difference(b.index))),'only_current_IDs':list(map(int,b.index.difference(a.index))),'common_numeric_columns':len(cols),'nonidentical_float_cells':len(diff),'maximum_absolute_difference':max((abs(r['difference']) for r in diff),default=0),'precision':'CSV cells preserve exact XLSX XML value text; differences use Pythonfloat64 parsing.','files':{p.name:sha(p) for p in out.iterdir()}}
 (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(manifest,indent=2))
if __name__=='__main__':main()
