"""Adversarial frozen-evaluation CSV validation; no production edits."""
from models_audit import *
import io


def main():
    source=M/'raw/synthetic_linear/predictions.csv'
    reader=csv.DictReader(io.StringIO(source.read_text()));fields=reader.fieldnames;records=list(reader)
    requests=[]
    for name,replacement in [('nan','NaN'),('infinity','inf'),('finite_mismatch','999')]:
        copy=[dict(r) for r in records];copy[0]['Predicted_ddG_kcal_mol']=replacement
        stream=io.StringIO();writer=csv.DictWriter(stream,fieldnames=fields,lineterminator='\n');writer.writeheader();writer.writerows(copy)
        path=M/'inputs'/f'evaluate_{name}_predictions.csv';freeze(path,stream.getvalue().encode())
        requests.append((name,path))
    phase={'source_prediction_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'inputs':{str(p.relative_to(A)):hashlib.sha256(p.read_bytes()).hexdigest() for _,p in requests},'binary_sha256':hashlib.sha256(BIN.read_bytes()).hexdigest(),'scientific_expectation':'Nonfinite predicted responses cannot be valid finite-model predictions; reject before computing accuracy metrics.'}
    freeze(M/'inputs/evaluate_edges_manifest_before_execution.json',(json.dumps(phase,indent=2)+'\n').encode())
    results=[]
    for name,path in requests:
        out=M/'raw'/f'evaluate_{name}';out.mkdir(exist_ok=True)
        if not(out/'command.json').exists():
            run(f'evaluate_{name}',['evaluate','--data',M/'inputs/synthetic_linear.sigpack','--metadata',M/'inputs/synthetic_linear_labels.csv','--model',M/'raw/synthetic_linear/report.json','--predictions',path,'--output',out/'evaluation.json'])
        cmd=json.loads((out/'command.json').read_text());results.append({'case':name,'input':str(path.relative_to(A)),'exit_code':cmd['exit_code'],'stdout':(out/'stdout').read_text(),'stderr':(out/'stderr').read_text(),'evaluation':json.loads((out/'evaluation.json').read_text()) if(out/'evaluation.json').exists() else None})
    save(M/'results/evaluate_edge_cases.json',results)


if __name__=='__main__':main()
