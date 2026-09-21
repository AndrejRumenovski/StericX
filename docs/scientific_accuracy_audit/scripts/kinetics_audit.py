"""Independent equation audit; observation and analysis are separate frozen phases.

Run with uv run --extra science python .../kinetics_audit.py freeze, then analyze.
Decimal references use SI constants, not StericX scientific kernels.
Frozen Python modules are used ONLY in observe_python(), explicitly as SUTs.
"""
from __future__ import annotations
import argparse, csv, hashlib, importlib.util, json, math, struct, subprocess, sys
from datetime import UTC, datetime
from decimal import Decimal as D, localcontext
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pandas as pd

AUDIT=Path(__file__).resolve().parents[1]
ROOT=AUDIT/'kinetics'
INP=ROOT/'inputs'
OUT=ROOT/'frozen_outputs'
RES=ROOT/'results'
FROZEN=AUDIT/'frozen/repository'
SUT=AUDIT/'frozen/bin/stericx'
OBSERVER=AUDIT/'frozen/bin/stericx-audit-observer'
# NIST 2022 constants table; thermochemical calorie is exactly 4.184 J.
KB=D('1.380649e-23'); H=D('6.62607015e-34'); NA=D('6.02214076e23')
R=KB*NA/D(4184)
HARTREE_KCAL=D('4.3597447222060e-18')*NA/D(4184)

def write_json(path, value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')

def clean(x):
    if isinstance(x,(float,np.floating)) and not math.isfinite(x):return str(x)
    if isinstance(x,dict):return {k:clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple,np.ndarray)):return [clean(v) for v in x]
    return x

def freeze_manifest(path, files, description):
    write_json(path,{'utc':datetime.now(UTC).isoformat(),'description':description,'files':[{'path':str(p.relative_to(AUDIT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size} for p in sorted(set(files))]})

def load_sut_module(name, filename):
    spec=importlib.util.spec_from_file_location(name,FROZEN/'scripts'/filename)
    m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m

def generate_inputs():
    if (ROOT/'inputs_manifest.json').exists():raise SystemExit('Inputs already frozen; use analyze')
    inputs=[]
    for t in (50,200,273.15,298.15,353.15,400,500,1000):
        for dg in (-100,-50,-20,-5,-2,-1,-1e-6,-0.0,0,1e-8,1e-6,0.1,1,1.82,2.17,5,10,20,30,40,50,100,200):
            inputs.append({'id':f'kinetics_{len(inputs):03}','op':'kinetics','temperature':t,'ddg':dg})
    for t,dg in [(0,1),(-1,1),('NaN',1),('Infinity',1),(298.15,'NaN'),(298.15,'Infinity'),(298.15,'-Infinity'),(1e-30,1),(1e-30,-1),(1e30,1),(298.15,-0.0)]:
        inputs.append({'id':f'edge_{len(inputs):03}','op':'kinetics','temperature':t,'ddg':dg})
    ensembles=[{'id':'single','energies':[4.0],'temperature':298.15},
       {'id':'equal_two','energies':[0,0],'temperature':298.15},
       {'id':'known_two','energies':[0,1],'temperature':298.15},
       {'id':'three','energies':[-2,-1,0],'temperature':298.15},
       {'id':'three_offset','energies':[999998,999999,1000000],'temperature':298.15},
       {'id':'extreme','energies':[0,1000],'temperature':298.15},
       {'id':'warmer','energies':[0,1],'temperature':500},
       {'id':'failed_one','energies':[0,1,2],'statuses':[-1,0,1],'temperature':298.15},
       {'id':'failed_all','energies':[0,1],'statuses':[-1,-1],'temperature':298.15},
       {'id':'nonfinite_one','energies':[0,'nan',1],'temperature':298.15},
       {'id':'window','energies':[0,1,8],'temperature':298.15,'energy_window':6}]
    write_json(INP/'kinetics.json',inputs);write_json(INP/'ensembles.json',ensembles)
    aggs=[]
    conformers=[{'buried_volume':x,'qvbur_min':x/10,'qvbur_max':x/2,'max_delta_qvbur':x/4,'near_vbur':x*.75,'far_vbur':x*.25} for x in (10,20,40)]
    for name,weights in [('equal',[1,1]),('nonunit',[1,3]),('three',[1,2,3]),('single',[1]),('extreme',[1,0]),('tiny',[1e-10,3e-10]),('zero',[0,0]),('negative',[-1,2]),('overflow',[3e38,3e38]),('permuted',[3,1])]:
        cs=conformers[:len(weights)]
        if name=='permuted':cs=list(reversed(cs))
        aggs.append({'id':name,'op':'aggregate','conformers':cs,'weights':weights})
    aggs.append({'id':'nonfinite_descriptor','op':'aggregate','conformers':[{'buried_volume':'NaN'}],'weights':[1]})
    write_json(INP/'aggregation.json',aggs)
    write_json(INP/'ee.json',{'values':[-101,-100,-99.999,-88,-4,0,4,88,99.999,100,101,'nan'],'temperatures':[298.15,353.15]})
    # CREST summaries with exact populations at the program's default temperature.
    w=1/(1+math.exp(-1/(float(R)*298.15)))
    (INP/'crest_default.log').write_text(f' Erel/kcal Etot weight/tot conformer set degen origin\n 1 0.0 -100.0 {w:.15g} {w:.15g} 1 1 input\n 2 1.0 -99.99840639856 {1-w:.15g} {1-w:.15g} 2 1 input\n T /K 298.15\n')
    (INP/'crest_negative.log').write_text(' Erel/kcal Etot weight/tot conformer set degen origin\n 1 0.0 -100.0 -.5 -.5 1 1 input\n 2 1.0 -99.99 1.5 1.5 2 1 input\n T /K 298.15\n')
    (INP/'crest_missing.log').write_text('No population table\n')
    (INP/'a.xyz').write_text('3\nAnalytic conformer A\nH 0 0 0\nC 0 0 2\nF 2 0 1\n')
    (INP/'b.xyz').write_text('3\nAnalytic conformer B\nH 0 0 0\nC 0 0 4\nF 4 0 2\n')
    (INP/'invalid.xyz').write_text('3\nInvalid missing atom\nH 0 0 0\nC 0 0 2\n')
    base={'Reaction_ID':'audit','Ligand_XYZ_Path':'a.xyz','Attach_Atom_Idx':0,'Primary_Bond_Vector_Idx':1,'NBO_Charge':0,'IR_Frequency':0,'Temp_K':298.15,'Exp_ddG_kcal_mol':0,'Conformer_XYZ_Paths':'a.xyz;b.xyz','Conformer_Relative_Energies_kcal_mol':'0;1','Conformer_Boltzmann_Weights':'1;3'}
    variants={'provided':{},'weights_missing':{'Conformer_Boltzmann_Weights':''},'tiny_weights':{'Conformer_Boltzmann_Weights':'1e-10;3e-10'},'negative_weights':{'Conformer_Boltzmann_Weights':'-1;2'},'missing_conformer':{'Conformer_XYZ_Paths':'a.xyz;missing.xyz'},'invalid_conformer':{'Conformer_XYZ_Paths':'a.xyz;invalid.xyz'},'negative_energy':{'Conformer_Relative_Energies_kcal_mol':'-1;0'},'relative_offset':{'Conformer_Relative_Energies_kcal_mol':'1;2'},'single':{'Conformer_XYZ_Paths':'a.xyz','Conformer_Boltzmann_Weights':'1','Conformer_Relative_Energies_kcal_mol':'0'},'temp_zero':{'Temp_K':0}}
    for name,changes in variants.items():
        with (INP/f'parse_{name}.csv').open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(base));writer.writeheader();writer.writerow({**base,**changes})
    freeze_manifest(ROOT/'inputs_manifest.json',[p for p in INP.iterdir() if p.is_file()],'All specified kinetic/conformer test inputs frozen before SUT execution')

def observe_python():
    prep=load_sut_module('frozen_audit_prepare_data','prepare_data.py')
    quantum=load_sut_module('frozen_audit_quantum','stericx_quantum.py')
    observations=[]
    for case in json.loads((INP/'ensembles.json').read_text()):
        es=[float(x) for x in case['energies']];n=len(es);statuses=case.get('statuses',[0]*n)
        try:
            # Isolate deterministic weighting/filtering from external embedding/force field.
            with patch.object(prep.AllChem,'EmbedMultipleConfs',return_value=tuple(range(n))),patch.object(prep.AllChem,'MMFFHasAllMoleculeParams',return_value=True),patch.object(prep.AllChem,'MMFFOptimizeMoleculeConfs',return_value=list(zip(statuses,es))):
                e=prep.embed_and_optimize('CP(C)C',20260919,n,0.1,case.get('energy_window',2000),case['temperature'])
            observations.append({'id':case['id'],'weights':e.boltzmann_weights,'relative_energies':e.relative_energies_kcal_mol,'ids':e.conformer_ids,'statuses':e.mmff_statuses})
        except Exception as exc:observations.append({'id':case['id'],'error':str(exc)})
    thermo=[]
    for t in (298.15,500):
        backend=object.__new__(quantum.QuantumBackend)
        backend.config=quantum.QuantumConfig(cache_dir=ROOT/'unused_cache',temperature_k=t)
        frames=[quantum.XyzFrame(('C',),np.zeros((1,3)),'synthetic',-100+x/float(HARTREE_KCAL)) for x in (0,1)]
        for name in ('default','missing','negative'):
            try:result=backend._conformer_thermodynamics(frames,INP/f'crest_{name}.log');thermo.append({'temperature':t,'summary':name,'result':result})
            except Exception as exc:thermo.append({'temperature':t,'summary':name,'error':str(exc)})
    ee=json.loads((INP/'ee.json').read_text());ees=[]
    for t in ee['temperatures']:
        inputs=pd.Series([float(x) for x in ee['values']]);out=prep.ee_to_ddg(inputs,t)
        ees.extend({'ee':x,'temperature':t,'ddg':y} for x,y in zip(inputs,out))
    write_json(OUT/'python.json',clean({'mmff':observations,'crest':thermo,'ee':ees,'scope':'Mocked energy inputs isolate actual frozen preparation logic; no claim to independent RDKit/MMFF chemical validity or external CREST execution.'}))

def freeze_outputs():
    if (ROOT/'sut_manifest.json').exists():raise SystemExit('SUT outputs already frozen; use analyze')
    generate_inputs()
    for name in ('kinetics','aggregation'):
        cases=json.loads((INP/f'{name}.json').read_text())
        result=subprocess.run([str(OBSERVER)],input=''.join(json.dumps(c)+'\n' for c in cases),text=True,capture_output=True,check=True)
        (OUT/f'{name}.jsonl').write_text(result.stdout);(OUT/f'{name}.stderr').write_text(result.stderr)
    observe_python()
    runs=[]
    for path in sorted(INP.glob('parse_*.csv')):
        name=path.stem;dest=OUT/f'{name}.sigpack'
        command=[str(SUT),'parse','--csv',str(path),'--xyz-dir',str(INP),'--output',str(dest)]
        r=subprocess.run(command,capture_output=True,text=True)
        runs.append({'id':name,'command':command,'returncode':r.returncode,'stdout':r.stdout,'stderr':r.stderr,'packed':str(dest.relative_to(AUDIT)) if dest.exists() else None})
    for dg in (0,1,20,-1):
        command=[str(SUT),'simulate','--ddg',str(dg),'--temp','298.15']
        r=subprocess.run(command,capture_output=True,text=True);runs.append({'id':f'simulate_{dg}','command':command,'returncode':r.returncode,'stdout':r.stdout,'stderr':r.stderr})
    write_json(OUT/'cli_runs.json',runs)
    freeze_manifest(ROOT/'sut_manifest.json',[p for p in OUT.iterdir() if p.is_file()]+[OBSERVER,SUT,Path(__file__)],'Raw SUT observations frozen BEFORE independent comparisons')

def dec(x):return D.from_float(float(x))

def boltzmann(energies,t,degeneracies=None):
    with localcontext() as ctx:
        ctx.prec=70
        es=[dec(x) for x in energies];e0=min(es);g=degeneracies or [1]*len(es)
        raw=[D(gi)*(-(e-e0)/(R*dec(t))).exp() for e,gi in zip(es,g)];z=sum(raw)
        return [float(x/z) for x in raw]

def kinetics_reference(dg,t):
    if not math.isfinite(dg) or not math.isfinite(t) or t<=0:return None
    if abs(dg/(float(R)*t))>1e6:return {'extreme_exponent':True}
    with localcontext() as ctx:
        ctx.prec=75
        g=dec(dg);temp=dec(t);x=g/(R*temp)
        rate=KB*temp/H*(-x).exp()
        # Independently derive two pathway probabilities via common energy offset.
        barriers=[D(0),g];low=min(barriers)
        rates=[(-(b-low)/(R*temp)).exp() for b in barriers]
        pr=100*rates[0]/sum(rates);ps=100*rates[1]/sum(rates)
        return {'rate':rate,'r_percent':pr,'s_percent':ps,'ee_percent':abs(pr-ps)}

def table(path,rows):
    pd.DataFrame(rows).to_csv(path,index=False,float_format='%.17g')

def analyze():
    for name in ('inputs_manifest.json','sut_manifest.json'):
        manifest=json.loads((ROOT/name).read_text())
        for f in manifest['files']:
            assert hashlib.sha256((AUDIT/f['path']).read_bytes()).hexdigest()==f['sha256'],f
    rows=[]
    for o in map(json.loads,(OUT/'kinetics.jsonl').read_text().splitlines()):
        dg=o['input_ddg'];t=o['input_temperature']
        if dg is None or t is None:rows.append({'id':o['id'],'status':'invalid_input','observed_nonfinite':json.dumps(o['_nonfinite'])});continue
        ref=kinetics_reference(dg,t)
        if ref is None or 'extreme_exponent' in ref:rows.append({'id':o['id'],'status':'invalid_temperature' if ref is None else 'extreme_f32_range','observed_nonfinite':json.dumps(o['_nonfinite'])});continue
        for k,v in ref.items():
            observed=o[k];expected=float(v)
            error=abs(observed-expected) if observed is not None and math.isfinite(expected) else None
            rows.append({'id':o['id'],'temperature':t,'ddg':dg,'descriptor':k,'sut':observed,'reference':expected,'reference_decimal':str(v),'abs_error':error,'relative_error':error/abs(expected) if error is not None and expected!=0 else None,'status':'nonfinite_output' if observed is None else 'underflow' if observed==0 and v!=0 and k=='rate' else 'finite'})
    table(RES/'kinetics.csv',rows)
    py=json.loads((OUT/'python.json').read_text());ensemble_rows=[]
    cases={c['id']:c for c in json.loads((INP/'ensembles.json').read_text())}
    for o in py['mmff']:
        c=cases[o['id']]
        if 'error' in o:ensemble_rows.append({'id':o['id'],'error':o['error']});continue
        es=[float(c['energies'][i]) for i in o['ids']];ref=boltzmann(es,c['temperature'])
        for i,s,r in zip(o['ids'],o['weights'],ref):ensemble_rows.append({'id':o['id'],'conformer':i,'sut_weight':s,'reference_weight':r,'abs_error':abs(s-r),'sut_status':c.get('statuses',[0]*len(c['energies']))[i]})
    table(RES/'mmff_weights.csv',ensemble_rows)
    crest_rows=[]
    for o in py['crest']:
        if 'error' in o:continue
        ref=boltzmann([0,1],o['temperature'])
        for i,(s,r) in enumerate(zip(o['result'],ref)):crest_rows.append({'summary':o['summary'],'temperature':o['temperature'],'conformer':i,'sut_weight':s['boltzmann_weight'],'reference_weight':r,'abs_error':abs(s['boltzmann_weight']-r),'expected_status':'invalid_negative_population' if o['summary']=='negative' else 'valid'})
    table(RES/'crest_weights.csv',crest_rows)
    agg_rows=[]
    aggs={x['id']:x for x in json.loads((INP/'aggregation.json').read_text())}
    mapping={'vbur_boltz':'buried_volume','qvbur_min_boltz':'qvbur_min','qvbur_max_boltz':'qvbur_max','max_delta_qvbur_boltz':'max_delta_qvbur','near_vbur_boltz':'near_vbur','far_vbur_boltz':'far_vbur'}
    for o in map(json.loads,(OUT/'aggregation.jsonl').read_text().splitlines()):
        c=aggs[o['id']]
        if 'error' in o:agg_rows.append({'id':o['id'],'error':o['error']});continue
        weights=[dec(np.float32(w)) for w in c['weights']];s=sum(weights)
        for key,ckey in mapping.items():
            vals=[float(x.get(ckey,0)) for x in c['conformers']]
            if not all(map(math.isfinite,vals)):agg_rows.append({'id':o['id'],'descriptor':key,'status':'nonfinite_descriptor_accepted','observed':o[key]});continue
            ref=float(sum(dec(np.float32(x))*w for x,w in zip(vals,weights))/s)
            agg_rows.append({'id':o['id'],'descriptor':key,'sut':o[key],'reference':ref,'abs_error':abs(o[key]-ref) if o[key] is not None else None})
    table(RES/'aggregation.csv',agg_rows)
    ee_rows=[]
    for o in py['ee']:
        ee=float(o['ee']);t=o['temperature']
        if not math.isfinite(ee):status='missing';ref=None
        elif abs(ee)>100:status='invalid_physical_percentage';ref=None
        elif abs(ee)==100:status='infinite_barrier_limit';ref=math.inf
        else:status='valid';ref=float(R)*t*math.log((100+abs(ee))/(100-abs(ee)))
        ee_rows.append({**o,'reference':ref,'status':status,'abs_error':abs(o['ddg']-ref) if ref is not None and math.isfinite(ref) else None})
    table(RES/'ee_to_ddg.csv',ee_rows)
    # Decode actual native packed files without StericX readers.
    packed=[]
    for path in OUT.glob('parse_*.sigpack'):
        numbers=struct.unpack('=16f',path.read_bytes())
        weight=.5 if path.stem=='parse_weights_missing' else 0 if path.stem=='parse_single' else .75
        expected=[3.7+2*weight,1.7,3.47+2*weight]
        packed.append({'case':path.stem,'l':numbers[0],'b1':numbers[1],'b5':numbers[2],'expected_l':expected[0],'expected_b1':expected[1],'expected_b5':expected[2],'max_abs_error':max(abs(a-b) for a,b in zip(numbers,expected)),'stored_energy_span':numbers[14],'conformer_count':numbers[13],'packed_hex':path.read_bytes().hex()})
    table(RES/'native_packed_ensembles.csv',packed)
    # Preserve source CSV conversion comparison for every measured row.
    source=pd.read_csv(FROZEN/'data/official/ni_hda_kraken.csv')
    targets=[]
    for _,row in source.iterrows():
        if pd.isna(row.get('ee')) or pd.isna(row.get('ddG_abs')):continue
        ee=abs(float(row['ee']));ddg=float(row['ddG_abs'])
        if not 0<ee<100:continue
        implied=ddg/(float(R)*math.log((100+ee)/(100-ee)))
        targets.append({'source_id':row.iloc[0],'ee':ee,'published_ddg':ddg,'implied_temperature':implied,'ddg_at_298_15':float(R)*298.15*math.log((100+ee)/(100-ee)),'ddg_at_353_15':float(R)*353.15*math.log((100+ee)/(100-ee))})
    table(RES/'ni_hda_target_temperature.csv',targets)
    # Same ΔΔG, arbitrarily different absolute barriers: analytic counterexample.
    pairs=[]
    for a,b in [(10,11),(20,21)]:
        rr=kinetics_reference(a,298.15)['rate'];ss=kinetics_reference(b,298.15)['rate'];pairs.append({'barrier_r':a,'barrier_s':b,'ddg_s_minus_r':b-a,'rate_r':str(rr),'rate_s':str(ss),'ratio_r_s':str(rr/ss)})
    write_json(RES/'same_selectivity_different_rates.json',pairs)
    frame=pd.DataFrame(rows);finite=frame[frame['abs_error'].notna()]
    summaries={k:{'n':len(g),'max_abs_error':float(g.abs_error.max()),'mae':float(g.abs_error.mean()),'max_relative_error':float(g.relative_error.max())} for k,g in finite.groupby('descriptor')}
    realistic=frame[(frame.descriptor=='rate') & (frame.temperature>=200) & (frame.ddg>=0) & (frame.ddg<=50)]
    summaries['rate_0_to_50_kcal_200_to_1000_K']={'n':len(realistic),'max_relative_error':float(realistic.relative_error.max())}
    summaries['physical_constants']={'R_kcal':str(R),'Hartree_kcal':str(HARTREE_KCAL),'kB':str(KB),'h':str(H),'Avogadro':str(NA)}
    summaries['invalid_inputs']=int((frame.status=='invalid_input').sum())
    summaries['f32_rate_underflows']=int((frame.status=='underflow').sum())
    summaries['f32_nonfinite_outputs']=int((frame.status=='nonfinite_output').sum())
    write_json(RES/'summary.json',clean(summaries))
    freeze_manifest(ROOT/'results_manifest.json',[p for p in RES.iterdir() if p.is_file()],'Independent results; references use equations and frozen primary constants')
    print(json.dumps(clean(summaries),indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('phase',choices=['freeze','analyze']);args=parser.parse_args()
    for p in (INP,OUT,RES):p.mkdir(parents=True,exist_ok=True)
    if args.phase=='freeze':freeze_outputs()
    else:analyze()
