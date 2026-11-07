"""One bounded numerical retry pass for point-invalid observed S+B fits.

The model and all original fit outputs are retained. A covariance-only flag is
not a trigger. Each mass receives at most two additional starts, with a common
EDM < 0.01 point-validity criterion. The already investigated m=12 is excluded.
Default action prints/saves the plan; --run starts the sequential pass only
after both original FitDiagnostics anchors have finished for all 61 masses.
"""
from pathlib import Path
import argparse, hashlib, json, math, shutil
import ROOT
from run_observed import execute, common, extract, C, R
from observed_helpers import read_stage, dump

ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
MASSES=[f'{n/10:.1f}' for n in range(60,121)]
EDM_MAX=.01

def finite(x):
    return isinstance(x,(int,float)) and math.isfinite(x)

def point_valid(q):
    return bool(q.get('status')==0 and finite(q.get('edm')) and 0<=q['edm']<EDM_MAX
                and finite(q.get('exact_deviance')) and finite(q.get('minimum')) and q['minimum']>0
                and finite(q.get('minimum_bin')) and q['minimum_bin']>0 and not q.get('r_boundary',False))

def ready_masses():
    return [m for m in MASSES if all((R/'observed'/f'm{m}'/'joint'/
            f'fitdiagnostics_anchor_{label}'/'combine.receipt.json').exists()
            for label in ['zero','expected'])]

def make_plan(masses=None):
    rows=[]
    for mass in masses or MASSES:
        if mass=='12.0':continue
        D=R/'observed'/f'm{mass}';path=D/'results.json'
        if not path.exists():continue
        a=json.loads(path.read_text());s=a.get('fit',{}).get('s',{})
        if not s or point_valid(s):continue
        b=a.get('fit',{}).get('b',{})
        rows.append(dict(mass=mass,initial_fit=s,background_seed_valid=point_valid(b),
                         initial_flags=a.get('flags',[])))
    return rows

def retry_one(record,timeout=1200):
    mass=record['mass'];D=R/'observed'/f'm{mass}';J=D/'joint'
    audit_path=D/'bounded_retry_pass.json'
    if audit_path.exists():
        old=json.loads(audit_path.read_text())
        if old.get('status')=='completed':return old
    assert record['background_seed_valid'],f'No valid interior B-only seed at {mass}'
    snapshot=J/'retry_pass_source_b.root';original=J/'workspace_postfit_b.root'
    if not snapshot.exists():shutil.copy2(original,snapshot)
    original_protocol=json.loads((D/'observed_protocol.json').read_text())
    lower,upper=original_protocol['initial_r_range']
    assert lower<0<upper
    prior=record['initial_fit'];rhat=prior.get('r')
    if not finite(rhat):rhat=0.
    # This is an initial numerical step only, never a prior or parameter bound.
    error=prior.get('error');step=max(error if finite(error) and error>0 else 0.,.05*abs(rhat),.001)
    trials=[]
    for label,start,strategy,tolerance in [('pass_zero',0.,1,.001),('pass_selected',rhat,2,.00001)]:
        directory=J/f'fitdiagnostics_anchor_retry_{label}';directory.mkdir(exist_ok=True)
        path=directory/'workspace.root'
        if not (directory/'combine.receipt.json').exists():
            f=ROOT.TFile.Open(str(snapshot));w=f.Get('w')
            w.var('r').setRange(lower,upper);w.var('r').setVal(start);w.var('r').setError(step)
            w.writeToFile(str(path));f.Close()
        command=[C,'-M','FitDiagnostics']+common(path,mass,strategy,tolerance)+[
            '--saveShapes','--saveNormalizations','--keepFailures','--robustFit','0',
            '--minos','none','--skipBOnlyFit','--rMin',str(lower),'--rMax',str(upper),
            '--preFitValue',str(start)]
        receipt=execute(command,directory,timeout=timeout)
        q,bins=read_stage(J/'workspace_seeded.root',directory/'fitDiagnostics_shared.root','s')
        q.update(directory=str(directory),receipt=receipt,initial_r=start,initial_r_error=step,
                 strategy=strategy,tolerance=tolerance,initial_snapshot=str(snapshot))
        q['uniform_point_valid']=point_valid(q) and receipt['exit_code']==0
        dump(directory/'retry_result.json',q);trials.append(q)
        print('RETRY',mass,label,q.get('r'),q.get('status'),q.get('covQual'),q.get('edm'),q.get('exact_deviance'),flush=True)
    updated=extract(mass)
    result=dict(status='completed',mass=float(mass),initial_fit=record['initial_fit'],
                initial_flags=record['initial_flags'],trials=trials,max_attempts=2,
                covariance_only_failure_not_trigger=True,EDM_max=EDM_MAX,
                unchanged_r_range=[lower,upper],original_outputs_preserved=True,
                source=str(original),immutable_seed=str(snapshot),
                seed_sha256=hashlib.sha256(snapshot.read_bytes()).hexdigest(),
                selected_point_valid=point_valid(updated['fit']['s']),selected_fit=updated['fit']['s'],
                remaining_flags=updated['flags'])
    dump(audit_path,result);return result

def main(args):
    masses=[f'{float(m):.1f}' for m in args.masses] if args.masses else MASSES
    assert all(m in MASSES for m in masses)
    ready=ready_masses();plan=make_plan(masses)
    payload=dict(ready_original_FD_masses=len(ready),required_original_FD_masses=61,
                 candidates=plan,excluded_masses=['12.0'],EDM_max=EDM_MAX,max_new_starts_per_mass=2,
                 unchanged_physical_model=True,worker_count=1)
    dump(R/'observed_retry_plan.json',payload)
    print('PLAN',len(ready),'of61FDready;candidates',[q['mass'] for q in plan],flush=True)
    if not args.run:return
    assert len(ready)==61,'Wait until both original FitDiagnostics anchors finish at all61 masses'
    pass_path=R/'observed_retry_pass.json'
    if pass_path.exists():
        previous=json.loads(pass_path.read_text())
        if previous.get('status')=='completed':
            print('Retry pass already completed; no additional fits.',flush=True)
            return
        payload=previous['plan'];plan=payload['candidates']
    results=[]
    dump(pass_path,dict(status='running',plan=payload,results=results))
    for record in plan:
        results.append(retry_one(record,args.timeout))
        dump(R/'observed_retry_pass.json',dict(status='running',plan=payload,results=results))
    dump(R/'observed_retry_pass.json',dict(status='completed',plan=payload,results=results))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--run',action='store_true')
    parser.add_argument('--masses',nargs='+');parser.add_argument('--timeout',type=int,default=1200)
    main(parser.parse_args())
