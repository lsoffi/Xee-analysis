"""Extract completed C0 FitDiagnostics trials without touching limit results."""
from pathlib import Path
import argparse,fcntl,json,math
import ROOT
from fit_helpers import dump,read_stage
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
P=Path(__file__).resolve().parents[2];R=P/'runs/20261031_separate_sr_c0_limits'


def extract(side,mass,directory=None):
    mass=f'{float(mass):.1f}';D=Path(directory or R/side/f'm{mass}');workspace=D/'workspace.root'
    with (D/'fit_extraction.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        result=dict(side=side,mass=float(mass),degree=6 if side=='low' else 4,CR_degree=4,envelope=False,fit={},flags=[],status='partial',limits_results_preserved=True)
        curves={};trials={}
        for stage in ['b','s']:
            values=[]
            for directory in sorted(D.glob('fitdiagnostics_anchor_*')):
                fd=directory/'fitDiagnostics_c0.root'
                receipts=[directory/'combine.receipt.json',directory/'receipt.json']
                receipt=next((p for p in receipts if p.exists()),None)
                if receipt is None:continue
                native=json.loads(receipt.read_text())
                if not fd.exists():
                    values.append(dict(readable=False,source=str(fd),failure='missing_fit_file',receipt=native));continue
                try:q,_=read_stage(workspace,fd,stage)
                except Exception as error:q=dict(readable=False,source=str(fd),failure=repr(error))
                q['receipt']=native;q['anchor']=directory.name;q['native_success']=native.get('exit_code')==0
                if native.get('exit_code')!=0:q['point_valid']=False
                values.append(q)
            trials[stage]=values
            usable=[q for q in values if q.get('readable') and q.get('point_valid')]
            finite=[q for q in values if q.get('readable') and q.get('exact_deviance') is not None and math.isfinite(q['exact_deviance'])]
            candidates=usable or finite
            if not candidates:continue
            chosen=min(candidates,key=lambda q:q['exact_deviance'])
            q,bins=read_stage(workspace,chosen['source'],stage,D/f'workspace_postfit_{stage}.root')
            q.update(receipt=chosen['receipt'],anchor=chosen['anchor'],native_success=chosen['native_success'])
            q['point_valid']=q['point_valid'] and q['native_success'];q['numerical_ok']=q['numerical_ok'] and q['native_success']
            result['fit'][stage]=q;curves[stage]=bins
            if not q['point_valid']:result['flags'].append('fit_'+stage+'_point_invalid')
            if q['covQual']!=3:result['flags'].append('fit_'+stage+'_covariance_flag')
            if not q['continuity_passed']:result['flags'].append('fit_'+stage+'_C0_density_mismatch')
            if not q['native_success']:result['flags'].append('fit_'+stage+'_native_exit_nonzero')
            if stage=='s':result.update(r=q['r'],r_error=q['error'])
        if len(result['fit'])==2:
            delta=result['fit']['b']['exact_deviance']-result['fit']['s']['exact_deviance'];result['delta_deviance_B_minus_SB']=delta
            if delta<-.01:result['flags'].append('S_minimum_worse_than_B')
            result['status']='complete'
        result['point_valid']=len(result['fit'])==2 and all(q['point_valid'] for q in result['fit'].values())
        result['numerical_ok']=len(result['fit'])==2 and all(q['numerical_ok'] for q in result['fit'].values())
        dump(D/'fit_results.json',result);dump(D/'bins.json',curves);dump(D/'fit_trials.json',trials)
        print(side,mass,result['status'],result.get('r'),result['flags'],flush=True)
        return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('side',choices=['low','high']);p.add_argument('mass');p.add_argument('--directory',type=Path);a=p.parse_args();extract(a.side,a.mass,a.directory)
