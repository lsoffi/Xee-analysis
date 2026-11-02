"""Read-only completion and fit-quality audit for the requested fit/limit scan."""
from pathlib import Path
import json,math,datetime
P=Path(__file__).resolve().parents[2];R=P/'runs/20261031_shared_sr7_floatmu_fitlimits'

def main():
    rows=[]
    for i in range(60,121):
        m=f'{i/10:.1f}';D=R/'observed'/('m'+m);p=D/'results.json'
        a=json.loads(p.read_text()) if p.exists() else {}
        receipts={x.name:json.loads((x/'combine.receipt.json').read_text()) for x in (D/'joint').glob('*') if x.is_dir() and (x/'combine.receipt.json').exists()}
        six=a.get('limits',{});fit=a.get('fit',{})
        fit_done=len(fit)==2 and all(receipts.get('fitdiagnostics_anchor_'+s,{}).get('exit_code')==0 for s in ['zero','expected'])
        lim_done=len(six)==6 and all(math.isfinite(v) and v>0 for v in six.values()) and receipts.get('limits',{}).get('exit_code')==0
        rows.append(dict(mass=float(m),fits_complete=fit_done,limits_complete=lim_done,
          fit_point_valid=bool(fit_done and all(v.get('point_valid',False) for v in fit.values())),
          fit_covariance_valid=bool(fit_done and all(v.get('numerical_ok',False) for v in fit.values())),
          observed_limit_crossing_verified=a.get('observed_limit_check',{}).get('verified',False),
          b_deviance=fit.get('b',{}).get('exact_deviance'),
          plots_complete=all((R/'fits'/('m'+m)/('fit_'+s+'.png')).exists() for s in ['SRL','CR','SRH','full_log']),
          flags=a.get('flags',[]),status=a.get('status','not_started')))
    bd=[v['b_deviance'] for v in rows if v['b_deviance'] is not None]
    a=dict(updated_at=datetime.datetime.now().astimezone().isoformat(),target=61,
       fits_complete=sum(v['fits_complete'] for v in rows),limits_complete=sum(v['limits_complete'] for v in rows),
       fit_point_valid=sum(v['fit_point_valid'] for v in rows),fit_covariance_valid=sum(v['fit_covariance_valid'] for v in rows),
       crossing_verified=sum(v['observed_limit_crossing_verified'] for v in rows),plots_complete=sum(v['plots_complete'] for v in rows),
       B_only_deviance_spread=max(bd)-min(bd) if bd else None,
       all_requested_methods_complete=all(v['fits_complete'] and v['limits_complete'] for v in rows),
       rows=rows,note='Method completion is distinct from numerical quality and scientific validation; native expected bands remain unvalidated.')
    p=R/'completion_audit.json';tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(a,indent=2));tmp.replace(p)
    print(json.dumps({k:v for k,v in a.items() if k!='rows'},indent=2))
if __name__=='__main__':main()
