"""Apply uniform EDM validity to existing records without changing raw fits.

Also preserve already-run control profiles separately and expose the same
immutable two-anchor point estimator used for every other toy.
"""
from pathlib import Path
import json, math, shutil, fcntl
from toys_helpers import atomic_json
from toys_worker import R, PROTOCOL


def refresh_fit(q):
    edm=q.get('edm')
    q['edm_valid']=bool(edm is not None and math.isfinite(edm) and 0<=edm<.01)
    q['fit_valid']=bool(q.get('fit_valid') and q['edm_valid'])
    q['numerical_validity_version']='status0_edm_lt_0p01_v1'
    return q


def _refresh_cr_locked(path):
    c=json.loads(path.read_text())
    for q in c['attempts']:refresh_fit(q)
    good=[q for q in c['attempts'] if q.get('fit_valid') and q['receipt']['exit_code']==0]
    if good:
        new=min(good,key=lambda q:q['deviance_with_aux'])
        same=all(abs(a['value']-b['value'])<1e-12 for a,b in zip(c['selected']['parameters'],new['parameters']))
        if not same:raise RuntimeError('CR selection would change; dependent anchors require bounded refit: '+str(path))
        c['selected']=new;c['valid']=True
    else:c['selected']=refresh_fit(c['selected']);c['valid']=False
    c['numerical_validity_version']='status0_edm_lt_0p01_v1';atomic_json(path,c)
    return c


def refresh_cr(path):
    with (path.parent / 'cr.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        return _refresh_cr_locked(path)


def main():
    rows=[];needs=[]
    for p in sorted((R/'fits').glob('m*/*/toy*/result.json')):
        q=json.loads(p.read_text())
        if q.get('numerical_validity_version')=='status0_edm_lt_0p01_v1' and q.get('protocol')==PROTOCOL:
            crpath=Path(q['CR_cache_directory'])/'cr_result.json' if q.get('CR_cache_directory') else None
            if crpath and crpath.exists():
                c=json.loads(crpath.read_text());edm=c['selected'].get('edm')
                new=dict(CR_edm=edm,CR_edm_valid=bool(edm is not None and math.isfinite(edm) and 0<=edm<.01))
                if any(q.get(k)!=v for k,v in new.items()):q.update(new);atomic_json(p,q)
            continue
        if q.get('exception'):continue
        backup=p.parent/'before_edm_queue_change';backup.mkdir(exist_ok=True)
        if not (backup/'result.json').exists():shutil.copy2(p,backup/'result.json')
        cr=refresh_cr(Path(q['CR_cache_directory'])/'cr_result.json')
        mp=p.parent/'multistart_v2.json';meta=json.loads(mp.read_text())
        for anchor in meta['anchors']:
            for a in anchor['attempts']:refresh_fit(a)
            good=[a for a in anchor['attempts'] if a.get('fit_valid') and a['receipt']['exit_code']==0]
            anchor['valid']=bool(good)
            if good:anchor['selected']=min(good,key=lambda a:a['deviance_with_aux'])
            else:refresh_fit(anchor['selected'])
        good=[a['selected'] for a in meta['anchors'] if a['valid']]
        if not good:
            needs.append(str(p));q['fit_valid']=False;q['valid']=False;q['edm_valid']=False
            q.setdefault('failure_flags',[]).append('joint_EDM_above_threshold_or_nonfinite')
            atomic_json(p,q);continue
        best=min(good,key=lambda a:a['deviance_with_aux'])
        meta.update(protocol=PROTOCOL,valid=True,both_anchors_valid=all(a['valid'] for a in meta['anchors']),selected_anchor=best['anchor_label'],selected_deviance=best['deviance_with_aux'],numerical_validity_version='status0_edm_lt_0p01_v1')
        atomic_json(mp,meta)
        if q.get('interval_method')=='profile68' and q['toy']==1:
            d=R/f"profiles/m{q['mass']:.1f}/{q['injection']}/toy001";d.mkdir(parents=True,exist_ok=True)
            pf=d/'result.json'
            if not pf.exists():
                profile=dict(q,point_result=str(p),point_r=best['r'],point_deviance=best['deviance_with_aux'],primary_point_estimate_unchanged=True,profile_valid=bool(q.get('interval_valid') and refresh_fit(dict(q))['fit_valid'] and cr['valid']),profile_found_lower_minimum=q['deviance_with_aux']<best['deviance_with_aux']-1e-4,profile_minus_point_r=q['r']-best['r'],coverage_evaluated=False,legacy_profile_completed_before_queue_change=True)
                atomic_json(pf,profile)
        # All masses and injections now use the same two-anchor point estimator.
        q.update(best)
        q.update(protocol=PROTOCOL,multistart=meta,multistart_valid=True,point_from_anchor=True,
                 numerical_validity_version='status0_edm_lt_0p01_v1',CR_valid=cr['valid'],CR_edm=cr['selected'].get('edm'),CR_edm_valid=cr['selected'].get('edm_valid',False),
                 profile_deferred=q['toy']==1,profile_computed=False,interval_method='hesse68',
                 fit_valid=bool(best['fit_valid'] and cr['valid']),joint_fit_valid=best['fit_valid'])
        q['bias']=q['r']-q['rtrue'];q['bias_sigma']=q['bias']/q['error'] if q['error']>0 else None
        q['hesse_valid']=q['covQual']==3 and q['error']>0 and math.isfinite(q['error'])
        q['interval_valid']=bool(best.get('interval_valid') and q['hesse_valid'] and q['fit_valid'])
        q['valid']=q['interval_valid'];q['pull']=q['bias_sigma'] if q['interval_valid'] else None
        q['covered68']=q['interval_low']<=q['rtrue']<=q['interval_high'] if q['interval_valid'] else None
        q['failure_flags']=[s for s in q.get('failure_flags',[]) if s not in ['interval_invalid','Hesse_covariance_not_full_quality','profile_failed_point_retained','joint_EDM_above_threshold_or_nonfinite']]
        if not q['hesse_valid']:q['failure_flags'].append('Hesse_covariance_not_full_quality')
        if not q['interval_valid']:q['failure_flags'].append('interval_invalid')
        if not q['edm_valid']:q['failure_flags'].append('joint_EDM_above_threshold_or_nonfinite')
        atomic_json(p,q)
        rows.append(dict(path=str(p),edm=q['edm'],point_valid=q['fit_valid'],point_r=q['r']))
    gof_updated=[];gof_cr_cache={}
    for path in sorted((R.parent/'gof').glob('toy*/m*/result.json')):
        q=json.loads(path.read_text());crpath=path.parent.parent/'cr_result.json'
        if not crpath.exists():continue
        if crpath not in gof_cr_cache:gof_cr_cache[crpath]=refresh_cr(crpath)
        cr=gof_cr_cache[crpath]
        new=dict(CR_valid=cr['valid'],CR_edm=cr['selected'].get('edm'),CR_edm_valid=cr['selected'].get('edm_valid',False),joint_EDM_available=bool(q.get('nominal_diagnostics',{}).get('available',False)),valid=bool(q.get('valid') and cr['valid']))
        if any(q.get(k)!=v for k,v in new.items()):
            q.update(new);atomic_json(path,q);gof_updated.append(str(path))
    out=dict(updated=rows,gof_CR_updated=gof_updated,requires_bounded_refit=needs,criterion='Finite EDM <0.01, status0 and existing positivity/boundary checks; no bias-dependent selection')
    atomic_json(R/'edm_refresh.json',out);print(json.dumps(out,indent=2))
    return out


if __name__=='__main__':main()
