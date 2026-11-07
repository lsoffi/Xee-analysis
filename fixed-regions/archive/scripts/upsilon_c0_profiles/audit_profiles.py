"""Independent read-only checks of saved conditional profiles and references."""
from pathlib import Path
import json, math, sys
from collections import Counter
import ROOT

P=Path(__file__).resolve().parents[2];S=P/'scripts/upsilon_c0_profiles';R=P/'runs/20261102_upsilon_c0_profiles'
sys.path.insert(0,str(S));sys.path.insert(0,str(P/'scripts/separate_sr_c0_limits'))
from extract import extract
from fit_helpers import read_stage,dump
from run import read,sha,stamp,MASSES
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)


def close(a,b,tol=1e-6):
    return a is not None and b is not None and math.isfinite(a) and math.isfinite(b) and abs(a-b)<=tol


def audit():
    errors=[];rows=[];points=[];receipts=[];pending_refinements=[]
    state=read(R/'runner.json')
    if state['status']!='completed':errors.append('Main runner not closed')
    rr=R/'reference_refinement_runner.json'
    if rr.exists() and read(rr)['status']!='completed':errors.append('Reference refinement runner not closed')
    reverse=R/'reverse_runner.json'
    if not reverse.exists() or read(reverse)['status']!='completed':errors.append('Reverse continuation not closed')
    branch=R/'branch_runner.json'
    if not branch.exists() or read(branch)['status']!='completed':errors.append('Bounded lower-candidate checks not closed')
    for side in ['low','high']:
        for mass in MASSES:
            D=R/side/f'm{mass:.1f}';label=f'{side}/{mass:.1f}'
            try:
                config=read(D/'configuration.json');profile=read(D/'profile.json');complete=read(D/'completion.json')
                assert profile['status']=='complete' and complete['completed'], 'Profile incomplete'
                assert sha(config['workspace'])==config['workspace_sha256'], 'Source workspace changed'
                assert sha(config['source_fit'])==config['source_fit_sha256'], 'Old FD result changed'
                assert (D/'profile_before_reverse.json').exists(), 'Primary profile was not preserved'
                assert (D/'profile_before_branch_refit.json').exists(), 'Profile before bounded branch refits not preserved'
                assert len(profile['points'])==len(config['grid']), 'Missing profile grid points'
                assert all(close(a,b,1e-8) for a,b in zip([q['r'] for q in profile['points']],config['grid'])), 'Grid differs from immutable plan'
                baseline=profile['reference_min_deviance']
                assert close(baseline,min(q['exact_deviance'] for q in profile['references'])), 'Reference is not lowest valid supplied candidate'
                saved_valid=0;usable=0;minimum_check_failures=0;cov3=0;maximum_c0=0.;minimum_q=math.inf
                for i,point in enumerate(profile['points']):
                    record=dict(side=side,mass=mass,r=point['r'],stored_point_valid=point.get('point_valid'),anchor=point.get('anchor'))
                    if not point.get('source'):
                        assert not point.get('point_valid'),'Missing native source for valid point'
                        points.append(record);continue
                    native=extract(config['workspace'],point['source'],point['r'])
                    assert native['readable'], f'Unreadable selected fit at r={point["r"]}: {native.get("failure")}'
                    assert native['status']==point['status'] and native['covQual']==point['covQual'], 'Stored fit status differs from native'
                    assert close(native['edm'],point['edm'],1e-12), 'Stored EDM differs from native'
                    assert close(native['exact_deviance'],point['exact_deviance']), 'Stored deviance differs from independent recomputation'
                    assert native['continuity_passed']==point['continuity_passed'], 'C0 result differs'
                    if point.get('point_valid'):assert native['point_valid'],'Stored valid point fails independent checks'
                    pd=D/('r_'+format(point['r'],'.12g').replace('.','p'))
                    alltrials=[read(p) for p in pd.glob('*/result.json')]
                    feasible=[q for q in alltrials if q.get('exact_deviance') is not None and math.isfinite(q['exact_deviance']) and q.get('minimum_bin',0)>0 and min(q.get('continuum_minima',{'x':-1}).values())>0 and q.get('continuity_passed')]
                    assert feasible,'No physically feasible candidate available'
                    lowest=min(feasible,key=lambda q:q['exact_deviance'])
                    if lowest['source']!=point['source']:
                        lowcheck=extract(config['workspace'],lowest['source'],point['r'])
                        assert lowcheck['readable'] and close(lowcheck['exact_deviance'],lowest['exact_deviance']),'Lower candidate deviance fails independent check'
                        assert lowcheck['minimum_bin']>0 and min(lowcheck['continuum_minima'].values())>0 and lowcheck['continuity_passed'],'Lower candidate not physically feasible'
                    gap=native['exact_deviance']-lowest['exact_deviance']
                    minimum_passed=gap<=.02
                    assert close(point.get('lower_feasible_candidate_gap'),gap),'Wrong lower-feasible candidate gap'
                    assert point.get('minimum_check_passed')==minimum_passed,'Wrong minimum-check decision'
                    assert point.get('profile_point_valid')==bool(point.get('point_valid') and minimum_passed),'Wrong profile-point usability decision'
                    if not minimum_passed:
                        assert 'unresolved_lower_deviance_candidate' in point['flags'],'Missing unresolved-candidate flag'
                        minimum_check_failures+=1
                    usable+=bool(point.get('profile_point_valid'))
                    assert close(point['q'],native['exact_deviance']-baseline), 'Profile q has wrong reference or factor of two'
                    if point.get('point_valid'):
                        saved_valid+=1;cov3+=native['covQual']==3;minimum_q=min(minimum_q,point['q'])
                        assert point['q']>=-1e-6,'Valid conditional point below chosen reference'
                    maximum_c0=max(maximum_c0,native['C0_max_relative_difference'])
                    for proc in ['dy','Zd','upsilon1s','upsilon2s']:
                        assert close(native['selected_yields'][proc],point['selected_yields'][proc],1e-5),'Selected analytic yield mismatch'
                    record.update(status=native['status'],covQual=native['covQual'],edm=native['edm'],
                        status_history=native['status_history'],num_invalid_NLL=native['num_invalid_NLL'],
                        exact_deviance=native['exact_deviance'],q=point['q'],C0=native['C0_max_relative_difference'],
                        lower_feasible_candidate_gap=gap,minimum_check_passed=minimum_passed,
                        profile_point_valid=point.get('profile_point_valid'),lowest_feasible_source=lowest['source'])
                    points.append(record)
                for ref in profile['references']:
                    if ref['kind'] in ['previous_SB','refined_unconditional']:
                        checked,_=read_stage(config['workspace'],ref['source'],'s')
                    else:checked=extract(config['workspace'],ref['source'],ref['r'])
                    assert checked['readable'] and close(checked['exact_deviance'],ref['exact_deviance']),'Reference deviance mismatch'
                    assert checked['status']==0 and math.isfinite(checked['edm']) and 0<=checked['edm']<.01,'Invalid reference minimizer'
                    assert checked['r']>=0 and min(checked['continuum_minima'].values())>0 and checked['minimum_bin']>0,'Invalid physical reference'
                improvement=min(q['exact_deviance'] for q in profile['references'] if q['kind'].startswith('previous_'))-min(q['exact_deviance'] for q in profile['points'] if q.get('point_valid'))
                if improvement>.001 and not (D/'refined_reference/completion.json').exists():pending_refinements.append(label)
                for path in sorted(D.glob('r_*/*/receipt.json')):
                    receipt=read(path);cmd=read(path.parent/'command.json')
                    assert '--customStartingPoint' in cmd and '--skipSBFit' in cmd, 'Not a conditional full-model native command'
                    assert '--freezeParameters' not in cmd, 'Unrequested parameter freezing'
                    assert cmd[cmd.index('-d')+1]==config['workspace'],'Command uses different workspace'
                    native_path=path.parent/'fitDiagnostics_conditional.root'
                    nf=ROOT.TFile.Open(str(native_path));fr=nf.Get('fit_b') if nf and not nf.IsZombie() else None
                    assert fr, 'Receipt has no readable native conditional result'
                    native_status=dict(status=int(fr.status()),covQual=int(fr.covQual()),edm=float(fr.edm()),
                        status_history=[dict(label=fr.statusLabelHistory(i),code=int(fr.statusCodeHistory(i))) for i in range(fr.numStatusHistory())])
                    nf.Close()
                    receipts.append(dict(side=side,mass=mass,source=str(path),exit_code=receipt['exit_code'],
                        kind='conditional',anchor=path.parent.name,tolerance=float(cmd[cmd.index('--cminDefaultMinimizerTolerance')+1]),**native_status))
                rp=D/'refined_reference/receipt.json'
                if rp.exists():
                    receipt=read(rp);cmd=read(rp.parent/'command.json');q=read(rp.parent/'result.json')
                    assert '--skipBOnlyFit' in cmd and '--customStartingPoint' not in cmd,'Wrong unconditional reference method'
                    assert cmd[cmd.index('--rMin')+1]=='0','Unconditional refinement permits negative r'
                    assert (D/'profile_before_reference_refinement.json').exists(),'Prior reference profile was not preserved'
                    receipts.append(dict(side=side,mass=mass,source=str(rp),exit_code=receipt['exit_code'],kind='reference_refinement',
                        anchor='refined_reference',status=q.get('status'),covQual=q.get('covQual'),edm=q.get('edm'),
                        status_history=q.get('status_history')))
                rows.append(dict(side=side,mass=mass,complete=True,planned=len(config['grid']),selected=len(profile['points']),
                    valid=saved_valid,profile_usable=usable,minimum_check_failures=minimum_check_failures,
                    cov3=cov3,invalid=len(profile['points'])-saved_valid,maximum_C0_difference=maximum_c0,
                    minimum_valid_q=minimum_q,reference_r=profile['reference_r'],reference_deviance=baseline,
                    reference_refinement=profile.get('reference_refinement'),q_at_zero=profile['points'][0]['q']))
            except Exception as error:
                errors.append(label+': '+type(error).__name__+': '+str(error))
    if pending_refinements:errors.append('Authorized reference refinements pending: '+','.join(pending_refinements))
    result=dict(audit_passed=not errors,read_only=True,fits_launched=0,created_at=stamp(),errors=errors,
        counts=dict(expected_profiles=8,complete_profiles=sum(q['complete'] for q in rows),selected_points=len(points),
            valid_points=sum(q.get('stored_point_valid',False) for q in points),native_receipts=len(receipts),
            profile_usable_points=sum(q.get('profile_point_valid',False) for q in points),
            failed_minimum_checks=sum(q.get('minimum_check_passed') is False for q in points),
            native_conditional_receipts=sum(q['kind']=='conditional' for q in receipts),
            native_reference_receipts=sum(q['kind']=='reference_refinement' for q in receipts),
            native_exit0=sum(q['exit_code']==0 for q in receipts),
            receipts_by_anchor=dict(Counter(q['anchor'] for q in receipts))),
        numerical_validity_separate_from_completeness=True,minimum_comparison_tolerance=.02,
        statistical_calibration_performed=False,
        global_minimum_guaranteed=False,rows=rows,selected_points=points,receipts=receipts)
    dump(R/'completion_audit.json',result)
    lines=['# Verifica dei profili sotto le Upsilon','',f'Integrità e completezza: {"PASS" if result["audit_passed"] else "FAIL"}.',
        '', 'Controlli di sola lettura: griglia, input invariati, parametro r fissato, parametri Υ liberi, devianza estesa con dieci vincoli ausiliari, positività, continuità e normalizzazione del profilo.',
        '', 'I punti con fallimenti del minimizzatore restano segnalati. Un punto convergente non è usabile nel profilo se esiste un altro candidato fisicamente ammissibile con devianza inferiore di oltre 0,02, anche quando questo candidato porta un flag del minimizzatore. Non sono state validate copertura o approssimazioni asintotiche; non si garantisce il minimo globale.',
        '', '| Analisi | Massa | Punti | Convergenti | Usabili nel profilo | cov3 | r riferimento | q(0) |','|---|---:|---:|---:|---:|---:|---:|---:|']
    for q in rows:lines.append(f'| {q["side"]} | {q["mass"]:.1f} | {q["selected"]} | {q["valid"]} | {q["profile_usable"]} | {q["cov3"]} | {q["reference_r"]:.6g} | {q["q_at_zero"]:.6g} |')
    if errors:lines.extend(['','Problemi:']+['- '+s for s in errors])
    (R/'completion_audit.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(dict(audit_passed=result['audit_passed'],counts=result['counts'],errors=errors)),flush=True)
    return result


if __name__=='__main__':
    q=audit();raise SystemExit(0 if q['audit_passed'] else 1)
