"""Read optional double-precision nominal-state GoF diagnostics.

The diagnostic runtime captures these values BEFORE saturated minimization.
Ordinary --trackParameters output is not used to infer the nominal state.
"""
from pathlib import Path
import math
import ROOT
from toys_helpers import continuum_positivity, exact_deviance

P=Path(__file__).resolve().parents[2]
DIAGNOSTIC_COMBINE=P/'work/shared_gof_diagnostic_runtime/v2/combine'
RTD_FLAG='GOF_TRACK_NOMINAL'
META={'minimize_ok','status','edm','covQual','fitresult_available','fitresult_state_max_abs_difference'}


def read_nominal_diagnostics(rootfile, workspacefile, expected_score=None):
    """Return nominal minimizer, positivity and saved-state checks; no fits run."""
    result={'available':False,'attempted':True,'outcome_complete':True,'runtime':str(DIAGNOSTIC_COMBINE),'capture':'before_saturated_minimization',
            'EDM_threshold':.01,'covQual_required':False,'valid':False,'flags':[]}
    tf=ROOT.TFile.Open(str(rootfile))
    if not tf or tf.IsZombie():
        result['flags']=['unreadable_GoF_ROOT'];return result
    tree=tf.Get('limit')
    if not tree or tree.GetEntries()!=1:
        tf.Close();result['flags']=['missing_or_nonunique_GoF_score'];return result
    tree.GetEntry(0);score=float(tree.limit)
    branches={b.GetName() for b in tree.GetListOfBranches()}
    captured={name[len('gof_nominal_'):]:float(getattr(tree,name)) for name in branches if name.startswith('gof_nominal_')}
    tf.Close()
    if not captured:
        result['flags']=['nominal_diagnostics_missing'];return result
    result.update(available=True,score=score,**{key:captured.get(key) for key in META})
    result['r']=captured.get('r')
    result['minimizer_ok']=captured.get('minimize_ok')==1 and captured.get('status')==0
    edm=captured.get('edm');result['edm_valid']=edm is not None and math.isfinite(edm) and 0<=edm<.01
    diff=captured.get('fitresult_state_max_abs_difference')
    result['fitresult_state_consistent']=captured.get('fitresult_available')==1 and diff is not None and math.isfinite(diff) and diff<=1e-7
    if not result['minimizer_ok']:result['flags'].append('nominal_minimizer_or_status_failure')
    if not result['edm_valid']:result['flags'].append('nominal_EDM_above_threshold_or_nonfinite')
    if not result['fitresult_state_consistent']:result['flags'].append('nominal_fitresult_state_mismatch')
    parameters={name:value for name,value in captured.items() if name not in META}
    result['parameters']=parameters
    if not all(math.isfinite(v) for v in parameters.values()):
        result['flags'].append('nonfinite_nominal_parameter');return result
    wf=ROOT.TFile.Open(str(workspacefile));w=wf.Get('w') if wf and not wf.IsZombie() else None
    if not w:
        result['flags'].append('unreadable_workspace_for_nominal_audit');return result
    for name,value in parameters.items():
        if w.cat(name):w.cat(name).setIndex(int(round(value)))
        elif w.var(name):w.var(name).setVal(value)
    obs=w.obj('ModelConfig').GetObservables()
    floating=[v.GetName() for v in w.obj('ModelConfig').GetPdf().getParameters(obs)
              if v.InheritsFrom('RooRealVar') and not v.isConstant()]
    missing=[name for name in floating if name not in parameters]
    result['missing_floating_parameters']=missing
    result['snapshot_complete']=not missing
    if missing:result['flags'].append('nominal_parameter_snapshot_incomplete')
    indices={name:int(round(value)) for name,value in parameters.items() if name.startswith('pdfindex_')}
    result['pdf_indices']=indices
    if len(indices)!=4:result['flags'].append('missing_shared_family_indices')
    try:
        minima=continuum_positivity(w);result['continuum_minima']=minima
        result['continuum_positive']=all(math.isfinite(v) and v>0 for v in minima.values())
    except Exception as error:
        result['continuum_positive']=False;result['positivity_exception']=repr(error)
    if not result['continuum_positive']:result['flags'].append('nonpositive_continuum')
    if result['snapshot_complete']:
        deviance,minbin=exact_deviance(w,True)
        result.update(deviance_with_aux=deviance,minimum_bin_expectation=minbin,total_bins_positive=math.isfinite(minbin) and minbin>0)
        correction=2*sum(p.getCorrection() for p in w.obj('ModelConfig').GetPdf().getComponents() if p.InheritsFrom('RooMultiPdf'))
        result['discrete_penalty']=correction
        result['score_reconstructed']=deviance+correction if deviance is not None else None
        result['score_state_difference']=result['score_reconstructed']-score if result['score_reconstructed'] is not None else None
        result['score_state_consistent']=result['score_state_difference'] is not None and abs(result['score_state_difference'])<1e-5
        if not result['total_bins_positive']:result['flags'].append('nonpositive_total_bin_expectation')
        if not result['score_state_consistent']:result['flags'].append('nominal_score_state_mismatch')
    if expected_score is not None:
        result['previous_score']=float(expected_score);result['previous_score_difference']=score-float(expected_score)
        result['previous_score_parity']=math.isfinite(score) and abs(result['previous_score_difference'])<1e-7
        if not result['previous_score_parity']:result['flags'].append('diagnostic_score_differs_from_original')
    result['valid']=not result['flags']
    wf.Close()
    return result
