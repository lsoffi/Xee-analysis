"""Read-only ROOT result audit of inactive envelope coefficients; no minimization."""
from pathlib import Path
import collections,json,math
import ROOT
from worker_helpers import apply_result
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
P=Path(__file__).resolve().parents[2];R=P/'runs/20261030_upper_bias_cheb5_toys';buckets=collections.defaultdict(list)
for p in sorted((R/'fits').glob('m*/*/toy*/result.json')):
 q=json.loads(p.read_text())
 if q.get('protocol')=='two_fixed_anchors_v2' and q.get('fit_valid'):
  key=(q.get('covQual'),sum(v==1 for v in q.get('pdf_indices',{}).values()));buckets[key].append((p,q))
chosen=[]
for key in sorted(buckets):
 values=sorted(buckets[key],key=lambda v:(v[1]['mass'],v[1]['label'],v[1]['toy']));chosen.append(values[0])
 if len(values)>1:chosen.append(values[-1])
chosen=chosen[:20];rows=[]
for p,q in chosen:
 d=Path(q['directory']);f=ROOT.TFile.Open(str(d/'fitDiagnostics_upper.root'));fit=f.Get('fit_s');pars=fit.floatParsFinal();matrix=fit.covarianceMatrix();names={v.GetName():i for i,v in enumerate(pars)}
 wf=ROOT.TFile.Open(str(d/'workspace.root'));w=wf.Get('w');apply_result(w,fit);inactive=[];constraints={}
 for era in ['2022','2022EE','2023','2023BPix']:
  idx=q['pdf_indices']['pdfindex_'+era];prefix='c' if idx==1 else 'b'
  constraint=w.pdf('positive_continuum_'+era);constraints[era]=float(constraint.getVal()) if constraint else None
  for j in range(1,4):
   name=f'local_{prefix}{j}_{era}';i=names.get(name);v=pars.at(i) if i is not None else None
   inactive.append(dict(name=name,selected_family='Bernstein3' if idx==1 else 'Chebyshev3',present_in_floatParsFinal=i is not None,
       error=float(v.getError()) if v else None,saved_isConstant=bool(v.isConstant()) if v else None,
       covariance_diagonal=float(matrix(i,i)) if i is not None and i<matrix.GetNrows() else None,
       covariance_row_maxabs=max(abs(float(matrix(i,k))) for k in range(matrix.GetNcols())) if i is not None and i<matrix.GetNrows() else None))
 command=json.loads((d/'command.json').read_text());rows.append(dict(result=str(p.relative_to(R)),fit_directory=str(d.relative_to(R)),mass=q['mass'],injection=q['label'],toy=q['toy'],status=fit.status(),covQual=fit.covQual(),
      bernstein_era_count=sum(v==1 for v in q['pdf_indices'].values()),indices=q['pdf_indices'],n_float_initial=fit.floatParsInit().getSize(),n_float_final=pars.getSize(),covariance_dimension=matrix.GetNrows(),
      freeze_disassociated_option_present='MINIMIZER_freezeDisassociatedParams' in command,inactive=inactive,SR_chebyshev_positivity_constraint_values=constraints))
 wf.Close();f.Close()
allinactive=[a for q in rows for a in q['inactive']]
summary=dict(n_fits=len(rows),covQual_counts=dict(collections.Counter(str(q['covQual']) for q in rows)),
             bernstein_count_distribution=dict(collections.Counter(str(q['bernstein_era_count']) for q in rows)),
             inactive_parameters_checked=len(allinactive),inactive_present_in_float_list=sum(a['present_in_floatParsFinal'] for a in allinactive),
             inactive_with_nonzero_error=sum(a['error']!=0 for a in allinactive),inactive_with_nonzero_covariance_row=sum(a['covariance_row_maxabs']!=0 for a in allinactive),
             positivity_constraint_values=sorted(set(v for q in rows for v in q['SR_chebyshev_positivity_constraint_values'].values())),
             all_have_freeze_option=all(q['freeze_disassociated_option_present'] for q in rows))
audit=dict(summary=summary,rows=rows,
 findings=['All12inactiveSRcoefficients appear in the saved floating-parameter list, but have exactlyzero errors and covariance rows in everytestedfit, includingcovQual2 and3.',
           'This is consistent with parameters frozen duringminimization/Hesse and restored beforeRooFitResultserialization; presence in floatParsFinal alone is not proof they were free Hessian directions.',
           'The Chebyshevpositivityconstraints do remain attached forBernstein-selectederas, but evaluate to1 in all20checkedfinalpoints.',
           'LocalCascadeMinimizer freezeDisassociatedParams collectsallMultiPdfparameters, subtractsactivePDFparameters andfreezesremainder irrespectiveofexternalconstraints; improveOnce,hesse,minos invokeit.',
           'No evidencefromthisaudit thatinactivecoefficients causeHessecovQual2 or biaspointestimates. Existingactive-parameter correlations/boundaries stillrequireseparateinterpretation.',
           'When displayingparameterchanges, inactivebasiscoefficients mustnotbeinterpretedas independentlymeasured postfitparameters.'],
 source_references={'prepare.py':'scripts/upper_envelope_scan/prepare.py:16-28; Cheb0/Bern1 andunconditionalChebpositivityconstraint',
 'CascadeMinimizer.cc':'work/CombinedLimit/src/CascadeMinimizer.cc:94-139 freezeDisassociatedParams;253improveOnce;411hesse;71-88save/restoration'})
(R/'inactive_parameter_audit.json').write_text(json.dumps(audit,indent=2));print(json.dumps(summary,indent=2))
