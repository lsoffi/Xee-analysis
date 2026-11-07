"""Read-only review of saved v2 toy results; no fits, no result changes."""
from pathlib import Path
import collections,datetime,json,math,statistics
from zoneinfo import ZoneInfo
P=Path(__file__).resolve().parents[2];R=P/'runs/20261030_upper_bias_cheb5_toys'
def finite(x):return isinstance(x,(int,float)) and math.isfinite(x)
def stats(v):
 v=[float(x) for x in v if finite(x)]
 return dict(n=len(v),mean=statistics.mean(v) if v else None,sd=statistics.stdev(v) if len(v)>1 else None)
def validpoint(q):
 return bool(q.get('readable') and q.get('fit_valid') and q.get('receipt',{}).get('exit_code')==0 and finite(q.get('deviance_with_aux')))
rows=[];violations=[];v2=[]
for p in sorted((R/'fits').glob('m*/*/toy*/result.json')):
 q=json.loads(p.read_text());q['_file']=str(p.relative_to(R));rows.append(q)
 if q.get('protocol')!='two_fixed_anchors_v2':continue
 v2.append(q);multi=q.get('multistart')
 if not multi and (p.parent/'multistart_v2.json').exists():multi=json.loads((p.parent/'multistart_v2.json').read_text())
 if multi:
  selecteds=[a.get('selected',{}) for a in multi.get('anchors',[])];valid=[a for a in selecteds if validpoint(a)]
  if valid:
   best=min(valid,key=lambda a:a['deviance_with_aux']);delta=multi.get('selected_deviance',best['deviance_with_aux'])-best['deviance_with_aux']
   if abs(delta)>1e-6:violations.append(dict(file=q['_file'],issue='metadata_does_not_choose_best_valid_anchor',delta=delta))
   if q.get('fit_valid') and finite(q.get('deviance_with_aux')) and q['deviance_with_aux']>best['deviance_with_aux']+1.001e-4:
    violations.append(dict(file=q['_file'],issue='valid_point_worse_than_best_anchor',delta=q['deviance_with_aux']-best['deviance_with_aux']))
  if q.get('fit_valid') and not valid:violations.append(dict(file=q['_file'],issue='valid_point_without_valid_anchor'))
  for a in multi.get('anchors',[]):
   candidates=[b for b in a.get('attempts',[]) if validpoint(b)]
   if candidates and a.get('valid') and a['selected']['deviance_with_aux']>min(b['deviance_with_aux'] for b in candidates)+1e-6:
    violations.append(dict(file=q['_file'],issue='anchor_not_minimum_valid_attempt'))
 if q.get('fit_valid'):
  for key,condition in [('CR_valid',q.get('CR_valid')),('status0',q.get('status')==0),('positive',q.get('positive')),('not_rboundary',not q.get('r_boundary')),('receipt0',q.get('receipt',{}).get('exit_code')==0),('finite_r',finite(q.get('r')))]:
   if not condition:violations.append(dict(file=q['_file'],issue='valid_point_failed_'+key))
 # Expected-median anchor and numerical range are independent of injection.
 if multi and finite(q.get('expected_median_anchor')):
  expected=q['expected_median_anchor'];values=sorted(a['initial_r'] for a in multi['anchors'])
  if values!=[0.,expected]:violations.append(dict(file=q['_file'],issue='anchors_not_fixed_zero_expected',values=values))
  wanted=max(100.,5*expected)
  if multi.get('r_range')!=[-wanted,wanted]:violations.append(dict(file=q['_file'],issue='r_range_depends_on_injection'))

cells=[]
for mass in [i/10 for i in range(105,121)]:
 for label in ['r0','rExpected']:
  rr=[q for q in v2 if q.get('mass')==mass and q.get('label')==label];good=[q for q in rr if q.get('fit_valid')];bad=[q for q in rr if not q.get('fit_valid')]
  bs=stats([q['r']-q['rtrue'] for q in good if finite(q.get('r')) and finite(q.get('rtrue'))]);alls=stats([q['r']-q['rtrue'] for q in rr if finite(q.get('r')) and finite(q.get('rtrue'))]);bads=stats([q['r']-q['rtrue'] for q in bad if finite(q.get('r')) and finite(q.get('rtrue'))]);delta=alls['mean']-bs['mean'] if finite(alls['mean']) and finite(bs['mean']) else None
  causes=collections.Counter()
  for q in bad:
   if q.get('CR_valid') is False:causes['CR_invalid']+=1
   if 'multistart_no_converged_anchor' in q.get('failure_flags',[]):causes['no_valid_anchor']+=1
   if q.get('r_boundary'):causes['r_boundary']+=1
   if q.get('status') is not None and q['status']!=0:causes['joint_status_nonzero']+=1
   if q.get('positive') is False:causes['nonpositive_pdf']+=1
   if q.get('exception'):causes['worker_exception']+=1
  cells.append(dict(mass=mass,injection=label,n=len(rr),valid=len(good),invalid=len(bad),finite_invalid=bads['n'],causes=dict(causes),bias_valid=bs,bias_all_finite=alls,bias_invalid_finite=bads,allfinite_minus_valid_mean=delta,mean_shift_over_valid_sd=delta/bs['sd'] if finite(delta) and bs['sd'] and bs['sd']>0 else None))
primary=[q for q in v2 if q.get('fit_valid')]
counts=dict(completed=len(rows),v2=len(v2),legacy=len(rows)-len(v2),point_valid=len(primary),point_invalid=len(v2)-len(primary),
 point_valid_covQual=dict(collections.Counter(str(q.get('covQual')) for q in primary)),
 point_valid_but_interval_invalid=sum(not q.get('interval_valid') for q in primary),
 point_valid_but_strict_invalid=sum(not q.get('valid') for q in primary),
 flags=dict(collections.Counter(f for q in v2 for f in q.get('failure_flags',[]))))
ordered=sorted([c for c in cells if c['n']>=10 and finite(c['mean_shift_over_valid_sd'])],key=lambda c:abs(c['mean_shift_over_valid_sd']),reverse=True)
q=dict(time=datetime.datetime.now(ZoneInfo('Europe/Rome')).isoformat(),snapshot_only=True,counts=counts,violations=violations,cells=cells,strongest_finite_invalid_selection_sensitivity=ordered[:8],
 interpretation=['Primary fit_valid does not require covQual3 or a valid uncertainty interval. Hesse covariance filtering is confined to secondary pull diagnostics.',
 'Anchor selection uses minimum absolute extended-bin deviance plus auxiliary penalties among converged positive internal-r candidates. Profile result may improve the minimum; no valid result may be worse by more than1e-4.',
 'Failure categories overlap. Finite failed estimates are not trustworthy physics fits; their difference from the valid-only mean is a sensitivity diagnostic for nonrandom numerical failures, not a correction.',
 'With50toys percell the MC precision on standardized mean bias is typically~0.14, so passing ±0.2 atcentralvaluealone is not a95%equivalence demonstration.',
 'Onlypreselectedtoy1profilecontrols exist bydesign; noensemblecoverage is estimated.'])
(R/'methodology_quality_review.json').write_text(json.dumps(q,indent=2));print(json.dumps(counts,indent=2));print('VIOLATIONS',len(violations));print(json.dumps(violations[:5]));print('LARGESTSELECTION',[{k:c[k] for k in ['mass','injection','n','valid','invalid','causes','mean_shift_over_valid_sd']} for c in ordered[:5]])
