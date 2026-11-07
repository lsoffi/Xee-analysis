#!/usr/bin/env python3
"""Prospective C0 toy validation runner. No campaign runs without --execute.

This is a new portable implementation for the CURRENT model, not a claim that
historical envelope toy results validate C0. A generator must be explicitly
provided for bias; use a fitted, vetted alternative with the same data layout.
"""
from pathlib import Path
import argparse, csv, hashlib, json, math, os, signal, subprocess, sys, time
P=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(P/'analysis/separate_sr_c0_limits'))

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,q):
 p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(q,indent=2,allow_nan=False)+'\n')
def native(command,directory):
 directory.mkdir(parents=True,exist_ok=True);receipt=directory/'receipt.json'
 if receipt.exists():
  assert json.loads((directory/'command.json').read_text())==list(map(str,command)),'Changed native command: use a new tag'
  return json.loads(receipt.read_text())
 if (directory/'process.json').exists():raise RuntimeError('Interrupted native attempt retained; inspect and use a fresh output tag: '+str(directory))
 save(directory/'command.json',list(map(str,command)));t=time.monotonic()
 with (directory/'combine.log').open('w') as f:
  result=subprocess.Popen(list(map(str,command)),cwd=directory,stdout=f,stderr=subprocess.STDOUT,start_new_session=True,
      env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1'))
  save(directory/'process.json',dict(pid=result.pid,started_unix=time.time(),timeout_seconds=1800))
  try:rc=result.wait(timeout=1800)
  except subprocess.TimeoutExpired:
   os.killpg(result.pid,signal.SIGTERM)
   try:result.wait(timeout=5)
   except subprocess.TimeoutExpired:os.killpg(result.pid,signal.SIGKILL);result.wait()
   rc=124
 q=dict(exit_code=rc,seconds=time.monotonic()-t,timeout_seconds=1800);save(receipt,q);return q

def stats(values):
 import numpy as np
 a=np.asarray(values,dtype=float);n=len(a)
 return dict(n=n,mean=float(a.mean()) if n else None,sd=float(a.std(ddof=1)) if n>1 else None,
   se=float(a.std(ddof=1)/math.sqrt(n)) if n>1 else None)

def summarize(rows,truth,observed=None):
 import numpy as np
 good=[q for q in rows if q.get('fit_valid')]
 allfinite=[q for q in rows if q.get('r') is not None and math.isfinite(q['r'])]
 d=np.array([q['r']-truth for q in good]);s=stats(d.tolist());normalized=None
 if len(d)>2 and d.std(ddof=1)>0:
  rng=np.random.default_rng(41022);samples=rng.choice(d,size=(2000,len(d)),replace=True)
  widths=samples.std(axis=1,ddof=1);samples=samples[widths>0];widths=widths[widths>0]
  b=samples.mean(axis=1)/widths
  normalized=dict(value=float(d.mean()/d.std(ddof=1)),bootstrap_se=float(b.std(ddof=1)),
                  bootstrap_ci95=np.quantile(b,[.025,.975]).tolist(),resamples=2000)
 q=dict(n_total=len(rows),n_valid=len(good),n_failed=len(rows)-len(good),bias=s,
   bias_all_finite=stats([x['r']-truth for x in allfinite]),normalized_bias=normalized,
   pull_hesse=stats([(x['r']-truth)/x['error'] for x in good if x['covQual']==3 and x['error']>0]),
   coverage_measured=False,conditional_on_numerical_success=True)
 if observed is not None:
  scores=[x['gof_score'] for x in rows if x.get('gof_score') is not None]
  n=len(scores);k=sum(x>observed for x in scores)
  q['gof']=dict(observed_score=observed,n_readable=n,n_requested=len(rows),n_gt_observed=k,
    tail_convention='strict T_toy > T_observed, matching the supplied analysis definition',
    p_raw=k/n if n else None,p_plus_one=(k+1)/(n+1) if n else None,
    binomial_mc_se=math.sqrt((k/n)*(1-k/n)/n) if n else None,
    status='all native scores readable' if n==len(rows) else 'incomplete; no unconditional p-value',
    count_policy='All finite native scores retained, including FD disagreement; inspect individual flags.',
    n_FD_disagreements=sum(x.get('gof_fd_agreement') is False for x in rows))
 return q

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('method',choices=['bias','gof'])
 ap.add_argument('--side',choices=['low','high'],required=True);ap.add_argument('--mass',type=float,required=True)
 ap.add_argument('--toys',type=int,default=100);ap.add_argument('--seed',type=int,default=512001)
 ap.add_argument('--generator',type=Path,help='Explicit vetted generator workspace (same SR/CR layout and signal mass). Required for bias.')
 ap.add_argument('--injection',default='0',help='0, expected, or a numeric r. Expected is native median limit, possibly provisional.')
 ap.add_argument('--tag',required=True,help='New output directory label; choose a new tag when changing protocol.')
 ap.add_argument('--execute',action='store_true');a=ap.parse_args();os.chdir(P)
 if a.toys<1:ap.error('--toys must be positive')
 if a.method=='bias' and a.generator is None:ap.error('Bias requires --generator. Same fitted model tests closure only.')
 if a.method=='gof' and a.injection!='0':ap.error('This GoF is generated under B-only: use injection 0.')
 d=P/'runs/20261031_separate_sr_c0_limits'/a.side/f'm{a.mass:.1f}'
 source=a.generator.resolve() if a.generator else d/'workspace_postfit_b.root'
 out=P/'runs/validation_c0'/a.tag/a.side/f'm{a.mass:.1f}'
 plan=dict(method=a.method,side=a.side,mass=a.mass,toys=a.toys,seed=a.seed,generator=str(source),
   injection=a.injection,fit_template=str(d/'workspace_postfit_b.root'),output=str(out),
   two_fit_anchors=['0','expected median'],r_floating=True,Ymeans_floating=True,
   generation='Conditional frequentist binned Poisson plus auxiliary measurements; no pre-generation refit.',
   validation_status='Prospective current-C0 validation; not a historical completed result.')
 print(json.dumps(plan,indent=2))
 if not a.execute:return
 import ROOT,numpy as np
 from fit_helpers import apply_fit,inspect_model,continuity,polynomial_minimum
 ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
 ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
 expected=json.loads((d/'results.json').read_text())['limits']['expected_50']
 rtrue=expected if a.injection=='expected' else float(a.injection)
 bounds=json.loads((d/'fit_protocol.json').read_text())['r_range'];bound=max(abs(x) for x in bounds)
 if abs(rtrue)>=bound:raise ValueError('Injected r outside fit range; define an expanded protocol first.')
 plan.update(r_true=rtrue,expected=expected,r_range=[-bound,bound],generator_sha256=sha(source),
             template_sha256=sha(d/'workspace_postfit_b.root'))
 if (out/'protocol.json').exists():assert json.loads((out/'protocol.json').read_text())==plan,'Changed protocol: use a new tag'
 else:save(out/'protocol.json',plan)
 # Freeze truth values only in a separate generator; fit template stays floating.
 gen=out/'generation';gen.mkdir(exist_ok=True);genws=gen/'workspace.root'
 if not genws.exists():
  gf=ROOT.TFile.Open(str(source));w=gf.Get('w');tf=ROOT.TFile.Open(str(d/'workspace_postfit_b.root'));tw=tf.Get('w')
  for name in ['mass_SR','mass_CR']:
   x=w.var(name);tx=tw.var(name)
   assert x and (x.getMin(),x.getMax(),x.getBins())==(tx.getMin(),tx.getMax(),tx.getBins()),'Generator layout mismatch'
  assert {x.GetName() for x in w.obj('ModelConfig').GetGlobalObservables()}=={x.GetName() for x in tw.obj('ModelConfig').GetGlobalObservables()}
  w.var('r').setRange(-1e6,1e6);w.var('r').setVal(rtrue)
  for v in w.allVars():
   if v.GetName().startswith(('local_','Ymean_','CMS_EXO25020_bkgScale_')):v.setConstant(True)
  w.saveSnapshot('generation',w.allVars());w.writeToFile(str(genws));gf.Close();tf.Close()
 channels=','.join(f'{region}_ch{i}' for region in ['SR','CR'] for i in range(1,5))
 command=[P/'bin/combine','-M','GenerateOnly',genws,'--snapshotName','generation','--saveToys',
  '--toysFrequentist','--bypassFrequentistFit','--expectSignal',str(rtrue),'--genBinnedChannels',channels,
  '-t',str(a.toys),'-s',str(a.seed),'-m',str(a.mass),'-n','_c0_generated']
 if native(command,gen)['exit_code']!=0:raise RuntimeError('Toy generation failed')
 files=list(gen.glob('higgsCombine*.root'));assert len(files)==1
 fh=ROOT.TFile.Open(str(files[0]));rows=[]
 common=['--cminDefaultMinimizerStrategy','2','--cminDefaultMinimizerTolerance','.001',
  '--cminDefaultMinimizerPrecision','1e-12','--X-rtd','MINIMIZER_freezeDisassociatedParams']
 def score(ws,folder,data='data_obs'):
  cmd=[P/'bin/combine','-M','GoodnessOfFit',ws,'-D',data,'--algo','saturated',
       '--rMin',str(-bound),'--rMax',str(bound),'-m',str(a.mass),'-n','_c0_gof']+common
  rc=native(cmd,folder);fs=list(folder.glob('higgsCombine*.root'))
  if rc['exit_code']!=0 or len(fs)!=1:return None
  f=ROOT.TFile.Open(str(fs[0]));t=f.Get('limit');value=None
  if t and t.GetEntries():t.GetEntry(0);value=float(t.limit)
  f.Close();return value if value is not None and math.isfinite(value) else None
 observed=score(d/'workspace_postfit_b.root',out/'gof_observed') if a.method=='gof' else None
 if a.method=='gof' and observed is None:raise RuntimeError('Observed GoF unreadable')
 for i in range(1,a.toys+1):
  td=out/f'toy{i:04d}';td.mkdir(exist_ok=True);result=td/'result.json'
  if result.exists():rows.append(json.loads(result.read_text()));continue
  toy=fh.Get(f'toys/toy_{i}');snapshot=fh.Get(f'toys/toy_{i}_snapshot');assert toy and snapshot
  ws=td/'workspace.root'
  f=ROOT.TFile.Open(str(d/'workspace_postfit_b.root'));w=f.Get('w')
  getattr(w,'import')(toy.Clone('toy_data'))
  for g in w.obj('ModelConfig').GetGlobalObservables():
   v=snapshot.find(g.GetName());assert v;g.setVal(v.getVal())
  assert all(not w.var(f'Ymean_{p}_{era}').isConstant() for p in ['upsilon1s','upsilon2s'] for era in ['2022','2022EE','2023','2023BPix'])
  w.writeToFile(str(ws));f.Close();attempts=[]
  for label,start in [('zero',0.),('expected',expected)]:
   ad=td/label
   cmd=[P/'bin/combine','-M','FitDiagnostics',ws,'-D','toy_data','--skipBOnlyFit','--keepFailures',
        '--robustFit','0','--minos','none','--rMin',str(-bound),'--rMax',str(bound),'--preFitValue',str(start),
        '-m',str(a.mass),'-n','_c0_toy']+common
   rc=native(cmd,ad);fr=ROOT.TFile.Open(str(ad/'fitDiagnostics_c0_toy.root'));wf=ROOT.TFile.Open(str(ws));w=wf.Get('w')
   fit=fr.Get('fit_s') if fr and not fr.IsZombie() else None
   if not fit:attempts.append(dict(anchor=label,fit_valid=False,reason='missing fit_s'));wf.Close();continue
   apply_fit(w,fit);metric,_=inspect_model(w,'s',data_name='toy_data');p=fit.floatParsFinal().find('r');cont=continuity(w,'s')
   row=dict(anchor=label,status=fit.status(),covQual=fit.covQual(),edm=fit.edm(),r=p.getVal(),error=p.getError(),
            deviance=metric['exact_deviance'],C0_passed=cont['maximum_relative_difference']<1e-8,
            r_boundary=abs(p.getVal())>.9999*bound,minimum_continuum=min(polynomial_minimum(w).values()))
   row['fit_valid']=bool(rc['exit_code']==0 and row['status']==0 and 0<=row['edm']<.01 and row['C0_passed']
          and not row['r_boundary'] and row['minimum_continuum']>0 and math.isfinite(row['deviance']))
   attempts.append(row);fr.Close();wf.Close()
  good=[x for x in attempts if x.get('fit_valid')];finite=[x for x in attempts if x.get('deviance') is not None and math.isfinite(x['deviance'])]
  selected=dict(min(good or finite,key=lambda x:x['deviance'])) if finite else dict(fit_valid=False,r=None)
  selected.update(toy=i,r_true=rtrue,attempts=attempts)
  if a.method=='gof':
   value=score(ws,td/'gof','toy_data');selected['gof_score']=value
   selected['gof_fd_agreement']=abs(value-selected['deviance'])<.02 if value is not None and selected.get('deviance') is not None else False
  # Keep failed attempts and all finite r values; never replace a toy.
  from fit_helpers import dump
  dump(result,selected);rows.append(json.loads(result.read_text()))
  summary=summarize(rows,rtrue,observed);save(out/'summary.json',dict(plan=plan,summary=summary,rows=rows))
  print(i,a.toys,selected.get('r'),selected['fit_valid'],flush=True)
 fh.Close();assert sha(source)==plan['generator_sha256'];assert sha(d/'workspace_postfit_b.root')==plan['template_sha256']
 print(out/'summary.json')
if __name__=='__main__':main()
