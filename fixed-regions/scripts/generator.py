#!/usr/bin/env python3
"""Prepare a candidate alternative Chebyshev SR generator, same CR/C0/DCBs.

This candidate must pass its own GoF and stability checks. Building/fitting it
is not evidence that it is an adequate generator. Historical envelope generator
scripts remain separate. Dry-run is the default; no minimization without --execute.
"""
from pathlib import Path
import argparse,hashlib,json,sys,subprocess
P=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(P/'analysis/separate_sr_c0_limits'))

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--side',required=True,choices=['low','high']);ap.add_argument('--mass',type=float,required=True)
 ap.add_argument('--degree',type=int,required=True);ap.add_argument('--tag',required=True)
 ap.add_argument('--execute',action='store_true');ap.add_argument('--build-only',action='store_true')
 a=ap.parse_args();directory=P/'runs/generators'/a.tag/a.side/f'm{a.mass:.1f}'
 if not 1<=a.degree<=10:ap.error('degree must be in 1..10')
 protocol=dict(side=a.side,mass=a.mass,degree=a.degree,CR_degree=4,C0=True,
       fit='B-only, eight Y means floating',starts=3,output=str(directory),candidate_not_validated=True)
 print(json.dumps(protocol,indent=2))
 if not(a.execute or a.build_only):return
 import ROOT
 import model
 from fit_helpers import read_stage,dump,polynomial_minimum
 from validation import native
 ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
 directory.mkdir(parents=True,exist_ok=True)
 if (directory/'protocol.json').exists():
  assert json.loads((directory/'protocol.json').read_text())==protocol,'Changed generator protocol: use a new tag'
 else:dump(directory/'protocol.json',protocol)
 attempts=[]
 for variant in range(3):
  wd=directory/f'start{variant}';path=wd/'workspace.root'
  if not path.exists():
   model.build(a.side,a.mass,path,degree_override=a.degree)
   if variant:
    f=ROOT.TFile.Open(str(path));w=f.Get('w')
    for era in model.ERAS:
     for k in [a.degree-1,a.degree]:
      if k>0:
       v=w.var(f'local_c{k}_{era}');v.setVal(v.getVal()+(1 if variant==1 else -1)*.0005)
    assert min(polynomial_minimum(w).values())>0
    copy=w.Clone('w');f.Close();copy.writeToFile(str(path))
    audit=path.with_suffix('.audit.json');meta=json.loads(audit.read_text())
    meta.update(workspace_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),generator_start_variant=variant)
    dump(audit,meta)
  if a.build_only:continue
  command=[P/'bin/combine','-M','FitDiagnostics',path,'-m',str(a.mass),'-n','_generator',
    '--skipSBFit','--keepFailures','--minos','none','--robustFit','0',
    '--cminDefaultMinimizerStrategy','2','--cminDefaultMinimizerTolerance','.0001',
    '--cminDefaultMinimizerPrecision','1e-12','--X-rtd','MINIMIZER_freezeDisassociatedParams']
  receipt=native(command,wd)
  q,_=read_stage(path,wd/'fitDiagnostics_generator.root','b');q.update(receipt=receipt,variant=variant)
  q['generator_valid']=bool(q.get('point_valid') and q.get('covQual')==3 and receipt['exit_code']==0)
  dump(wd/'result.json',q);attempts.append(q)
 if a.build_only:return
 valid=[q for q in attempts if q['generator_valid']]
 dump(directory/'attempts.json',attempts)
 if len(valid)<2:raise RuntimeError('Need at least two converged generator starts with covariance quality 3. Inspect retained attempts.')
 best=min(valid,key=lambda q:q['exact_deviance'])
 ws=Path(best['source']).parent/'workspace.root'
 q,_=read_stage(ws,best['source'],'b',directory/'workspace_postfit_b.root')
 q.update(degree=a.degree,CR_degree=4,candidate_not_validated=True,
          start_deviance_spread=max(x['exact_deviance'] for x in valid)-min(x['exact_deviance'] for x in valid))
 dump(directory/'result.json',q);print(directory/'workspace_postfit_b.root')
if __name__=='__main__':main()
