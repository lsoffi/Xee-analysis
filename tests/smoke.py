#!/usr/bin/env python3
"""Structural/numerical smoke tests; builds/inspects models, performs no fit."""
from pathlib import Path
import argparse,importlib.util,json,sys,hashlib
import numpy as np
P=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(P/'analysis/separate_sr_c0_limits'))
import ROOT
from fit_helpers import inspect_model,continuity,polynomial_minimum
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--reference-project',type=Path);a=ap.parse_args()
 report=dict(scope='No minimization. Syntax, imported input hashes, parameter structure, physical C0, prediction closure and optional original-model parity.',models=[])
 for folder in ['analysis','scripts']:
  for path in (P/folder).rglob('*.py'):compile(path.read_text(),str(path),'exec')
 report['syntax']=True
 for side,expected,npar in [('low',314857,40),('high',52105,32)]:
  path=P/'runs/20261031_separate_sr_c0_limits'/side/'m9.5/workspace.root'
  if not path.exists():raise RuntimeError('Build both m9.5 workspaces first.')
  f=ROOT.TFile.Open(str(path));w=f.Get('w');audit=json.loads(path.with_suffix('.audit.json').read_text())
  assert hashlib.sha256(Path(audit['source']).read_bytes()).hexdigest()==audit['source_sha256']
  assert w.data('data_obs').sumEntries()==expected
  coefficients=[v for v in w.allVars() if v.GetName().startswith('local_c')]
  assert len(coefficients)==npar and all(not v.isConstant() for v in coefficients)
  assert len([v for v in w.allVars() if v.GetName().startswith('Ymean_') and not v.isConstant()])==8
  assert all(w.function('local_yield_'+era) and w.var('local_yield_'+era+'_CR') for era in ['2022','2022EE','2023','2023BPix'])
  assert not any(p.InheritsFrom('RooMultiPdf') for p in w.allPdfs())
  cont=continuity(w,'b');assert cont['maximum_relative_difference']<1e-8
  metric,bins=inspect_model(w,'b');assert min(polynomial_minimum(w).values())>0
  record=dict(side=side,data_events=expected,coefficients=len(coefficients),C0max=cont['maximum_relative_difference'])
  if a.reference_project:
   rf=ROOT.TFile.Open(str(a.reference_project/'runs/20261031_separate_sr_c0_limits'/side/'m9.5/workspace.root'));rw=rf.Get('w')
   oldmetric,oldbins=inspect_model(rw,'b')
   differences=[float(np.max(np.abs(np.array(v['total'])-oldbins[k]['total']))) for k,v in bins.items()]
   assert max(differences)<1e-8
   record['original_prediction_max_absolute_difference']=max(differences);rf.Close()
  f.Close();report['models'].append(record)
 # Test the exact sample-statistic convention and failed-toy preservation.
 spec=importlib.util.spec_from_file_location('validation',P/'scripts/validation.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
 rows=[dict(r=x,fit_valid=True,error=1,covQual=3) for x in [0.,1.,2.]]+[dict(r=100.,fit_valid=False,error=1,covQual=2)]
 q=m.summarize(rows,1.)
 assert q['bias']['mean']==0 and q['bias']['sd']==1 and q['n_failed']==1
 assert q['bias_all_finite']['n']==4 and q['coverage_measured'] is False
 report['summary_conventions']=True;report['passed']=True
 out=P/'docs/release_smoke_test.json';out.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
