from pathlib import Path
import json,hashlib,numpy as np,ROOT
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
P=Path.cwd();R=P/'runs/20261030_gof_thesis48';reference={};reference_parameters={};out=[]
for t in json.loads((R/'tasks.json').read_text()):
 d=Path(t['path']);a=json.loads((d/'audit.json').read_text());f=ROOT.TFile.Open(str(d/'workspace.root'));w=f.Get('w');w.var('r').setVal(0);vec=[]
 # Compare model definitions at shared background parameter values, not different fit seeds.
 pars={v.GetName():v.getVal() for v in w.allVars() if not v.isConstant() and v.GetName() not in ['mass','mass_CR','r','evt_weight']}
 shared=reference_parameters.setdefault(t['era'],pars)
 for name,value in shared.items():
  if w.var(name):w.var(name).setVal(value)
 w.cat('pdfindex_'+t['era']).setIndex(0)
 for idx in [0,1]:
  w.cat('pdfindex_'+t['era']).setIndex(idx)
  for ch in a['channels']:
   name='mass_CR' if ch.startswith('CR_') else 'mass';x=w.var(name);obs=ROOT.RooArgSet(x);b=w.obj('ModelConfig_bonly').GetPdf().getPdf(ch);s=w.obj('ModelConfig').GetPdf().getPdf(ch)
   for j in range(x.getBins()):
    x.setBin(j);vb=b.getVal(obs)*b.expectedEvents(obs);vs=s.getVal(obs)*s.expectedEvents(obs);assert abs(vb-vs)<1e-7*max(1,abs(vb));assert vb>0;vec.append(vb)
   for pr in ['upsilon1s','upsilon2s']:assert b.dependsOn(w.var(f'Ymean_{pr}_{t["era"]}'))
 ref=reference.setdefault(t['era'],np.array(vec));delta=float(np.max(np.abs(np.array(vec)-ref)/np.maximum(1,np.abs(ref))));
 assert delta<1e-8,(t,delta)
 # Only one discrete envelope index, shared between SR+CR through physical Y parameters.
 assert w.genobj('discreteParams').getSize()==1
 assert hashlib.sha256(Path(a['source']).read_bytes()).hexdigest()==a['source_sha256']
 out.append(dict(era=t['era'],mass=t['mass'],max_relative_background_difference_at_shared_parameters=delta,SB_at_r0_equals_B=True,positive_bin_means=True,source_unchanged=True));f.Close()
(R/'validation.json').write_text(json.dumps(out,indent=2));print('audited',len(out),'max relative model difference',max(v['max_relative_background_difference_at_shared_parameters'] for v in out))
