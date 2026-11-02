from pathlib import Path
import ROOT,json,hashlib,numpy as np
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=Path.cwd();R=P/'runs/20261030_sr6_cr4_scan';rows=[];reference_data=None
for result in sorted(R.glob('m*/results.json')):
 d=result.parent;a=json.loads(result.read_text());audit=json.loads((d/'joint/audit.json').read_text());inp=json.loads((d/'signal/audit.json').read_text());f=ROOT.TFile.Open(str(d/'joint/workspace.root'));w=f.Get('w');sourcefile=ROOT.TFile.Open(inp['source']);sourcew=sourcefile.Get('w');eras=['2022','2022EE','2023','2023BPix'];checks=[]
 for i,e in enumerate(eras,1):
  for pr in ['upsilon1s','upsilon2s']:
   v=w.var(f'Ymean_{pr}_{e}');sg=w.function(f'Ysigma_{pr}_{e}');assert v and not v.isConstant() and sg.isConstant();original=sourcew.pdf(f'shapeBkg_{pr}_ch{i}');assert abs(sg.getVal()-original.sigmaL().getVal())<1e-12
   for j,tail in enumerate([original.alphaL(),original.nL(),original.alphaR(),original.nR()]):
    val=w.function(f'Ytail_{pr}_{e}_{j}');assert val.isConstant() and abs(val.getVal()-tail.getVal())<1e-12
   assert w.obj('ModelConfig').GetPdf().getPdf('SR_ch'+str(i)).dependsOn(v) and w.obj('ModelConfig').GetPdf().getPdf('CR_ch'+str(i)).dependsOn(v)
  mp=w.pdf('local_envelope_'+e);idx=w.cat('pdfindex_'+e);vals=[]
  for j in [0,1]:idx.setIndex(j);vals.append(mp.getCorrection())
  assert abs(vals[0]-vals[1])<1e-12;checks.append(dict(era=e,equal_penalties=vals))
 ds=w.genobj('discreteParams');assert ds and ds.getSize()==4;assert w.var('mass').getMax()<=w.var('mass_CR').getMin()+1e-8
 dat=w.data('data_obs');content=[]
 for k in range(dat.numEntries()):
  row=dat.get(k);cat=row.find('CMS_channel').getLabel();obsname='mass_CR' if cat.startswith('CR_') else 'mass';content.append((cat,round(row.getRealValue(obsname),12),round(dat.weight(),10)))
 dhash=hashlib.sha256(json.dumps(sorted(content)).encode()).hexdigest()
 if reference_data is None:reference_data=dhash
 assert dhash==reference_data,(a['mass'],'different data')
 unchanged=hashlib.sha256(Path(inp['source']).read_bytes()).hexdigest()==inp['source_sha256'];assert unchanged
 pars=json.loads((d/'joint/quality.json').read_text());assert all(p['continuum_min']>0 for p in pars.values())
 keys=['expected_2p5','expected_16','expected_50','expected_84','expected_97p5'];ordered=all(a['limits'][u]<a['limits'][v] for u,v in zip(keys,keys[1:]));rows.append(dict(mass=a['mass'],source_unchanged=unchanged,data_hash=dhash,identical_data_across_masses=True,disjoint=True,shared_mu=True,fixed_sigmas=True,discrete_indices=4,penalties=checks,ordered_limits=ordered,numerical_ok=a['numerical_ok'],CLs_verified=a['limit_crossing']['verified']));sourcefile.Close();f.Close()
(R/'validation.json').write_text(json.dumps(rows,indent=2));print('Audited',len(rows),'points;',sum(a['numerical_ok'] for a in rows),'numerically checked')
