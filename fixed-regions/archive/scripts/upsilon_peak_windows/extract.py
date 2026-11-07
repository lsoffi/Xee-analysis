from pathlib import Path
import ROOT,json,hashlib
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=Path.cwd();R=P/'runs/20261029_upsilon_peak_windows';base=P/'data/mass_scan/261028/cards_region2_data_envelope_allCorrections_binned/ee';rows=[];ys=[];hashes={};eras=['2022','2022EE','2023','2023BPix'];ref={};maxdiff=0
for d in sorted([d for d in base.iterdir() if d.name.replace('.','').isdigit()],key=lambda d:float(d.name)):
 m=float(d.name)
 if not 6.5<=m<=11:continue
 src=d/'Xee_ee_4_allYears.root';hashes[str(src)]=hashlib.sha256(src.read_bytes()).hexdigest();f=ROOT.TFile.Open(str(src));w=f.Get('w');x=w.var('mass');obs=ROOT.RooArgSet(x)
 for i,era in enumerate(eras,1):
  for proc in ['Zd','upsilon1s','upsilon2s']:
   pdf=w.pdf(('shapeSig_' if proc=='Zd' else 'shapeBkg_')+proc+f'_ch{i}');v=dict(mass=m,era=era,process=proc,low=x.getMin(),high=x.getMax(),class_name=pdf.ClassName())
   for key,method in [('mean','x0'),('sigmaL','sigmaL'),('sigmaR','sigmaR'),('alphaL','alphaL'),('alphaR','alphaR'),('nL','nL'),('nR','nR')]:v[key]=getattr(pdf,method)().getVal()
   v['checks']=[]
   for xx in [6.01,8.5,9.2,9.45,9.8,10.02,10.5,11.99]:x.setVal(xx);v['checks'].append([xx,pdf.getVal(obs)])
   if proc=='Zd':rows.append(v)
   else:
    key=(proc,era)
    if key not in ref:ref[key]=v;ys.append(v)
    else:
     diff=max(abs(v[k]-ref[key][k]) for k in ['mean','sigmaL','sigmaR','alphaL','alphaR','nL','nR']);maxdiff=max(diff,maxdiff);assert diff<1e-10
 f.Close()
(R/'parameters.json').write_text(json.dumps(dict(signal=rows,upsilon=ys),indent=2));(R/'provenance.json').write_text(json.dumps(dict(hashes=hashes,signal_shapes=len(rows),upsilon_shapes=len(ys),max_upsilon_parameter_change_across_mass_inputs=maxdiff),indent=2));print(len(rows),len(ys),maxdiff)
