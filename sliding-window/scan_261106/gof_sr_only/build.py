import ROOT as R,json,sys
from pathlib import Path
R.gSystem.Load('libHiggsAnalysisCombinedLimit')
d=Path(__file__).resolve().parent;m=sys.argv[1]
src=d.parent/'coverage_all_env'/('m'+m)/'envelope.root'
if m=='1.4':src=d.parent/'coverage_m14_env/envelope.root'
f=R.TFile.Open(str(src));w=f.Get('w');cat=w.cat('CMS_channel');masks=R.RooArgList();fit=[];evaluate=[];labels=[]
for i in range(cat.numTypes()):
 cat.setBin(i);label=cat.getLabel();name='mask_'+label
 v=R.RooRealVar(name,name,0,0,1);v.setConstant(True);getattr(w,'import')(v)
 masks.add(w.var(name));fit.append(name+'=0');evaluate.append(name+('=1' if label.startswith('CR') else '=0'));labels.append(label)
for n in ('ModelConfig','ModelConfig_bonly'):
 pdf=w.obj(n).GetPdf();pdf.addChannelMasks(masks);assert pdf.channelMasks().getSize()==len(labels)
assert sum(x.startswith('SR') for x in labels)==4
p=d/('m'+m);p.mkdir(exist_ok=True)
(p/'fit.txt').write_text(','.join(fit));(p/'eval.txt').write_text(','.join(evaluate))
w.writeToFile(str(p/'workspace.root'))
(p/'audit.json').write_text(json.dumps({'source':str(src),'channels':labels,'fit':fit,'eval':evaluate},indent=2))
