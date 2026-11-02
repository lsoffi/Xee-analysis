"""Independent fixed-r MultiDimFit checks of one native robust68 interval."""
from pathlib import Path
import argparse, json, ROOT
from worker_helpers import apply_result, execute, atomic_json
ROOT.gROOT.SetBatch(True); ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
P=Path(__file__).resolve().parents[2]
p=argparse.ArgumentParser();p.add_argument('result');a=p.parse_args()
rp=Path(a.result).resolve();q=json.loads(rp.read_text());d=rp.parent/'profile_endpoint_audit';d.mkdir(exist_ok=True)
f=ROOT.TFile.Open(str(Path(q['directory'])/'workspace.root'));w=f.Get('w');ff=ROOT.TFile.Open(str(Path(q['directory'])/'fitDiagnostics_upper.root'));fit=ff.Get('fit_s');apply_result(w,fit);w.var('mass').setVal(11.2);w.var('mass_CR').setVal(9.8);w.writeToFile(str(d/'workspace.root'));ff.Close();f.Close()
results=[]
for label,value in [('low',q['interval_low']),('high',q['interval_high'])]:
 wd=d/label;wd.mkdir(exist_ok=True)
 cmd=[str(P/'bin/combine'),'-M','MultiDimFit',str(d/'workspace.root'),'-D','toy_data','-m',str(q['mass']),'-n','_endpoint','--algo','fixed','--fixedPointPOIs','r='+str(value),'--saveFitResult','--cminDefaultMinimizerStrategy','2','--cminDefaultMinimizerTolerance','0.00001','--cminDefaultMinimizerPrecision','1e-12','--cminRunAllDiscreteCombinations','--X-rtd','MINIMIZER_freezeDisassociatedParams','-v','0']
 receipt=execute(cmd,wd,180)
 result=dict(label=label,fixed_r=value,receipt=receipt)
 fs=list(wd.glob('higgsCombine*.root'))
 if len(fs)==1:
  h=ROOT.TFile.Open(str(fs[0]));tr=h.Get('limit');result['points']=[dict(r=float(row.r),deltaNLL=float(row.deltaNLL),quantile=float(row.quantileExpected)) for row in tr];h.Close()
 results.append(result)
 print(label,result,flush=True)
threshold=float(.5*ROOT.Math.chisquared_quantile_c(1-.68,1))
atomic_json(d/'audit.json',dict(target_deltaNLL=threshold,source=str(rp),source_r=q['r'],native_interval=[q['interval_low'],q['interval_high']],results=results))
