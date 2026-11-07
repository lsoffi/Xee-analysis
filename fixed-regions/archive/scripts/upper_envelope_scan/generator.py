from pathlib import Path
import ROOT,json,sys,subprocess
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
P=Path.cwd();R=P/'runs/20261030_upper_envelope_scan';mode=sys.argv[1];degree=int(sys.argv[2]);mass=sys.argv[3] if len(sys.argv)>3 else '10.8';rtrue=float(sys.argv[4]) if len(sys.argv)>4 else 0
D=R/f'gof/degree{degree}/generation' if mode=='gof' else R/f'bias/degree{degree}/m{mass}/{sys.argv[5]}/generation';D.mkdir(parents=True,exist_ok=True)
source=R/f'degree{degree}/m{mass}/joint'/('workspace_postfit_b.root' if mode=='gof' else 'workspace_seeded.root');f=ROOT.TFile.Open(str(source));w=f.Get('w');keep=[]
if mode=='bias':
    pars=json.loads((P/'runs/20261030_upper_sr_cr_mu_fixed/degree3/results.json').read_text())['parameters'];vals={a['name']:a['value'] for a in pars};reps=[]
    for era in ['2022','2022EE','2023','2023BPix']:
        coeff=[ROOT.RooConstVar(f'generator_c{j}_{era}','',vals[f'local_c{j}_{era}']) for j in range(1,4)]
        pdf=ROOT.RooChebychev('generator_cheb3_'+era,'',w.var('mass'),ROOT.RooArgList(*coeff));keep.extend(coeff+[pdf]);getattr(w,'import')(pdf);reps.append(f'local_envelope_{era}='+pdf.GetName())
    for k,v in vals.items():
        if w.var(k) and (k.endswith('_CR') or k.startswith('local_yield_') or 'bkgScale_upsilon' in k):w.var(k).setVal(v)
    for mcname in ['ModelConfig','ModelConfig_bonly']:
        mc=w.obj(mcname);new=w.factory(f'EDIT::{mcname}_generator('+mc.GetPdf().GetName()+','+','.join(reps)+')');assert new;mc.SetPdf(new)
for c in w.allCats():
    if c.GetName().startswith('pdfindex'):c.setConstant(True)
w.genobj('discreteParams').removeAll()
for v in w.allVars():
    if v.GetName().startswith(('local_','Ymean_','CMS_EXO25020_bkgScale_')):v.setConstant(True)
w.var('r').setVal(rtrue);w.var('r').setRange(-10000,10000);w.var('mass').setVal(11.2);w.saveSnapshot('generation',w.allVars());w.writeToFile(str(D/'generator.root'));f.Close()
seed=731000+degree*1000+round(float(mass)*10)+round(rtrue*13);n='100' if mode=='gof' else '-1'
cmd=[str(P/'bin/combine'),'-M','GenerateOnly',str(D/'generator.root'),'--snapshotName','generation','--saveToys','--toysFrequentist','--bypassFrequentistFit','--expectSignal',str(rtrue),'--genBinnedChannels',','.join(f'{reg}_ch{i}' for reg in ['SR','CR'] for i in range(1,5)),'-t',n,'-s',str(seed),'-m',mass,'-n','_generated']
(D/'command.json').write_text(json.dumps(cmd,indent=2))
with (D/'generate.log').open('w') as h:rc=subprocess.run(cmd,cwd=D,stdout=h,stderr=subprocess.STDOUT,timeout=180).returncode
assert rc==0;(D/'truth.json').write_text(json.dumps(dict(mode=mode,degree=degree,mass=float(mass),rtrue=rtrue,seed=seed,generator='matched B fit' if mode=='gof' else 'alternative Chebyshev3 fitted to data')))
print(D,flush=True)
