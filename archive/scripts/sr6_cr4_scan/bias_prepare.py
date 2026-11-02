from pathlib import Path
import ROOT,json,sys,numpy as np
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=Path.cwd();R=P/'runs/20261030_sr6_cr4_scan';mass=sys.argv[1];family=sys.argv[2];d=R/'bias'/family/('m'+mass);d.mkdir(parents=True,exist_ok=True);source=R/('m'+mass)/'joint/workspace.root';f=ROOT.TFile.Open(str(source));w=f.Get('w');keep=[];rep=[]
folder='20261030_sr_bernstein7_bias' if family=='bern7' else '20261030_sr_cheb7_bern_bias';rows=json.loads((P/'runs'/folder/'generators.json').read_text())
for era in ['2022','2022EE','2023','2023BPix']:
 a=next(a for a in rows if a['degree']==family and a['era']==era);p=np.array(a['p']);nt=8;nc=13
 if family=='bern7':
  coeff=[ROOT.RooConstVar(f'generator_{family}_{era}_{j}','',float(p[j])) for j in range(8)];pdf=ROOT.RooBernstein('generator_'+era,'',w.var('mass'),ROOT.RooArgList(*coeff));yieldSR=sum(p[:8])/8
 else:
  coeff=[ROOT.RooConstVar(f'generator_{family}_{era}_{j}','',float(p[j]/p[0])) for j in range(1,8)];pdf=ROOT.RooChebychev('generator_'+era,'',w.var('mass'),ROOT.RooArgList(*coeff));ci=np.polynomial.chebyshev.chebint(p[:8]);yieldSR=(np.polynomial.chebyshev.chebval(1,ci)-np.polynomial.chebyshev.chebval(-1,ci))/2
 keep.extend(coeff+[pdf]);getattr(w,'import')(pdf);rep.append('local_envelope_'+era+'='+pdf.GetName());w.var('local_yield_'+era).setVal(yieldSR)
 co=p[8:13];ci=np.polynomial.chebyshev.chebint(co);w.var('local_yield_'+era+'_CR').setVal((np.polynomial.chebyshev.chebval(1,ci)-np.polynomial.chebyshev.chebval(-1,ci))/2)
 for j in range(1,5):w.var(f'local_c{j}_{era}_CR').setVal(co[j]/co[0])
 for j,pr in enumerate(['upsilon1s','upsilon2s']):w.var(f'CMS_EXO25020_bkgScale_{pr}_{era}').setVal(p[nc+j]);w.var(f'Ymean_{pr}_{era}').setVal(p[nc+2+j])
for mcname in ['ModelConfig','ModelConfig_bonly']:
 mc=w.obj(mcname);new=w.factory(f'EDIT::{mcname}_generator('+mc.GetPdf().GetName()+','+','.join(rep)+')');assert new;mc.SetPdf(new)
for c in w.allCats():
 if c.GetName().startswith('pdfindex'):c.setConstant(True)
w.genobj('discreteParams').removeAll()
for v in w.allVars():
 if v.GetName().startswith(('local_','Ymean_','CMS_EXO25020_bkgScale_')):v.setConstant(True)
w.var('r').setVal(0);w.saveSnapshot('generation',w.allVars());w.writeToFile(str(d/'generator.root'));print('generator',mass,family,flush=True)
