from pathlib import Path
import json,time,traceback
import numpy as np
import ROOT,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path('/eos/user/s/soffi/www/Xee/261106/scan_SR_CR/low_logpoly2');out.mkdir(exist_ok=True)
records=json.loads(Path('records.json').read_text());eras=['2022','2022EE','2023','2023BPix']
def draw(rec,path,dest):
 f=ROOT.TFile.Open(str(path));fit=f.Get('fit_s');regs=['SR0']+['CR_'+g for g in rec['selected_CRs']];fig,axs=plt.subplots(4,len(regs),figsize=(6*len(regs),13),squeeze=False);fig.suptitle(f"S+B · m = {rec['mass']:.1f} GeV · SR exp(pol. grado 2) · C0 ai bordi comuni · status {fit.status()} · covQual {fit.covQual()}",fontsize=16)
 for i,era in enumerate(eras):
  for j,reg in enumerate(regs):
   lo,hi=rec['sr_segments'][0] if j==0 else {'omega_phi':(.650,1.155),'charmonium':(2.8,4.),'upsilon':(9.095,10.505)}[reg[3:]];n=round((hi-lo)*200);x=lo+(np.arange(n)+.5)/200;ax=axs[i,j];d=f.Get('shapes_fit_s/'+reg+'_'+era);g=d.Get('data');y=np.array([g.GetPointY(k) for k in range(n)]);ax.errorbar(x,y,yerr=np.sqrt(y),fmt='.',ms=2,color='black',label='Dati');
   for name in ['continuum','omega','phi','jpsi','psi2s','upsilon1s','upsilon2s','Zd','total']+ [key.GetName() for key in d.GetListOfKeys() if 'ups' in key.GetName().lower() and key.GetName() not in ['upsilon1s','upsilon2s']]:
    h=d.Get(name)
    if not h:continue
    v=np.array([h.GetBinContent(k+1) for k in range(n)]);ax.step(x,v,where='mid',lw=1.5,label={'total':'Totale','Zd':'Segnale','continuum':'Continuo'}.get(name,name));
    if name=='total':
     chi=np.sum((y-v)**2/v) if np.all(v>0) else float('nan');ax.set_title(f'{era} · {reg} · χ² Pearson = {chi:.1f} ({n} bin)')
   ax.set(xlabel='Massa [GeV]',ylabel='Eventi / 5 MeV',xlim=(lo,hi));ax.legend(fontsize=7,ncol=3)
 fig.tight_layout(rect=[0,0,1,.96]);tmp=dest.with_suffix('.tmp.png');fig.savefig(tmp,dpi=120);tmp.replace(dest);plt.close(fig);f.Close()

import sys,shutil,os
m=sys.argv[1];d=Path(sys.argv[2]);dest=out/('m'+m+'.png')
if not dest.exists():
 rec=next(r for r in records if r['mass']==float(m))
 local=Path(os.environ['_CONDOR_SCRATCH_DIR'])/'early.png'
 draw(rec,next(d.glob('fitDiagnostics*.root')),local)
 tmp=out/('.early_m'+m+'.png');shutil.copyfile(local,tmp)
 if not dest.exists():tmp.replace(dest)
 else:tmp.unlink()
