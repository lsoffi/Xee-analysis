import csv,json,math,time
from pathlib import Path
import ROOT
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT.gROOT.SetBatch(True)
run=Path.cwd()
web=Path('/eos/user/s/soffi/www/Xee/261106/scan_SR_CR/bias_low20')
web.mkdir(exist_ok=True)
jobs=[line.split() for line in (run/'points.txt').read_text().splitlines()]
rows=[]
for m,g,inj,seed in jobs:
 for model in ('cheb','log'):
  p=run/'results'/seed/(model+'.root')
  f=ROOT.TFile.Open(str(p));t=f.Get('tree_fit_sb');assert t and t.GetEntries()==1,p
  t.GetEntry(0);r=float(t.r)
  fit=f.Get('fit_s')
  cq=int(fit.covQual()) if fit else -1
  rows.append(dict(m=float(m),generator=g,model=model,inj=float(inj),seed=int(seed),r=r,err=float(t.rErr),bias=r-float(inj),status=int(t.fit_status),covQual=cq,good=int(t.fit_status)==0 and math.isfinite(r) and abs(r)<79.999))
  f.Close()
def save(name,data):
 with (web/name).open('w') as o:
  writer=csv.DictWriter(o,fieldnames=list(data[0]));writer.writeheader();writer.writerows(data)
save('toys.csv',rows)
groups=[]
for key in sorted({(r['m'],r['inj'],r['generator'],r['model']) for r in rows}):
 m,inj,g,model=key
 sample=[r for r in rows if (r['m'],r['inj'],r['generator'],r['model'])==key];assert len(sample)==20
 good=[r for r in sample if r['good']];b=np.array([r['bias'] for r in good])
 sd=b.std(ddof=1);groups.append(dict(m=m,inj=inj,generator=g,model=model,n=len(b),cov3=sum(r['covQual']==3 for r in good),mean_bias=b.mean(),sd=sd,mc_error=sd/np.sqrt(len(b)),normalized_bias=b.mean()/sd))
save('bias.csv',groups)
fig,axs=plt.subplots(1,2,figsize=(10,4.5),sharey=True)
for ax,signal in zip(axs,[False,True]):
 for k,(g,model,col) in enumerate([('cheb','cheb','#0072b2'),('cheb','log','#d55e00'),('log','cheb','#009e73'),('log','log','#a05195')]):
  subset=[r for r in groups if bool(r['inj'])==signal and r['generator']==g and r['model']==model]
  ax.errorbar([r['m']+(k-1.5)*.006 for r in subset],[r['mean_bias'] for r in subset],yerr=[r['mc_error'] for r in subset],fmt='o',color=col,label=f'{g} → {model}')
 ax.axhline(0,color='0.4',lw=1);ax.set(xticks=[.4,.5,.6],xlabel='Massa [GeV]',title='Iniezione al limite mediano atteso' if signal else 'Solo fondo')
axs[0].set_ylabel('Bias medio su r ± errore MC');axs[0].legend(fontsize=8);fig.tight_layout();fig.savefig(web/'bias.png',dpi=150)
for r in groups:print(r['m'],r['inj'],r['generator'],r['model'],r['n'],r['cov3'],round(r['mean_bias'],2),round(r['mc_error'],2),round(r['normalized_bias'],2))
