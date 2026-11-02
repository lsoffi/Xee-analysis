import json,csv,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
p=Path('outputs');bounds=json.load(open(p/'degree7_native_bounds.json'));tags=['degree7_native_free','degree7_native_limit1','degree7_native_limit2'];labels=['Coefficienti liberi','Limiti rigidi ±1σ pre-fit','Limiti rigidi ±2σ pre-fit'];summary=[];params=[]
fig,axs=plt.subplots(3,2,figsize=(13,12),layout='constrained',sharex=True,sharey=True)
for i,(tag,title) in enumerate(zip(tags,labels)):
 m=json.load(open(p/f'combine_2022_{tag}_pre_post.json'));r=m['parameters'][0]
 row=dict(case=title,r=r['s'],r_lo=r['s_lo'],r_hi=r['s_hi'],status_b=m['b']['status'],covQual_b=m['b']['covQual'],status_s=m['s']['status'],covQual_s=m['s']['covQual'])
 for a in m['parameters']:
  a=dict(a);a['case']=title
  if a['parameter'].startswith('t'):
   a['shift_b_in_saved_sigma']=(a['b']-a['prefit'])/a['prefit_saved_error'];a['shift_s_in_saved_sigma']=(a['s']-a['prefit'])/a['prefit_saved_error']
  params.append(a)
 for j,(stage,stitle) in enumerate([('fit_b','Solo fondo'),('fit_s','Segnale + fondo')]):
  ax=axs[i,j];a=np.genfromtxt(p/f'combine_2022_{tag}_{stage}_shapes.csv',delimiter=',',names=True);x=a['mass'];n=a['data'];mu=a['total'];ed=np.r_[a['low'],a['high'][-1]];scale=(1/75)/a['width'];row[stage+'_deviance']=float(2*np.sum(mu-n+n*np.log(n/mu)))
  ax.errorbar(x,n*scale,yerr=np.sqrt(n)*scale,fmt='.',color='black',ms=4,lw=.7,label='Dati 2022')
  for key,c,lab,ls in [('total','#c62828','Totale','-'),('dy','#1565c0','Chebyshev grado 7','--'),('upsilon1s','#ef8a00','Υ(1S), norma fissa','-.'),('upsilon2s','#9c27b0','Υ(2S)',':'),('signal','#188341','Segnale 9 GeV','-')]:
   y=a[key]*scale
   if key=='signal' and np.any(y<0):y=np.abs(y);lab='|Segnale|, r negativo';ls='--'
   if np.any(y>0):ax.stairs(np.where(y>0,y,np.nan),ed,color=c,lw=1.5,ls=ls,label=lab,baseline=None)
  ax.set_yscale('log');ax.set_xlim(6,12);ax.set_ylim(.001,10000);ax.grid(alpha=.12,which='both');ax.set_title(title+' — '+stitle,fontsize=10)
  ax.legend(fontsize=7,loc='lower left',ncol=2)
  if stage=='fit_s':ax.text(.04,.94,f'r = {r["s"]:.3f} −{abs(r["s_lo"]):.3f}/+{r["s_hi"]:.3f}',transform=ax.transAxes,va='top',fontsize=9,bbox=dict(facecolor='white',alpha=.9,edgecolor='none'))
  if j==0:ax.set_ylabel('Eventi / 0.01333 GeV')
  if i==2:ax.set_xlabel('Massa [GeV]')
 summary.append(row)
fig.suptitle('2022 — fit 6–12 GeV, 450 bin originali, Chebyshev grado 7\nConfronto libero / limiti pre-fit ±1σ / ±2σ; stessi vincoli e opzioni Combine',fontsize=13)
for ext in ['png','pdf']:fig.savefig(p/f'degree7_native_bounded_fits.{ext}',dpi=170)
