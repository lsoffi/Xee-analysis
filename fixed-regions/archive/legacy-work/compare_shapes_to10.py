import numpy as np,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
p=Path('outputs');cases=[('free','Liberi'),('limit1','Limiti ±1σ'),('limit2','Limiti ±2σ')];components=[('dy','Chebyshev'),('upsilon1s','Υ(1S)'),('upsilon2s','Υ(2S)'),('signal','Segnale a 9 GeV')]
fig,axs=plt.subplots(4,3,figsize=(16,12),layout='constrained',sharex=True);fr,ars=plt.subplots(4,3,figsize=(16,11),layout='constrained',sharex=True);summ=[]
for j,(case,title) in enumerate(cases):
 pref='combine_2022_degree7_to10_'+case;ds={s:np.genfromtxt(p/f'{pref}_{s}_shapes.csv',delimiter=',',names=True) for s in ['prefit','fit_b','fit_s']};a=ds['prefit'];x=a['mass'];ed=np.r_[a['low'],a['high'][-1]]
 for i,(key,label) in enumerate(components):
  ax=axs[i,j];ar=ars[i,j];positive=[];sums={s:float(ds[s][key].sum()) for s in ds};summ.append(dict(case=case,component=key,**sums))
  for st,c,ls,lab in [('prefit','#424242','--','Pre-fit'),('fit_b','#1766ba',':','Post-fit B'),('fit_s','#cc2d2d','-','Post-fit S+B')]:
   y=ds[st][key];positive.extend(y[y>0]);ax.stairs(np.where(y>0,y,np.nan),ed,color=c,ls=ls,lw=1.6,label=lab,baseline=None)
   if y.sum()>0:
    base=ds['prefit'][key];ratio=(y/y.sum())/(base/base.sum());ar.plot(x,ratio,color=c,ls=ls,lw=1.5,label=lab)
  ax.set_yscale('log');ax.set_ylim(min(positive)*.6,max(positive)*1.8);ax.set_title(label+' — '+title,fontsize=10);ax.grid(alpha=.13,which='both');ax.legend(fontsize=7)
  ar.set_title(label+' — '+title,fontsize=10);ar.grid(alpha=.2);ar.legend(fontsize=7)
  if key in ['upsilon1s','upsilon2s']:ar.set_ylim(.98,1.02)
  if key=='signal':ax.text(.03,.05,'B: segnale nullo',transform=ax.transAxes,color='#1766ba',fontsize=8)
  if j==0:ax.set_ylabel('Eventi / 0.01333 GeV');ar.set_ylabel('Shape post/pre\n(integrali a 1)')
  if i==3:ax.set_xlabel('Massa [GeV]');ar.set_xlabel('Massa [GeV]')
  ax.set_xlim(6,10);ar.set_xlim(6,10)
fig.suptitle('Fit 6–10 GeV — ogni componente pre/post, Chebyshev grado 7, binning originale\nCurve con normalizzazioni del fit; varie configurazioni hanno covQual 2',fontsize=13)
fr.suptitle('Fit 6–10 GeV — cambiamento della sola forma\nRapporti post/pre dopo normalizzazione a integrale 1; non includono bande di incertezza',fontsize=13)
for ext in ['png','pdf']:
 fig.savefig(p/f'degree7_to10_componenti_pre_post.{ext}',dpi=160);fr.savefig(p/f'degree7_to10_componenti_shape_ratio.{ext}',dpi=160)
(p/'degree7_to10_componenti_integrali.json').write_text(json.dumps(summ,indent=2));print(json.dumps(summ,indent=2))
# Chebyshev comparison between ranges for the same bounded/free case, common region.
f2,ax2=plt.subplots(1,3,figsize=(15,4.4),layout='constrained')
for ax,(case,title) in zip(ax2,cases):
 for rng,color,ls in [('native','#1766ba','--'),('to10','#cc2d2d','-')]:
  a=np.genfromtxt(p/f'combine_2022_degree7_{rng}_{case}_fit_s_shapes.csv',delimiter=',',names=True);ax.plot(a['mass'],a['dy'],color=color,ls=ls,lw=1.6,label='Fit 6–12' if rng=='native' else 'Fit 6–10')
 ax.set_yscale('log');ax.set_xlim(6,10);ax.set_title('Chebyshev S+B — '+title);ax.set_xlabel('Massa [GeV]');ax.set_ylabel('Eventi / 0.01333 GeV');ax.legend(fontsize=8);ax.grid(alpha=.15,which='both')
f2.suptitle('Confronto del Chebyshev post-fit tra i due range nella regione comune 6–10',fontsize=13)
for ext in ['png','pdf']:f2.savefig(p/f'degree7_chebyshev_confronto_range_6_10.{ext}',dpi=170)
