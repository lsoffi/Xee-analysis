import numpy as np,json,csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
p=Path('outputs');refwidth=8/75
stages=['prefit','fit_b','fit_s'];sets=[];allshifts=[]
for tag,upper in [('rebin8',12),('to11_rebin8',11)]:
 pref='combine_2022_'+tag
 ds={stage:np.genfromtxt(p/f'{pref}_{stage}_shapes.csv',delimiter=',',names=True) for stage in stages};sets.append((tag,upper,ds))
 meta=json.loads((p/f'{pref}_pre_post.json').read_text())
 fig,axs=plt.subplots(1,2,figsize=(13,5.2),layout='constrained')
 for ax,stage,title in zip(axs,['fit_b','fit_s'],['Solo fondo','Segnale + fondo']):
  a=ds[stage];ed=np.r_[a['low'],a['high'][-1]];factor=refwidth/a['width'];n=a['data']
  ax.errorbar(a['mass'],n*factor,yerr=np.sqrt(n)*factor,fmt='.',color='black',ms=5,lw=.9,label='Dati 2022')
  for key,c,lab,ls in [('total','#c62828','Totale post-fit','-'),('dy','#1565c0','Chebyshev','--'),('upsilon1s','#ef8a00','Υ(1S), norma fissa','-.'),('upsilon2s','#9c27b0','Υ(2S)',':'),('signal','#188341','Segnale 9 GeV' if stage=='fit_s' else 'Segnale nullo (r = 0)','-')]:
   y=a[key]*factor;ax.stairs(np.where(y>0,y,np.nan),ed,color=c,lw=1.7,label=lab,ls=ls,baseline=None)
  ax.set_yscale('log');ax.set_xlim(8,10.5);ax.set_ylim(.001,2000);ax.set_xlabel('Massa [GeV]');ax.set_ylabel('Eventi / 0.10667 GeV');ax.grid(alpha=.15,which='both');ax.legend(fontsize=8,loc='lower left');q=meta['b' if stage=='fit_b' else 's'];ax.set_title(title+f" — status {q['status']}, covQual {q['covQual']}")
 fig.suptitle(f'Combine FitDiagnostics — range 6–{upper} GeV, rebin ×8\nChebyshev + Υ(1S) + Υ(2S); strategia 0, robustFit 1')
 for ext in ['png','pdf']:fig.savefig(p/f'{pref}_fitdiagnostics_log.{ext}',dpi=170)
 for row in meta['parameters']:
  sigma=row['prefit_saved_error'];constrained=row['parameter'].startswith('CMS_eff_') or row['parameter'] in ['CMS_scale_e','lumi_1','CMS_EXO25020_signalModelSigmaNuisance_2022']
  for stage in ['b','s']:
   shift=(row[stage]-row['prefit'])/sigma if sigma>0 and row[stage+'_err'] is not None else None
   allshifts.append(dict(range=f'6-{upper}',parameter=row['parameter'],fit=stage,shift=shift,reference='Gaussian prior' if constrained else 'saved error, not prior' if sigma>0 else 'undefined',beyond_2=abs(shift)>2 if shift is not None else None,covQual=meta[stage]['covQual']))
with (p/'rebin8_parameter_shifts.csv').open('w') as f:
 wr=csv.DictWriter(f,fieldnames=allshifts[0].keys());wr.writeheader();wr.writerows(allshifts)
# Component comparison and ratios
fig,axs=plt.subplots(4,2,figsize=(13,13),layout='constrained',sharex='col');fr,ars=plt.subplots(4,2,figsize=(13,11),layout='constrained',sharex='col')
summary=[]
for j,(tag,upper,ds) in enumerate(sets):
 a=ds['prefit'];x=a['mass'];ed=np.r_[a['low'],a['high'][-1]];factor=refwidth/a['width']
 for i,(key,label) in enumerate([('dy','Chebyshev'),('upsilon1s','Υ(1S)'),('upsilon2s','Υ(2S)'),('signal','Segnale 9 GeV')]):
  ax=axs[i,j];ar=ars[i,j];positive=[]
  summary.append(dict(range=f'6-{upper}',component=key,**{s:float(ds[s][key].sum()) for s in stages}))
  for stage,c,ls,lab in [('prefit','#424242','--','Pre-fit'),('fit_b','#1766ba',':','Post-fit B'),('fit_s','#cc2d2d','-','Post-fit S+B')]:
   yy=ds[stage][key];y=yy*factor;positive.extend(y[y>0]);ax.stairs(np.where(y>0,y,np.nan),ed,color=c,lw=1.8,ls=ls,label=lab,baseline=None)
   if yy.sum()>0:
    pre=ds['prefit'][key];ratio=(yy/yy.sum())/(pre/pre.sum());ar.plot(x,ratio,color=c,ls=ls,lw=1.7,label=lab)
  ax.set_yscale('log');ax.set_ylim(min(positive)*.6,max(positive)*1.7);ax.set_title(f'{label} — 6–{upper} GeV');ax.set_ylabel('Eventi / 0.10667 GeV');ax.grid(alpha=.15,which='both');ax.legend(fontsize=8)
  ar.set_title(f'{label} — 6–{upper} GeV');ar.set_ylabel('Post/pre, integrali a 1');ar.grid(alpha=.15);ar.legend(fontsize=8)
  if key in ['upsilon1s','upsilon2s']:ar.set_ylim(.98,1.02)
  if key=='signal':ax.text(.02,.04,'B: segnale nullo',transform=ax.transAxes,color='#1766ba',fontsize=9)
 axs[-1,j].set_xlim(6,upper);ars[-1,j].set_xlim(6,upper);axs[-1,j].set_xlabel('Massa [GeV]');ars[-1,j].set_xlabel('Massa [GeV]')
fig.suptitle('Rebin ×8 — componenti pre-fit e post-fit\nAltezze riportate alla larghezza nominale; ultimo bin residuo incluso')
fr.suptitle('Rebin ×8 — variazione della sola forma\nOgni componente normalizzata a integrale 1 nel proprio range')
for ext in ['png','pdf']:
 fig.savefig(p/f'rebin8_componenti_pre_post.{ext}',dpi=160);fr.savefig(p/f'rebin8_componenti_shape_ratio.{ext}',dpi=160)
with (p/'rebin8_componenti_integrali.csv').open('w') as f:
 wr=csv.DictWriter(f,fieldnames=summary[0].keys());wr.writeheader();wr.writerows(summary)
print(json.dumps(summary,indent=2))
