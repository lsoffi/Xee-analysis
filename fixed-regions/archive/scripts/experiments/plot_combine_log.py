import numpy as np,json,csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
p=Path('outputs');meta=json.loads((p/'combine_2022_pre_post.json').read_text());rows=[]
for a in meta['parameters']:
 sigma=a['prefit_saved_error'];constrained=a['parameter'] in ['CMS_scale_e','CMS_EXO25020_signalModelSigmaNuisance_2022','CMS_eff_e_id','CMS_eff_e_reco','CMS_eff_trigger_e','lumi_1']
 for tag in ['b','s']:
  shift=(a[tag]-a['prefit'])/sigma if sigma>0 and a[tag+'_err'] is not None else None
  rows.append(dict(parameter=a['parameter'],fit=tag,prefit=a['prefit'],postfit=a[tag],reference_sigma=sigma if sigma>0 else None,reference_type='Gaussian constraint' if constrained else 'saved workspace error; not a prior' if sigma>0 else 'no prefit sigma',shift_in_reference_sigma=shift,beyond_2sigma=abs(shift)>2 if shift is not None else None))
with (p/'combine_2022_parameter_shifts.csv').open('w') as f:
 wr=csv.DictWriter(f,fieldnames=rows[0].keys());wr.writeheader();wr.writerows(rows)
for lim,suffix in [((8,10.5),''),((6,12),'_fullrange')]:
 fig,axs=plt.subplots(1,2,figsize=(14,5.6),layout='constrained')
 for ax,tag,title in zip(axs,['fit_b','fit_s'],['Solo fondo','Segnale + fondo']):
  a=np.genfromtxt(p/('combine_2022_'+tag+'_shapes.csv'),delimiter=',',names=True);x=a['mass'];ed=np.r_[x-1/150,x[-1]+1/150]
  ax.errorbar(x,a['data'],yerr=np.sqrt(a['data']),fmt='.',ms=3,color='black',elinewidth=.6,label='Dati 2022',zorder=5)
  for k,c,l,ls in [('total','#c62828','Totale post-fit','-'),('dy','#1565c0','Chebyshev','--'),('upsilon1s','#ef8a00','Υ(1S), normalizzazione fissa','-.'),('upsilon2s','#9c27b0','Υ(2S)',':'),('signal','#188341','Segnale a 9 GeV' if tag=='fit_s' else 'Segnale: nullo (r = 0)','-')]:
   y=np.where(a[k]>0,a[k],np.nan)
   ax.stairs(y,ed,color=c,lw=1.8,label=l,linestyle=ls)
  ax.set_yscale('log');ax.set_xlim(*lim);ax.set_ylim(.001,700 if lim[0]==6 else 200);ax.set_xlabel('Massa [GeV]');ax.set_ylabel('Eventi / 0.01333 GeV');ax.set_title(title);ax.grid(alpha=.15,which='both');ax.legend(fontsize=8,loc='lower left',framealpha=.95)
 fig.suptitle('Combine FitDiagnostics — 2022, Chebyshev + Υ(1S) + Υ(2S)\nShape post-fit; strategia 0, robustFit 1; asse verticale logaritmico',fontsize=12)
 for ext in ['png','pdf']:fig.savefig(p/('combine_2022_fitdiagnostics_log'+suffix+'.'+ext),dpi=180)
print('Beyond 2:',[r for r in rows if r['beyond_2sigma']]);print('largest shifts:',sorted([r for r in rows if r['shift_in_reference_sigma'] is not None],key=lambda r:abs(r['shift_in_reference_sigma']),reverse=True)[:4])
