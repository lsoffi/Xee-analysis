import json, ROOT, numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path('/eos/user/s/soffi/www/Xee/261106/scan_SR_CR/confronto_4_modelli_low')
for stage in ['b','s']:
 fig,axes=plt.subplots(2,4,figsize=(14,5.5),sharex='col',gridspec_kw={'height_ratios':[1.55,1]})
 allpull=[]
 for j,(name,folder,color,marker) in enumerate([('Cheb 2','low_cheb2_totalC0','#116598','o'),('Log 2','low_logpoly2','#bf5026','x')]):
  p=Path(folder)/'results/m0.5';meta=json.loads((p/'model.json').read_text());wf=ROOT.TFile.Open(str(p/'workspace.root'));w=wf.Get('w');f=ROOT.TFile.Open(str(next(p.glob('fitDiagnostics*.root'))));fit=f.Get('fit_'+stage);pars=fit.floatParsFinal()
  for par in pars:w.var(par.GetName()).setVal(par.getVal())
  aux=0
  for constraint in meta['constraints']:
   v,mu,sigma=list(w.pdf(constraint).servers());aux+=((v.getVal()-mu.getVal())/sigma.getVal())**2
  total=aux;n=len(meta['constraints'])-len(pars);summary=[];k=0
  for ch in meta['channels']:
   h=f.Get('shapes_fit_'+stage+'/'+ch['label']+'/total');y=np.array(ch['counts']);v=np.array([h.GetBinContent(i+1) for i in range(len(y))]);assert np.all(v>0)
   z=(y-v)/np.sqrt(v);total+=z@z;n+=len(z)
   if ch['region']!='SR0':continue
   edges=np.array(ch['edges']);x=(edges[:-1]+edges[1:])/2;ax,bot=axes[:,k];k+=1
   if j==0:
    ax.errorbar(x,y,np.sqrt(y),fmt='.',color='black',ms=4,label='Dati',lw=.8);ax.set_title(ch['era'])
    bot.axhspan(-1,1,color='0.94');bot.axhline(0,color='.5',lw=.8)
    for l in [-2,2]:bot.axhline(l,color='.7',ls=':',lw=.8)
    bot.set_xlabel('Massa [GeV]')
   ax.stairs(v,edges,color=color,label=name,lw=1.5,linestyle='--' if j else '-')
   bot.plot(x+(.00032 if j else -.00032),z,marker,ms=4,color=color,mfc='none')
   ax.text(.02,.96-.09*j,f'{name}  χ² / bin = {z@z/len(z):.2f}',transform=ax.transAxes,va='top',fontsize=8,color=color)
   ax.margins(y=.23);summary.append([ch['era'],round(z@z,2),round(max(abs(z)),2),int(sum(abs(z)>2))]);allpull.extend(z)
  print(stage,name,round(total,2),n,round(total/n,4),summary,flush=True)
  f.Close();wf.Close()
 lim=max(3,np.ceil(max(abs(np.array(allpull)))*2)/2)
 for a in axes[1]:a.set_ylim(-lim,lim)
 axes[0,0].set_ylabel('Eventi / 5 MeV');axes[1,0].set_ylabel('Pull di Pearson');fig.legend(*axes[0,0].get_legend_handles_labels(),loc='upper right',ncol=3,frameon=False)
 fig.suptitle('0,5 GeV · SR · '+('solo fondo' if stage=='b' else 'segnale + fondo'),x=.06,ha='left');fig.tight_layout(rect=[0,0,1,.94]);fig.savefig(out/f'diagnostica_05_grado2_{stage}.png',dpi=150);plt.close(fig)
