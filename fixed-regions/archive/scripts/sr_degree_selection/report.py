from fit import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import beta
import csv
O=json.loads((R/'observed.json').read_text());T=[json.loads(s) for s in (R/'toys.jsonl').read_text().splitlines()];A=[json.loads(s) for s in (R/'signal_asimov.jsonl').read_text().splitlines()] if (R/'signal_asimov.jsonl').exists() else [];summary=[]
def combine_rows(rows):
 out=[]
 for k in sorted(set(x['toy'] for x in rows)):
  fs=[x for x in rows if x['toy']==k]
  if len(fs)==4:out.append(dict(toy=k,valid=all(f['valid'] for f in fs),D=sum(f['D'] for f in fs),D_SR=sum(f['D_SR'] for f in fs),D_CR=sum(f['D_CR'] for f in fs)))
 return out
for d in range(3,8):
 obs=[f for f in O if f['degree']==d];D=sum(f['D'] for f in obs);toys=combine_rows([f for f in T if f['degree']==d and f['generator_degree']==d]);good=[x for x in toys if x['valid']];k=sum(x['D']>=D for x in good);s=dict(degree=d,D=D,D_SR=sum(f['D_SR'] for f in obs),D_CR=sum(f['D_CR'] for f in obs),observed_valid=all(f['valid'] for f in obs),gof_total=len(toys),gof_valid=len(good),gof_exceedances=k,gof_p_conditional=k/len(good) if good else None,gof_failure_bounds=[k/len(toys),(k+len(toys)-len(good))/len(toys)] if toys else None)
 if d<7:
  nxt=combine_rows([f for f in T if f['degree']==d+1 and f['generator_degree']==d]);pairs=[(a,b) for a,b in zip(toys,nxt) if a['toy']==b['toy'] and a['valid'] and b['valid']];delta=D-sum(f['D'] for f in O if f['degree']==d+1);kk=sum(a['D']-b['D']>=delta for a,b in pairs);s.update(LR_delta=delta,LR_valid=len(pairs),LR_total=len(toys),LR_exceedances=kk,LR_conditional_p=kk/len(pairs) if pairs else None,LR_failure_bounds=[kk/len(toys),(kk+len(toys)-len(pairs))/len(toys)] if toys else None,negative_LR=sum(a['D']<b['D']-.001 for a,b in pairs))
 if d>=4:
  ss=[a for a in A if a['degree']==d];s['signal_tests']=len(ss);s['signal_valid']=sum(a['valid'] for a in ss);s['signal_max_abs_bias_sigma']=max([abs(a['bias_sigma']) for a in ss if a['valid']]+[0])
 summary.append(s)
(R/'summary.json').write_text(json.dumps(summary,indent=2))
for d in range(3,8):
 for region,sl in [('SR',slice(0,NS)),('CR',slice(NS,338))]:
  fig,axs=plt.subplots(2,2,figsize=(11,7));xx=((edges[:-1]+edges[1:])/2)[sl]
  for ei,ax in enumerate(axs.flat):
   f=next(f for f in O if f['degree']==d and f['era']==eras[ei]);nn=np.array(f['counts'])[sl];ax.errorbar(xx,nn,yerr=np.sqrt(nn),fmt='k.',ms=2.5,label='Dati');ax.plot(xx,np.array(f['mu'])[sl],color='#0072B2',label='Totale');ax.plot(xx,np.array(f['continuum'])[sl],'--',color='#E69F00',label='Continuo');ax.plot(xx,np.array(f['Y'][0])[sl],color='#009E73',label='Υ(1S)');ax.plot(xx,np.array(f['Y'][1])[sl],color='#D55E00',label='Υ(2S)');ax.set_title(eras[ei]+(' — NON VALIDATO' if not f['valid'] else ''));ax.set_xlabel('Massa [GeV]');ax.set_ylabel('Eventi / bin originale');ax.grid(alpha=.2)
  axs[0,0].legend(fontsize=8);fig.suptitle(f'Fit simultaneo solo fondo: SR grado {d} / CR grado 4 — {region}\nμ e normalizzazioni Υ condivisi SR–CR; σ e code MC fisse');fig.tight_layout();fig.savefig(R/f'fit_degree{d}_{region}.png',dpi=140);plt.close(fig)
# Per-era SR residuals across degrees
fig,axes=plt.subplots(4,5,figsize=(16,9),sharex=True,sharey=True)
for ei,e in enumerate(eras):
 for j,d in enumerate(range(3,8)):
  f=next(f for f in O if f['degree']==d and f['era']==e);nn=np.array(f['counts'])[:NS];mu=np.array(f['mu'])[:NS];res=np.sign(nn-mu)*np.sqrt(np.maximum(0,2*(mu-nn+xlogy(nn,nn/mu))));ax=axes[ei,j];ax.axhline(0,color='k',lw=.5);ax.plot((edges[:NS]+edges[1:NS+1])/2,res,'.',ms=2);ax.set_title(f'{e} · grado {d}',fontsize=9);ax.set_ylim(-4.5,4.5);ax.grid(alpha=.2)
fig.suptitle('Residui di devianza nella SR — fit simultaneo con CR invariata');fig.tight_layout();fig.savefig(R/'residuals.png',dpi=150);plt.close(fig)
if A:
 fig,axes=plt.subplots(2,2,figsize=(12,8),sharex=True,sharey=True)
 for row,gen in enumerate([7,8]):
  for col,inj in enumerate(['zero','expected']):
   ax=axes[row,col]
   for d in [4,5,6,7]:
    aa=sorted([a for a in A if a['generator']==gen and a['degree']==d and a['injection']==inj],key=lambda a:a['mass']);ax.plot([a['mass'] for a in aa],[a['bias_sigma'] if a['valid'] else np.nan for a in aa],'o-',ms=3,label=f'Fit SR grado {d}')
    bad=[a for a in aa if not a['valid']];ax.plot([a['mass'] for a in bad],[a['bias_sigma'] for a in bad],'rx',ms=6)
   ax.axhspan(-.2,.2,color='grey',alpha=.2);ax.axhline(0,color='k',lw=.7);ax.set_title(f'Generatore SR {gen}, '+('r=0' if inj=='zero' else 'r=limite atteso di riferimento'));ax.set_xlabel('Massa ipotizzata [GeV]');ax.set_ylabel('(r fit − r vero) / errore Hesse');ax.grid(alpha=.2)
 axes[0,0].legend(fontsize=8);fig.suptitle('Segnale spurio e recupero DCB: test Asimov simultaneo SR+CR\nNessuna fluttuazione statistica; forma segnale MC nominale, non test di copertura');fig.tight_layout();fig.savefig(R/'signal_bias.png',dpi=150);plt.close(fig)
html='<meta charset="utf-8"><style>body{font:17px system-ui;max-width:1200px;margin:auto}img{max-width:100%}</style><h1>Scelta grado SR — CR definitiva invariata</h1><p><a href="REPORT.md">Report</a> · <a href="summary.json">Risultati numerici</a></p><img src="signal_bias.png"><img src="residuals.png">'
for d in range(3,8):html+=f'<details><summary>SR grado {d} / CR grado 4</summary><img src="fit_degree{d}_SR.png"><img src="fit_degree{d}_CR.png"></details>'
(R/'index.html').write_text(html);print(json.dumps(summary,indent=2))
