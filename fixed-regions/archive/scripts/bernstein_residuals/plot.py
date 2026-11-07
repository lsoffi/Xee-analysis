from pathlib import Path
import json,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');SRC=P/'runs/20261029_alternative_background_scan/fits/bernstein7';R=P/'runs/20261029_bernstein_residuals';ERAS=['2022','2022EE','2023','2023BPix'];tables={};stats=[];checks=[]
for era in ERAS:
 a=np.genfromtxt(SRC/f'm7.5/fitdiagnostics/{era}_fit_b_shapes.csv',delimiter=',',names=True);D=a['data'];B=a['total'];pull=(D-B)/np.sqrt(B);term=B-D;positive=D>0;term[positive]+=D[positive]*np.log(D[positive]/B[positive]);dev=np.sign(D-B)*np.sqrt(np.maximum(2*term,0));tables[era]=dict(a=a,pull=pull,deviance=dev)
 for d in sorted(SRC.glob('m*')):
  b=np.genfromtxt(d/f'fitdiagnostics/{era}_fit_b_shapes.csv',delimiter=',',names=True);checks.append(dict(era=era,mass=float(d.name[1:]),max_background_difference=float(np.max(np.abs(b['total']-B))),identical=bool(np.array_equal(b['total'],B))))
 stats.append(dict(era=era,events=float(D.sum()),background=float(B.sum()),pearson_mean=float(pull.mean()),pearson_std=float(pull.std(ddof=1)),max_abs_pearson=float(np.max(np.abs(pull))),poisson_deviance=float(np.sum(dev**2))))
 out=np.c_[a['mass_GeV'],a['low'],a['high'],D,B,a['dy'],a['upsilon1s'],a['upsilon2s'],D-B,pull,dev]
 np.savetxt(R/f'residuals_{era}.csv',out,delimiter=',',header='mass_GeV,low,high,data,background,continuum,upsilon1s,upsilon2s,data_minus_background,pearson_residual,poisson_deviance_residual',comments='')
assert all(c['identical'] for c in checks)
(R/'background_mass_independence.json').write_text(json.dumps(dict(workspaces=51,eras=4,all_identical=True,checks=checks),indent=2))
with (R/'residual_summary.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(stats[0]));w.writeheader();w.writerows(stats)
# Four separate eras, unchanged original binning.
for label,lo,hi in [('6_12',6,12),('6_8',6,8),('8_10p5',8,10.5),('10_12',10,12)]:
 fig,axs=plt.subplots(4,1,figsize=(13,9),sharex=True,layout='constrained');extent=max(4.5,max(np.max(np.abs(t['pull'][(t['a']['mass_GeV']>=lo)&(t['a']['mass_GeV']<=hi)])) for t in tables.values())+.4)
 for ax,era in zip(axs,ERAS):
  t=tables[era];a=t['a'];mask=(a['mass_GeV']>=lo)&(a['mass_GeV']<=hi);ax.axhspan(-2,2,color='#e8edf1');ax.axhspan(-1,1,color='#d4dfe7');ax.axhline(0,color='#c8382b',lw=1);ax.plot(a['mass_GeV'][mask],t['pull'][mask],'.',color='#1f608c',ms=3.3);ax.set_ylabel(f'{era}\n(D−B)/√B');ax.set_ylim(-extent,extent);ax.set_xlim(lo,hi);ax.grid(alpha=.15)
 axs[-1].set_xlabel('Massa dielettronica [GeV]');fig.suptitle('Bernstein grado 7 identificabile — residui del fit solo fondo\nShape Υ(1S), Υ(2S) assunte valide; fit 6–12 GeV, bin originali 13,33 MeV',fontsize=13)
 for ext in ['png','pdf']:fig.savefig(R/f'pearson_four_eras_{label}.{ext}',dpi=180)
 plt.close(fig)
# Overview with data and the B-only model.
fig=plt.figure(figsize=(14,10));outer=fig.add_gridspec(2,2,hspace=.3,wspace=.2)
for cell,era in zip(outer,ERAS):
 gs=cell.subgridspec(2,1,height_ratios=[3,1],hspace=.05);ax=fig.add_subplot(gs[0]);res=fig.add_subplot(gs[1],sharex=ax);t=tables[era];a=t['a'];e=np.r_[a['low'],a['high'][-1]]
 ax.errorbar(a['mass_GeV'],a['data'],yerr=np.sqrt(a['data']),fmt='.',color='black',ms=2,lw=.4,label='Dati');ax.stairs(a['total'],e,baseline=None,color='#c8382b',lw=1.5,label='Fondo totale B-only');ax.stairs(a['dy'],e,baseline=None,color='#246eae',lw=1,ls='--',label='Continuo Bernstein');ax.set_yscale('log');ax.set_ylim(.5,max(a['data'])*2);ax.set_xlim(6,12);ax.set_title(era);ax.set_ylabel('Eventi / 13,33 MeV');ax.tick_params(labelbottom=False);ax.grid(alpha=.15,which='both');ax.legend(fontsize=8)
 res.axhspan(-2,2,color='#e8edf1');res.axhline(0,color='#c8382b',lw=1);res.plot(a['mass_GeV'],t['pull'],'.',ms=2.5,color='#1f608c');res.set_ylim(-5.5,5.5);res.set_ylabel('(D−B)/√B',fontsize=9);res.set_xlabel('Massa [GeV]');res.grid(alpha=.15)
fig.suptitle('Dati e fit solo fondo — Bernstein grado 7 identificabile\nShape Υ assunte valide; normalizzazioni e nuisance originali del fit preservate',fontsize=13,y=.985);fig.subplots_adjust(top=.91,bottom=.07,left=.07,right=.98)
for ext in ['png','pdf']:fig.savefig(R/f'data_and_Bonly_four_eras.{ext}',dpi=180)
plt.close(fig)
# Sum is a display diagnostic, not a new combined fit or significance.
a=tables[ERAS[0]]['a'];D=sum(tables[e]['a']['data'] for e in ERAS);B=sum(tables[e]['a']['total'] for e in ERAS);pull=(D-B)/np.sqrt(B)
fig,ax=plt.subplots(figsize=(13,4),layout='constrained');ax.axhspan(-2,2,color='#e8edf1');ax.axhline(0,color='#c8382b',lw=1);ax.plot(a['mass_GeV'],pull,'.',ms=3.5,color='#1f608c');ax.set_xlim(6,12);lim=max(5.5,float(np.max(np.abs(pull)))+.4);ax.set_ylim(-lim,lim);ax.set_xlabel('Massa dielettronica [GeV]');ax.set_ylabel('(ΣD−ΣB)/√ΣB');ax.grid(alpha=.15);ax.set_title('Somma delle quattro ere — residui B-only, binning originale\nDiagnostica descrittiva: nessuna incertezza del modello inclusa')
for ext in ['png','pdf']:fig.savefig(R/f'pearson_sum_eras_6_12.{ext}',dpi=180)
plt.close(fig)
np.savetxt(R/'residuals_sum_eras.csv',np.c_[a['mass_GeV'],D,B,pull],delimiter=',',header='mass_GeV,data_sum,background_sum,pearson_residual',comments='')
(R/'README.md').write_text('''# Residui del solo Bernstein

Modello attivo: Bernstein grado 7 identificabile + Υ(1S) + Υ(2S), fit B-only simultaneo delle quattro ere, range 6–12 GeV, 450 bin originali per era.

**Le shape delle Υ sono assunte valide su indicazione dell’utente.** Nessuna modifica delle loro forme e nessun nuovo fit sono stati eseguiti per questi grafici. Sono preservate le normalizzazioni e le nuisance del fit già effettuato; l’analisi dei residui è concentrata sul continuo, assumendo corrette le componenti risonanti.

Il fondo B-only e i dati sono esattamente identici in tutti i 51 workspace dello scan, verificati bin per bin nelle quattro ere. È usato il risultato archiviato a massa ipotizzata 7,5 GeV come riferimento univoco: non si tratta di un fit con segnale a 7,5. Nessuna selezione della massa basata sui residui.

## Grafici

- [Dati, fondo e residui](data_and_Bonly_four_eras.png)
- [Residui 6–12, quattro ere](pearson_four_eras_6_12.png)
- [Zoom 6–8](pearson_four_eras_6_8.png)
- [Zoom 8–10,5](pearson_four_eras_8_10p5.png)
- [Zoom 10–12](pearson_four_eras_10_12.png)
- [Somma delle ere, diagnostica](pearson_sum_eras_6_12.png)

Tutti i grafici sono disponibili anche in PDF. I CSV per era contengono dati, fondo, componenti, differenza D−B, residuo Pearson e residuo di devianza Poisson firmato.

I residui mostrati sono (D−B)/√B. Le bande ±1 e ±2 sono riferimenti di scala, non intervalli di confidenza del modello. Non includono l’incertezza della previsione né le correlazioni indotte dal fit agli stessi dati. Non sono significanze locali del segnale. Anche il grafico sommato è una diagnostica descrittiva.

Gli altri tre modelli sono archiviati nello studio precedente, con tutti i file conservati e compatibilità dei collegamenti storici. Vedere archive_manifest.json.
''')
print(json.dumps(dict(eras=4,bins_per_era=450,background_identical_across_51_masses=True,summary=stats),indent=2))
