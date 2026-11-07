from pathlib import Path
import json,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path.cwd();R=P/'runs/20261030_sr6_cr4_scan';A=json.loads((R/'bias_summary.json').read_text());fig,axs=plt.subplots(2,2,figsize=(13,8),sharex=True);grid=np.arange(60,91)/10
for j,rt in enumerate([0,1]):
 for family,col in [('bern7','tab:blue'),('cheb7','tab:orange')]:
  rows=[a for a in A if a['family']==family and a['rtrue']==rt and 'bias' in a];by={a['mass']:a for a in rows};valid=[a for a in rows if a.get('valid')];bad=[a for a in rows if not a.get('valid')]
  for ax,key in [(axs[j,0],'bias_sigma'),(axs[j,1],'bias')]:
   ax.plot(grid,[by[m][key] if m in by and by[m]['valid'] else np.nan for m in grid],color=col,label=f'Generatore {family}')
   ax.scatter([a['mass'] for a in valid],[a[key] for a in valid],s=18,color=col)
   ax.scatter([a['mass'] for a in bad if a[key] is not None],[a[key] for a in bad if a[key] is not None],marker='x',s=35,color=col)
 axs[j,0].axhspan(-.2,.2,color='gray',alpha=.15);axs[j,0].set_ylabel(f'r iniettato = {rt}\nBias / σ(r)');axs[j,1].set_ylabel('r fitted − r iniettato')
 for ax in axs[j]:ax.axhline(0,color='gray',lw=.7);ax.grid(alpha=.2);ax.legend(fontsize=9);ax.set_xlim(5.95,9.05)
for ax in axs[-1]:ax.set_xlabel('Massa ipotizzata [GeV]')
fig.suptitle(f'Bias Asimov — envelope Chebyshev6/Bernstein6, SR+CR — {len(A)}/124 fit\nNuisance originali profilate; croci: verifica numerica aperta; r=1 non è il limite atteso')
fig.tight_layout()
for ext in ['png','pdf']:fig.savefig(R/f'bias.{ext}',dpi=150)
plt.close(fig)
