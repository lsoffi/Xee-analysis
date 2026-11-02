from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import norm
P=Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_common95_cheb2_no_upsilon';fig,axs=plt.subplots(1,2,figsize=(13,5.2))
for dirname,label,color in [('20261029_common95_cheb2','Con Υ(1S), Υ(2S)','#0072B2'),('20261029_common95_cheb2_no_upsilon','Senza Υ','#D55E00')]:
 rows=json.loads((P/'runs'/dirname/'summary.json').read_text())
 for ax,key in zip(axs,['observed','p0']):
  def val(r):return r.get('limits',{}).get('observed',np.nan) if key=='observed' else float(norm.logsf(r['Z'])/np.log(10))
  a=[r for r in rows if key!='observed' or r.get('limit_crossing',{}).get('verified')]
  ax.plot([r['mass'] for r in a],[val(r) for r in a],color=color,lw=1,label=label)
  for ok,marker in [(True,'o'),(False,'x')]:
   b=[r for r in a if r.get('numerical_ok',r['fit_ok'])==ok];ax.scatter([r['mass'] for r in b],[val(r) for r in b],color=color,marker=marker,s=18)
axs[0].set_yscale('log')
for ax in axs:ax.set_xlabel('Massa ipotizzata [GeV]');ax.grid(alpha=.2);ax.legend(fontsize=9)
axs[0].set_ylabel('Limite osservato 95% CLs su r');axs[1].set_ylabel('log₁₀(p-value locale)');fig.suptitle('Stesse finestre comuni 95% — Chebyshev 2 + DCB, con e senza Υ');fig.text(.5,.01,'Croci: controlli numerici non superati. Confronto diagnostico, nessuna validazione toy.',ha='center');fig.tight_layout(rect=[0,.04,1,.95])
for ext in ['png','pdf']:fig.savefig(R/f'comparison.{ext}',dpi=170)
