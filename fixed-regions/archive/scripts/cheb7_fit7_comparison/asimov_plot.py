from pathlib import Path
import json,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path('runs/20261030_cheb7_fit7_comparison');OLD=Path('runs/20261030_cheb7_bias_scan');old=json.loads((OLD/'asimov_summary.json').read_text());new=json.loads((R/'asimov_summary.json').read_text());rows=[dict(r,model='fixed3' if r['model']=='fixed' else 'sliding3',injection=r['injected'],valid_point=r['b_status']==0 and r['s_status']==0 and r['delta_nll']<=.01) for r in old]+new;styles=[('fixed3','tab:red','Grado3 fisso'),('sliding3','tab:blue','Grado3 sliding'),('fixed7','tab:green','Grado7 fisso')]
fig,axes=plt.subplots(1,3,figsize=(15,4.5))
for model,color,label in styles:
 for ax,scenario in zip(axes[:2],['zero','injected']):
  a=sorted([r for r in rows if r['model']==model and r['scenario']==scenario],key=lambda r:r['mass']);ax.plot([r['mass'] for r in a],[r['r']-r['injection'] if r['valid_point'] else np.nan for r in a],'o-',ms=3,color=color,label=label)
 a=sorted([r for r in rows if r['model']==model and r['scenario']=='zero'],key=lambda r:r['mass']);axes[2].plot([r['mass'] for r in a],[r['Z'] if r['valid_point'] else np.nan for r in a],'o-',ms=3,color=color,label=label)
for ax,title in zip(axes,['Asimov senza segnale: bias','Asimov con segnale: bias','Asimov senza segnale: Z nominale']):ax.set_title(title,fontsize=10);ax.set_xlabel('Massa [GeV]');ax.grid(alpha=.2);ax.legend(fontsize=8);ax.axhline(0,color='gray')
axes[0].set_ylabel('r fitted - r vero');axes[1].set_ylabel('r fitted - r vero');axes[2].set_ylabel('Z');fig.suptitle('Stesso generatore Chebyshev7 + Y, senza fluttuazioni statistiche');fig.tight_layout();fig.savefig(R/'asimov_overlay.png',dpi=150);fig.savefig(R/'asimov_overlay.pdf');plt.close(fig)
