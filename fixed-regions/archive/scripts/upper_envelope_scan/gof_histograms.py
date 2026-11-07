from pathlib import Path
import json,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path.cwd()/'runs/20261030_upper_envelope_scan';fig,axes=plt.subplots(2,4,figsize=(16,7));colors={1:'#D55E00',2:'#0072B2'}
for di,d in enumerate([1,2]):
    for mi,m in enumerate([10.5,11.,11.5,12.]):
        ax=axes[di,mi];q=json.loads((R/f'degree{d}/m{m:.1f}/results.json').read_text());toys=[json.loads(p.read_text()) for p in (R/f'gof/degree{d}').glob(f'toy*/m{m:.1f}/result.json')];v=np.array([t['score'] for t in toys]);obs=q['gof_observed'];k=int(sum(v>obs))
        edges=np.linspace(min(v.min(),obs)-5,max(v.max(),obs)+5,18);ax.hist(v,bins=edges,color=colors[d],alpha=.65,label='Toy B-only, fit S+B');ax.axvline(obs,color='black',ls='--',label='Dati');ax.set_title(f'Grado {d}, m={m:.1f}: p={k}/{len(v)}');ax.set_xlabel('Statistica saturated Combine');ax.set_ylabel('Toy');ax.grid(alpha=.15)
        if mi==0:ax.legend(fontsize=8)
fig.suptitle('Distribuzioni GoF: esempio di quattro ipotesi, quattro ere combinate\nPrimo stadio CR ripetuto in ogni toy; flag CR inclusi e documentati nel report',fontsize=13);fig.tight_layout(rect=(0,0,1,.91));fig.savefig(R/'gof_distributions.png',dpi=150);fig.savefig(R/'gof_distributions.svg');plt.close(fig)
