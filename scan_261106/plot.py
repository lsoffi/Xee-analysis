import json
import math
import os
import time
from pathlib import Path
import ROOT
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

web=Path('/eos/user/s/soffi/www/Xee/261106/scan_SR_CR')
web.mkdir(exist_ok=True)
masses=[float(m) for m in Path('masses.txt').read_text().split()]
def extract(m,method):
    d=Path('results')/f'm{m:.1f}'
    if not (d/(method+'.done')).exists():return None
    files=list(d.glob(f'higgsCombine*.{method}.*.root'))
    if len(files)!=1:return None
    f=ROOT.TFile.Open(str(files[0]))
    if not f or f.IsZombie():return None
    t=f.Get('limit')
    vals=sorted((float(x.quantileExpected),float(x.limit)) for x in t) if t else []
    f.Close()
    if method=='Significance' and len(vals)==1 and math.isfinite(vals[0][1]):
        return .5*math.erfc(max(0,vals[0][1])/math.sqrt(2))
    if method=='AsymptoticLimits' and len(vals)==6:
        q,v=zip(*vals)
        if np.allclose(q,[-1,.025,.16,.5,.84,.975],atol=.002) and all(math.isfinite(x) and x>0 for x in v) and sorted(v[1:])==list(v[1:]):return v
def save(fig,name):
    tmp=web/(name+'.tmp.png')
    fig.savefig(tmp,dpi=135,bbox_inches='tight');tmp.replace(web/(name+'.png'));plt.close(fig)
def render():
    L=[extract(m,'AsymptoticLimits') for m in masses]; P=[extract(m,'Significance') for m in masses]
    groups=np.split(np.arange(len(masses)),np.where(np.diff(masses)>.101)[0]+1)
    for label,data in [('limite',L),('pvalue',P)]:
        count=sum(v is not None for v in data)
        if not count:continue
        fig,ax=plt.subplots(figsize=(9,5))
        for group in groups:
            x=np.array(masses)[group]
            if label=='limite':
                a=np.array([L[i] if L[i] else [np.nan]*6 for i in group]).T
                for low,hi,color in [(1,5,'#f5cf52'),(2,4,'#85ce70')]:
                    ax.fill_between(x,a[low],a[hi],color=color);ax.vlines(x,a[low],a[hi],color=color)
                ax.plot(x,a[3],'k--',label='Mediana attesa');ax.plot(x,a[0],'ko-',ms=3,label='Osservato')
            else:ax.plot(x,[P[i] if P[i] is not None else np.nan for i in group],'o-',ms=3)
        if label=='limite':
            from matplotlib.patches import Patch
            h,l=ax.get_legend_handles_labels();ax.legend(handles=h[:2]+[Patch(color='#85ce70',label='±1σ'),Patch(color='#f5cf52',label='±2σ')])
        ax.set(yscale='log',xlabel='Massa [GeV]',ylabel='Limite al 95% CL su r' if label=='limite' else 'p-value locale osservato',title=f'SR + CR · 4 ere · Chebyshev 6 · {count}/{len(masses)} punti')
        if label=='pvalue':
            low,high=ax.get_ylim()
            high=max(high,.55)
            ax.set_ylim(low,high)
            right=ax.twinx()
            right.set_yscale('log')
            right.set_ylim(low,high)
            ticks=[(z,.5*math.erfc(z/math.sqrt(2))) for z in range(39)]
            visible=[(z,p) for z,p in ticks if low<=p<=high]
            right.set_yticks([p for z,p in visible])
            right.set_yticklabels([str(z) for z,p in visible])
            right.minorticks_off()
            right.set_ylabel('Significanza locale Z [σ]')
        ax.grid(alpha=.2);save(fig,label)
    summary=dict(masses=masses,limits=L,p_local=P)
    tmp=web/'results.tmp';tmp.write_text(json.dumps(summary));tmp.replace(web/'results.json')
    return all(x is not None for x in L+P)
for _ in range(720):
    if render():break
    time.sleep(60)
