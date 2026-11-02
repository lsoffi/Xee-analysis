from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import chi2

P=Path(__file__).resolve().parents[1]
O=P/'runs/20261030_data_above_10p5'
d=json.loads((O/'spectrum.json').read_text())
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
colors=['#0072B2','#D55E00','#009E73','#CC79A7']
edges=np.linspace(*d['observable_range'],d['bins_per_era']+1)
centers=(edges[1:]+edges[:-1])/2
mask=centers>d['cr_end']
counts=[]
def errors(n):
    alpha=1-0.682689492137
    low=np.zeros_like(n,dtype=float)
    sel=n>0
    low[sel]=chi2.ppf(alpha/2,2*n[sel])/2
    high=chi2.ppf(1-alpha/2,2*(n+1))/2
    return np.array([n-low,high-n])
def draw(ax,n,color):
    ax.stairs(n[mask],edges[338:],color=color,lw=1,alpha=.65)
    ax.errorbar(centers[mask],n[mask],yerr=errors(n[mask]),fmt='o',ms=3,
                color=color,elinewidth=.7,capsize=0)
    ax.axvline(d['cr_end'],color='#b07800',ls='--',lw=1.2)
    ax.axvline(12,color='#555555',ls=':',lw=1.2)
    ax.set_xlim(10.49,12.025);ax.set_ylim(bottom=0)
    ax.set_xticks([10.5,10.75,11,11.25,11.5,11.75,12])
    ax.set_xlabel(r'$m_{ee}$ [GeV]');ax.set_ylabel('Eventi / 13,33 MeV')
    ax.grid(alpha=.2)
fig,axes=plt.subplots(2,2,figsize=(13,8.4))
for ax,s,color in zip(axes.flat,d['summary'],colors):
    n=np.array([r['count'] for r in sorted(d['rows'],key=lambda r:r['bin']) if r['era']==s['era']]);counts.append(n)
    draw(ax,n,color)
    ax.set_title(f"{s['era']}  |  {s['events_above_cr']:,} eventi".replace(',',' '),fontweight='bold')
fig.suptitle('Dati sopra la CR: 10,5067–12 GeV, quattro ere\nBinning originale, nessun fit',fontsize=16,y=.98)
fig.text(.5,.015,'Tratteggio: fine CR (10,5067 GeV)   |   Puntinato: limite del workspace (12 GeV)   |   Errori Poisson 68,3%',ha='center',fontsize=10)
fig.tight_layout(rect=(0,.045,1,.92))
fig.savefig(O/'data_above_10p5_eras.png',dpi=170);fig.savefig(O/'data_above_10p5_eras.svg');plt.close(fig)
n=np.sum(counts,axis=0)
fig,ax=plt.subplots(figsize=(12,5.4));draw(ax,n,'#12314b')
ax.set_title(f'Somma delle quattro ere: {n[mask].sum():,} eventi sopra la CR'.replace(',',' '),fontsize=16)
fig.text(.5,.025,'Dati binned originali; errori Poisson 68,3%. Il bordo a 12 GeV è il limite degli input, non una fine fisica dello spettro.',ha='center',fontsize=10)
fig.tight_layout(rect=(0,.065,1,1));fig.savefig(O/'data_above_10p5_sum.png',dpi=170);fig.savefig(O/'data_above_10p5_sum.svg');plt.close(fig)
print(O)
