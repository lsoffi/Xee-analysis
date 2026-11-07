exec(open('xee_width_plot.py').read().split('fig,ax=')[0])
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import MultipleLocator
fig,(ax,b)=plt.subplots(2,1,figsize=(10,8),sharex=True,gridspec_kw={'height_ratios':[1,1.25]})
groups=np.split(np.arange(len(m)),np.where(np.diff(m)>.11)[0]+1)
x=np.array(m)
for era,c in zip(['2022','2022EE','2023','2023BPix'],['#0072b2','#d55e00','#009e73','#a4518c']):
 for j,g in enumerate(groups):
    ax.plot(x[g],[d[era,k] for k in x[g]],color=c,label=era if j==0 else None,marker='.',ms=3)
ax.legend(ncol=4,loc='upper left');ax.set_ylabel('Larghezza [GeV]');ax.set_title('Quanto è larga la finestra al 95% del segnale?',loc='left')
for lo,hi,name in [(.650,1.155,'ω / φ'),(2.8,4,'J/ψ, ψ(2S)'),(9.095,10.505,'Υ')]:
 b.axhspan(lo,hi,color='#e7dfd1',zorder=0);b.text(9.12,(lo+hi)/2,name,va='center',ha='right',fontsize=9,color='#624e33')
lo,hi=np.array([row['sr_segments'][0] for row in r]).T
for j,g in enumerate(groups):
 b.fill_between(x[g],lo[g],hi[g],color='#82b8d5',alpha=.8,label='SR comune alle 4 ere' if j==0 else None)
 b.plot(x[g],lo[g],color='#246485',lw=.85);b.plot(x[g],hi[g],color='#246485',lw=.85)
 b.plot(x[g],x[g],':',color='#173c50',lw=1.2,label='Massa ipotizzata' if j==0 else None)
for k in [7,9]:
 i=m.index(k); b.plot([k,k],[lo[i],hi[i]],color='#123651',lw=2,marker='_',ms=8)
 b.annotate(f'{hi[i]-lo[i]:.2f} GeV'.replace('.',','),(k,(lo[i]+hi[i])/2),xytext=(-7,0),textcoords='offset points',ha='right',va='center',fontsize=9,bbox=dict(fc='white',ec='none',alpha=.9,pad=2))
b.set(xlabel='Ipotesi di massa [GeV]',ylabel='Massa dielettronica nel fit [GeV]',ylim=(0,10.65),xlim=(.2,9.27))
b.set_title('Quali masse dielettroniche entrano nella SR?',loc='left');b.legend(loc='upper left',ncol=2)
for a in [ax,b]: a.grid(axis='both',alpha=.15);a.xaxis.set_major_locator(MultipleLocator(1))
fig.tight_layout()
fig.savefig('/eos/user/s/soffi/www/Xee/261106/scan_SR_CR/finestre_SR_per_era.png',dpi=170)
