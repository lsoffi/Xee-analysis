import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import chi2
from pathlib import Path
a=np.genfromtxt('outputs/xee_2022_bonly_fit_bins.csv',delimiter=',',names=True);x=(a['low']+a['high'])/2;ed=np.r_[a['low'],a['high'][-1]];n=a['data'];mu=a['total']
fig,axs=plt.subplots(2,2,figsize=(13,6.4),gridspec_kw={'height_ratios':[3,1]},sharex='col',layout='constrained')
for j in range(2):
 ax=axs[0,j];ax.errorbar(x,n,yerr=np.sqrt(n),fmt='.',color='black',ms=3,elinewidth=.7,label='Dati 2022')
 for k,c,l in [('total','#d33f36','Fondo totale, post-fit'),('dy','#2171b5','Chebyshev, grado 7'),('upsilon1s','#f39c12','Υ(1S)'),('upsilon2s','#ad579f','Υ(2S)')]:ax.stairs(a[k],ed,color=c,lw=1.6,label=l)
 ax.set_ylabel('Eventi / 0.01333 GeV');ax.set_ylim(bottom=0);ax.grid(alpha=.15)
 z=(n-mu)/np.sqrt(mu);axs[1,j].plot(x,z,'.',color='black',ms=3);axs[1,j].axhline(0,color='red',lw=1);axs[1,j].set_ylabel('(Dati−fit)/√fit');axs[1,j].set_xlabel('Massa [GeV]');axs[1,j].grid(alpha=.15);axs[1,j].set_ylim(-4,4)
axs[0,0].legend(fontsize=8);axs[0,0].set_xlim(6,12);axs[0,1].set_xlim(8,10.5);axs[0,1].set_ylim(0,120)
fig.suptitle('2022 — fit solo fondo: Chebyshev + Υ(1S) + Υ(2S)\n10 parametri liberi; status = 0, covQual = 3; devianza Poisson = 424.19 / 440',fontsize=12)
for ext in ['png','pdf']:fig.savefig('outputs/xee_2022_bonly_fit.'+ext,dpi=180)
print('Nominal asymptotic goodness-of-fit p:',chi2.sf(424.19034241,440));print('sum data/model',n.sum(),mu.sum())
