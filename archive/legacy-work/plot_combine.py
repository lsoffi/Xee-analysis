import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig,axs=plt.subplots(1,2,figsize=(13,4.8),layout='constrained')
for ax,tag,title in zip(axs,['fit_b','fit_s'],['Solo fondo','Segnale + fondo']):
 a=np.genfromtxt('outputs/combine_2022_'+tag+'_shapes.csv',delimiter=',',names=True);x=a['mass'];ed=np.r_[x-1/150,x[-1]+1/150]
 ax.errorbar(x,a['data'],yerr=np.sqrt(a['data']),fmt='.',ms=3,color='black',elinewidth=.6,label='Dati 2022')
 for k,c,l in [('total','#d33f36','Totale post-fit'),('dy','#2171b5','Chebyshev'),('upsilon1s','#f39c12','Υ(1S), norma fissa'),('upsilon2s','#ad579f','Υ(2S)'),('signal','#29934b','Segnale a 9 GeV')]:ax.stairs(a[k],ed,color=c,lw=1.5,label=l)
 ax.set_xlim(8,10.5);ax.set_ylim(0,120);ax.set_xlabel('Massa [GeV]');ax.set_ylabel('Eventi / 0.01333 GeV');ax.set_title(title);ax.grid(alpha=.15)
axs[0].legend(fontsize=8);fig.suptitle('Combine v11.1.0 — FitDiagnostics, 2022, solo Chebyshev\nStrategia 0, robustFit 1; entrambi i fit: status 0, covQual 3',fontsize=12)
for ext in ['png','pdf']:fig.savefig('outputs/combine_2022_fitdiagnostics.'+ext,dpi=180)
