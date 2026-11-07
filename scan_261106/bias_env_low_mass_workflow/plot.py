import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
web=Path('/eos/user/s/soffi/www/Xee/261106/scan_SR_CR')
fig,ax=plt.subplots()
rows=[]
for g in ('C2','L2','C3','L3'):
 vals=[]
 for i in (4,5,6):
  rec=next(r for r in csv.DictReader((web/f'bias_env_0{i}_200/bias.csv').open()) if r['generator']==g)
  rows.append(dict(mass=i/10,**rec))
  vals.append([float(rec[k]) for k in ('bias_norm','CI_low','CI_high')])
 y,l,h=zip(*vals)
 off={'C2':-.003,'L2':-.001,'C3':.001,'L3':.003}[g]
 ax.errorbar([x+off for x in (.4,.5,.6)],y,yerr=[[a-b for a,b in zip(y,l)],[b-a for a,b in zip(y,h)]],marker='o',capsize=3,label=g)
ax.axhspan(-.2,.2,color='green',alpha=.09)
ax.axhline(0,color='gray',lw=.7)
ax.set(xticks=[.4,.5,.6],xlabel='Mass hypothesis [GeV]',ylabel='Mean fitted r / toy standard deviation',title='Cheb 2 + LogPol 3 envelope, injected r = 0')
ax.legend(title='Toy generator');fig.tight_layout()
dest=web/'bias_env_low_mass';dest.mkdir(exist_ok=True)
for ext in ('pdf','png'):fig.savefig(dest/('bias_vs_mass.'+ext))
with (dest/'bias_vs_mass.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
