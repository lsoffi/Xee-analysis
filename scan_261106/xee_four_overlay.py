from pathlib import Path
import json
import numpy as np
import ROOT,matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
web=Path('/eos/user/s/soffi/www/Xee/261106/scan_SR_CR/confronto_4_modelli_low')
for m in ('0.4','0.5','0.6'):
 md=json.loads((Path('low_logpoly3/results')/('m'+m)/'model.json').read_text())
 figs,axes=plt.subplots(2,4,figsize=(15,6),sharex='col',gridspec_kw={'height_ratios':[2,1]})
 for n,ch in enumerate(c for c in md['channels'] if c['region']=='SR0'):
  edges=np.array(ch['edges']);x=(edges[1:]+edges[:-1])/2;y=np.array(ch['counts'])
  ax,res=axes[:,n];ax.errorbar(x,y,np.sqrt(y),fmt='k.',ms=3,label='Dati')
  for label,name in [('Cheb2','low_cheb2_totalC0'),('Cheb3','low_cheb3_totalC0'),('Cheb6','low_cheb6_totalC0' if m=='0.6' else ''),('Log3','low_logpoly3_totalC0' if m=='0.6' else 'low_logpoly3')]:
   f=ROOT.TFile.Open(str(next((Path(name)/'results'/('m'+m)).glob('fitDiagnostics*.root'))));h=f.Get('shapes_fit_b/'+ch['label']+'/total')
   v=np.array([h.GetBinContent(i+1) for i in range(len(y))]);f.Close()
   ax.step(x,v,where='mid',label=label,lw=1.4);res.plot(x,(y-v)/np.sqrt(v),'.-',ms=2,lw=.9)
  ax.set_title(ch['era']);res.axhline(0,color='.5',lw=.7);res.set_xlabel('Massa [GeV]')
 axes[0,0].set_ylabel('Eventi / 5 MeV');axes[1,0].set_ylabel('(Dati − fit) / √fit');axes[0,0].legend(fontsize=8)
 figs.suptitle('Fit simultaneo di solo fondo · SR · ipotesi '+m+' GeV');figs.tight_layout();figs.savefig(web/('modelli_m'+m+'.png'),dpi=135);plt.close(figs)
