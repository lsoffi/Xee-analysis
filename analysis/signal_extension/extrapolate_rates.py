from pathlib import Path
import json,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path.cwd();R=P/'runs/20261030_upper_envelope_scan';A=json.loads((R/'signal_rates_existing.json').read_text());x=np.array([a['mass'] for a in A]);y=np.array([a['rates'] for a in A]);eras=['2022','2022EE','2023','2023BPix'];out=[]
def fit(e,end,degree):
    sel=(x>=10)&(x<=end);t=x[sel]-end;anchor=y[np.argmin(abs(x-end)),e];B=np.stack([t**k for k in range(1,degree+1)],axis=1);coef=np.linalg.lstsq(B,np.log(y[sel,e]/anchor),rcond=None)[0]
    return anchor,coef
fig,axes=plt.subplots(2,2,figsize=(12,8))
for e,(era,ax) in enumerate(zip(eras,axes.flat)):
    trials=[]
    for degree in [1,2]:
        for end in [10.5,10.8]:
            anchor,c=fit(e,end,degree);te=x>end;pred=anchor*np.exp(sum(c[k-1]*(x[te]-end)**k for k in range(1,degree+1)));trials.append(dict(degree=degree,training=[10,end],validation=[float(x[te][0]),11.],relative_rms=float(np.sqrt(np.mean((pred/y[te,e]-1)**2))),maximum_relative_error=float(np.max(abs(pred/y[te,e]-1)))))
    anchor,c=fit(e,11.,2);assert c[0]+2*c[1]*1<0 and c[0]<0
    altanchor,altc=fit(e,11.,1);m=np.arange(111,121)/10;pred=anchor*np.exp(c[0]*(m-11)+c[1]*(m-11)**2);alt=altanchor*np.exp(altc[0]*(m-11))
    out.append(dict(era=era,fit_range=[10.,11.],formula='N(m)=N(11)*exp(b*(m-11)+c*(m-11)^2)',N11=float(anchor),b=float(c[0]),c=float(c[1]),holdout=trials,predictions=[dict(mass=float(mm),rate=float(v),loglinear_alternative=float(vv)) for mm,v,vv in zip(m,pred,alt)]))
    sel=x>=9.5;ax.plot(x[sel],y[sel,e],'k.',label='Rese originali');xx=np.linspace(10,12,200);ax.plot(xx,anchor*np.exp(c[0]*(xx-11)+c[1]*(xx-11)**2),color='#0072B2',label='Fit log-quadratico');ax.plot(m,alt,'--',color='#E69F00',label='Alternativa log-lineare');ax.axvline(11,color='gray',ls=':');ax.axvspan(11,12,color='#fff0d0',alpha=.6);ax.set_title(era);ax.set_xlabel('Massa ipotizzata [GeV]');ax.set_ylabel('Resa nominale a r=1 (dominio 6–12)');ax.grid(alpha=.2);ax.legend(fontsize=8)
fig.suptitle('Estrapolazione delle rese del segnale: fit 10–11 GeV, ancorato al valore a 11\nOltre 11 GeV: ipotesi di modello; incertezze relative di normalizzazione mantenute a quelle di 11',fontsize=12);fig.tight_layout(rect=(0,0,1,.92));fig.savefig(R/'signal_rate_extrapolation.png',dpi=160);fig.savefig(R/'signal_rate_extrapolation.svg');plt.close(fig)
(R/'signal_rate_extrapolation.json').write_text(json.dumps(out,indent=2));print(json.dumps([dict(era=a['era'],prediction12=a['predictions'][-1],validation=a['holdout']) for a in out],indent=2))
