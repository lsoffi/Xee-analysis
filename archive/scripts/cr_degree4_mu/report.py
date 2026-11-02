from run import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import csv
F=json.loads((R/'fits.json').read_text());summary=[];pars=[]
for idx,cfg in enumerate(S):
 fs=[next(f for f in F if f['idx']==idx and f['era']==e) for e in eras];rec=dict(idx=idx,lo=cfg['lo'],hi=cfg['hi'],degree=cfg['degree'],D=sum(f['D'] for f in fs),valid_eras=sum(f['valid'] for f in fs),mu_boundary_eras=sum(any(f['boundary'][2:]) for f in fs),prediction=np.sum([f['prediction'] for f in fs],axis=0).tolist(),error=np.sqrt(np.sum(np.array([f['prediction_error'] for f in fs])**2,axis=0)).tolist());summary.append(rec)
 fig,axes=plt.subplots(2,2,figsize=(11,7))
 for ei,(f,ax) in enumerate(zip(fs,axes.flat)):
  xx=f['x'];ax.errorbar(xx,f['counts'],yerr=np.sqrt(f['counts']),fmt='k.',ms=3,label='Dati');f3=next(z for z in F if z['idx']==2*(idx//2) and z['era']==f['era']);ax.plot(xx,f3['mu'],':',color='grey',label='Totale grado 3');ax.plot(xx,f['mu'],color='#0072B2',label=f"Totale grado {cfg['degree']}");ax.plot(xx,f['continuum'],'--',color='#E69F00',label='Continuo');ax.plot(xx,f['Y'][0],color='#009E73',label='Υ(1S)');ax.plot(xx,f['Y'][1],color='#D55E00',label='Υ(2S)');stat='OK' if f['valid'] else 'NON VALIDATO';bd='; μ al bordo' if any(f['boundary'][2:]) else '';ax.set_title(f"{f['era']}: {stat}{bd}",fontsize=10,color='black' if f['valid'] and not bd else '#B33');ax.grid(alpha=.2);ax.set_xlabel('Massa [GeV]');ax.set_ylabel('Eventi / bin originale')
  for j,pr in enumerate(['Υ1S','Υ2S']):
   nt=cfg['degree']+1;pars.append(dict(idx=idx,era=f['era'],resonance=pr,mu=f['p'][nt+2+2*j],mu_error=f['errors'][nt+2+2*j],mu_nom=f['nominal_shapes'][j][0],sigma_fixed=f['p'][nt+3+2*j],scale=f['p'][nt+j],fit_valid=f['valid'],mu_boundary=f['boundary'][2+2*j]))
 axes[0,0].legend(fontsize=7);fig.suptitle(f"CR {cfg['lo']:.3f}–{cfg['hi']:.3f} GeV · Chebyshev {cfg['degree']}\nμ delle due Υ liberi; σ e code MC fisse");fig.tight_layout();fig.savefig(R/f'fit_{idx}.png',dpi=130);plt.close(fig)
(R/'summary.json').write_text(json.dumps(summary,indent=2))
with (R/'parameters.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(pars[0]));w.writeheader();w.writerows(pars)
fig,axes=plt.subplots(1,2,figsize=(12,4.5));ids=[0,5,10];labels=['9,093–10,507','9,200–10,400','9,293–10,307']
for d,col in [(3,'#0072B2'),(4,'#D55E00')]:
 a=[summary[2*wi+d-3] for wi in ids]
 for j,ax in enumerate(axes):ax.errorbar(np.arange(3)+(d-3.5)*.08,[s['prediction'][j] for s in a],yerr=[s['error'][j] for s in a],fmt='o-',color=col,label=f'Grado {d}');ax.set_xticks(range(3));ax.set_xticklabels(labels,rotation=10);ax.set_xlabel('CR [GeV]');ax.set_ylabel('Eventi Υ previsti');ax.set_title(['Intervallo fisso 6–9,093 GeV','Intervallo fisso 10,507–12 GeV'][j]);ax.grid(alpha=.2)
axes[0].legend();fig.suptitle('Stabilità delle code Υ: μ liberi, σ e code MC fisse\nErrori Hesse diagnostici; non intervalli profilati');fig.tight_layout();fig.savefig(R/'tail_comparison.png',dpi=150);plt.close(fig)
html='<meta charset="utf-8"><style>body{font:17px system-ui;max-width:1200px;margin:auto}img{max-width:100%}</style><h1>CR: Chebyshev 3 contro 4, soli μ liberi</h1><p><a href="REPORT.md">Report</a> · <a href="parameters.csv">Parametri</a> · <a href="pilot_summary.json">Toy</a></p><img src="tail_comparison.png">'
for s in summary:html+=f'<details><summary>CR{s["lo"]:.3f}–{s["hi"]:.3f}, grado{s["degree"]}; {s["valid_eras"]}/4 ere valide, {s["mu_boundary_eras"]} al bordo</summary><img src="fit_{s["idx"]}.png"></details>'
(R/'index.html').write_text(html);print('REPORT PLOTS DONE')
