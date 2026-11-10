import json, pathlib, csv, math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=pathlib.Path(__file__).resolve().parent
out=pathlib.Path('/eos/user/s/soffi/www/Xee/261106/scan_SR_CR/limits_recommended_envelope'); out.mkdir(parents=True,exist_ok=True)
rows=[]
for f in (p/'results').glob('*.json'):
 d=json.load(open(f)); q=d['quantiles']; levels=[.025,.16,.5,.84,.975,-1]
 assert len(q)==6
 vals=[]
 for level in levels:
  hits=[v['limit'] for v in q if abs(v['quantile']-level)<1e-5]
  assert len(hits)==1 and math.isfinite(hits[0]) and hits[0]>=0
  vals.append(hits[0])
 assert vals[:5]==sorted(vals[:5])
 rows.append([d['mass']]+vals)
rows.sort(); assert len(rows)==19
with open(out/'limits.csv','w') as f:
 w=csv.writer(f);w.writerow(['mass_GeV','expected_minus2sigma','expected_minus1sigma','expected_median','expected_plus1sigma','expected_plus2sigma','observed']);w.writerows(rows)
fig,ax=plt.subplots(figsize=(10,6))
for j,group in enumerate([[r for r in rows if r[0]<1],[r for r in rows if r[0]>1]]):
 x=[r[0] for r in group]
 ax.fill_between(x,[r[1] for r in group],[r[5] for r in group],color='#f6cc43',label='Expected ±2σ' if j==0 else None)
 ax.fill_between(x,[r[2] for r in group],[r[4] for r in group],color='#69bb63',label='Expected ±1σ' if j==0 else None)
 ax.plot(x,[r[3] for r in group],'k--',lw=1.7,label='Median expected' if j==0 else None)
 ax.plot(x,[r[6] for r in group],'ko-',ms=4,lw=1.5,label='Observed' if j==0 else None)
ax.axvspan(.65,1.154,color='gray',alpha=.13)
ax.text(.90,.96,'ω/φ CR',transform=ax.get_xaxis_transform(),ha='center',va='top',color='dimgray')
ax.set_yscale('log');ax.set_xlabel('Signal mass hypothesis [GeV]');ax.set_ylabel('95% CL upper limit on signal strength r')
ax.set_title('Simultaneous SR + selected CRs · four eras',loc='left',fontsize=14,pad=30)
ax.text(0,1.01,'C3 + L3 below φ; C4 + L4 between φ and J/ψ',transform=ax.transAxes,fontsize=10)
ax.grid(axis='y',which='both',alpha=.18);ax.legend(loc='upper right',frameon=False);ax.set_xlim(.35,2.75)
fig.text(.12,.015,'Native Combine AsymptoticLimits · trigger profiled · 19 mass hypotheses',fontsize=9)
fig.tight_layout(rect=[0,.04,1,1])
for ext in ['png','pdf']:fig.savefig(out/('limits.'+ext),dpi=180)
(out/'index.html').write_text('<html><body><h1>Sliding-window envelope limits</h1><p>19 masses; C3+L3 below phi, C4+L4 between phi and J/psi. Four eras, simultaneous SR+selected CRs, trigger profiled. Native AsymptoticLimits, strategy 0, MINIMIZER_freezeDisassociatedParams. Six finite quantiles and ordered expected bands checked; this is not an independent minimum audit.</p><img src="limits.png" style="max-width:100%"><p><a href="limits.pdf">PDF</a> · <a href="limits.csv">CSV</a></p></body></html>')
print(out)
