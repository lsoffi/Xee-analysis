#!/usr/bin/env python3
"""Aggregate completed prospective C0 validation cells versus mass."""
from pathlib import Path
import argparse,csv,json,os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parents[1]
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--tag',required=True);a=ap.parse_args()
 base=P/'runs/validation_c0'/a.tag;rows=[]
 for p in sorted(base.glob('*/m*/summary.json')):
  q=json.loads(p.read_text());s=q['summary'];b=s.get('normalized_bias') or {};g=s.get('gof') or {};plan=q['plan']
  rows.append(dict(side=plan['side'],mass=plan['mass'],injection=plan['injection'],method=plan['method'],
   n_total=s['n_total'],n_valid=s['n_valid'],n_failed=s['n_failed'],bias=s['bias']['mean'],sd=s['bias']['sd'],
   normalized_bias=b.get('value'),bias_error=b.get('bootstrap_se'),gof_p=g.get('p_raw'),
   gof_n=g.get('n_readable'),gof_fd_disagreements=g.get('n_FD_disagreements'),source=str(p)))
 if not rows:raise SystemExit('No summary cells found; execute validation first.')
 with (base/'summary_vs_mass.csv').open('w',newline='') as f:
  wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
 for method in ['bias','gof']:
  selected=[q for q in rows if q['method']==method]
  if not selected:continue
  fig,axes=plt.subplots(1,2,figsize=(12,4.5),sharey=True)
  for side,ax in zip(['low','high'],axes):
   rr=sorted([q for q in selected if q['side']==side],key=lambda q:q['mass'])
   key='normalized_bias' if method=='bias' else 'gof_p'
   good=[q for q in rr if q[key] is not None]
   ax.errorbar([q['mass'] for q in good],[q[key] for q in good],
      yerr=[q['bias_error'] for q in good] if method=='bias' else None,fmt='o-',color='#007E87')
   flagged=[q for q in good if q['n_failed'] or (q.get('gof_fd_disagreements') or 0)]
   ax.scatter([q['mass'] for q in flagged],[q[key] for q in flagged],marker='x',color='#B22222',s=80,label='Flag numerici')
   ax.axhline(0 if method=='bias' else .05,color='gray',ls='--')
   if method=='gof':ax.set_ylim(0,1)
   ax.set_title('SR bassa + CR' if side=='low' else 'SR alta + CR');ax.set_xlabel('Massa [GeV]');ax.grid(alpha=.2)
  axes[0].set_ylabel('Media bias / SD dei toy' if method=='bias' else 'p-value GoF saturated')
  fig.suptitle('C0 - '+method+' - leggere insieme ai fallimenti; nessuna copertura misurata')
  fig.tight_layout()
  for ext in ['png','pdf']:fig.savefig(base/f'{method}_vs_mass.{ext}',dpi=160)
  plt.close(fig)
 print(base/'summary_vs_mass.csv')
if __name__=='__main__':main()
