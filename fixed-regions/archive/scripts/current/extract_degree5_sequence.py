import ROOT,csv,json
from pathlib import Path
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
p=Path('outputs');eras=[('2022','ch1'),('2022EE','ch2'),('2023','ch3'),('2023BPix','ch4')];summ=[];rows=[]
for era in [e for e,_ in eras]+['combined']:
 fb=ROOT.TFile.Open(str(p/(f'fitDiagnostics_degree5_stageB_{era}.root' if era!='combined' else 'fitDiagnostics_degree5_Binitialized_combined_check.root')))
 fs=ROOT.TFile.Open(str(p/(f'fitDiagnostics_degree5_Binitialized_{era}.root' if era!='combined' else 'fitDiagnostics_degree5_Binitialized_combined_check.root')))
 wf=ROOT.TFile.Open(str(p/f'Xee_degree5_eraStart_{era}.root'));w=wf.Get('w');sf=ROOT.TFile.Open(str(p/f'Xee_degree5_Binitialized_{era}.root'));sw=sf.Get('w');b=fb.Get('fit_b');s=fs.Get('fit_s');r=s.floatParsFinal().find('r')
 summ.append(dict(era=era,events=w.data('data_obs').sumEntries(),r=r.getVal(),r_lo=r.getAsymErrorLo(),r_hi=r.getAsymErrorHi(),status_b=b.status(),covQual_b=b.covQual(),status_s=s.status(),covQual_s=s.covQual(),strategy=1 if era=='combined' else 0))
 names=list(dict.fromkeys([a.GetName() for a in b.floatParsFinal()]+[a.GetName() for a in s.floatParsFinal()]))
 for name in names:
  vb=b.floatParsFinal().find(name);vs=s.floatParsFinal().find(name);rows.append(dict(era=era,parameter=name,initial_own_era=w.var(name).getVal(),B_fit=vb.getVal() if vb else 0,SB_initial=sw.var(name).getVal(),SB_post=vs.getVal() if vs else None,SB_error=vs.getError() if vs else None,SB_error_lo=vs.getAsymErrorLo() if vs else None,SB_error_hi=vs.getAsymErrorHi() if vs else None))
 for e,ch in (eras if era=='combined' else [next(t for t in eras if t[0]==era)]):
  for stage,src in [('fit_b',fb),('fit_s',fs)]:
   dr=src.Get('shapes_'+stage+'/'+ch);h=dr.Get('total');g=dr.Get('data')
   with open(p/f'degree5_sequence_{era}_{ch}_{stage}_shapes.csv','w') as o:
    wr=csv.writer(o);wr.writerow(['mass','data','total','signal','dy','upsilon1s','upsilon2s'])
    for i in range(1,h.GetNbinsX()+1):
     width=h.GetBinWidth(i);wr.writerow([h.GetBinCenter(i),g.GetPointY(i-1)*width]+[dr.Get(k).GetBinContent(i)*width for k in ['total','Zd','dy','upsilon1s','upsilon2s']])
for name,arr in [('degree5_sequence_all_eras_summary',summ),('degree5_sequence_parameters',rows)]:
 (p/(name+'.json')).write_text(json.dumps(arr,indent=2))
 with open(p/(name+'.csv'),'w') as f:
  wr=csv.DictWriter(f,fieldnames=arr[0].keys());wr.writeheader();wr.writerows(arr)
lines=['# Grado5 — inizializzazione separata B → S+B per era','','## Procedura','','Range6–12GeV,450 bin originali/era. Coefficienti Chebyshev, normalizzazioni Chebyshev/Υ(1S)/Υ(2S) e nuisance floating con i vincoli originali. Coefficienti MC delle DCB costanti. Nessun parametro del fondo congelato e nessun vincolo aggiuntivo ricavato dai dati.','', 'Per ogni era: (1) costruzione di un grado5 positivo dalla proiezione del continuum originale di quella stessa era; (2) fit ai dati della sola somma Chebyshev+Υ(1S)+Υ(2S), mediante FitDiagnostics --skipSBFit; (3) salvataggio dei coefficienti e delle3 normalizzazioni del risultato fit_b nel workspace di partenza del fit S+B; (4) FitDiagnostics --skipBOnlyFit --preFitValue0, lasciando liberi tutti i parametri del fondo. Il punto iniziale S+B è quindi un fit B ai dati della stessa era, non il seed2022 copiato.','', 'I parametri del fit B sono usati solo come inizializzazione: non vengono aggiunti come prior e non si contano due volte i dati. Il file degree5_B_to_SB_initialization_audit.json documenta la corrispondenza esatta dei valori B/S+B iniziali.','', '## Risultati','','| Era | r (68% robustFit) | status B/S+B | covQual B/S+B |','|---|---:|---|---|']
for a in summ:lines.append(f"| {a['era']} | {a['r']:.5f} −{abs(a['r_lo']):.5f}/+{a['r_hi']:.5f} | {a['status_b']}/{a['status_s']} | {a['covQual_b']}/{a['covQual_s']} |")
lines+=['','Il combinato è inizializzato dai4 fit B separati. Con strategia0 il suo fit S+B ha covQual2; il controllo con strategia1/tolleranza0.001 produce covQual3 ed è quello riportato nella tabella e nei plot. Per i4 fit individuali si mantiene strategia0.','', 'Significanza locale combinata con strategia0: Z=0, p0=0.5. Nessun eccesso positivo a massa9 per questo modello. La scelta del grado non è validata dalla sola convergenza.','', '## Coefficienti e normalizzazioni: partenza B, risultato B, partenza S+B, risultato S+B','']
for era in [e for e,_ in eras]+['combined']:
 lines+=['### '+era,'','| Parametro | Seed era | Fit B | Inizio S+B | Post S+B |','|---|---:|---:|---:|---:|']
 for a in rows:
  if a['era']!=era:continue
  if not(a['parameter'].startswith('t') or 'bkgScale' in a['parameter'] or a['parameter']=='r'):continue
  val=f"{a['SB_post']:.7g} ± {a['SB_error']:.3g}" if a['SB_post'] is not None else '—'
  if a['parameter']=='r':val=f"{a['SB_post']:.6g} −{abs(a['SB_error_lo']):.4g}/+{a['SB_error_hi']:.4g}"
  lines.append(f"| `{a['parameter']}` | {a['initial_own_era']:.7g} | {a['B_fit']:.7g} | {a['SB_initial']:.7g} | {val} |")
 lines+=['']
(p/'degree5_sequence_fit_report.md').write_text('\n'.join(lines)+'\n');print(json.dumps(summ,indent=2))
