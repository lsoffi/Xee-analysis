import ROOT,csv,json
from pathlib import Path
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');ROOT.gROOT.SetBatch(True)
ERAS=[('2022','ch1'),('2022EE','ch2'),('2023','ch3'),('2023BPix','ch4')];summary=[];parameters=[]
for era in [e for e,_ in ERAS]+['combined']:
 f=ROOT.TFile.Open(f'outputs/fitDiagnostics_degree5_native_frozenTail_{era}.root');wf=ROOT.TFile.Open(f'outputs/Xee_degree5_native_frozenTail_{era}.root');w=wf.Get('w');b=f.Get('fit_b');s=f.Get('fit_s');r=s.floatParsFinal().find('r');summary.append(dict(era=era,events=w.data('data_obs').sumEntries(),r=r.getVal(),r_lo=r.getAsymErrorLo(),r_hi=r.getAsymErrorHi(),status_b=b.status(),covQual_b=b.covQual(),edm_b=b.edm(),status_s=s.status(),covQual_s=s.covQual(),edm_s=s.edm()))
 for stage,fit in [('fit_b',b),('fit_s',s)]:
  for v in fit.floatParsFinal():
   initial=w.var(v.GetName());parameters.append(dict(era=era,fit=stage,parameter=v.GetName(),prefit=initial.getVal(),postfit=v.getVal(),error=v.getError(),error_lo=v.getAsymErrorLo(),error_hi=v.getAsymErrorHi()))
 for _,ch in (ERAS if era=='combined' else [next(pair for pair in ERAS if pair[0]==era)]):
  for stage in ['prefit','fit_b','fit_s']:
   dr=f.Get('shapes_'+stage+'/'+ch);h=dr.Get('total');g=dr.Get('data');n=0
   with open(f'outputs/degree5_frozenTail_{era}_{ch}_{stage}_shapes.csv','w') as o:
    wr=csv.writer(o);wr.writerow(['mass','data','total','signal','dy','upsilon1s','upsilon2s'])
    for i in range(1,h.GetNbinsX()+1):
     width=h.GetBinWidth(i);count=g.GetPointY(i-1)*width;n+=count;wr.writerow([h.GetBinCenter(i),count]+[dr.Get(k).GetBinContent(i)*width for k in ['total','Zd','dy','upsilon1s','upsilon2s']])
   assert h.GetNbinsX()==450
for name,rows in [('degree5_frozenTail_all_eras_summary',summary),('degree5_frozenTail_all_eras_parameters',parameters)]:
 Path('outputs/'+name+'.json').write_text(json.dumps(rows,indent=2))
 with open('outputs/'+name+'.csv','w') as out:
  wr=csv.DictWriter(out,fieldnames=rows[0].keys());wr.writeheader();wr.writerows(rows)
lines=['# Chebyshev grado5 — tutte le ere, binning originale','','Range6–12GeV,450 bin per era; segnale DCB a massa9; Υ(1S) congelata; t4 congelato per2022, t3 congelato perlealtreere; altri coefficienti liberi; r tra−30 e30; strategia0,robustFit1. Fit individuali e combinato con r e nuisance comuni condivisi secondo il workspace originale.','', '| Era | Eventi | r (68% robustFit) | status B/S+B | covQual B/S+B |','|---|---:|---:|---|---|']
for a in summary:lines.append(f"| {a['era']} | {a['events']:.0f} | {a['r']:.5f} −{abs(a['r_lo']):.5f}/+{a['r_hi']:.5f} | {a['status_b']}/{a['status_s']} | {a['covQual_b']}/{a['covQual_s']} |")
lines+=['','Per inizializzare i nuovi Chebyshev di grado5 in modo positivo, in tutte le ere sono stati usati i primi5 coefficienti originali2022 come valori iniziali; non sono stati applicati vincoli aggiuntivi. Le normalizzazioni, le shape delle Υ e del segnale, e le nuisance sono quelle della rispettiva era. Gli errori dei coefficienti iniziali non costituiscono prior.','', 'La dipendenza dal grado è un controllo di sensibilità del modello, non una dimostrazione che un eccesso sia spurio o che il grado5 sia validato. Gli intervalli di r non sono significatività del segnale.']
Path('outputs/degree5_frozenTail_all_eras_report.md').write_text('\n'.join(lines)+'\n');print(json.dumps(summary,indent=2))
