import ROOT,csv,json
from pathlib import Path
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
f=ROOT.TFile.Open('outputs/fitDiagnostics_inclusive_2022_chebyshev_to11.root');wf=ROOT.TFile.Open('outputs/Xee_ee_4_2022_chebyshev_to11.root');w=wf.Get('w');b=f.Get('fit_b');s=f.Get('fit_s')
names=['r']+['t%d_2022_envelope'%i for i in range(7)]+['CMS_EXO25020_bkgScale_'+p+'_2022' for p in ['dy','upsilon1s','upsilon2s']]+['CMS_scale_e','CMS_EXO25020_signalModelSigmaNuisance_2022','CMS_eff_e_id','CMS_eff_e_reco','CMS_eff_trigger_e','lumi_1']
rows=[]
for name in names:
 v=w.var(name);row={'parameter':name,'prefit':v.getVal(),'prefit_saved_error':v.getError()}
 for label,fit in [('b',b),('s',s)]:
  q=fit.floatParsFinal().find(name)
  if q:row.update({label:q.getVal(),label+'_err':q.getError(),label+'_lo':q.getAsymErrorLo(),label+'_hi':q.getAsymErrorHi()})
  else:row.update({label:0 if name=='r' else v.getVal(),label+'_err':None,label+'_lo':None,label+'_hi':None})
 rows.append(row)
with open('outputs/combine_2022_to11_pre_post.csv','w') as o:
 wr=csv.DictWriter(o,fieldnames=rows[0].keys());wr.writeheader();wr.writerows(rows)
Path('outputs/combine_2022_to11_pre_post.json').write_text(json.dumps({'b':{'status':b.status(),'covQual':b.covQual(),'edm':b.edm()},'s':{'status':s.status(),'covQual':s.covQual(),'edm':s.edm()},'parameters':rows},indent=2))
def fmt(row,stage):
 if row[stage+'_err'] is None:return f"{row[stage]:.5g} (fisso)"
 if row['parameter']=='r' and stage=='s':return f"{row[stage]:.5g} −{abs(row[stage+'_lo']):.4g}/+{row[stage+'_hi']:.4g}"
 return f"{row[stage]:.6g} ± {row[stage+'_err']:.3g}"
lines=['# Combine FitDiagnostics — 2022, solo Chebyshev','', 'Eseguito con Combine v11.1.0, commit 6c543b6. Range 6–11 GeV; segnale DCB a massa 9 GeV. Υ(1S) congelata alla normalizzazione originale.','', '| Parametro | Pre-fit | Post-fit B | Post-fit S+B |','|---|---:|---:|---:|']
for row in rows:lines.append(f"| `{row['parameter']}` | {row['prefit']:.7g} | {fmt(row,'b')} | {fmt(row,'s')} |")
lines+=['','Le colonne pre-fit riportano i valori iniziali del workspace, non prior per i parametri liberi. Le nuisance di scala/risoluzione/efficienza/luminosità hanno vincoli gaussiani standard; i coefficienti Chebyshev e le normalizzazioni libere non hanno questi vincoli.','',f"Fit B: status {b.status()}, covQual {b.covQual()}, EDM {b.edm():.5g}.",f"Fit S+B: status {s.status()}, covQual {s.covQual()}, EDM {s.edm():.5g}.",'','Le shape delle Υ rimangono fisse. Le nuisance CMS_scale_e e signalModelSigmaNuisance agiscono sul segnale.','', 'Il workspace derivato contiene solo ch1/2022 e il Chebyshev estratto dal RooMultiPdf originale; le altre ere e il Bernstein sono esclusi. I dati conservano i 375 conteggi selezionati, somma 33056. Il nome mH120 del file ausiliario è il default di Combine e non cambia la massa 9 codificata nelle formule del segnale.']
Path('outputs/combine_2022_to11_pre_post.md').write_text('\n'.join(lines)+'\n')
for obj in ['shapes_prefit','shapes_fit_b','shapes_fit_s']:
 dr=f.Get(obj+'/ch1');print(obj,[(k.GetName(),k.GetClassName()) for k in dr.GetListOfKeys()] if dr else 'missing')
for row in rows:print(row['parameter'],row['prefit'],fmt(row,'b'),fmt(row,'s'))
