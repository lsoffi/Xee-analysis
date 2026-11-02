from pathlib import Path
import json,math,datetime,re
P=Path.cwd();R=P/'runs/20261030_upper_envelope_scan'
assert len(list(R.glob('degree*/m*/results.json')))==32
assert len(list(R.glob('bias/degree*/m*/*/m*/result.json')))==96
for degree in [1,2]:
    for k in range(105,121):
        rows=[json.loads(p.read_text()) for p in (R/f'gof/degree{degree}').glob(f'toy*/m{k/10:.1f}/result.json')]
        assert len(rows)==100 and set(q['toy'] for q in rows)==set(range(1,101)),(degree,k,len(rows))
        assert all(math.isfinite(q['score']) for q in rows)
import audit
import gof_histograms
import report
report.report()
a=json.loads((R/'validation.json').read_text());s=json.loads((R/'summary.json').read_text());now=datetime.datetime.now().astimezone().isoformat(timespec='seconds')
cfg=json.loads((R/'configuration.json').read_text());cfg.update(status='completed_with_documented_flags',completed_at=now);(R/'configuration.json').write_text(json.dumps(cfg,indent=2))
counts=s['counts'];bad_CR=[str(p.parent.relative_to(R)) for p in R.glob('gof/degree*/toy*/cr_result.json') if not json.loads(p.read_text())['valid']]
stats={str(d):dict(p_min=min(q['p'] for q in s['gof'] if q['degree']==d),p_max=max(q['p'] for q in s['gof'] if q['degree']==d)) for d in [1,2]}
meta=dict(status='completed_with_documented_flags',completed_at=now,counts=counts,gof_per_point=100,gof_CR_failed_datasets=bad_CR,gof_summary=stats,observed_FD_status0_cov3=a['all_observed_FD_status0_cov3'],CLs_observed_crossings=a['all_CLs_crossings_verified'],expected_bands='provisional',bias='Asimov only; neither degree uniformly passes 0.2 sigma',above11='parametric shape and extrapolated yields; no extrapolation uncertainty nuisance',source_audit='validation.json')
(R/'completion.json').write_text(json.dumps(meta,indent=2))
html=(R/'index.html').read_text()
for link in re.findall(r'(?:href|src)="([^"]+)"',html):
    if not link.startswith(('http','#')):assert (R/link).exists(),link
assert len(list(R.glob('degree*/m*/fit_SR.png')))==32 and len(list(R.glob('degree*/m*/fit_CR.png')))==32
readme=P/'README.md';text=readme.read_text();marker='## Scan SR alta 10,5–12: envelope gradi 1 e 2'
if marker not in text:
    section=f'''{marker}

[Report con tutti i fit e le curve](runs/20261030_upper_envelope_scan/index.html) · [Pacchetto segnale 11,1–12,0](runs/20261030_upper_envelope_scan/signal_extension_11p1_12p0.zip). Completati 32 punti osservati, 96 bias Asimov e 3200 score GoF (100 per massa/grado, toy riutilizzati fra masse). SR 10,506667–12, CR 9,093333–10,506667 Cheb4; μ Υ dalla sola CR poi fisse, rese libere condivise. Envelope Cheb1/Bern1 e Cheb2/Bern2 separati. Nessuno supera uniformemente il criterio di bias ±0,2σ. Bande attese provvisorie, 26 bias con flag; GoF p grado1 {stats['1']['p_min']:.2f}–{stats['1']['p_max']:.2f}, grado2 {stats['2']['p_min']:.2f}–{stats['2']['p_max']:.2f}, con {len(bad_CR)} dataset CR non convergenti conservati e segnalati. Sopra11: DCB parametriche originali, rese log-quadratiche 10–11 estrapolate, lnN relative a11; incertezza di estrapolazione non inclusa. Studio esplorativo, non validazione dei due gradi.

'''
    readme.write_text(text.replace('# Analisi statistica Xee\n\n','# Analisi statistica Xee\n\n'+section,1))
with (P/'AGENTS.md').open('a') as h:
    h.write(f'''\n## COMPLETATO scan SR alta envelope1/2 fino12 — {now}\n\nOutput runs/20261030_upper_envelope_scan/index.html, completion.json e validation.json. Tutti32FD/limits/significance,96Asimov,3200GoF=32×100 finiti;64FDstatus0cov3,32CLsobservedcrossverificati,32SR+32CRplot. Rese11.1–12 logquad10–11 ancorate11, DCBparamoriginali, originaliimmutati. Bias70/96numericalvalid,26flag; a r0 nessun punto valido passa.2σ perentrambigradi. GoF p1={stats['1']['p_min']:.2f}–{stats['1']['p_max']:.2f},p2={stats['2']['p_min']:.2f}–{stats['2']['p_max']:.2f}; CRfailed={bad_CR},mantenutidenominatore,failureboundsJSON. Toy031grado2ritentatoCRstatus0cov3,16scorerifatti,oldarchiviati. Nessunensemblebiasrichiesto(Asimovsolo). BandeexpectednativePROVVISORIE,nonreinterpretarecomevalidate; pilotminosnewfallitoNaNnonusato. 2flagFDvsZsono multimodalità grado1m11 egrado2m12, ramo fisicopositivoverificatoconstrainr0start10;negativoFDconservato. Risultati nominali/condizionali alleμfissate;nonraccomandaregrado1o2comevalidatidalbias. Reportinclude8GoFistogrammi,rese/shapeplot,CSV,ROOT/logcomandi,zipsegnali. Nonriavviarecampagna névecchirunner.\n''')
print(json.dumps(meta,indent=2))
