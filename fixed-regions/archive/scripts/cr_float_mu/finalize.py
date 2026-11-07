from run import *
import hashlib
F=json.loads((R/'fits.json').read_text());rows=[json.loads(l) for l in (R/'pilot.jsonl').read_text().splitlines()];truth=json.loads((old/'checks.json').read_text())['truth_alternative'];out=[]
for idx in [2,32]:
 for kind in ['gof','bias','asimov']:
  a=[r for r in rows if r['idx']==idx and r['kind']==kind];toys=[]
  for k in sorted(set(r['toy'] for r in a)):
   fs=[r for r in a if r['toy']==k];assert len(fs)==4;toys.append(dict(toy=k,valid=all(f['valid'] for f in fs),boundary=any(any(f['boundary'][2:]) for f in fs),D=sum(f['D'] for f in fs),prediction=np.sum([f['prediction'] for f in fs],axis=0).tolist(),error=np.sqrt(np.sum(np.array([f['prediction_error'] for f in fs])**2,axis=0)).tolist()))
  good=[t for t in toys if t['valid']];rec=dict(idx=idx,kind=kind,total=len(toys),valid=len(good),shape_boundary=sum(t['boundary'] for t in toys),toys=toys)
  if kind=='gof':
   obs=[f for f in F if f['idx']==idx];rec.update(observed_valid=all(f['valid'] for f in obs),observed_D=sum(f['D'] for f in obs),p_reported=None,reason='Pilot only; validity assessed below')
  elif kind=='bias' and good:
   p=np.array([t['prediction'] for t in good]);err=np.mean([t['error'] for t in good],axis=0);rec.update(conditional_bias_sigma=((p.mean(axis=0)-truth)/err).tolist(),truth=truth,interpretation='Conditional on successful fits; not validated bias/coverage')
  out.append(rec)
(R/'pilot_summary.json').write_text(json.dumps(out,indent=2));assert hashlib.sha256(Path(I['source']).read_bytes()).hexdigest()==I['sha256']
with (R/'REPORT.md').open('a') as f:
 f.write('\n\n## Pilot rapido dei toy\nPer rispettare la richiesta di rapidità: 50 toy GoF e 50 toy dal generatore alternativo Bernstein4 per ciascuna delle due finestre, grado3. Non ripetuta la campagna200toy×15configurazioni precedente. Stessi primi50 pseudo-dataset alternativi precedenti, ora rifittati con μ liberi. Due seed per ogni fit; inclusi due Asimov.\n')
 for r in out:f.write(f"\n- Configurazione {r['idx']} ({S[r['idx']]['lo']:.3f}–{S[r['idx']]['hi']:.3f}), {r['kind']}: {r['valid']}/{r['total']} toy con tutte4ere valide; {r['shape_boundary']} con almeno μ al bordo.")
 f.write('\n\nNon si riportano p GoF calibrati né si certifica bias/copertura: fallimenti e bordi impediscono una conclusione affidabile. Le medie condizionate salvate nel JSON sono diagnostiche, non criteri di scelta della CR. La maggiore flessibilità non risolve automaticamente la degenerazione. Per misurare entrambe le shape servono verifiche di identificabilità, profili e/o informazioni esterne sui parametri scarsamente determinati; non fissare i nuovi valori instabili nel fit di ricerca.\n')
completion=dict(observed_era_fits=192,observed_valid=sum(f['valid'] for f in F),observed_shape_boundary=sum(any(f['boundary'][2:]) for f in F),pilot_era_fits=len(rows),input_unchanged=True,baseline_unchanged=True,scope='48 observed configurations plus 50-toy pilot in two windows degree3; no full repeat of previous200toy campaign',decision='No validated narrower CR or new Upsilon shape measurement');(R/'completion.json').write_text(json.dumps(completion,indent=2))

for r in out:print({k:v for k,v in r.items() if k!='toys'})
