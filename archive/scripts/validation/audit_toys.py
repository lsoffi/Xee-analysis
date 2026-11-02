import ROOT,pathlib,json,math
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee/runs/20261029_background_validation');rows=[]
for meta in sorted(P.glob('toys/m*/*/*/batch_*/scenario.json')):
 m=json.loads(meta.read_text());p=next((meta.parent/'generate').glob('higgsCombine*.root'));f=ROOT.TFile.Open(str(p))
 for i in range(1,11):
  d=f.Get(f'toys/toy_{i}');g=f.Get(f'toys/toy_{i}_snapshot');assert d and g;assert d.numEntries()==1800,(str(p),i,d.numEntries());assert g.getSize()==10
  for j in range(d.numEntries()):d.get(j);assert d.weight()>=0 and abs(d.weight()-round(d.weight()))<1e-9
  rows.append(dict(mass=m['mass'],generator=m['generator'],injection=m['injection'],seed=m['seed'],toy=i,events=d.sumEntries(),bins=d.numEntries(),global_observables=g.getSize()))
assert len(rows)==1000,len(rows)
(P/'toy_integrity_audit.json').write_text(json.dumps(dict(datasets_verified=len(rows),bins_each=1800,global_observables_each=10,all_counts_nonnegative_integers=True,min_events=min(r['events'] for r in rows),max_events=max(r['events'] for r in rows),datasets=rows),indent=2));print('Verified',len(rows),'toy datasets')
