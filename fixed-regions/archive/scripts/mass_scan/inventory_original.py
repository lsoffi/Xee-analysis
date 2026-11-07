import ROOT,pathlib,json,hashlib,sys
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');ROOT.gROOT.SetBatch(True)
base=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee/data/mass_scan/261028/cards_region2_data_envelope_allCorrections_binned/ee');out=pathlib.Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True);items=[]
for d in sorted([d for d in base.iterdir() if d.is_dir() and d.name.replace('.','',1).isdigit()],key=lambda d:float(d.name)):
 for era in ['2022','2022EE','2023','2023BPix','allYears']:
  path=d/f'Xee_ee_4_{era}.root';f=ROOT.TFile.Open(str(path));w=f.Get('w');x=w.var('mass');data=w.data('data_obs');cats=w.allCats();discrete=[a.GetName() for a in cats if 'Envelope' in a.GetName()];entry=dict(mass=float(d.name),mass_label=d.name,era=era,path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),fit_range=[x.getMin(),x.getMax()],bins=x.getBins(),entries=data.numEntries(),sum_weights=data.sumEntries(),model_class=w.pdf('model_s').ClassName(),discrete_parameters=discrete)
  assert entry['fit_range']==[6.,12.],entry
  items.append(entry);f.Close()
(out/'input_inventory.json').write_text(json.dumps(items,indent=2));print('Audited',len(items),'original workspaces;',len(items)//5,'mass points. All mass ranges 6–12.')
