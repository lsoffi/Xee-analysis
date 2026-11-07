import ROOT,sys
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
f=ROOT.TFile.Open(sys.argv[1]);idx=int(sys.argv[2]);out=ROOT.TFile.Open(sys.argv[3],'RECREATE');out.mkdir('toys');out.cd('toys')
for suffix in ['', '_snapshot']:
 obj=f.Get(f'toys/toy_{idx}{suffix}');assert obj;obj.Write(f'toy_1{suffix}')
out.Close()
