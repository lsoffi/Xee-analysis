import csv,json,os,sys,time
from pathlib import Path
import ROOT
ROOT.gROOT.SetBatch(True)
base=Path.cwd()
web=Path('/eos/user/s/soffi/www/Xee/261106/scan_SR_CR/low_logpoly3')
web.mkdir(exist_ok=True)
eras=['2022','2022EE','2023','2023BPix']
def extract(path,metadata):
    f=ROOT.TFile.Open(str(path))
    records=[]
    for stage in ['b','s']:
        fit=f.Get('fit_'+stage)
        chisq=0.
        n=0
        for ch in metadata['channels']:
            if ch['region']!='SR0':continue
            h=f.Get('shapes_fit_'+stage+'/'+ch['label']+'/total')
            expected=[h.GetBinContent(i+1) for i in range(len(ch['counts']))]
            assert min(expected)>0
            chisq+=sum((y-v)**2/v for y,v in zip(ch['counts'],expected))
            n+=len(expected)
        records.append((stage,fit.status(),fit.covQual(),round(chisq,2),n))
    f.Close()
    return records
for attempt in range(720):
    rows=[]
    for m in ['0.4','0.5','0.6']:
        p=base/'low_logpoly3/results'/('m'+m)
        if not (p/'model.json').exists() or not (p/'FitDiagnostics.done').exists():continue
        md=json.loads((p/'model.json').read_text())
        for name,d in [('Cheb6',base/'results'/('m'+m)),('logpoly3',p)]:
            for values in extract(next(d.glob('fitDiagnostics*.root')),md):rows.append((m,name,*values))
    if rows:
        tmp=web/'confronto.tmp'
        with tmp.open('w') as output:
            w=csv.writer(output);w.writerow(['mass_GeV','model','fit','status','covQual','chi2_Pearson_SR_4eras','n_bins_SR']);w.writerows(rows)
        tmp.replace(web/'confronto.csv')
    if len(rows)==12:break
    time.sleep(30)
