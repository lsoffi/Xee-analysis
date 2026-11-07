"""Read original observed bins; no fits and no source modifications."""
from pathlib import Path
import ROOT
import json
import hashlib
import csv

ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
P = Path(__file__).resolve().parents[1]
out = P / 'runs/20261030_data_above_10p5'
out.mkdir(parents=True, exist_ok=True)
src = P / 'data/mass_scan/261028/cards_region2_data_envelope_allCorrections_binned/ee/7.5/Xee_ee_4_allYears.root'
sha = hashlib.sha256(src.read_bytes()).hexdigest()
f = ROOT.TFile.Open(str(src))
w = f.Get('w')
x = w.var('mass')
data = w.data('data_obs')
lo, hi, nb = x.getMin(), x.getMax(), x.getBins()
bw = (hi-lo)/nb
eras = ['2022', '2022EE', '2023', '2023BPix']
rows = []
for k in range(data.numEntries()):
    row = data.get(k)
    ch = row.find('CMS_channel').getLabel()
    m, n = row.getRealValue('mass'), data.weight()
    assert abs(n-round(n)) < 1e-8 and n >= 0
    i = round((m-lo)/bw-.5)
    assert abs(m-(lo+(i+.5)*bw)) < 1e-8
    rows.append(dict(era=eras[int(ch[2:])-1], bin=i,
                     low=lo+i*bw, high=lo+(i+1)*bw, center=m, count=int(round(n))))
cr_end = 10.506666666666668
summary = []
for era in eras:
    a = sorted([r for r in rows if r['era']==era], key=lambda r:r['bin'])
    assert len(a)==nb and len({r['bin'] for r in a})==nb
    populated = [r for r in a if r['count']>0]
    upper = [r for r in a if r['low']>=cr_end-1e-9]
    summary.append(dict(era=era,total=sum(r['count'] for r in a),
                        events_above_cr=sum(r['count'] for r in upper),
                        upper_bins=len(upper), empty_upper_bins=sum(r['count']==0 for r in upper),
                        last_populated_bin=populated[-1]))
info = dict(source=str(src), sha256=sha, observable_range=[lo,hi], bins_per_era=nb,
            bin_width=bw, cr_end=cr_end, dataset_class=data.ClassName(),
            total_events=data.sumEntries(), rows=rows, summary=summary,
            note='Binned observed counts; no event-by-event maximum is recoverable. The workspace ends at 12 GeV.')
f.Close()
assert hashlib.sha256(src.read_bytes()).hexdigest()==sha
(out/'spectrum.json').write_text(json.dumps(info, indent=2))
with (out/'spectrum.csv').open('w') as stream:
    writer=csv.DictWriter(stream, fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
print(json.dumps(summary, indent=2))
