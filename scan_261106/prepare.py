import re
import csv
import json
import math
import shutil
from pathlib import Path
import ROOT

base = Path.cwd()
web = Path('/eos/user/s/soffi/www/Xee/261106')
source = Path('/eos/user/n/npalmeri/shared')
eras = ['2022','2022EE','2023','2023BPix']
records = []
with (web/'tabella_SR_CR_5MeV.csv').open() as f:
    for row in csv.DictReader(f):
        mass = float(row['Massa [GeV]'])
        low, high = map(float,re.findall(r'[0-9]+(?:\.[0-9]+)?',row['SR effettiva [GeV]']))
        assert round(200*(high-low)) == int(row['Bin SR'])
        records.append(dict(eligible=True,mass=mass,sr_segments=[[low,high]],selected_CRs=[s.strip() for s in row['CR scelte'].split(',')]))
assert len(records)==69
(base/'records.json').write_text(json.dumps(records))
with (base/'bins.csv').open('w') as f:
    writer=csv.writer(f);writer.writerow(['era','bin_low_GeV','bin_high_GeV','events'])
    for era in eras:
        path,=(source/'datasets_261106/261106'/era).glob('*.root')
        fin=ROOT.TFile.Open(str(path));data=fin.Get('w').data('data_obs')
        for i in range(data.numEntries()):
            center=data.get(i).getRealValue('mass')
            edge=round(200*center-.5)
            assert abs(center-(edge+.5)/200)<1e-6
            writer.writerow([era,edge/200,(edge+1)/200,data.weight()])
        fin.Close()
for src in (base.parent/'cr_261106').glob('*/*/display_fit.json'):
    q=json.loads(src.read_text());p=q['fit']['values'];q['physical_parameters']=[]
    for k,(name,_) in enumerate(q['cr']['resonances']):
        nominal=q['signal_parameters'][k]
        shape={key:nominal[key] for key in ('alphaL','alphaR','nL','nR')}
        if q['tail_count']:
            assert q['cr']['fractional_n'] and q['tail_count']==8
            for j,key in enumerate(shape):shape[key]*=math.exp(p[7+4*k+j])
        shape.update(mean=nominal['mean']+nominal['sigmaL']*p[3+k],sigma=nominal['sigmaL']*math.exp(p[5+k]))
        q['physical_parameters'].append(dict(name=name,fitted=shape,events_in_CR=math.exp(p[1+k])))
    dst=base/'cr'/src.parent.parent.name/src.parent.name/src.name
    dst.parent.mkdir(parents=True,exist_ok=True);dst.write_text(json.dumps(q))
assert len(list((base/'cr').glob('*/*/display_fit.json')))==12
for rec in records:
    m=rec['mass']; reg=0 if m<2.4 else 1 if m<7 else 2
    rel=Path(f'cards_region{reg}_data_envelope_allCorrections_binned/ee/{m:.1f}')
    dst=base/'signal'/rel;dst.mkdir(parents=True,exist_ok=True)
    for era in eras:
        name=f'Xee_ee_4_{era}.root'
        shutil.copyfile(source/'datacards_261021/261021'/rel/name,dst/name)
(base/'masses.txt').write_text(''.join(f"{r['mass']:.1f}\n" for r in records))
