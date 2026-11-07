import json,csv,matplotlib,sys
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
base=Path('/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee/scan_261106')
r=json.loads((base/'records.json').read_text()); m=[x['mass'] for x in r]
d={(row['era'],round(float(row['mass']),1)):float(row['width']) for row in csv.DictReader(open(sys.argv[1]))}
fig,ax=plt.subplots(figsize=(9,5))
def line(y,label,**kw):
    for i in range(len(m)):
        if i==0 or m[i]-m[i-1]>.11: start=i
        if i+1==len(m) or m[i+1]-m[i]>.11: ax.plot(m[start:i+1],y[start:i+1],label=label if start==0 else None,**kw)
for era,c in zip(['2022','2022EE','2023','2023BPix'],['#0072b2','#d55e00','#009e73','#ad5092']):
    line([d[era,x] for x in m],era,color=c,marker='.',ms=3)
line([sum(b-a for a,b in row['sr_segments']) for row in r],'SR effettiva comune',color='black',ls='--',lw=1.4)
ax.set(xlabel='Ipotesi di massa [GeV]',ylabel='Larghezza della finestra [GeV]',title='Finestre al 95% della DCB per era e SR usata nel fit')
ax.legend();ax.grid(alpha=.2);fig.tight_layout()
fig.savefig('/eos/user/s/soffi/www/Xee/261106/scan_SR_CR/larghezze_finestre_95.png',dpi=160)
