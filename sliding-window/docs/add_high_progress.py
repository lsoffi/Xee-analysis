"""Rebuild section 6 from the checked-in dated partial-result snapshot."""
from pathlib import Path
from io import BytesIO
import csv,json,math
from datetime import datetime
from zoneinfo import ZoneInfo
from pypdf import PdfReader,PdfWriter
from reportlab.pdfgen.canvas import Canvas
from reportlab.lib.utils import simpleSplit,ImageReader
base=Path(__file__).resolve().parent;data=base/'data/high_partial'
status=json.loads((data/'status.json').read_text());snap=json.loads((data/'high_note_snapshot.json').read_text())
rows=list(csv.DictReader((data/'bias.csv').open()));lim=[{k:float(v) for k,v in r.items()} for r in csv.DictReader((data/'limits.csv').open()) if float(r['mass'])>=4.1]
stamp=datetime.strptime(status['updated'],'%Y-%m-%d %H:%M UTC').replace(tzinfo=ZoneInfo('UTC')).astimezone(ZoneInfo('Europe/Madrid')).strftime('%d November %Y, %H:%M CET')
# Run plot_high_limits.py first to regenerate the limit figure.
assert (data/'limits_high.png').exists()
buf=BytesIO();c=Canvas(buf,pagesize=(842,650));page=14

def title(t):c.setFillColorRGB(.09,.22,.31);c.setFont('Helvetica-Bold',21);c.drawString(36,608,t)
def para(y,t):
 c.setFont('Helvetica',10.5)
 for line in simpleSplit(t,'Helvetica',10.5,765):c.drawString(38,y,line);y-=14
 return y-10

def end():
 global page
 c.setFont('Helvetica',9);c.drawString(38,20,'PRELIMINARY | Snapshot: '+stamp);c.drawRightString(804,20,str(page));c.showPage();page+=1

def plot(name):
 im=ImageReader(str(data/name));w,h=im.getSize();s=min(768/w,375/h);c.drawImage(im,(842-w*s)/2,212,w*s,h*s)

title('6   High SR: ongoing studies above the J/psi')
y=para(571,'Scope: 50 hypotheses from 4.1 to 9.0 GeV, in 0.1 GeV steps. C4+L4 is a candidate envelope, not yet a final validated choice. The four-era simultaneous SR+selected-CR likelihood follows Sections 1-3, with signal in all channels, profiled trigger uncertainty and total-density C0 at actual common boundaries.')
y=para(y,'Completion snapshot: '+stamp+'. Counts refer to saved outputs, not a claim that every fit or limit is statistically validated. Figures on the following pages use the same partial-results snapshot.')
for t in [f"FitDiagnostics: {50-len(snap['missing_fits'])}/50 reports. Missing masses: "+', '.join(map(str,snap['missing_fits']))+' GeV.',f"Asymptotic limits: {len(lim)}/50 high-region points. Missing masses: "+', '.join(map(str,snap['missing_limits']))+' GeV.',f"Bias: {status['bias_received']}/40000 toy results received; {status['bias_valid']} valid, {status['bias_invalid']} invalid, {40000-status['bias_received']} still missing. Each of 50 masses has two generators (C4, L4), two injections (zero, median expected) and 200 toys per configuration.",f"Upsilon-boundary diagnosis: {snap['upsilon_results']}/108 outputs saved; {108-snap['upsilon_results']} pending. Conditional likelihood profiles and single-point CLs checks target 8.6 and 9.0 GeV."]:
 y=para(y,t)
y=para(y,'Recovery is in progress for the missing fits, limits and generator prefits at 8.5/8.8 GeV. Toy recovery reuses the original seeds and excludes active jobs and saved results. One Condor job per toy is retained. Numerical retries do not relax the model or validity criteria.')
y=para(y,'No high-region GoF or coverage conclusion is available. The pending programme is to finish and validate the fits and injections, resolve the observed-limit anomalies and then determine whether the candidate envelope is suitable.')
assert y>35;y=None;end()
title('6.1   Preliminary high-region limits');plot('limits_high.png')
y=para(192,'48/50 high-region limits are available; 7.7 and 8.5 GeV are missing and left as gaps. Expected bands are the native AsymptoticLimits quantiles. All included outputs contain the six finite quantiles with ordered expected bands; these checks do not independently validate the profile minima or CLs crossings.')
y=para(y,'The unusually high observed limits at 8.6 and 9.0 GeV remain under investigation. They are shown without clipping or replacement. These values alone establish neither a physical excess nor a trigger-related explanation. The high-region envelope is still undergoing bias and numerical validation.');assert y>35;end()
for tag,img,label in [('r0','bias_r0.png','zero signal'),('expected','bias_expected.png','median expected signal')]:
 rr=[r for r in rows if r['tag']==tag];received=sum(int(r['received']) for r in rr);valid=sum(int(r['valid']) for r in rr)
 complete=sum(int(r['received'])==200 for r in rr);plotted=sum(int(r['valid'])>=20 for r in rr)
 fullm=[i/10 for i in range(41,91) if sum(int(r['received']) for r in rr if float(r['mass'])==i/10)==400]
 empty=[i/10 for i in range(41,91) if not any(int(r['received']) for r in rr if float(r['mass'])==i/10)]
 title('6.'+('2' if tag=='r0' else '3')+'   Bias with '+label);plot(img)
 y=para(192,f'{received}/20000 toys received; {valid} valid and {received-valid} invalid; {20000-received} pending. Complete toy counts at {len(fullm)}/50 masses (both generators), leaving {50-len(fullm)} masses incomplete. {complete}/100 generator-mass samples have all 200 receipts; {plotted}/100 have at least 20 valid toys for plotting. Masses with no receipts: '+str(len(empty))+' (listed by mass in bias.csv).' )
 y=para(y,'Bias = mean(rhat - rinj) / sample SD. Both fitted branches must be valid; the winner minimises D+16. Bars are pointwise 95% bootstrap intervals. Open symbols indicate incomplete samples; full symbols indicate all 200 receipts, not necessarily 200 valid fits. The +/-0.2 lines are diagnostic targets, not confidence limits.')
 y=para(y,'This partial sample does not establish closure. Completion may change the estimates, and invalid fits remain reported separately. The measured trigger nuisance stays profiled, including in the injected-signal study. See bias.csv for per-mass counts and exclusions.');assert y>35;end()
c.save();w=PdfWriter();r=PdfReader(base/'Xee_analysis_note.pdf')
for p in r.pages[:13]:w.add_page(p)
for p in PdfReader(buf).pages:w.add_page(p)
for t,p in [('Common methodology',0),('4 Low SR',3),('5 Middle SR',5),('6 High SR - preliminary status',13),('6.1 High limits',14),('6.2 High bias: r=0',15),('6.3 High bias: injected signal',16)]:w.add_outline_item(t,p)
w.add_metadata({'/Title':'Xee analysis: low, middle and ongoing high SR studies'})
w.write(base/'Xee_analysis_note.pdf');print('Saved 17 pages:',stamp)
