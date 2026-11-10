"""Insert/update low validation pages. Run last, after legacy section builders and add_todos.py."""
from pathlib import Path
from io import BytesIO
import json
from pypdf import PdfReader,PdfWriter
from reportlab.pdfgen.canvas import Canvas
from reportlab.lib.utils import simpleSplit,ImageReader
base=Path(__file__).resolve().parent;data=base/'data/low_validation';d=json.load(open(data/'results.json'))
buf=BytesIO();c=Canvas(buf,pagesize=(842,650))
def title(t):
 c.setFillColorRGB(.09,.22,.31);c.setFont('Helvetica-Bold',20);c.drawString(36,610,t)
def para(y,t):
 c.setFont('Helvetica',10.5)
 for line in simpleSplit(t,'Helvetica',10.5,765):c.drawString(38,y,line);y-=14
 return y-9
def image(name,top,height):
 im=ImageReader(str(data/name));w,h=im.getSize();s=min(768/w,height/h);c.drawImage(im,(842-w*s)/2,top-h*s,w*s,h*s)
def end():
 c.setFont('Helvetica',9);c.drawString(38,20,'Low SR validation | Complete toy samples | 10 November 2026');c.showPage()
title('4.2   Low SR: coverage with the trigger profiled')
y=para(578,'C3+L3 envelope; four-era simultaneous SR+selected-CR likelihood. At each mass, inject the median expected signal strength and generate 200 toys with C3 and 200 with L3 (1,200 total). The measured trigger uncertainty remains profiled in the inference.')
image('coverage_low.png',535,310)
y=para(214,'For each toy, evaluate CLs at its injected signal strength. Coverage here is the fraction with CLs >= 0.05: toys that do not exclude the injected signal. Error bars are exact two-sided 95% binomial intervals, not uncertainty bands on the signal limit.')
y=para(y,'C3 / L3 non-exclusion fractions: 0.4 GeV: 189/200 / 189/200 (94.5% / 94.5%); 0.5 GeV: 189/200 / 192/200 (94.5% / 96.0%); 0.6 GeV: 193/200 / 188/200 (96.5% / 94.0%). All six intervals contain 95%.')
y=para(y,'No undercoverage is resolved at this toy precision for these generators and injections. This does not establish coverage for arbitrary background shapes, other signal strengths or model misspecification. The coverage result does not override the GoF tension on the following page.')
assert y>35;end()
title('4.3   Low SR: SR-only goodness of fit')
y=para(578,'Native Combine saturated GoF, with r free in the simultaneous SR+CR fit. All channels enter the fit; the statistic is evaluated in the SR only using the fit/evaluation masks. Use 500 frequentist toys per mass (1,500 total), one Condor job per toy, and the C3+L3 envelope.')
image('gof_low.png',520,255)
y=para(259,'The red line is the observed SR statistic; histograms show the toy distribution. The local p-value is the fraction of toys with a statistic at least as large as observed. Orange bins mark this tail; intervals quantify finite-toy precision.')
for r in d['gof']:
 lo,hi=r['p_CI95'];y=para(y,f"{r['mass']:.1f} GeV: observed statistic {r['observed']:.3f}; p = {r['tail']}/500 = {r['p']:.3f}; exact 95% binomial interval [{lo:.3f}, {hi:.3f}].")
y=para(y,'The 0.4 and 0.5 GeV fits show no GoF tension. At 0.6 GeV, p=0.006 flags poor agreement with the tested model: low-SR validation is not uniformly satisfactory. This is a local goodness-of-fit discrepancy, not a signal significance; overlapping windows are correlated. Residual localisation and numerical checks are needed before assigning its cause.')
assert y>35;end()
title('4.4   Low SR: goodness of fit versus mass')
image('gof_low_vs_mass.png',575,395)
y=para(166,'Summary of the same 500-toy samples per mass shown on the preceding page. Points are local SR-only GoF p-values, with exact 95% binomial intervals. The dashed line marks p=0.05; it is a reference threshold, not a global significance criterion.')
y=para(y,'The p-values are 0.458, 0.534 and 0.006 at 0.4, 0.5 and 0.6 GeV, respectively. The 0.6 GeV interval [0.001, 0.017] lies below 0.05. This flags local model disagreement and does not identify its cause or establish a signal excess. The mass windows are correlated.')
assert y>35;end();c.save()
r=PdfReader(base/'Xee_analysis_note.pdf');old=[p for p in r.pages if not any(t in (p.extract_text() or '') for t in ['4.2   Low SR: coverage with the trigger profiled','4.3   Low SR: SR-only goodness of fit','4.4   Low SR: goodness of fit versus mass'])]
new=list(PdfReader(buf).pages);pages=old[:5]+new+old[5:];w=PdfWriter()
for i,p in enumerate(pages):
 width=float(p.mediabox.width);ov=BytesIO();cc=Canvas(ov,pagesize=(width,float(p.mediabox.height)));cc.setFillColorRGB(1,1,1);cc.rect(width-53,10,45,23,fill=1,stroke=0);cc.setFillColorRGB(.09,.22,.31);cc.setFont('Helvetica',9);cc.drawRightString(width-38,20,str(i+1));cc.save();p.merge_page(PdfReader(ov).pages[0]);w.add_page(p)
for t,i in [('Common methodology',0),('4 Low SR',3),('4.2 Low coverage',5),('4.3 Low goodness of fit',6),('4.4 Low GoF versus mass',7),('5 Middle SR',8),('6 High SR - dated preliminary status',16),('7 TODO',20)]:w.add_outline_item(t,i)
w.add_metadata({'/Title':'Xee analysis: low, middle and high SR validation'})
w.write(base/'Xee_analysis_note.pdf');print('Saved',len(pages),'pages')
