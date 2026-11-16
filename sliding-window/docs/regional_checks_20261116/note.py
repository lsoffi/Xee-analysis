from pathlib import Path
from io import BytesIO
import sys
from pypdf import PdfReader,PdfWriter
from reportlab.pdfgen.canvas import Canvas
from reportlab.lib.utils import simpleSplit
base=Path(__file__).resolve().parent
p=Path(sys.argv[1]);reader=PdfReader(p);assert len(reader.pages)==35
buf=BytesIO();c=Canvas(buf,pagesize=(842,650));titles=[]
def page(title,sub,img,paras):
 titles.append(title);c.setFillColorRGB(.10,.32,.39);c.setFont('Helvetica-Bold',11);c.drawString(40,614,sub.upper());c.setFillColorRGB(.07,.12,.16);c.setFont('Helvetica-Bold',21);c.drawString(40,578,title)
 if img:c.drawImage(str(base/(img+'.png')),40,305,760,255,preserveAspectRatio=True)
 y=275 if img else 523
 for label,body in paras:
  c.setFont('Helvetica-Bold',12);c.drawString(42,y,label);y-=18;c.setFont('Helvetica',12)
  for line in simpleSplit(body,'Helvetica',12,750):c.drawString(42,y,line);y-=16
  y-=15
 assert y>33
 c.setFont('Helvetica',9);c.drawString(40,24,'Xee  /  Validation inventory  /  16 November 2026');c.drawRightString(800,24,str(35+len(titles)));c.showPage()
page('Where each regional validation stands','Reader guide to Sections 4-8',None,[
('Low SR | 0.4-0.6 GeV | C3+L3 used for limits','Bias at zero and median expected injection: 3/3 masses each (p. 37). SR-only GoF: 3/3 (pp. 7-8). Coverage: 3/3 (p. 6). The local discrepancy at 0.6 GeV remains open.'),
('Middle SR | 1.2-2.7 GeV | C4+L4','Bias at both injections: 16/16 masses each (p. 38). SR-only GoF: 16/16 (pp. 15-16). Coverage: 16/16 (pp. 13-14). The trigger comparison uses paired seeds (p. 31).'),
('High SR | 4.1-9.0 GeV | provisional C4+L4','Bias at both injections: 50/50 masses each (p. 39). SR-only GoF: 33/50 masses have estimates (p. 33). Coverage: 50/50 (p. 21). The separate C3+L3 comparison covers positive injection only (p. 32).'),
('How to read the counts','Bias is the mean fitted-minus-injected signal strength divided by its sample standard deviation on accepted toys. The green band marks the working target. GoF tests SR bins after fitting all selected channels. Coverage counts how often the CLs upper limit retains the injected strength. Counts refer to available mass points; fit validity and outcomes are given separately.')])
page('Low SR | the two signal injections','Section 4 validation detail','low',[
('Bias','The C3+L3 envelope is evaluated on toys valid for all three candidate envelopes in the paired comparison. At zero injection, 0/12 central estimates exceed |0.2| (maximum 0.161). At positive injection, 8/12 exceed it (maximum 0.271). Each count spans three masses and four generators.'),
('GoF and coverage','With 500 GoF toys per mass, the local p-values are 0.458, 0.534 and 0.006 at 0.4, 0.5 and 0.6 GeV. All six exact coverage intervals include 95% at the mass-specific median expected injections. The spectral tension at 0.6 GeV remains relevant even with satisfactory coverage.')])
page('Middle SR | the two signal injections','Section 5 validation detail','middle',[
('Bias','The quartic background-only study has 6,359/6,400 valid toys and 0/32 central bias estimates outside the target. At positive injection 6,277/6,400 toys pass and 28/32 estimates exceed the target. At 1.4 GeV the paired fixed-trigger comparison identifies a sizeable profiling effect (p. 31).'),
('GoF and coverage','Local GoF p-values fall below 0.05 at 2.2 (0.024), 2.6 (0.034) and 2.7 GeV (0.028). All 32 pointwise coverage intervals include 95% at the tested injections with the measured trigger uncertainty profiled. C4+L4 is the working model with these qualifications.')])
page('High SR | recovered zero-injection sample','Section 6 validation detail','high',[
('Bias','All 20,000 quartic toy receipts are present for each injection. At zero, 19,997 toys have two valid fits and 3/100 mass-generator central estimates exceed |0.2|. With signal injected, 19,661 toys pass and 97/100 estimates exceed the target. These counts supersede the interim samples on pp. 19-20.'),
('GoF and coverage','GoF uses 15,853 saved toys, with estimates at 33/50 masses and different toy counts per mass. Four local p-values are below 0.05 (p. 33). Coverage includes all 50 masses and both generators; one of the 100 pointwise intervals is wholly below 95%. Quartic injected-signal bias remains widespread.')])
page('9 GeV | interpreting the limit comparisons','Section 8 diagnostic conclusions',None,[
('What the three fits change','C0 matches the value at the physical junction. Adding C1 with a flat continuum in the Upsilon CR also forces zero continuum slope at the SR edge. Floating the CR slope lets data select a slope common to the two sides. The trigger uncertainty is profiled in all three models.'),
('Observed limit relative to its own expected sensitivity (see p. 34)','Observed / expected median: 113.95 / 1.40 (C0), 14.50 / 1.20 (C0+C1, fixed slope), 8.66 / 2.02 (C0+C1, free slope). Floating the slope lowers the observed limit and weakens expected sensitivity. Even then, 8.66 exceeds the upper 95% expected band edge of 6.11.'),
('What the component fits show','For quartic C0+C1 with fixed slope, moving from the background-only to the signal-plus-background fit adds 6,156 fitted SR signal events; the continuum falls by 4,036 and the Upsilon components by 838. This demonstrates compensation at the best fits. It does not quantify compensation along the full upper-limit profile.'),
('Additional diagnostics','The cubic free-slope signal-plus-background fit gives r = -2.94; the quartic masked-sideband fit gives r = 2.42. Both report status 0 and covariance quality 3. Different fit ranges and models make these diagnostic best fits unsuitable as substitutes for a limit or a validated prediction in the masked bins.'),
('Crossing and interpretation','Two initialisations agree at the tested fixed strengths. At 9 GeV the C0 checks at r = 20, 30 and 40 do not bracket the 0.05 crossing. HybridNew pilots at 8.6 and 9 GeV check sampling; a validated full crossing is still outstanding. The cause of the extreme C0 result remains unresolved.')])
c.save();w=PdfWriter(clone_from=reader);w.append(PdfReader(buf))
for i,t in enumerate(titles):w.add_outline_item(t,35+i)
tmp=p.with_suffix('.new.pdf');w.write(tmp);assert len(PdfReader(tmp).pages)==40;tmp.replace(p)
