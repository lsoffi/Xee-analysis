"""Apply regional headings, navigation and a high-SR work-in-progress section.
Run after the historical note-update scripts; retains the first 13 scientific pages.
"""
from pathlib import Path
from io import BytesIO
from pypdf import PdfReader, PdfWriter
from reportlab.pdfgen.canvas import Canvas
from reportlab.lib.colors import HexColor
from reportlab.lib.utils import simpleSplit
base=Path(__file__).resolve().parent
path=base/'Xee_analysis_note.pdf'
r=PdfReader(path); w=PdfWriter()
headings={3:(76,'4   Low SR: below the phi (0.4-0.6 GeV)'),4:(76,'4.1   Low-SR envelope selection and bias studies'),5:(60,'5   Middle SR: phi to J/psi (1.2-2.7 GeV)')}
for i,p in enumerate(r.pages[:13]):
 if i in headings:
  top,title=headings[i]; b=BytesIO(); c=Canvas(b,pagesize=(842,650))
  c.setFillColorRGB(1,1,1);c.rect(32,650-top,780,35,fill=1,stroke=0)
  c.setFillColor(HexColor('#17374f'));c.setFont('Helvetica-Bold',22);c.drawString(36,650-top+8,title);c.save()
  p.merge_page(PdfReader(b).pages[0])
 w.add_page(p)
b=BytesIO();c=Canvas(b,pagesize=(842,650))
c.setFillColor(HexColor('#17374f'));c.setFont('Helvetica-Bold',22)
c.drawString(36,603,'6   High SR: J/psi to Upsilon (4.1-9.0 GeV)')
y=562
sections=[
('6.1   Scope and candidate model','This section records the programme in progress, not a validated high-region recommendation. The 50 hypotheses follow the common SR windows and selected CRs in the sliding-window table. Hypothesis centres remain outside the CRs. C4+L4 is the candidate envelope; the simultaneous four-era likelihood, signal in SR and CR, nuisance constraints and total-density C0 conditions are those described in Sections 1-3.'),
('6.2   Fit diagnostics and limits','FitDiagnostics checks the fitted components, residuals, fit status and covariance quality. Native AsymptoticLimits provides the expected bands and observed limits. Exhaustive discrete-branch minimisation is used in the current high-region campaign after diagnosing a minimiser crash. Remaining missing or flagged fits must be resolved before presenting the campaign as complete.'),
('6.3   Bias studies and the Upsilon boundary','Generate with C4 and L4 and refit with the C4+L4 envelope at r = 0 and at the mass-specific median expected limit. Each toy has its own reproducible seed and Condor job. Both branches must satisfy the recorded validity criteria before envelope selection. The anomalous observed limits at 8.6 and 9.0 GeV are being investigated using conditional likelihood profiles and single-point CLs checks; neither a trigger explanation nor a physical excess is established by the available limit values alone.'),
('6.4   Pending conclusion','The final high-region section will collect the completed fit, limit and bias figures, numerical exclusions and the boundary diagnosis. No final high-region envelope choice is made here. Low- and middle-region results remain in Sections 4 and 5, including the middle-region GoF qualifications; this structural update does not change any statistical result.')]
for title,text in sections:
 c.setFont('Helvetica-Bold',14);c.drawString(38,y,title);y-=23
 c.setFont('Helvetica',11.5)
 for line in simpleSplit(text,'Helvetica',11.5,765):c.drawString(38,y,line);y-=16
 y-=19
assert y>35,y
c.setFont('Helvetica',9);c.drawString(38,20,'Work in progress | Updated 10 November 2026');c.drawRightString(804,20,'14');c.save()
w.add_page(PdfReader(b).pages[0])
common=w.add_outline_item('Common methodology: windows, CRs and simultaneous likelihood',0)
for title,page in [('Signal windows',0),('Control regions',1),('Simultaneous likelihood',2)]:w.add_outline_item(title,page,parent=common)
low=w.add_outline_item('4 Low SR: below phi',3);w.add_outline_item('Envelope selection and bias',4,parent=low)
mid=w.add_outline_item('5 Middle SR: phi to J/psi',5)
for title,page in [('Trigger at 1.4 GeV',6),('Fixed-trigger bias',8),('Coverage',9),('SR-only goodness of fit',11)]:w.add_outline_item(title,page,parent=mid)
w.add_outline_item('6 High SR: J/psi to Upsilon - work in progress',13)
w.add_metadata({'/Title':'Xee analysis: common methodology, low, middle and high SR'})
tmp=path.with_suffix('.tmp.pdf');w.write(tmp);tmp.replace(path)
