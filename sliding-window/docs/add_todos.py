"""Append or replace the planned-validation TODO page; run after other note builders."""
from pathlib import Path
from io import BytesIO
from pypdf import PdfReader, PdfWriter
from reportlab.pdfgen.canvas import Canvas
from reportlab.lib.utils import simpleSplit
base=Path(__file__).resolve().parent
path=base/'Xee_analysis_note.pdf'
r=PdfReader(path)
pages=[p for p in r.pages if '7   TODO: planned validation checks' not in (p.extract_text() or '')]
b=BytesIO();c=Canvas(b,pagesize=(842,650))
c.setFillColorRGB(.09,.22,.31);c.setFont('Helvetica-Bold',21)
c.drawString(36,608,'7   TODO: planned validation checks')
y=572
items=[
('Status: planned, not launched (10 November 2026).','The checks below are deferred work. No new toy campaign is authorised by this documentation update.'),
('1. Validate the C0-only continuum junction at a common SR/CR boundary.','C0 enforces equal physical event density at the boundary; C1 would also enforce equal slope. For a common event selection and no physical threshold at that boundary, a kink is a modelling approximation rather than an expected physical feature. Fixing the CR continuum shape makes C1 a potentially strong transfer of its local slope into the SR.'),
('2. Generate independent smooth-background alternatives.','Use plausible, positive continuum models that are smooth across the actual common boundary, including alternatives outside the fitted C0 family. Preserve the shared SM-resonance and signal construction, the simultaneous SR+CR geometry and all measured nuisance uncertainties. Prioritise mass hypotheses close to the Upsilon CR boundary.'),
('3. Compare C0-only and C1-constrained fits.','Use paired toy samples with zero and median-expected signal injections. Check convergence, signal bias, interval or CLs-limit coverage and sensitivity to the boundary placement. Retain invalid toys and numerical failures in the reporting; use one Condor job per toy with reproducible seeds.'),
('4. State the validation limits and decision before adoption.','Closure when both generation and fitting use the same C0 construction is not sufficient to establish robustness to a smooth true background. Define acceptance criteria before inspecting this comparison. Keep C0 as the current candidate pending validation; do not select the junction prescription to move observed limits into expected bands. The measured trigger uncertainty remains unchanged and profiled in nominal inference.')]
for title,text in items:
 c.setFont('Helvetica-Bold',11)
 for line in simpleSplit(title,'Helvetica-Bold',11,765):c.drawString(38,y,line);y-=15
 c.setFont('Helvetica',10.5)
 for line in simpleSplit(text,'Helvetica',10.5,765):c.drawString(38,y,line);y-=14
 y-=15
assert y>40,y
c.setFont('Helvetica',9);c.drawString(38,20,'Planned checks | No new campaign launched');c.drawRightString(804,20,str(len(pages)+1));c.save()
w=PdfWriter()
for p in pages:w.add_page(p)
w.add_page(PdfReader(b).pages[0])
for title,i in [('Common methodology',0),('4 Low SR',3),('5 Middle SR',5),('6 High SR - preliminary status',13),('7 TODO: planned validation checks',len(pages))]:w.add_outline_item(title,i)
w.add_metadata({'/Title':'Xee analysis: low, middle, high SR and planned validation'})
w.write(path)
print('Pages:',len(w.pages),'remaining y:',y)
