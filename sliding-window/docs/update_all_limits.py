"""Replace the limits page by heading, preserving all other pages and outlines."""
from pathlib import Path
from io import BytesIO
import json
from pypdf import PdfReader,PdfWriter
from reportlab.pdfgen.canvas import Canvas
from reportlab.lib.utils import simpleSplit,ImageReader
base=Path(__file__).resolve().parent
d=base/'data/limits_all_sr';s=json.loads((d/'limits_all_sr.json').read_text())
assert s['points']==69 and not s['missing']
p=base/'Xee_analysis_note.pdf';r=PdfReader(p)
i=next(i for i,p in enumerate(r.pages) if '6.1   Preliminary high-region limits' in p.extract_text() or '6.1   Sliding-window limits: all regions' in p.extract_text())
b=BytesIO();c=Canvas(b,pagesize=(842,650))
c.setFillColorRGB(.09,.22,.31);c.setFont('Helvetica-Bold',21);c.drawString(36,608,'6.1   Sliding-window limits: all regions')
im=ImageReader(str(d/'limits_all_sr.png'));iw,ih=im.getSize();scale=min(768/iw,405/ih);c.drawImage(im,(842-iw*scale)/2,190,iw*scale,ih*scale)
y=170
for text in [f"Updated {s['updated']}: all 69 hypotheses are available (3 low, 16 middle, 50 high). The newly completed 8.5 GeV point has observed r95 = 0.287766 and median expected r95 = 0.75. Expected bands are the native 68% and 95% quantiles.","All six quantiles are finite and the expected bands are ordered at every point. This is an output-integrity check, not validation of profile minima or CLs crossings. The high-region observed limits at 8.6 and 9.0 GeV remain under numerical investigation and are shown unchanged.","Low SR uses C3+L3; middle and high use C4+L4. The high envelope remains provisional. The low-SR GoF tension at 0.6 GeV remains unresolved. This limits update supersedes the earlier limit counts on the preceding historical status page; the following bias figures retain their own dated snapshot."]:
 c.setFont('Helvetica',10.5)
 for line in simpleSplit(text,'Helvetica',10.5,765):c.drawString(38,y,line);y-=14
 y-=8
assert y>32
c.setFont('Helvetica',9);c.drawString(38,20,'PRELIMINARY | Limits updated: '+s['updated']);c.drawRightString(804,20,str(i+1));c.save()
w=PdfWriter();newreader=PdfReader(b)
for j,page in enumerate(r.pages):w.add_page(newreader.pages[0] if j==i else page)
for title,j in [('Common methodology',0),('4 Low SR',3),('5 Middle SR',8),('6 High SR',16),('6.1 All-region limits',17),('6.2 Bias r=0',18),('6.3 Injected bias',19),('7 TODO',20)]:w.add_outline_item(title,j)
w.write(p)
print('Updated page',i+1,'of',len(w.pages))
