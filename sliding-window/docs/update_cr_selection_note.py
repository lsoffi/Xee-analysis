"""Add the CR-selection validation rationale to section 3 of the analysis note."""
from pathlib import Path
from io import BytesIO
from pypdf import PdfReader, PdfWriter
from reportlab.pdfgen.canvas import Canvas
from reportlab.lib.utils import simpleSplit
from reportlab.lib.colors import HexColor

path = Path(__file__).resolve().parent / 'Xee_analysis_note.pdf'
reader = PdfReader(path)
text = reader.pages[2].extract_text()
starts = ['For each candidate CR,', 'Each era contributes', 'A single signal strength',
          'Within each era,', 'The baseline limit scan', 'Baseline Combine settings:']
paragraphs = []
for i, start in enumerate(starts):
    begin = text.index(start)
    end = text.index(starts[i+1]) if i+1 < len(starts) else len(text)
    paragraph = ' '.join(text[begin:end].split())
    if i == 0:
        paragraph = paragraph.split('Validation of the CR-selection threshold.')[0].strip()
    paragraphs.append(paragraph)
paragraphs.insert(1, 'Validation of the CR-selection threshold. The 95% requirement is an operational choice intended to retain nearly all predicted SM resonance tails while limiting the number of control regions. It has no confidence-level interpretation and has not been statistically optimised. The remaining fraction is relative to the resonant contribution alone and does not directly bound the impact on the fitted signal or limits. Validate this choice by adding the next-ranked CR and comparing the fitted signal strength and limits, especially at masses where the selected CR set changes. This stability test remains to be performed.')
buffer = BytesIO()
c = Canvas(buffer, pagesize=(842,650))
c.setFillColor(HexColor('#15586b'))
c.setFont('Helvetica-Bold',10)
c.drawString(36,621,'Xee   /   DIELECTRON SEARCH                         ANALYSIS NOTE   /   9 NOVEMBER 2026')
c.setFillColor(HexColor('#162d39'))
c.setFont('Helvetica-Bold',23)
c.drawString(36,585,'3   Simultaneous likelihood for each mass')
y=550
for paragraph in paragraphs:
    c.setFont('Helvetica',11)
    for line in simpleSplit(paragraph,'Helvetica',11,766):
        c.drawString(36,y,line)
        y-=14
    y-=9
assert y>35, y
c.setFont('Helvetica',10)
c.drawRightString(806,18,'3')
c.save()
writer=PdfWriter()
replacement=PdfReader(buffer).pages[0]
for i,page in enumerate(reader.pages):
    writer.add_page(replacement if i==2 else page)
writer.add_metadata({'/Title':'Xee dielectron search | Analysis note'})
output=path.with_suffix('.updated.pdf')
writer.write(output)
output.replace(path)
