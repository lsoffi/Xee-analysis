"""Build the window table with region-wide merged recommendation cells."""
import csv
import html
from pathlib import Path
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer, PageBreak
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib import colors
from reportlab.lib.pagesizes import A3, landscape

base = Path(__file__).resolve().parent
source = base / 'tabella_SR_CR_5MeV.csv'
with source.open() as f:
    reader = csv.DictReader(f)
    fields = reader.fieldnames
    rows = list(reader)
column = 'Envelope raccomandato in SR'
for row in rows:
    mass = float(row[fields[0]])
    if 0.4 <= mass <= 0.6:
        row[column] = 'C2+L3 o C3+L3'
    elif 1.2 <= mass <= 2.7:
        row[column] = 'C4+L4 (*)'
with source.open('w') as f:
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)
notes = ('(*) C4+L4 is adopted for the simultaneous SR+CR fits between phi and J/psi (1.2-2.7 GeV). '
         'With median expected signal injected and only CMS_eff_trigger_e fixed to its injected value of zero, '
         '6400/6400 toys pass the existing two-branch validity criteria. 30/32 central bias estimates lie within +/-0.2; '
         'C4 at 1.9 GeV gives +0.203 and L4 at 2.6 GeV gives +0.214. Pointwise 95% bootstrap intervals retain finite-toy uncertainty. '
         'This is a diagnostic with known trigger efficiency. The nominal fits retain the measured trigger uncertainty and profile it. '
         'Coverage at 1.4 GeV with the trigger profiled is still being checked, so adoption of this envelope does not establish coverage.')
styles = getSampleStyleSheet()
styles['Normal'].fontSize = 9
styles['Normal'].leading = 12
styles['BodyText'].fontSize = 8
styles['BodyText'].leading = 10

def spans(chunk):
    result = []
    for value in ('C2+L3 o C3+L3', 'C4+L4 (*)'):
        indices = [i for i,r in enumerate(chunk) if r[column] == value]
        if indices:
            result.append((indices[0], indices[-1]))
    return result

def para(value):
    return Paragraph(html.escape(value).replace('∪',' U ').replace('—','-').replace('–','-'), styles['BodyText'])

story = []
for start in range(0, len(rows), 25):
    chunk = rows[start:start+25]
    story.extend([Paragraph('Sliding-window SR and CR configuration', styles['Heading1']),
                  Paragraph('Updated 9 November 2026. Recommendations are shared across each indicated mass range.',styles['Normal']),Spacer(1,12)])
    data = [[para(x) for x in fields]] + [[para(r[k]) for k in fields] for r in chunk]
    commands = [('BACKGROUND',(0,0),(-1,0),colors.HexColor('#dceaf0')),('VALIGN',(0,0),(-1,-1),'MIDDLE'),
                ('GRID',(0,0),(-1,-1),.3,colors.grey),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5)]
    for first,last in spans(chunk):
        for i in range(first+1,last+1):data[i+1][-1]=''
        commands.append(('SPAN',(-1,first+1),(-1,last+1)))
    table = Table(data,colWidths=[60,125,125,150,42,135,185,40,140],repeatRows=1)
    table.setStyle(TableStyle(commands))
    story.extend([table,Spacer(1,14),Paragraph(notes,styles['Normal'])])
    if start+25<len(rows):story.append(PageBreak())
SimpleDocTemplate(str(base/'tabella_SR_CR_5MeV.pdf'),pagesize=landscape(A3),leftMargin=35,rightMargin=35,topMargin=28,bottomMargin=28).build(story)
merged = spans(rows)
skip = {i for first,last in merged for i in range(first+1,last+1)}
starts = {first:last-first+1 for first,last in merged}
body=[]
for i,row in enumerate(rows):
    cells=[]
    for k in fields:
        if k==column and i in skip:continue
        attr=f' rowspan="{starts[i]}"' if k==column and i in starts else ''
        cells.append(f'<td{attr}>{html.escape(row[k])}</td>')
    body.append('<tr>'+''.join(cells)+'</tr>')
(base/'tabella_SR_CR.html').write_text('<!doctype html><html lang="it"><meta charset="utf-8"><title>Sliding-window table</title><style>body{font:15px Arial;margin:24px}table{border-collapse:collapse}td,th{border:1px solid #aaa;padding:8px}th{background:#dceaf0;position:sticky;top:0}td[rowspan]{background:#eef6f8;font-weight:bold}p{max-width:1100px;line-height:1.5}</style><h1>Sliding-window SR and CR configuration</h1><p>'+html.escape(notes)+'</p><p><a href="tabella_SR_CR_5MeV.pdf">PDF</a> | <a href="tabella_SR_CR_5MeV.csv">CSV</a></p><table><thead><tr>'+''.join('<th>'+html.escape(k)+'</th>' for k in fields)+'</tr></thead><tbody>'+''.join(body)+'</tbody></table></html>')
assert len([r for r in rows if r[column]=='C4+L4 (*)'])==16
assert len([r for r in rows if r[column]=='C2+L3 o C3+L3'])==3
print(len(rows),'rows; merged ranges',merged)
