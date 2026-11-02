"""Render the Italian methods guide as a compact, checked PDF.

Requires reportlab. The source text stays editable in docs/METODO.md.
"""
from pathlib import Path
import argparse
import html
import os
import re
import textwrap
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, PageBreak, Table, TableStyle,
    KeepTogether, Preformatted, Flowable, Image,
)

REPO = Path(__file__).resolve().parents[1]
NAVY = colors.HexColor('#153447')
TEAL = colors.HexColor('#007B83')
INK = colors.HexColor('#213746')
MUTED = colors.HexColor('#526877')
PALE = colors.HexColor('#EFF6F7')
ORANGE = colors.HexColor('#A55915')
AMBER = colors.HexColor('#FFF3E3')
GREY = colors.HexColor('#D8E3E7')


def load_fonts():
    roots = [Path(os.environ.get('XEE_FONT_DIR', '/nonexistent')),
             Path('/usr/share/fonts/truetype/dejavu'),
             Path.home() / '.cache/codex-runtimes/codex-primary-runtime/dependencies/native/libreoffice-headless/libreoffice/LibreOfficeDev.app/Contents/Resources/fonts/truetype']
    try:
        import matplotlib
        roots.append(Path(matplotlib.get_data_path())/'fonts/ttf')
    except ImportError:
        pass
    root = next((p for p in roots if (p/'DejaVuSans.ttf').exists()), None)
    if root is None:
        raise RuntimeError('Set XEE_FONT_DIR to a directory containing DejaVuSans fonts.')
    for name, file in [('Body','DejaVuSans.ttf'), ('BodyBold','DejaVuSans-Bold.ttf'),
                       ('BodyItalic','DejaVuSans-Oblique.ttf'), ('Mono','DejaVuSansMono.ttf')]:
        pdfmetrics.registerFont(TTFont(name,str(root/file)))
    pdfmetrics.registerFontFamily('Body',normal='Body',bold='BodyBold',italic='BodyItalic',boldItalic='BodyBold')


def styles():
    return {
        'body': ParagraphStyle('body',fontName='Body',fontSize=9.2,leading=13.4,textColor=INK,spaceAfter=7),
        'small': ParagraphStyle('small',fontName='Body',fontSize=7.6,leading=10.6,textColor=MUTED,spaceAfter=5),
        'h1': ParagraphStyle('h1',fontName='BodyBold',fontSize=22,leading=27,textColor=NAVY,spaceBefore=4,spaceAfter=15,keepWithNext=True),
        'h2': ParagraphStyle('h2',fontName='BodyBold',fontSize=13.2,leading=17,textColor=TEAL,spaceBefore=12,spaceAfter=7,keepWithNext=True),
        'h3': ParagraphStyle('h3',fontName='BodyBold',fontSize=10.3,leading=14,textColor=NAVY,spaceBefore=7,spaceAfter=5,keepWithNext=True),
        'code': ParagraphStyle('code',fontName='Mono',fontSize=6.9,leading=9.8,textColor=INK,backColor=PALE,borderPadding=7,spaceBefore=4,spaceAfter=9),
        'bullet': ParagraphStyle('bullet',fontName='Body',fontSize=9.1,leading=13.2,textColor=INK,leftIndent=12,firstLineIndent=-9,spaceAfter=5),
        'callout': ParagraphStyle('callout',fontName='Body',fontSize=9,leading=13,textColor=INK,backColor=AMBER,borderColor=colors.HexColor('#E8C69C'),borderWidth=.5,borderPadding=8,spaceBefore=8,spaceAfter=10),
        'cell': ParagraphStyle('cell',fontName='Body',fontSize=7.5,leading=10.5,textColor=INK),
        'cellhead': ParagraphStyle('cellhead',fontName='BodyBold',fontSize=7.5,leading=10.5,textColor=colors.white),
    }


def inline(text):
    text = text.replace('\u2011','-').replace('\u2013','-').replace('\u2014','-')
    text = html.escape(text)
    text = re.sub(r'\[([^\]]+)\]\(([^)]+)\)',r'<a href="\2" color="#007B83">\1</a>',text)
    text = re.sub(r'`([^`]+)`',r'<font name="Mono" size="8">\1</font>',text)
    text = re.sub(r'\*\*([^*]+)\*\*',r'<b>\1</b>',text)
    text = re.sub(r'(?<!\*)\*([^*]+)\*(?!\*)',r'<i>\1</i>',text)
    return text


class Configuration(Flowable):
    """Explicit diagram of the two alternative eight-channel likelihoods."""
    def __init__(self,width):
        super().__init__(); self.width=width; self.height=430

    def draw(self):
        c=self.canv; w=self.width
        def label(x,y,s,size=9,bold=False,color=INK):
            c.setFont('BodyBold' if bold else 'Body',size);c.setFillColor(color);c.drawString(x,y,s)
        def box(x,y,bw,bh,title,lines,bg=PALE):
            c.setFillColor(bg);c.setStrokeColor(GREY);c.roundRect(x,y,bw,bh,7,fill=1,stroke=1)
            label(x+12,y+bh-20,title,11,True,TEAL)
            for i,line in enumerate(lines):label(x+12,y+bh-39-15*i,line,8.6)
        def arrow(x1,y1,x2,y2):
            c.setStrokeColor(TEAL);c.setLineWidth(1.5);c.line(x1,y1,x2,y2)
            c.line(x2,y2,x2-5,y2+3);c.line(x2,y2,x2-5,y2-3)
        label(0,416,'4 ERE: 2022  |  2022EE  |  2023  |  2023BPix',10,True)
        label(0,398,'Stessi bin originali: 13,333 MeV. Un fit simultaneo per ciascuna analisi.',8.4,color=MUTED)
        gap=30; bw=(w-gap)/2
        box(0,281,bw,100,'ANALISI BASSA',[
            'SR: 6 - 9,093333 GeV',
            'Continuo Chebyshev 6',
            '232 bin per era',
            'Masse testate: 6,0 - 10,5 GeV'])
        box(bw+gap,281,bw,100,'CR NELLO STESSO FIT',[
            'CR: 9,093333 - 10,506667 GeV',
            'Continuo Chebyshev 4',
            '106 bin per era',
            'Segnale + Υ(1S) + Υ(2S) presenti'])
        arrow(bw+3,330,bw+gap-3,330)
        label(6,262,'Raccordo C0 a 9,093333 GeV, separatamente in ogni era.',8.6,True)
        box(0,145,bw,100,'CR NELLO STESSO FIT',[
            'CR: 9,093333 - 10,506667 GeV',
            'Continuo Chebyshev 4',
            '106 bin per era',
            'Stessi dati CR dell’analisi bassa'])
        box(bw+gap,145,bw,100,'ANALISI ALTA',[
            'SR: 10,506667 - 12 GeV',
            'Continuo Chebyshev 4',
            '112 bin per era',
            'Masse testate: 9,1 - 12,0 GeV'])
        arrow(bw+3,194,bw+gap-3,194)
        label(6,126,'Raccordo C0 a 10,506667 GeV, separatamente in ogni era.',8.6,True)
        c.setFillColor(NAVY);c.roundRect(0,38,w,72,7,fill=1,stroke=0)
        label(12,91,'IN CIASCUN FIT: 8 CANALI = 4 ERE × (SR + CR)',10,True,colors.white)
        label(12,72,'Unico r del segnale; Υ condivise fra le due regioni della stessa era.',8.7,color=colors.white)
        label(12,55,'Coefficienti dei due continui distinti e liberi; rese collegate dal C0.',8.7,color=colors.white)
        label(0,20,'Le due analisi condividono la CR: risultati correlati, non da moltiplicare.',8.3,True,ORANGE)
        label(0,4,'Nessun envelope nella configurazione C0. Passo delle ipotesi: 0,1 GeV.',8.3,color=MUTED)


def table(lines,sty,width):
    rows=[]
    for line in lines:
        cells=[x.strip() for x in line.strip().strip('|').split('|')]
        if all(re.fullmatch(r'[:\- ]+',x or '-') for x in cells):continue
        rows.append(cells)
    n=max(map(len,rows)); widths=[width/n]*n
    if n==2:widths=[width*.27,width*.73]
    elif n==3:widths=[width*.23,width*.37,width*.40]
    content=[[Paragraph(inline(s),sty['cellhead'] if i==0 else sty['cell']) for s in r] for i,r in enumerate(rows)]
    t=Table(content,colWidths=widths,repeatRows=1,hAlign='LEFT')
    t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),NAVY),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,PALE]),('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),7),('RIGHTPADDING',(0,0),(-1,-1),7),('TOPPADDING',(0,0),(-1,-1),7),('BOTTOMPADDING',(0,0),(-1,-1),7),('LINEBELOW',(0,0),(-1,0),.5,NAVY)]))
    return [t,Spacer(1,10)]


def parse(source,sty,width):
    lines=source.splitlines(); story=[];i=0
    while i<len(lines):
        line=lines[i].strip()
        if not line:i+=1;continue
        if line=='<!-- PAGE -->':story.append(PageBreak());i+=1;continue
        if line=='[SCHEMA_CONFIGURAZIONE]':story.append(Configuration(width));i+=1;continue
        if line.startswith('<!--'):i+=1;continue
        figure = re.fullmatch(r'!\[([^\]]*)\]\(([^)]+)\)',line)
        if figure:
            picture=Image(str(REPO/'docs'/figure.group(2)))
            max_height=235 if figure.group(2).endswith('gof_historical.png') else 285
            scale=min(width/picture.imageWidth,max_height/picture.imageHeight)
            picture.drawWidth=picture.imageWidth*scale
            picture.drawHeight=picture.imageHeight*scale
            picture.hAlign='CENTER'
            story.append(KeepTogether([Spacer(1,5),picture,Paragraph(inline(figure.group(1)),sty['small'])]))
            i+=1;continue
        if line.startswith('```'):
            block=[];i+=1
            while i<len(lines) and not lines[i].strip().startswith('```'):
                block.extend(textwrap.wrap(lines[i],width=99,replace_whitespace=False,drop_whitespace=False,subsequent_indent='  ') or ['']);i+=1
            story.append(Preformatted('\n'.join(block),sty['code']));i+=1;continue
        if line.startswith('|'):
            block=[]
            while i<len(lines) and lines[i].strip().startswith('|'):block.append(lines[i]);i+=1
            story.extend(table(block,sty,width));continue
        if line.startswith('#'):
            level=len(line)-len(line.lstrip('#')); text=line[level:].strip()
            story.append(Paragraph(inline(text),sty['h1' if level==1 else 'h2' if level==2 else 'h3']));i+=1;continue
        if re.match(r'^([-*]|\d+\.)\s',line):
            story.append(Paragraph(inline(re.sub(r'^[-*]\s','• ',line)),sty['bullet']));i+=1;continue
        paragraph=[line];i+=1
        while i<len(lines) and lines[i].strip() and not re.match(r'^(#|\||```|<!--|\[SCHEMA|[-*] |\d+\. )',lines[i].strip()):paragraph.append(lines[i].strip());i+=1
        text=' '.join(paragraph)
        style=sty['callout'] if 'Manca / da rifare:' in text else sty['small'] if text.startswith('Fonte:') or text.startswith('Fonti:') else sty['body']
        story.append(Paragraph(inline(text),style))
    return story


def footer(c,doc):
    w,h=A4;c.saveState();c.setStrokeColor(GREY);c.setLineWidth(.5);c.line(42,39,w-42,39)
    c.setFont('Body',7);c.setFillColor(MUTED)
    c.drawString(42,26,'Xee | Configurazione, risultati e riproduzione | 2 novembre 2026')
    c.drawRightString(w-42,26,str(doc.page))
    if doc.page>1:
        c.setFont('BodyBold',7.2);c.setFillColor(TEAL);c.drawString(42,h-29,'MANUALE OPERATIVO / SR + CR')
    c.restoreState()


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path,default=REPO/'docs/METODO.md');ap.add_argument('--output',type=Path,default=REPO/'docs/Xee_guida_SR_CR.pdf');args=ap.parse_args()
    load_fonts();sty=styles();args.output.parent.mkdir(parents=True,exist_ok=True)
    doc=SimpleDocTemplate(str(args.output),pagesize=A4,rightMargin=42,leftMargin=42,topMargin=45,bottomMargin=52,title='Xee - Configurazione SR+CR e guida alla riproduzione',author='Analisi Xee',subject='Configurazione corrente, studi disponibili e verifiche da completare')
    story=parse(args.source.read_text(),sty,A4[0]-84)
    doc.build(story,onFirstPage=footer,onLaterPages=footer)
    print(args.output)


if __name__=='__main__':main()
