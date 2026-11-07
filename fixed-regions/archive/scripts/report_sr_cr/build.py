from pathlib import Path
import io,json,csv,math,hashlib,copy,re
from xml.sax.saxutils import escape
from datetime import datetime
from zoneinfo import ZoneInfo
from reportlab.pdfgen import canvas
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import Paragraph,Table,TableStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.utils import ImageReader
from pypdf import PdfReader,PdfWriter,Transformation
P=Path.cwd();TMP=P/'tmp/pdfs/sr_cr_report';OUT=P/'output/pdf';R=P/'runs/20261030_sr6_cr4_scan';G=P/'runs/20261030_gof_profiled_combined'
FONT=Path('/Users/imb/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/libreoffice-headless/libreoffice/LibreOfficeDev.app/Contents/Resources/fonts/truetype')
for name,file in [('DV','DejaVuSans.ttf'),('DV-B','DejaVuSans-Bold.ttf'),('DV-I','DejaVuSans-Oblique.ttf'),('Mono','DejaVuSansMono.ttf')]:pdfmetrics.registerFont(TTFont(name,str(FONT/file)))
pdfmetrics.registerFontFamily('DV',normal='DV',bold='DV-B',italic='DV-I',boldItalic='DV-B')
W,H=841.8898,595.2756;M=38;CW=W-2*M;NAVY=colors.HexColor('#12314b');BLUE=colors.HexColor('#0072B2');TEAL=colors.HexColor('#007c80');ORANGE=colors.HexColor('#b06a00');RED=colors.HexColor('#a12d29');GRAY=colors.HexColor('#526474');PALE=colors.HexColor('#edf4f8')
c=canvas.Canvas(str(TMP/'main.pdf'),pagesize=(W,H));c.setTitle('Analisi Xee - Strategia simultanea SR+CR e risultati');c.setAuthor('Analisi statistica Xee');c.setSubject('Configurazione, motivazioni, fit, bias, goodness of fit, limiti e significanze; 30 ottobre 2026')
A=json.loads((R/'summary.json').read_text());B=json.loads((R/'bias_summary.json').read_text());GS=json.loads((G/'summary.json').read_text());gb={a['mass']:a for a in GS};S=json.loads((TMP/'signal_acceptance.json').read_text());eras=['2022','2022EE','2023','2023BPix'];pages=[];y=0;source_files=set()
def cleantext(text):
 text=str(text).replace('−','-').replace('–','-').replace('—','-')
 text=re.sub(r'(?<=\d)(?=(?:toy|fit|limiti|masse|punti|score|GeV|MeV|campioni|era|combinazioni|ID|volte|superamenti|fallimenti|dataset|seed|parametri|soluzioni|coordinate|shape|sigma|σ))',' ',text)
 text=re.sub(r'\b(grado|su|con|da|in|fra|tra|tutti|tutte|al)(?=\d)', r'\1 ',text)
 for a,b in [('SRenvelope','SR envelope'),('SR/CR','SR/CR'),('rtrue','rtrue'),('fitGoF','fit GoF'),('scoreGoF','score GoF'),('rtoy','r toy'),('Υliberi','Υ liberi'),('Υfisse','Υ fisse'),('ΥMCfisse','Υ MC fisse'),('σfloating','σ floating'),('SRunica','SR unica'),('SR+CR','SR+CR'),('fit4','fit 4'),('CR4','CR4'),('rfree','r libero')]:text=text.replace(a,b)
 return text
def fmt(v,n=3):return f'{v:.{n}g}'.replace('.',',')
def par(text,x,width,top,size=10,leading=None,color=NAVY,bold=False):
 st=ParagraphStyle('p',fontName='DV-B' if bold else 'DV',fontSize=size,leading=leading or size*1.4,textColor=color,spaceAfter=0,allowWidows=0,allowOrphans=0)
 p=Paragraph(cleantext(text),st);ww,hh=p.wrap(width,1000);p.drawOn(c,x,top-hh);return hh

def page(title,section='STRATEGIA E RISULTATI',sub=None):
 global y
 if pages:c.showPage()
 pages.append(dict(title=title,section=section))
 c.setFillColor(NAVY);c.rect(0,H-9,W,9,fill=1,stroke=0)
 c.setFont('DV-B',8);c.setFillColor(TEAL);c.drawString(M,H-30,'Xee  |  '+cleantext(section))
 h=par(title,M,CW,H-44,21,25,NAVY,True);y=H-44-h-12
 if sub:y-=par(sub,M,CW,y,9.3,13,GRAY)+9
 c.setStrokeColor(colors.HexColor('#c9d6df'));c.line(M,31,W-M,31);c.setFont('DV',7.3);c.setFillColor(GRAY);c.drawString(M,18,'30 ottobre 2026  |  SR+CR, quattro ere  |  Risultati e limiti di validazione')
 c.drawRightString(W-M,18,f'{len(pages)}')

def p(t,size=10,space=8,color=NAVY):
 global y
 h=par(t,M,CW,y,size,None,color);y-=h+space
 if y<38:raise RuntimeError(f'Overflow page {len(pages)} {pages[-1]["title"]} y={y}')

def box(t,title=None,color=TEAL):
 global y
 text=(f'<b>{title}</b><br/>' if title else '')+t
 st=ParagraphStyle('box',fontName='DV',fontSize=10,leading=14,textColor=NAVY);pp=Paragraph(cleantext(text),st);ww,hh=pp.wrap(CW-26,1000)
 c.setFillColor(PALE if color!=RED else colors.HexColor('#fff0e9'));c.roundRect(M,y-hh-19,CW,hh+19,5,fill=1,stroke=0);c.setFillColor(color);c.rect(M,y-hh-19,4,hh+19,fill=1,stroke=0);pp.drawOn(c,M+13,y-hh-9);y-=hh+29
 if y<38:raise RuntimeError(f'Box overflow {len(pages)}')

def table(headers,rows,widths=None,size=9,pad=5):
 global y
 widths=widths or [CW/len(headers)]*len(headers)
 if abs(sum(widths)-1)<.01:widths=[CW*a for a in widths]
 def cell(v,head=False):return Paragraph(cleantext(v),ParagraphStyle('cell',fontName='DV-B' if head else 'DV',fontSize=size,leading=size*1.28,textColor=colors.white if head else NAVY))
 data=[[cell(v,True) for v in headers]]+[[cell(v) for v in row] for row in rows];tt=Table(data,colWidths=widths,repeatRows=1,hAlign='LEFT');tt.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),NAVY),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,PALE]),('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),pad),('RIGHTPADDING',(0,0),(-1,-1),pad),('TOPPADDING',(0,0),(-1,-1),pad),('BOTTOMPADDING',(0,0),(-1,-1),pad),('LINEBELOW',(0,0),(-1,0),.5,NAVY)]));ww,hh=tt.wrap(CW,1000)
 if y-hh<38:raise RuntimeError(f'Table overflow page {len(pages)}: {hh} > {y-38}')
 tt.drawOn(c,M,y-hh);y-=hh+10

def fig(path,height=290,caption=None):
 global y
 path=Path(path);source_files.add(path);im=ImageReader(str(path));iw,ih=im.getSize();ww=min(CW,height*iw/ih);hh=ww*ih/iw
 if y-hh<38:raise RuntimeError(f'Figure overflow {len(pages)}')
 c.drawImage(im,M+(CW-ww)/2,y-hh,width=ww,height=hh,mask='auto');y-=hh+6
 if caption:p(caption,8.5,8,GRAY)

def twofig(left,right,height=260,caption=None):
 global y
 for i,path in enumerate([left,right]):
  path=Path(path);source_files.add(path);im=ImageReader(str(path));iw,ih=im.getSize();ww=min((CW-14)/2,height*iw/ih);hh=ww*ih/iw;c.drawImage(im,M+i*(CW+14)/2+((CW-14)/2-ww)/2,y-hh,width=ww,height=hh,mask='auto')
 y-=height+8
 if caption:p(caption,8.5,8,GRAY)

def source(path):source_files.add(P/path);return P/path
exec((P/'scripts/report_sr_cr/content.py').read_text())
c.save()
# Merge the original vector fit pages, retaining their native size and adding navigation.
main=PdfReader(TMP/'main.pdf');writer=PdfWriter();writer.append(main,import_outline=False);core_count=len(main.pages)
for i,meta in enumerate(pages):writer.add_outline_item(meta['title'],i)
appendix_start=len(writer.pages);parent=writer.add_outline_item('Atlante completo: 31 masse, SR e CR, quattro ere',appendix_start)
fit_index=[]
for a in A:
 m=a['mass'];mass_parent=None
 for region in ['SR','CR']:
  path=R/f'm{m:.1f}/fit_{region}.pdf';source_files.add(path);rr=PdfReader(path);assert len(rr.pages)==1
  original=rr.pages[0];ow=float(original.mediabox.width);oh=float(original.mediabox.height);pagew=ow;pageh=oh+55
  sheet=writer.add_blank_page(pagew,pageh);sheet.merge_transformed_page(original,Transformation().translate(0,32))
  bio=io.BytesIO();overlay=canvas.Canvas(bio,pagesize=(pagew,pageh));overlay.setFillColor(NAVY);overlay.setFont('DV-B',9);overlay.drawString(25,pageh-16,f'ATLANTE FIT  |  mX = {m:.1f} GeV  |  {region}');overlay.setFont('DV',8)
  overlay.drawString(25,18,f'r fitted = {a["r"]:.5g}   |   Z locale = {a["Z"]:.3f}   |   pGoF S+B = {gb[m]["p"]:.2f}   |   FD B/S: status 0, covQual 3')
  overlay.drawRightString(pagew-25,18,str(len(writer.pages)));overlay.save();sheet.merge_page(PdfReader(bio).pages[0]);pi=len(writer.pages)-1
  if region=='SR':mass_parent=writer.add_outline_item(f'Massa {m:.1f} GeV',pi,parent=parent)
  writer.add_outline_item(f'{region}: quattro ere e componenti',pi,parent=mass_parent)
  fit_index.append(dict(mass=m,region=region,page=pi+1,source=str(path)))
writer.add_metadata({'/Title':'Xee - Strategia simultanea SR+CR, scelte e risultati','/Author':'Analisi statistica Xee','/Subject':'Report dettagliato; scan 6-9 GeV, CR Upsilon, envelope Chebyshev6/Bernstein6','/Keywords':'Xee, CMS Combine, SR, CR, DCB, Upsilon, bias, GoF, CLs'})
out=OUT/'Xee_report_strategia_SR_CR_2026-10-30.pdf'
for pg in writer.pages:pg.compress_content_streams()
# Preserve each imported plot's resources independently to avoid resource aliasing.
with out.open('wb') as f:writer.write(f)
manifest=dict(pdf=str(out),core_pages=core_count,fit_pages=len(fit_index),total_pages=len(writer.pages),sections=pages,fit_index=fit_index,sources=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(source_files) if p.exists()])
(TMP/'manifest.json').write_text(json.dumps(manifest,indent=2));(OUT/'Xee_report_strategia_SR_CR_manifest.json').write_text(json.dumps(manifest,indent=2));print(json.dumps({k:manifest[k] for k in ['pdf','core_pages','fit_pages','total_pages']},indent=2))
