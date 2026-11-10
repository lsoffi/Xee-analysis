"""Append/rebuild the middle-region GoF section after the 11 existing pages."""
from pathlib import Path
from io import BytesIO
import csv
from pypdf import PdfReader,PdfWriter
from reportlab.pdfgen.canvas import Canvas
from reportlab.lib.utils import simpleSplit,ImageReader
from reportlab.lib.colors import HexColor
base=Path(__file__).resolve().parent
note=base/'Xee_analysis_note.pdf';rows=list(csv.DictReader((base/'data/gof_middle_sr.csv').open()))
assert len(rows)==16 and all(int(r['n'])==500 for r in rows)
buf=BytesIO();c=Canvas(buf,pagesize=(842,650))
def title(s):
 c.setFillColor(HexColor('#17374f'));c.setFont('Helvetica-Bold',21);c.drawString(36,603,s)
def para(y,s):
 c.setFillColor(HexColor('#223744'));c.setFont('Helvetica',11.5)
 for line in simpleSplit(s,'Helvetica',11.5,765):c.drawString(38,y,line);y-=16
 return y-10
def footer(n):
 c.setFillColor(HexColor('#223744'));c.setFont('Helvetica',9);c.drawString(38,20,'Updated 10 November 2026 | Xee analysis');c.drawRightString(804,20,str(n));c.showPage()
title('5.4   Goodness of fit in the phi-to-J/psi region')
y=para(568,'Scope: 16 mass hypotheses, 1.2-2.7 GeV in 0.1 GeV steps, using the adopted C4+L4 envelope and four eras. The signal windows and selected control regions follow the sliding-window table. The fit includes SR and CR data, signal in both, the measured nuisance constraints (including trigger efficiency), and the existing total-density C0 relations at actual common boundaries.')
y=para(y,'The signal strength r is free in the fit. A shared discrete index chooses C4 or L4, and the continuous nuisance parameters are profiled. This study uses the existing workspace parameter ranges; it does not impose an additional r range on the GoF command.')
y=para(y,'Why evaluate the SR separately? A global SR+CR goodness-of-fit statistic also tests the CR shapes, including the charmonium description. To examine the SR while retaining the constraints provided by the CR data, all channel masks are set to zero during minimisation, and only the CR masks are set to one when Combine evaluates the statistic. The same mask protocol is used for the observed data and for every toy.')
y=para(y,'Implementation: native Combine GoodnessOfFit with --algo saturated and --cminDefaultMinimizerStrategy 0, using --setParametersForFit for all masks = 0 and --setParametersForEval for CR masks = 1, SR masks = 0. These settings apply to both nominal and saturated-model stages. This is the native masked saturated-likelihood statistic, not a manually calculated Pearson chi-square; shared auxiliary-constraint handling remains that of Combine.')
y=para(y,'Calibration: 500 frequentist toys per mass, generated and analysed with -t 1 --toysFreq -s SEED. Each toy runs as its own Condor job and immediately saves a result and exit receipt. The seeds are unique: 91000000 + 1000*i + j, with i = 12,...,27 and j = 0,...,499. All 8000 results are present and finite; the observed statistic is calculated once per mass without the toy-generation options.')
y=para(y,'The p-value is k/500, where k counts toys with a statistic greater than or equal to the observed value. Error bars are exact two-sided Clopper-Pearson 95% binomial intervals, describing finite-toy precision. The 0.05 reference is a local diagnostic threshold. Overlapping windows make the 16 tests correlated; no global multiple-testing calibration is claimed.')
y=para(y,'Reproduction: scan_261106/gof_sr_only builds the masked workspaces and evaluates observed statistics; gof_sr_single_toy contains the per-toy runner and full seed manifest. The standalone report script rebuilds the plot and CSV. This GoF tests SR agreement under a simultaneous fit; it does not establish CR adequacy, coverage, or absence of signal-estimator bias. Output completeness alone is not an independent global-minimum audit.')
assert y>35,y
footer(12)
title('Middle-region GoF results and operational status')
im=ImageReader(str(base/'figures/gof_middle_sr.png'));w,h=im.getSize();sc=min(768/w,315/h);c.drawImage(im,(842-w*sc)/2,265,w*sc,h*sc)
y=para(246,'Thirteen of the sixteen masses have p >= 0.05. The three values below the local threshold are: 2.2 GeV, 12/500 = 0.024 (95% interval 0.0125-0.0415); 2.6 GeV, 17/500 = 0.034 (0.0199-0.0539); and 2.7 GeV, 14/500 = 0.028 (0.0154-0.0465). The complete numerical table is stored in docs/data/gof_middle_sr.csv.')
y=para(y,'The calculation is complete, but the outcome is not a uniform GoF pass across the middle region. C4+L4 remains the adopted working configuration; the local tensions at 2.2, 2.6 and 2.7 GeV are retained as qualifications. The mask procedure does not remove the CR data from parameter estimation, and a favourable SR-only p-value must not be interpreted as validation of the CR model.')
y=para(y,'Operational handover: no further middle-region production is queued by this update. The active model-selection and diagnostic programme moves to the high region, 4.1-9.0 GeV: complete FitDiagnostics, limits and C4/L4 signal-injection studies, and investigate the anomalous observed limits near the Upsilon CR. The recorded middle-region GoF qualifications remain part of the analysis record rather than being reclassified as resolved.')
assert y>35,y
footer(13);c.save();reader=PdfReader(note);assert len(reader.pages)>=11
writer=PdfWriter()
for page in reader.pages[:11]:writer.add_page(page)
writer.append(PdfReader(buf));writer.add_metadata({'/Title':'Xee dielectron search | Analysis note'})
tmp=note.with_suffix('.updated.pdf');writer.write(tmp);tmp.replace(note)
