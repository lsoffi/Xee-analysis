"""Append the completed coverage study to the first nine analysis-note pages."""
from pathlib import Path
from io import BytesIO
from pypdf import PdfReader, PdfWriter
from reportlab.pdfgen.canvas import Canvas
from reportlab.lib.utils import simpleSplit, ImageReader
from reportlab.lib.colors import HexColor

base=Path(__file__).resolve().parent
note=base/'Xee_analysis_note.pdf'
buffer=BytesIO()
c=Canvas(buffer,pagesize=(842,650))
def title(text):
    c.setFillColor(HexColor('#17374f'));c.setFont('Helvetica-Bold',22);c.drawString(36,600,text)
def para(y,text):
    c.setFillColor(HexColor('#223744'));c.setFont('Helvetica',12)
    for line in simpleSplit(text,'Helvetica',12,762):
        c.drawString(38,y,line);y-=17
    return y-12
def footer(number):
    c.setFont('Helvetica',10);c.drawString(38,20,'Updated 9 November 2026');c.drawRightString(804,20,str(number));c.showPage()
title('5.3   Coverage with the trigger uncertainty profiled')
y=para(563,'The injected-signal bias study measures the average displacement of the fitted signal strength. Coverage tests a different property of the inference procedure: how often it retains the true injected signal among the values allowed by a nominal 95% CLs exclusion test. A biased point estimate can coexist with adequate coverage.')
y=para(y,'At each of the 16 masses from 1.2 to 2.7 GeV, we inject a known signal strength equal to the previously calculated median expected limit at that mass. This probes an intensity near the analysis sensitivity. We generate 200 toys with C4 and 200 with L4, for 6400 toys in total. The existing 400 toys at 1.4 GeV are reused.')
y=para(y,'The toys use the saved generator snapshots and the GenerateOnly frequentist-toy configuration. Each toy is analysed with the simultaneous SR plus selected CR model and the shared C4+L4 envelope index profiled across the four eras. The trigger-efficiency nuisance is profiled with its measured constraint. Both branches have the same parameter penalty. Common-boundary C0 relations follow the selected branch, including the CR continuum normalisation at 1.2 and 1.3 GeV.')
y=para(y,'Combine AsymptoticLimits --singlePoint evaluates CLs at the injected r, using minimizer strategy 0 and MINIMIZER_freezeDisassociatedParams, with rMax = 150. A toy excludes its true injection if CLs < 0.05. The estimated coverage is 1 minus the number of such exclusions divided by 200, calculated separately for every mass and generator.')
y=para(y,'For example, the C4 sample at 1.4 GeV excludes the true injection in 8 of 200 toys. Its coverage estimate is therefore 192/200 = 96%. The bars on the next page are exact two-sided Clopper-Pearson 95% binomial intervals. They quantify finite-toy uncertainty on the coverage estimate. These intervals replace the approximate Wilson intervals used in the initial chat update.')
y=para(y,'All 6400 unique seed results are present, with finite CLs in [0,1]. The log scan found no exception, failed-fit or NaN messages. These output and log checks do not prove that every conditional minimisation found its global minimum. The two boundary workspaces passed branch-by-branch expected-count and bin-shape comparisons after including the CR normalisation switch.')
assert y>40
footer(10)
title('Coverage results and implications for the envelope')
im=ImageReader(str(base/'figures/coverage_all_env.png'));w,h=im.getSize();scale=min(768/w,315/h);c.drawImage(im,(842-w*scale)/2,245,w*scale,h*scale)
y=para(225,'Coverage estimates range from 92.5% to 98.0%. Every one of the 32 pointwise 95% intervals contains 95%. The lowest estimate is C4 at 1.7 GeV: 15/200 exclusions, or 92.5% coverage, with an exact interval of 87.9-95.7%. Thus no individual configuration shows statistically resolved undercoverage at this precision.')
y=para(y,'This supports using C4+L4 with the measured trigger uncertainty profiled for the tested region and injections. The positive bias of the fitted-r mean has not translated into resolved undercoverage in this study. An extra systematic uncertainty is not motivated by that bias alone.')
y=para(y,'The intervals are pointwise, and 200 toys cannot establish exact 95% coverage or rule out modest deficits. The conclusion applies to the C4/L4 generators, saved nuisance truth settings and one injected strength per mass. It does not validate all possible background shapes, signal strengths or simultaneous coverage across the mass scan.')
assert y>30
footer(11)
c.save();r=PdfReader(note);writer=PdfWriter()
for page in r.pages[:9]:writer.add_page(page)
writer.append(PdfReader(buffer));writer.add_metadata({'/Title':'Xee dielectron search | Analysis note'})
tmp=note.with_suffix('.updated.pdf');writer.write(tmp);tmp.replace(note)
