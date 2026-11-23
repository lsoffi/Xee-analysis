from pathlib import Path
from io import BytesIO
import argparse
from pypdf import PdfReader, PdfWriter
from reportlab.pdfgen.canvas import Canvas
from reportlab.lib.utils import simpleSplit
p=argparse.ArgumentParser();p.add_argument('--base',required=True);p.add_argument('--output',required=True);a=p.parse_args()
here=Path(__file__).resolve().parent
r=PdfReader(a.base);assert len(r.pages)==40
w=PdfWriter(clone_from=r);w.append(PdfReader(here/'current_configuration.pdf'))
buf=BytesIO();c=Canvas(buf,pagesize=(842,650));page=43

def header(title):
 c.setFillColorRGB(.1,.35,.42);c.setFont('Helvetica-Bold',11);c.drawString(40,615,'SECTION 8 / 9 GEV UPDATE / 23 NOVEMBER 2026')
 c.setFillColorRGB(.08,.12,.18);c.setFont('Helvetica-Bold',20);c.drawString(40,579,title)
def paras(blocks,y):
 for s in blocks:
  c.setFont('Helvetica',11)
  for line in simpleSplit(s,'Helvetica',11,755):c.drawString(42,y,line);y-=15
  y-=10
 assert y>35,y
def end():
 global page
 c.setFont('Helvetica',9);c.drawString(40,23,'Xee analysis note / supersedes the three-configuration comparison on pp. 34 and 40');c.drawRightString(800,23,str(page));page+=1;c.showPage()
header('Completed limits with floating LogPol2 and LogPol3 CR continua')
c.drawImage(str(here/'figures/limits_m9_configurations.png'),45,260,750,300,preserveAspectRatio=True,anchor='c')
paras([
'The SR envelope remains C4+L4, with C0+C1 at its common boundary with the Upsilon CR. The two added configurations float the LogPol2 or LogPol3 continuum coefficients in that CR. The measured trigger systematic remains profiled. These are changes to the background model, not changes to the trigger measurement.',
'LogPol2: expected median 2.14453, observed 2.13458; expected 68% band [1.03071, 4.31577], 95% band [0.61153, 6.50398]. LogPol3: median 3.07031, observed 5.25228; 68% band [1.33425, 6.17887], 95% band [0.76158, 9.31170]. Both observed values lie within their own expected 68% bands.',
'This demonstrates strong model dependence at 9 GeV. It does not prove the earlier excess was numerical or that either new model has correct coverage. The original C0 crossing remains unverified. Both new expected calculations report approximate covariance warnings; retain these qualifications alongside the completed quantiles.'
],240);end()
header('FitDiagnostics, crossing checks and implementation recovery')
paras([
'FitDiagnostics uses the same simultaneous SR+CR likelihood with r in [-30,30], robustFit=1, strategy 0, saved shapes and normalizations. LogPol2 B-only and S+B fits have status 0 and covariance quality 3. The S+B estimate is r = -0.07813 with asymmetric errors -0.65090/+0.54393. LogPol3 has status 0 but covariance quality 2 in both fits; r = 0.42321 with errors -0.51447/+1.54517. Covariance quality 2 is a retained limitation, not a fully accurate covariance.',
'The following two pages show saved S+B total prefit/postfit shapes, individual postfit components and Pearson residuals for all four eras and both regions. B-only figures are included in the accompanying figure directory. Pearson chi2/Nbin on these figures is not chi2/ndf and is not a calibrated GoF p-value.',
'All four final expected/observed jobs exited zero on 19 November. The saved results contain five ordered expected quantiles and one observed value per model. Observed CLs values bracket 0.05: LogPol2, r=2.126664 gives 0.05027 and r=2.134579 gives 0.04994; LogPol3, r=5.252283 gives 0.05010 and r=5.263534 gives 0.04990. This checks the local crossing, not the global minimum at every scan point.',
'The long runtime was traced to repeated RooFit dependency-dirty propagation during numerical gradients. A quadrature integral represented by 282 intermediate bin nodes created redundant graph paths. An equivalent compiled quadrature node, cancellation of common normalization factors within normalized shapes, flattened identical positivity checks, and fresh-workspace import reduced this overhead. A larger stack limit addressed separately observed stack exhaustion.',
'These changes preserve the likelihood rather than alter the statistical model. Saved audits compare all bin expectations across branches, channels and S+B/B-only models. The compiled comparison agrees below 1e-14; positivity and the recorded coefficient-variation comparisons are exactly equal. Representative update time decreased from 6.46 s to 0.033 s, with subsequent optimized updates around 0.003-0.007 s. These microbenchmarks do not guarantee a full-limit runtime.',
'Documentation audit: the previous 40-page repository note already contains low/middle/high bias, GoF and coverage summaries, trigger-freeze diagnostics, C0/C1, initialization and HybridNew discussions. This update adds the previously missing final CR LogPol2/3 results, fit figures, numerical recovery and the current full-spectrum configuration. Historical partial studies remain labelled as such; no missing result is inferred from a submission or elapsed time.'
],540);end()
for d in [2,3]:
 header(f'9 GeV shapes and residuals: floating CR LogPol{d}')
 c.drawImage(str(here/f'figures/fit_L{d}_s.png'),160,48,520,510,preserveAspectRatio=True,anchor='c');end()
c.save();w.append(PdfReader(buf))
w.add_outline_item('Current simultaneous-fit configuration',40)
w.add_outline_item('9 GeV: final LogPol2/3 comparisons and numerical checks',42)
w.write(a.output)
assert len(PdfReader(a.output).pages)==46
