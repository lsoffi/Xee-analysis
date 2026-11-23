from pathlib import Path
from io import BytesIO
import shutil
from pypdf import PdfReader, PdfWriter
from reportlab.pdfgen.canvas import Canvas
from reportlab.lib.utils import simpleSplit
import pypdfium2

here = Path(__file__).resolve().parent
reader = PdfReader(BytesIO()) if False else None
pages = [
('Current simultaneous-fit configuration', [
('Scope and geometry', 'Configuration authorised on 23 November 2026 for 69 hypotheses: 0.4-0.6, 1.2-2.7 and 4.1-9.0 GeV, in 0.1 GeV steps. The existing sliding-window table defines the common SR for all four eras. It encloses the era-specific shortest intervals containing 95% of the full signal DCB, with outward bin alignment. CR bins are removed from the SR to avoid double counting.'),
('Control-region selection', 'Use the selected CRs and their bin-aligned ranges from the table. The minimal set of one or two CRs must account for at least 95% of the predicted SM-resonance tails in the actual SR. This selection uses resonance contributions, not the non-resonant continuum. It is a geometric selection rule, not a statistical confidence level.'),
('Signal regions', 'The low SR uses an envelope of Chebyshev degree 3 and LogPol degree 3 (C3+L3). The middle and high SRs use C4+L4. The alternative background branch is profiled with a shared discrete index across channels. The paired branches have equal degree; the implementation uses zero additional discrete correction.'),
('Control-region continua', 'Omega/phi: exponential, with floating slope and normalization. Charmonium: decreasing exponential, with floating negative slope and normalization. Upsilon: LogPol2, exp(c1*t + c2*t*t), with both coefficients and normalization floating; t maps the CR range onto [-1,1]. The continuum is normalized over its region. Boundary constraints can make regional yields derived quantities rather than independent parameters.'),
('Shared signal and SM resonances', 'Fit SR and selected CR channels simultaneously in 2022, 2022EE, 2023 and 2023BPix, with one common signal strength r. Include signal in the SR and every selected CR. Use shared DCB shapes and regional integrals to transfer component yields. SM-resonance means and normalizations float; adopted widths and tail parameters stay fixed. Retain the existing signal-rate and regional-acceptance convention.')]),
('Boundary conditions and limit calculation', [
('C0 and C1 in every SR', 'Impose equality of the total background density and its first mass derivative only at a genuinely shared SR-CR boundary. This applies to low, middle and high SRs. If only one of two selected CRs touches the SR, constrain only that junction. Impose no matching across an excluded interval or gap.'),
('What is matched', 'The constraints act on physical event densities, not separately unit-normalized PDFs. Shared resonance contributions cancel in the boundary equations. Tails of resonances associated with another selected CR are explicitly included when solving for the continuum value and derivative. The shared signal is also continuous and cancels in the matching equations. Consequently, a simple equality of the two continuum PDFs alone is not the general condition.'),
('Nuisance parameters', 'Retain the original ten Gaussian auxiliary constraints and signal systematic uncertainties. In particular, the measured trigger-efficiency uncertainty is profiled in the nominal result. Freezing the trigger in earlier toy studies was a diagnostic and is not part of this configuration.'),
('Combine and execution', 'Run native AsymptoticLimits with --cminDefaultMinimizerStrategy 0 and --X-rtd MINIMIZER_freezeDisassociatedParams, adding --cminRunAllDiscreteCombinations for the envelope. Expected and observed limits run as separate jobs on the same workspace. Combine the five expected quantiles with the observed result to form the limit plot.'),
('Implementation checks and interpretation', 'Before limits, check nuisance inclusion, frozen SM parameters, event totals, finite nonnegative nominal bin predictions, C0/C1 under parameter perturbations, and per-bin equivalence of the standalone branches, envelope and optimized graph. Validate result quantiles and retain numerical warnings. These structural checks do not establish statistical coverage or guarantee a global minimum.'),
('Relation to earlier studies', 'This is the current production configuration, not a claim that the new scan has finished. Earlier bias, GoF and coverage results in these notes refer to their documented models. They do not automatically validate the combination of floating CR continua, Upsilon LogPol2 and C0+C1 in every SR. The measured observed limit is reported with its own expected bands; agreement between them is not a model-selection criterion.')])]
buf = BytesIO()
c = Canvas(buf, pagesize=(842,650))
for i, (title, blocks) in enumerate(pages):
    c.setFillColorRGB(.1,.35,.42); c.setFont('Helvetica-Bold',11)
    c.drawString(40,615,'CURRENT FIT / 23 NOVEMBER 2026')
    c.setFillColorRGB(.08,.12,.18); c.setFont('Helvetica-Bold',22)
    c.drawString(40,579,title)
    y=544
    for heading, paragraph in blocks:
        c.setFont('Helvetica-Bold',12); c.drawString(42,y,heading); y-=18
        c.setFont('Helvetica',11)
        for line in simpleSplit(paragraph,'Helvetica',11,752):
            c.drawString(42,y,line); y-=15
        y-=12
    assert y>40, y
    c.setFont('Helvetica',9); c.drawString(40,23,'Xee analysis note / current configuration')
    c.drawRightString(800,23,str(40+i+1)); c.showPage()
c.save()
w=PdfWriter(); w.append(PdfReader(buf)); w.write(here/'current_configuration.pdf')
