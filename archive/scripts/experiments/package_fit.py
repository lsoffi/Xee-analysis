import json, pathlib, shutil
r=json.load(open('outputs/dcb_fit_results.json'))
pars=r['parameters']; names=['mu','sigma','alphaL','nL','alphaR','nR']
s='''#include <TFile.h>
#include <TMatrixDSym.h>
#include <TNamed.h>
#include <TParameter.h>
#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooCrystalBall.h>
#include <RooExtendPdf.h>
#include <RooAbsData.h>
int main(){
TFile input("/Users/imb/Downloads/dataset_region2_envelope_2023.root");
auto original=input.Get<RooWorkspace>("w");
RooWorkspace out("dcb_workspace");
out.import(*original->data("data_obs_resonant"));
auto mass=out.var("mass");
'''
for n in names+['yield_']:
 p=pars[n];s+=f'RooRealVar {n}("{n}","{n}",{p["value"]:.17g}); {n}.setError({p["error_hesse"]:.17g});\n'
s+='RooCrystalBall dcb("dcb","Double-sided Crystal Ball",*mass,mu,sigma,alphaL,nL,alphaR,nR);\nRooExtendPdf model("dcb_extended","Extended DCB",dcb,yield_); out.import(model);\nTFile output("outputs/dcb_fit.root","RECREATE"); out.Write();\nTMatrixDSym covariance(7);\n'
for i,row in enumerate(r['covariance']):
 for j,v in enumerate(row):s+=f'covariance({i},{j})={v:.17g};\n'
s+='covariance.Write("covariance_mu_sigma_alphaL_nL_alphaR_nR_yield");\n'
s+=f'TParameter<double>("chi2",{r["chi2"]:.17g}).Write(); TParameter<int>("ndof",{r["ndof"]}).Write();\n'
s+='TNamed("method","Chi-square with stored SumW2 errors; 129 zero-variance bins omitted; bin-integrated expectations. nR at lower bound; tail error unreliable.").Write();\noutput.Close();return 0;}\n'
pathlib.Path('work/save_model.cpp').write_text(s)
# Standalone numerical reproducer, with exported histogram included.
prefix=pathlib.Path('work/fit_dcb.py').read_text().split('qs=np.interp')[0]
prefix=prefix.replace("'work/data_obs_resonant.csv'","str(BASE/'data_obs_resonant.csv')").replace("'work/data_metadata.txt'","str(BASE/'data_metadata.txt')")
prefix="from pathlib import Path\nBASE=Path(__file__).resolve().parent\n"+prefix
suffix=pathlib.Path('work/fit_weighted.py').read_text().split('\n',1)[1]
combined=(prefix+'\n'+suffix).replace("'outputs/", "str(BASE/'").replace(".json'", ".json')").replace(".csv'", ".csv')").replace(".png'", ".png')").replace(".pdf'", ".pdf')")
# Explicit path substitutions avoid changing unrelated strings.
combined=prefix+'\n'+suffix
for stem in ['dcb_fit_results.json','dcb_covariance.csv','dcb_fit.png','dcb_fit.pdf','dcb_binned_fit.csv']:
 combined=combined.replace("'outputs/"+stem+"'", "str(BASE/'"+stem+"')")
pathlib.Path('outputs/refit_dcb.py').write_text(combined)
shutil.copy('work/data_obs_resonant.csv','outputs/data_obs_resonant.csv');shutil.copy('work/data_metadata.txt','outputs/data_metadata.txt')
report='''Fit Double Crystal Ball — data_obs_resonant

Input: dataset_region2_envelope_2023.root, workspace w.
Observable mass, range [6, 12], 450 native bins. Unit not specified in ROOT.
Weighted RooDataHist: sum(weights) = 1967902.5082.
Method: chi-square using the stored SumW2 bin errors, expected counts integrated
in each bin, seven free parameters including normalization. The 129 bins with
zero stored variance are excluded from chi-square; 321 bins are used. No rebinning.
This assumes independent bins and the stored errors; any inter-bin covariance
from data construction is unavailable. Empty zero-error bins do not constrain
this chi-square, so this is not a Poisson fit of all 450 bins.

'''
for n,p in pars.items(): report+=f'{n}: {p["value"]:.9g} +/- {p["error_hesse"]:.6g} (HESSE)\n'
report+=f'\nchi-square / ndof = {r["chi2"]:.6g} / {r["ndof"]}; ratio = {r["chi2"]/r["ndof"]:.5g}\nAsymptotic p-value = {r["pvalue_asymptotic"]:.4g} (approximate, not calibrated with toys).\n'
report+='''
MINOS intervals for mu and sigma are in the JSON file. Nine starting points
were tested. Minimum valid, covariance accurate locally; nR reaches its lower
bound 1.01 and its symmetric error must not be interpreted as a measurement.
The fit is not satisfactory over the full interval: a single DCB misses the
broad low-mass component and does not support a reliable global shape model.
The fitted yield is NOT the histogram integral and must not be interpreted as
an unbiased physical event yield with this misspecified model.
Consider a separate continuum component or a predeclared narrower fit range
before using the parameters for physics. Neither was silently added here.

DCB convention: t=(mass-mu)/sigma; Gaussian core for -alphaL <= t <= alphaR;
left/right power-law tails with exponents nL/nR, continuous in value and slope.
Normalization is over the finite fit range [6,12].

ROOT output: workspace dcb_workspace with data_obs_resonant, dcb and
 dcb_extended; parameter errors and covariance saved. No RooFitResult object
is fabricated: minimization was performed with iminuit.
Reproduce: install numpy scipy iminuit matplotlib, then run refit_dcb.py.
All the exported data and scripts must remain in the same directory.
Original input file was not modified.
'''
pathlib.Path('outputs/dcb_fit_report.txt').write_text(report)
