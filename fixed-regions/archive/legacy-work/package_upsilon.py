import json,pathlib
r=json.load(open('outputs/upsilon_dcb_rebin2_results.json'));base=pathlib.Path('outputs')
prefix=pathlib.Path('work/fit_dcb.py').read_text().split('qs=np.interp')[0]
prefix=prefix.replace("'work/data_obs_resonant.csv'","str(BASE/'data_obs_resonant.csv')").replace("'work/data_metadata.txt'","str(BASE/'data_metadata.txt')")
source="from pathlib import Path\nBASE=Path(__file__).resolve().parent\n"+prefix+'\n'+pathlib.Path('work/fit_upsilon.py').read_text().split('\n',1)[1]
for n in ['upsilon_dcb_rebin2_results.json','upsilon_dcb_rebin2_bins.csv','upsilon_dcb_rebin2.png','upsilon_dcb_rebin2.pdf']:source=source.replace("'outputs/"+n+"'","str(BASE/'"+n+"')")
(base/'refit_upsilon_rebin2.py').write_text(source)
report=f'''Double-sided Crystal Ball — Υ(1S), signal only
Dataset data_obs_resonant from workspace w. No background component.
Nominal range 8.8–9.8; actual range {r['actual_range_aligned_to_bins']}
chosen to include complete bins after rebinning adjacent native bins by two.
Rebinning preserves sum(w) and sum(w^2). Native width 0.0133333,
rebinned width 0.0266667; 38 bins, none empty.

Fit: iminuit ExtendedBinnedNLL weighted Böhm–Zech approximation, exact bin
integrals of the DCB normalized over the fit window; floating normalization.
This approximation assumes appropriate independent-event weights and no
unmodeled inter-bin correlations. The construction of the weights was not
provided. It is not an ordinary Poisson likelihood for the bin contents.

'''
for p,v in r['parameters'].items():report+=f"{p} = {v['value']:.9g} +/- {v['error_hesse']:.6g} (local HESSE)\n"
report+='\nMINOS profile intervals:\n'
for p,v in r['minos'].items():report+=f"{p}: {v['lower']:+.7g}, {v['upper']:+.7g}; valid={v['valid']}\n"
report+=f'''\nWeighted deviance / nominal ndof = {r['weighted_deviance']:.5f}/{r['ndof_nominal']}
Nominal asymptotic p-value = {r['asymptotic_pvalue_nominal']:.4f}.
This is not calibrated for parameter boundaries; it is only indicative.
SumW2 residual chi-square (diagnostic, not minimized) = {r['chi2_sumw2_diagnostic']:.5f}.

alphaR reaches the lower boundary and nR the upper boundary. Tail parameters
and Gaussian-core sigma are poorly identified. The broad asymmetric MINOS
interval for sigma is the relevant uncertainty; the small symmetric HESSE
error is not sufficient. Do not interpret sigma alone as a robust physical
resolution. The plot describes the signal shape in the window reasonably,
but parameter identification is not resolved by the improved fit statistic.
Nine starting points were tested. Units were not specified in the input.
The old full-range chi-square 560/314 used another objective/range/binning
and cannot be treated as a nested-model likelihood-ratio comparison.

ROOT output: workspace upsilon_dcb_workspace with rebinned data, dcb and
dcb_extended; covariance, weighted deviance and nominal ndof stored separately.
No RooFitResult is fabricated (minimization used iminuit).
Reproduce with refit_upsilon_rebin2.py and the supplied data_obs_resonant.csv
and data_metadata.txt in the same directory; requires numpy scipy iminuit matplotlib.
'''
(base/'upsilon_dcb_rebin2_report.txt').write_text(report)
# Reuse the verified ROOT writer, but replace dataset, ranges and fit values.
old=json.load(open('outputs/dcb_fit_results.json'))
s=pathlib.Path('work/save_model.cpp').read_text()
a=s.index('TFile input(');b=s.index('RooRealVar mu(')
import numpy as np
bins=np.genfromtxt('outputs/upsilon_dcb_rebin2_bins.csv',delimiter=',',names=True)
initial='RooWorkspace out("upsilon_dcb_workspace");\nRooRealVar massvar("mass","mass",'+str(r['actual_range_aligned_to_bins'][0])+','+str(r['actual_range_aligned_to_bins'][1])+');\nmassvar.setBins(38);\nTH1D histogram("rebinned","rebinned",38,8.8,'+str(r['actual_range_aligned_to_bins'][1])+');histogram.Sumw2();\n'
for i,row in enumerate(bins):initial+=f'histogram.SetBinContent({i+1},{row["sumw"]:.17g});histogram.SetBinError({i+1},{np.sqrt(row["sumw2"]):.17g});\n'
initial+='RooDataHist data("data_obs_resonant","signal, rebin 2",RooArgList(massvar),RooFit::Import(histogram));out.import(data);auto mass=out.var("mass");\n'
s=s[:a]+initial+s[b:]
s=s.replace('#include <fstream>','#include <fstream>\n#include <TH1D.h>\n#include <RooDataHist.h>\n#include <RooArgList.h>')
import re
for n,v in r['parameters'].items():
 lo,hi=v['limits'];s=re.sub(r'RooRealVar '+n+r'\([^;]+; '+n+r'\.setError\([^;]+;',f'RooRealVar {n}("{n}","{n}",{v["value"]:.17g},{lo},{hi}); {n}.setError({v["error_hesse"]:.17g});',s)
for i,row in enumerate(r['covariance']):
 for j,v in enumerate(row):s=re.sub(r'covariance\('+str(i)+','+str(j)+r'\)=[^;]+;',f'covariance({i},{j})={v:.17g};',s)
s=s.replace('outputs/dcb_fit.root','outputs/upsilon_dcb_rebin2.root')
s=re.sub(r'TParameter<double>\("chi2",[^;]+;',f'TParameter<double>("weighted_deviance",{r["weighted_deviance"]}).Write();',s)
s=re.sub(r'TParameter<int>\("ndof",[^;]+;',f'TParameter<int>("nominal_ndof",{r["ndof_nominal"]}).Write();',s)
s=re.sub(r'TNamed\("method",[^;]+;', 'TNamed("method","Weighted Bohm-Zech extended binned likelihood. Signal only, rebin 2, 8.8-9.813333. alphaR/nR at limits; use MINOS errors from JSON.").Write();',s)
s=s[:s.index('TFile check(')]+'return 0;}\n'
pathlib.Path('work/save_upsilon.cpp').write_text(s)
