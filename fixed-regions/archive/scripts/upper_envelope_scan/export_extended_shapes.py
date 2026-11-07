from pathlib import Path
import ROOT,json
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
R=Path.cwd()/'runs/20261030_upper_envelope_scan';f=ROOT.TFile.Open(str(R/'extended_signal_shapes/parametric_signal_shapes.root'));w=f.Get('w');x=w.var('mass');obs=ROOT.RooArgSet(x);grid=[9+i/400 for i in range(1201)];rows=[]
for k in range(110,121):
    mass=k/10;w.var('signal_mass').setVal(mass)
    for era in ['2022','2022EE','2023','2023BPix']:
        pdf=w.pdf('signal_'+era);vals=[]
        for value in grid:x.setVal(value);vals.append(pdf.getVal(obs))
        rows.append(dict(mass=mass,era=era,density=vals))
(R/'extended_shapes_plot_data.json').write_text(json.dumps(dict(x=grid,rows=rows,normalization_range=[6,12])))
