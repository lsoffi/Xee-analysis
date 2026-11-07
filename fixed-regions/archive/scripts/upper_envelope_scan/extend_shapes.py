from pathlib import Path
import ROOT,json,hashlib
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
P=Path.cwd();O=P/'runs/20261030_upper_envelope_scan/extended_signal_shapes';O.mkdir(exist_ok=True)
src=P/'data/mass_scan/261028/cards_region2_data_envelope_allCorrections_binned/ee/11.0/Xee_ee_4_allYears.root'
f=ROOT.TFile.Open(str(src));w=f.Get('w');nw=ROOT.RooWorkspace('w');mh=ROOT.RooRealVar('signal_mass','signal mass hypothesis',11.,6.,12.);mh.setConstant(True);keep=[mh]
eras=['2022','2022EE','2023','2023BPix'];oldmass=w.function('11');assert oldmass and oldmass.getVal()==11
for i,era in enumerate(eras,1):
    old=w.pdf(f'shapeSig_Zd_ch{i}');custom=ROOT.RooCustomizer(old,'parametric');custom.replaceArg(oldmass,mh);pdf=custom.build();pdf.SetName('signal_'+era);getattr(nw,'import')(pdf,ROOT.RooFit.RecycleConflictNodes());keep.extend([custom,pdf])
nw.var('mass').setVal(11.2);nw.var('signal_mass').setVal(11.)
checks=[]
for mass in [10.5,10.8,11.]:
    rf=ROOT.TFile.Open(str(src.parent.parent/f'{mass:.1f}'/src.name));rw=rf.Get('w');nw.var('signal_mass').setVal(mass)
    for i,era in enumerate(eras,1):
        for nuisance in ['CMS_scale_e',f'CMS_EXO25020_signalModelSigmaNuisance_{era}']:
            for shift in [-1,0,1]:
                nw.var(nuisance).setVal(shift);rw.var(nuisance).setVal(shift)
                a=nw.pdf('signal_'+era);b=rw.pdf(f'shapeSig_Zd_ch{i}')
                pairs=[(a.x0(),b.x0()),(a.sigmaL(),b.sigmaL()),(a.alphaL(),b.alphaL()),(a.nL(),b.nL()),(a.alphaR(),b.alphaR()),(a.nR(),b.nR())]
                delta=max(abs(x.getVal()-y.getVal()) for x,y in pairs);assert delta<1e-11,(mass,era,nuisance,shift,delta)
                checks.append(dict(mass=mass,era=era,nuisance=nuisance,shift=shift,max_parameter_difference=delta))
            nw.var(nuisance).setVal(0);rw.var(nuisance).setVal(0)
    rf.Close()
rows=[]
for k in range(111,121):
    mass=k/10;nw.var('signal_mass').setVal(mass);nw.var('mass').setVal(mass)
    for era in eras:
        p=nw.pdf('signal_'+era);a=dict(mass=mass,era=era,mean=p.x0().getVal(),sigma=p.sigmaL().getVal(),alphaL=p.alphaL().getVal(),nL=p.nL().getVal(),alphaR=p.alphaR().getVal(),nR=p.nR().getVal());assert all(a[q]>0 for q in ['sigma','alphaL','alphaR','nL','nR']);rows.append(a)
    nw.writeToFile(str(O/f'signal_shapes_m{mass:.1f}.root'))
nw.writeToFile(str(O/'parametric_signal_shapes.root'))
(O/'audit.json').write_text(json.dumps(dict(source=str(src),source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),shape_only=True,normalizations_available=False,definition='Original parametric RooFormulaVar graph; replace mass constant 11 with signal_mass, preserving MC coefficients and nuisance dependence. No yield assumptions.',checks=checks,parameters=rows),indent=2))
print('Created 10 mass points, 4 eras, shape only; validation checks',len(checks))
