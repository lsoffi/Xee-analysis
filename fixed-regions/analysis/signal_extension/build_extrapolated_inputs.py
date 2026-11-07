from pathlib import Path
import ROOT,json,hashlib
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
P=Path.cwd();R=P/'runs/20261030_upper_envelope_scan';source=P/'data/mass_scan/261028/cards_region2_data_envelope_allCorrections_binned/ee/11.0/Xee_ee_4_allYears.root';rates=json.loads((R/'signal_rate_extrapolation.json').read_text());audit=[]
for k in range(111,121):
    mass=k/10;D=R/f'extrapolated_inputs/m{mass:.1f}';D.mkdir(parents=True,exist_ok=True);f=ROOT.TFile.Open(str(source));w=f.Get('w');sf=ROOT.TFile.Open(str(R/f'extended_signal_shapes/signal_shapes_m{mass:.1f}.root'));sw=sf.Get('w');nw=ROOT.RooWorkspace('w');keep=[]
    custom=ROOT.RooCustomizer(w.pdf('model_s'),'extended')
    for i,a in enumerate(rates,1):
        era=a['era'];newshape=sw.pdf('signal_'+era).clone(f'shapeSig_Zd_ch{i}');oldrate=w.function(f'n_exp_binch{i}_proc_Zd');w.var('r').setVal(1);oldval=oldrate.getVal();pred=next(p['rate'] for p in a['predictions'] if abs(p['mass']-mass)<1e-8)
        newrate=ROOT.RooFormulaVar(f'extrapolated_rate_Zd_{era}',f'@0*{pred/oldval:.17g}',ROOT.RooArgList(oldrate));newrate.setAttribute('combine.signal',True);newrate.setStringAttribute('combine.channel',f'ch{i}');newrate.setStringAttribute('combine.process','Zd')
        custom.replaceArg(w.pdf(f'shapeSig_Zd_ch{i}'),newshape);custom.replaceArg(oldrate,newrate);keep.extend([newshape,newrate]);audit.append(dict(mass=mass,era=era,rate=pred,reference_rate=oldval))
    sp=custom.build();sp.SetName('model_s');getattr(nw,'import')(sp,ROOT.RooFit.RecycleConflictNodes());getattr(nw,'import')(w.pdf('model_b'),ROOT.RooFit.RecycleConflictNodes())
    for mcname,pdfname in [('ModelConfig','model_s'),('ModelConfig_bonly','model_b')]:
        old=w.obj(mcname);mc=ROOT.RooStats.ModelConfig(mcname,nw);mc.SetPdf(nw.pdf(pdfname))
        for get,setter in [('GetParametersOfInterest','SetParametersOfInterest'),('GetObservables','SetObservables'),('GetNuisanceParameters','SetNuisanceParameters'),('GetGlobalObservables','SetGlobalObservables')]:
            args=ROOT.RooArgSet()
            for v in getattr(old,get)():
                if nw.arg(v.GetName()):args.add(nw.arg(v.GetName()))
            getattr(mc,setter)(args)
        getattr(nw,'import')(mc);keep.append(mc)
    getattr(nw,'import')(w.data('data_obs'));nw.var('r').setVal(1)
    # Regional builder expects this conventional node name; point it to the scaled rate.
    for i,a in enumerate(rates,1):
        node=nw.function(f'extrapolated_rate_Zd_{a["era"]}');assert abs(node.getVal()-next(p['rate'] for p in a['predictions'] if p['mass']==mass))<1e-6
    nw.writeToFile(str(D/'Xee_ee_4_allYears.root'));f.Close();sf.Close()
(R/'extrapolated_inputs/audit.json').write_text(json.dumps(dict(original_source=str(source),original_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),rates=audit,relative_lognormal_uncertainties='held at m11 values',parametric_shape=True,original_data_preserved=True),indent=2))
