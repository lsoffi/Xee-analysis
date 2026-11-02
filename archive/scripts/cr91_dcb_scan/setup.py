from pathlib import Path
P=Path.cwd();S=P/'scripts/cr91_dcb_scan'
s=(P/'scripts/common95_cheb2/prepare.py').read_text();s=s.replace("MASS=sys.argv[1];MODEL='chebyshev2'","MASS=sys.argv[1];REGION=sys.argv[2];MODEL='chebyshev3'").replace("P/'runs/20261029_common95_cheb2'/('m'+MASS)","P/'runs/20261029_cr91_dcb_scan'/('m'+MASS)/REGION")
a=s.index(' c1=');b=s.index(' yn=ROOT.RooRealVar',a);src=(P/'scripts/upsilon_control_cheb3/prepare_signal.py').read_text();aa=src.index(' c1=');bb=src.index(' yn=ROOT.RooRealVar',aa);s=s[:a]+src[aa:bb]+s[b:]
s=s.replace("f=ROOT.TFile.Open(str(src));", "if REGION=='control':LO=9.093333333333334;HI=10.506666666666668;NB=106\nf=ROOT.TFile.Open(str(src));")
s=s.replace("if proc=='Zd':assert frac>=.95-1e-7", "if proc=='Zd' and REGION=='signal':assert frac>=.95-1e-7")
s=s.replace('x.setVal(float(MASS))','x.setVal((LO+HI)/2)').replace("nw.var('mass').setVal(float(MASS))","nw.var('mass').setVal((LO+HI)/2)").replace('Degree2','Degree3');(S/'prepare.py').write_text(s)
s=(P/'scripts/cr_windows/cr91_105/degree3/build_joint.py').read_text().replace('import ROOT,pathlib,json','import ROOT,pathlib,json,sys\nMASS=sys.argv[1]');s=s.replace("P/'runs/20261029_cr_windows/cr91_105/degree3'","P/'runs/20261029_cr91_dcb_scan'/('m'+MASS)").replace("P/'runs/20261029_upsilon_control_cheb3/signal/workspace.root'","R/'signal/workspace.root'")
a=s.index("custom=ROOT.RooCustomizer");b=s.index("cat=ROOT.RooCategory",a);block=s[a:b];block=block.replace("custom=ROOT.RooCustomizer(cr.obj('ModelConfig_bonly').GetPdf(),'CR')","custom=ROOT.RooCustomizer(cr.obj(mcname).GetPdf(),'CR')").replace("['upsilon1s','upsilon2s']","['Zd','upsilon1s','upsilon2s']")
# Keep both S+B and B-only CR models, with the physical signal intensity shared.
block='crpdfs={}\nfor mcname in [\'ModelConfig\',\'ModelConfig_bonly\']:\n'+''.join(' '+line+'\n' for line in block.splitlines())+" crpdfs[mcname]=cp\n";s=s[:a]+block+s[b:];s=s.replace("source=sr.obj(mcname);sim=", "cp=crpdfs[mcname];source=sr.obj(mcname);sim=")
s=s.replace("nw.var('mass').setVal(7.5)","nw.var('mass').setVal(float(MASS))").replace('CR_signal_omitted_negligible_gaussian_acceptance=True','CR_signal_included=True').replace("for proc in ['upsilon1s','upsilon2s']:","for proc in ['upsilon1s','upsilon2s']:")
s=s.replace("# Customize CR graph:","assert sr.var('mass').getMax() <= cr.var('mass').getMin()+1e-8\n# Customize CR graph:")
(S/'build_joint.py').write_text(s)
s=(P/'scripts/upsilon_control_cheb3/run_joint.py').read_text().replace("P/'runs/20261029_upsilon_control_cheb3'/'joint'","P/'runs/20261029_cr91_dcb_scan'/('m'+MASS)/'joint'");(S/'run.py').write_text(s)
s=(P/'scripts/common95_cheb2_gaussian/quality.py').read_text().replace("P/'runs/20261029_common95_cheb2_gaussian'/('m'+MASS)","P/'runs/20261029_cr91_dcb_scan'/('m'+MASS)/'joint'");(S/'quality.py').write_text(s)
