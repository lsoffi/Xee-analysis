import ROOT,json
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit');ROOT.gROOT.SetBatch(True)
for case in ['free','limit1','limit2']:
 f=ROOT.TFile.Open('outputs/Xee_ee_4_2022_chebyshev_degree7_to10_'+case+'.root');w=f.Get('w');x=w.var('mass');obs=ROOT.RooArgSet(x);pdf=w.pdf('dy_component_1_2022');m=json.load(open('outputs/combine_2022_degree7_to10_'+case+'_pre_post.json'));params={a['parameter']:a for a in m['parameters']}
 for a in m['parameters']:w.var(a['parameter']).setVal(a['s'])
 def values():
  out=[]
  for mass in [9.5,9.75,10]:x.setVal(mass);out.append(pdf.getVal(obs))
  return out
 post=values();base=post[2]-2*post[1]+post[0];tests=[]
 for i in range(7):
  name='t%d_2022_envelope'%i;v=w.var(name);old=v.getVal();v.setVal(params[name]['prefit']);yy=values();v.setVal(old)
  change=(base-(yy[2]-2*yy[1]+yy[0]))/post[1]
  tests.append((name,change))
 print(case,'relative curvature',base/post[1],'leave-one-pre contributions',sorted(tests,key=lambda q:abs(q[1]),reverse=True))
