from pathlib import Path
import ROOT,json,sys,math,re,hashlib,csv,numpy as np
ROOT.gROOT.SetBatch(True);ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
P=Path.cwd();m=sys.argv[1];degree=int(sys.argv[2]);D=P/f'runs/20261030_upper_envelope_scan/degree{degree}/m{m}';J=D/'joint'
f=ROOT.TFile.Open(str(J/'fitdiagnostics/fitDiagnostics_upper.root'));fw=ROOT.TFile.Open(str(J/'workspace_seeded.root'));w=fw.Get('w');eras=['2022','2022EE','2023','2023BPix'];result=dict(mass=float(m),degree=degree,fit={});parameter_rows=[]
for stage in ['b','s']:
    fit=f.Get('fit_'+stage);pars=[]
    for v in fit.floatParsFinal():
        w.var(v.GetName()).setVal(v.getVal());a=dict(name=v.GetName(),value=v.getVal(),error=v.getError(),low=v.getMin(),high=v.getMax());pars.append(a);parameter_rows.append(dict(stage=stage,**a))
    w.var('r').setVal(0 if stage=='b' else fit.floatParsFinal().find('r').getVal())
    indices={};mins=[]
    for era in eras:
        bp=fit.floatParsFinal().find('local_b1_'+era);cp=fit.floatParsFinal().find('local_c1_'+era)
        idx=1 if bp and bp.getError()>0 and (not cp or cp.getError()==0) else 0;w.cat('pdfindex_'+era).setIndex(idx);indices[era]=idx
        for suffix,deg in [('',degree),('_CR',4)]:
            if suffix=='' and idx==1:mins.append(min([1]+[w.var(f'local_b{k}_{era}').getVal() for k in range(1,deg+1)]));continue
            coeff=[1]+[w.var(f'local_c{k}_{era}{suffix}').getVal() for k in range(1,deg+1)];roots=np.polynomial.chebyshev.chebroots(np.polynomial.chebyshev.chebder(coeff));points=[-1,1]+[q.real for q in roots if abs(q.imag)<1e-8 and -1<q.real<1];mins.append(float(min(np.polynomial.chebyshev.chebval(points,coeff))))
    for mu in [v for v in w.allVars() if v.GetName().startswith('Ymean_')]:assert mu.isConstant() and not fit.floatParsFinal().find(mu.GetName())
    a=dict(status=fit.status(),covQual=fit.covQual(),edm=fit.edm(),parameters=pars,indices=indices,minimum=min(mins),r=w.var('r').getVal(),error=fit.floatParsFinal().find('r').getError() if stage=='s' else None)
    result['fit'][stage]=a
    w.saveSnapshot('postfit_'+stage,w.allVars());w.writeToFile(str(J/f'workspace_postfit_{stage}.root'))
    if stage=='s':result.update(r=a['r'],r_error=a['error'])
for role,key in [('significance','Z'),('gof_observed','gof_observed')]:
    rf=ROOT.TFile.Open(str(next((J/role).glob('higgsCombine*.root'))));tr=rf.Get('limit');tr.GetEntry(0);result[key]=float(tr.limit);rf.Close()
result['p0']=.5*math.erfc(result['Z']/math.sqrt(2))
rf=ROOT.TFile.Open(str(next((J/'limits').glob('higgsCombine*.root'))));t=rf.Get('limit');vals=[]
for e in t:vals.append((float(e.quantileExpected),float(e.limit)))
assert len(vals)==6;result['limits']={k:min(vals,key=lambda a:abs(a[0]-q))[1] for k,q in [('observed',-1),('expected_2p5',.025),('expected_16',.16),('expected_50',.5),('expected_84',.84),('expected_97p5',.975)]};rf.Close()
log=(J/'limits/combine.log').read_text();pairs=re.findall(r'At r =\s*([\d.eE+-]+):[^\n]*CLs\s*=\s*([\d.eE+-]+)',log)
rr,cl=min(pairs,key=lambda a:abs(float(a[0])-result['limits']['observed']));result['observed_limit_check']=dict(r=float(rr),CLs=float(cl),verified=abs(float(cl)-.05)<.002)
t=f.Get('tree_fit_sb');t.GetEntry(0);result['Z_FD']=math.sqrt(max(0,-2*float(t.nll_nll0))) if result['r']>0 else 0
result['significance_consistent']=abs(result['Z_FD']-result['Z'])<.05
result['numerical_ok']=all(a['status']==0 and a['covQual']==3 and a['minimum']>0 for a in result['fit'].values()) and result['significance_consistent'] and result['observed_limit_check']['verified'] and abs(result['r'])<99
for region in ['signal','control']:
    audit=json.loads((D/region/'audit.json').read_text());assert hashlib.sha256(Path(audit['source']).read_bytes()).hexdigest()==audit['source_sha256']
(D/'results.json').write_text(json.dumps(result,indent=2))
with (D/'parameters.csv').open('w') as h:
    writer=csv.DictWriter(h,fieldnames=parameter_rows[0].keys());writer.writeheader();writer.writerows(parameter_rows)
print(m,degree,result['numerical_ok'],result['r'],result['Z'],result['limits'],flush=True)
