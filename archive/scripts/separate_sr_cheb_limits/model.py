"""Separate low or high SR plus the common CR; preparation only, no fits.
Single local Chebyshev6 (low) or Chebyshev4 (high); CR Chebyshev4.
All eight Upsilon means and normalizations float in both-region likelihood.
Original full-[6,12] DCB rates are multiplied by dynamic regional acceptances.
"""
from pathlib import Path
import argparse, json, hashlib, math
import numpy as np
from numpy.polynomial import Chebyshev, Polynomial
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
P=Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee')
R=P/'runs/20261031_separate_sr_cheb_limits'
ERAS=['2022','2022EE','2023','2023BPix']
REGIONS={'SRL':(6.,9.093333333333334,232),
         'CR':(9.093333333333334,10.506666666666668,106),
         'SRH':(10.506666666666668,12.,112)}
CR_REFERENCE=P/'runs/20261030_upper_sr_cr_mu_fixed/cr_only/results.json'

def source_path(mass):
    mass=f'{float(mass):.1f}'
    if float(mass)<=11:
        return P/'data/mass_scan/261028/cards_region2_data_envelope_allCorrections_binned/ee'/mass/'Xee_ee_4_allYears.root'
    return P/'runs/20261030_upper_envelope_scan/extrapolated_inputs'/('m'+mass)/'Xee_ee_4_allYears.root'

def _import(w,obj):
    getattr(w,'import')(obj,ROOT.RooFit.RecycleConflictNodes(),ROOT.RooFit.Silence())
    return w.arg(obj.GetName()) or w.obj(obj.GetName())

def _poly_integrals(degree,family,lo,hi):
    if family=='cheb':
        ans=[]
        for k in range(degree+1):
            a=Chebyshev.basis(k).integ()
            ans.append(3*(a((hi-9)/3)-a((lo-9)/3)))
        return ans
    ans=[]
    for k in range(degree+1):
        p=math.comb(degree,k)*Polynomial([0,1])**k*Polynomial([1,-1])**(degree-k)
        a=p.integ();ans.append(6*(a((hi-6)/6)-a((lo-6)/6)))
    return ans

def _poly_integral_node(w,name,coeff,weights,constant_one=False):
    if constant_one:
        expr=f'({weights[0]:.17g})'+''.join(f'+@{k-1}*({weights[k]:.17g})' for k in range(1,len(weights)))
    else:expr='+'.join(f'@{k}*({v:.17g})' for k,v in enumerate(weights))
    return _import(w,ROOT.RooFormulaVar(name,expr,ROOT.RooArgList(*coeff)))

def _minimum_cheb(coeff):
    rr=np.polynomial.chebyshev.chebroots(np.polynomial.chebyshev.chebder(coeff))
    pp=[-1,1]+[q.real for q in rr if abs(q.imag)<1e-8 and -1<q.real<1]
    return float(min(np.polynomial.chebyshev.chebval(pp,coeff)))

def _positive_constraint(w,name,coeff):
    checks=[]
    for t in np.linspace(-1,1,257):
        values=[Chebyshev.basis(k)(t) for k in range(1,len(coeff)+1)]
        checks.append('(1'+''.join(f'+@{j}*({v:.17g})' for j,v in enumerate(values))+'>1e-9)')
    return _import(w,ROOT.RooGenericPdf(name,'&&'.join(checks),ROOT.RooArgList(*coeff)))

def _dcb_acceptance(w,proc,era,region,pdf,regions):
    lo,hi,_=regions[region]
    params=ROOT.RooArgList(pdf.x0(),pdf.sigmaL(),pdf.alphaL(),pdf.nL(),pdf.alphaR(),pdf.nR())
    primitives=[]
    for edge in [6,lo,hi,12]:
        name=f'dcb_primitive_{str(edge).replace(".","p")}_{proc}_{era}'
        node=w.function(name)
        if not node:
            t=f'(({edge:.17g}-@0)/@1)';al='abs(@2)';ar='abs(@4)';nl='@3';nr='@5'
            gl=f'sqrt(TMath::Pi()/2)*TMath::Erf(-{al}/sqrt(2))';gr=f'sqrt(TMath::Pi()/2)*TMath::Erf({ar}/sqrt(2))'
            bl=f'({nl}/{al}-{al})';br=f'({nr}/{ar}-{ar})'
            left=f'pow({nl}/{al},{nl})*exp(-0.5*{al}*{al})*((abs({nl}-1)<1e-8)?log(({bl}+{al})/({bl}-{t})):(pow({bl}-{t},1-{nl})-pow({bl}+{al},1-{nl}))/({nl}-1))'
            right=f'pow({nr}/{ar},{nr})*exp(-0.5*{ar}*{ar})*((abs({nr}-1)<1e-8)?log(({br}+{t})/({br}+{ar})):(pow({br}+{t},1-{nr})-pow({br}+{ar},1-{nr}))/(1-{nr}))'
            expr=f'({t}<-{al})?({gl}+{left}):(({t}>{ar})?({gr}+{right}):(sqrt(TMath::Pi()/2)*TMath::Erf({t}/sqrt(2))))'
            node=_import(w,ROOT.RooFormulaVar(name,expr,params))
        primitives.append(node)
    return _import(w,ROOT.RooFormulaVar(f'local_acceptance_{proc}_{era}_{region}','(@2-@1)/(@3-@0)',ROOT.RooArgList(*primitives)))

SEEDS={
    'low':P/'runs/20261030_sr6_cr4_scan/m7.5/joint/fitdiagnostics/fitDiagnostics_local.root',
    'high':P/'runs/20261030_upper_sr_cr/degree3/joint/fitdiagnostics/fitDiagnostics_upper.root'}

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def regions_for(side):
    return {'SR':REGIONS['SRL' if side=='low' else 'SRH'],'CR':REGIONS['CR']}

def initial_values(side):
    """Previously converged B fit values are numerical starts, never priors."""
    source=SEEDS[side];handle=ROOT.TFile.Open(str(source));fit=handle.Get('fit_b')
    assert fit and fit.status()==0 and fit.edm()<.01
    vals={v.GetName():float(v.getVal()) for v in fit.floatParsFinal()}
    indices={v.GetName():v.getIndex() for v in fit.constPars() if v.InheritsFrom('RooAbsCategory')}
    provenance=dict(path=str(source),sha256=sha(source),fit='fit_b',status=fit.status(),covQual=fit.covQual(),EDM=fit.edm(),era_transformations={})
    for era in ERAS:
        if side=='high':
            vals[f'local_c4_{era}']=0.
            provenance['era_transformations'][era]='Exact degree3-to-degree4 embedding, c4=0; same local SR domain'
        elif indices.get('pdfindex_'+era,0)==1:
            # The old envelope selected Bernstein6 for this era. Convert its
            # fitted local polynomial exactly, then fix Chebyshev c0=1.
            b=[1.]+[vals[f'local_b{k}_{era}'] for k in range(1,7)]
            u=Polynomial([.5,.5]);poly=sum((math.comb(6,k)*b[k]*u**k*(1-u)**(6-k) for k in range(7)),Polynomial([0.]))
            c=poly.convert(kind=Chebyshev).coef;c=c/c[0]
            for k in range(1,7):vals[f'local_c{k}_{era}']=float(c[k])
            provenance['era_transformations'][era]='Exact fitted Bernstein6-to-Chebyshev6 conversion on same local domain, c0 normalized to 1'
        else:
            provenance['era_transformations'][era]='Fitted Chebyshev6 coefficients on identical local SR domain'
        degree=6 if side=='low' else 4
        assert _minimum_cheb([1]+[vals[f'local_c{k}_{era}'] for k in range(1,degree+1)])>1e-8
        assert _minimum_cheb([1]+[vals[f'local_c{k}_{era}_CR'] for k in range(1,5)])>1e-8
    handle.Close();return vals,provenance

def build(side,mass,output=None):
    assert side in ['low','high']
    mass=f'{float(mass):.1f}';degree=6 if side=='low' else 4;regions=regions_for(side)
    output=Path(output or R/side/f'm{mass}'/'workspace.root')
    if output.exists():raise FileExistsError(f'Preserving existing workspace: {output}')
    output.parent.mkdir(parents=True,exist_ok=True)
    source=source_path(mass);source_hash=sha(source)
    sf=ROOT.TFile.Open(str(source));old=sf.Get('w');assert old
    w=ROOT.RooWorkspace('w');keep=[];seed,seed_provenance=initial_values(side)
    for reg,(lo,hi,nb) in regions.items():
        x=ROOT.RooRealVar('mass_'+reg,'mass [GeV]',(lo+hi)/2,lo,hi);x.setBins(nb);_import(w,x)
    cat=ROOT.RooCategory('CMS_channel','CMS_channel')
    for reg in regions:
        for i in range(1,5):cat.defineType(f'{reg}_ch{i}')
    _import(w,cat);cat=w.cat('CMS_channel');old.var('r').setVal(1)
    mcold=old.obj('ModelConfig');nuis=list(mcold.GetNuisanceParameters());globs=list(mcold.GetGlobalObservables())
    constraints=list(mcold.GetPdf().getAllConstraints(mcold.GetObservables(),mcold.GetNuisanceParameters(),False));assert len(constraints)==10
    for con in constraints:_import(w,con)
    counts={f'{reg}_ch{i}':0. for reg in regions for i in range(1,5)};rows=[]
    original=old.data('data_obs')
    for j in range(original.numEntries()):
        row=original.get(j);xx=row.getRealValue('mass');ch=row.find('CMS_channel').getLabel();weight=original.weight()
        selected=[reg for reg,(lo,hi,_) in regions.items() if lo<=xx<hi]
        assert len(selected)<=1
        if selected:
            reg=selected[0];counts[reg+'_'+ch]+=weight;rows.append((reg,ch,xx,weight))
    positivity=[];era_audits=[]
    for i,era in enumerate(ERAS,1):
        ch=f'ch{i}';coeffs={};dy_rates={}
        for reg in regions:
            suffix='' if reg=='SR' else '_CR';deg=degree if reg=='SR' else 4
            coeffs[reg]=[]
            for k in range(1,deg+1):
                name=f'local_c{k}_{era}{suffix}';v=seed[name];bound=2 if k==1 else 1
                assert -bound<v<bound,(name,v)
                coeffs[reg].append(_import(w,ROOT.RooRealVar(name,'',v,-bound,bound)))
            positivity.append(_positive_constraint(w,'positive_continuum_'+era+suffix,coeffs[reg]))
            name='local_yield_'+era+suffix;n=counts[reg+'_'+ch]
            dy_rates[reg]=_import(w,ROOT.RooRealVar(name,'regional continuum yield',seed[name],0,max(10.,5*n)))
        fullrates={}
        for proc in ['Zd','upsilon1s','upsilon2s']:
            rate=old.function(f'extrapolated_rate_Zd_{era}') if proc=='Zd' and float(mass)>11 else old.function(f'n_exp_bin{ch}_proc_{proc}')
            assert rate;fullrates[proc]=_import(w,rate)
            if proc!='Zd':
                scale=w.var(f'CMS_EXO25020_bkgScale_{proc}_{era}');scale.setVal(seed[scale.GetName()]);scale.setConstant(False)
                op=old.pdf(f'shapeBkg_{proc}_{ch}');name=f'Ymean_{proc}_{era}'
                mu=ROOT.RooRealVar(name,'',seed[name],op.x0().getVal()-.2,op.x0().getVal()+.2);mu.setConstant(False);_import(w,mu)
                _import(w,ROOT.RooConstVar(f'Ysigma_{proc}_{era}','',op.sigmaL().getVal()))
                for k,arg in enumerate([op.alphaL(),op.nL(),op.alphaR(),op.nR()]):_import(w,ROOT.RooConstVar(f'Ytail_{proc}_{era}_{k}','',arg.getVal()))
        record=dict(era=era,processes={})
        for reg,(lo,hi,nb) in regions.items():
            x=w.var('mass_'+reg);cont=_import(w,ROOT.RooChebychev(f'local_continuum_{era}_{reg}','',x,ROOT.RooArgList(*coeffs[reg])))
            dy_rate=dy_rates[reg];dy_rate.setStringAttribute('combine.channel',reg+'_'+ch);dy_rate.setStringAttribute('combine.process','dy')
            pdfs={};rates={}
            for proc in ['Zd','upsilon1s','upsilon2s']:
                if proc=='Zd':
                    op=old.pdf('shapeSig_Zd_'+ch);cu=ROOT.RooCustomizer(op,reg);cu.replaceArg(old.var('mass'),x);pdf=cu.build();pdf.SetName(f'shape_{proc}_{era}_{reg}');keep.extend([cu,pdf]);pdf=_import(w,pdf)
                else:
                    sg=w.function(f'Ysigma_{proc}_{era}');tails=[w.function(f'Ytail_{proc}_{era}_{k}') for k in range(4)]
                    pdf=_import(w,ROOT.RooCrystalBall(f'shape_{proc}_{era}_{reg}','',x,w.var(f'Ymean_{proc}_{era}'),sg,sg,*tails))
                acceptance=_dcb_acceptance(w,proc,era,reg,pdf,regions)
                rate=_import(w,ROOT.RooFormulaVar(f'local_rate_{proc}_{era}_{reg}','@0*@1',ROOT.RooArgList(fullrates[proc],acceptance)))
                rate.setAttribute('combine.signal',proc=='Zd');rate.setStringAttribute('combine.channel',reg+'_'+ch);rate.setStringAttribute('combine.process',proc)
                pdfs[proc]=pdf;rates[proc]=rate
                record['processes'].setdefault(proc,{})[reg]=dict(acceptance=acceptance.getVal(),rate_r1=rate.getVal(),full_rate_r1=fullrates[proc].getVal())
            for stage in ['s','b']:
                pp=[cont,pdfs['upsilon1s'],pdfs['upsilon2s']];rr=[dy_rate,rates['upsilon1s'],rates['upsilon2s']]
                if stage=='s':pp=[pdfs['Zd']]+pp;rr=[rates['Zd']]+rr
                _import(w,ROOT.RooAddPdf(f'pdf_{stage}_{reg}_{ch}','',ROOT.RooArgList(*pp),ROOT.RooArgList(*rr)))
        for proc in ['Zd','upsilon1s','upsilon2s']:
            assert 0<sum(q['acceptance'] for q in record['processes'][proc].values())<=1+1e-10
        era_audits.append(record)
    for stage,mcname in [('s','ModelConfig'),('b','ModelConfig_bonly')]:
        sim=ROOT.RooSimultaneousOpt('joint_'+mcname,'',cat)
        for reg in regions:
            for i in range(1,5):sim.addPdf(w.pdf(f'pdf_{stage}_{reg}_ch{i}'),f'{reg}_ch{i}')
        for con in constraints+positivity:sim.addExtraConstraint(w.pdf(con.GetName()))
        _import(w,sim);mc=ROOT.RooStats.ModelConfig(mcname,w);mc.SetPdf(w.pdf(sim.GetName()))
        mc.SetObservables(ROOT.RooArgSet(*[w.var('mass_'+reg) for reg in regions],cat));mc.SetParametersOfInterest(ROOT.RooArgSet(w.var('r')))
        mc.SetNuisanceParameters(ROOT.RooArgSet(*[w.var(v.GetName()) for v in nuis]));mc.SetGlobalObservables(ROOT.RooArgSet(*[w.var(v.GetName()) for v in globs]))
        getattr(w,'import')(mc);keep.append(mc)
    weight=ROOT.RooRealVar('evt_weight','',0,1e9);obs=ROOT.RooArgSet(*[w.var('mass_'+reg) for reg in regions],cat,weight)
    data=ROOT.RooDataSet('data_obs','',obs,ROOT.RooFit.WeightVar(weight))
    for reg,ch,xx,n in rows:
        cat.setLabel(reg+'_'+ch);w.var('mass_'+reg).setVal(xx);weight.setVal(n);data.add(obs,n)
    getattr(w,'import')(data)
    for reg,(lo,hi,_) in regions.items():w.var('mass_'+reg).setVal((lo+hi)/2)
    w.var('r').setVal(0);w.var('r').setRange(-1000,1000);w.var('r').setConstant(False)
    assert data.numEntries()==4*sum(nb for lo,hi,nb in regions.values())
    assert abs(data.sumEntries()-sum(counts.values()))<1e-8
    assert len([v for v in w.allVars() if v.GetName().startswith('Ymean_') and not v.isConstant()])==8
    assert not any(p.InheritsFrom('RooMultiPdf') for p in w.allPdfs())
    w.writeToFile(str(output));sf.Close();assert sha(source)==source_hash
    audit=dict(side=side,mass=float(mass),degree=degree,CR_degree=4,family='cheb',source=str(source),source_sha256=source_hash,
               regions=regions,counts=counts,entries=data.numEntries(),events=data.sumEntries(),eras=era_audits,
               initial_values=seed_provenance,initial_values_only=True,means_fixed=False,mean_constraints_from_CR=False,
               simultaneous_floating_means=[f'Ymean_{proc}_{era}' for era in ERAS for proc in ['upsilon1s','upsilon2s']],
               independent_SR_CR_continua=True,local_polynomial_domains=True,discrete_categories=0,envelope_penalty=0,
               nuisance_constraints=[a.GetName() for a in constraints],original_signal_nuisances=True,signal_rate_reference=[6,12],
               signal_rate_extrapolated=float(mass)>11,workspace_sha256=sha(output))
    output.with_suffix('.audit.json').write_text(json.dumps(audit,indent=2)+'\n');return output

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('side',choices=['low','high']);parser.add_argument('mass');parser.add_argument('--output');args=parser.parse_args()
    print(build(args.side,args.mass,args.output),flush=True)

