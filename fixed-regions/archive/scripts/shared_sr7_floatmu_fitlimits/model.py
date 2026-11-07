"""Shared physical-mass polynomial across two SRs, independent CR continuum.

Native RooFit PDFs with a fixed polynomial reference range [6,12]. The SR
continuum normalization is its count in SRL union SRH; yields are linked by
analytic integrals. Original signal graph and normalization uncertainties are
preserved. Both Upsilon means in every era float jointly across all regions;
previous CR-only estimates are starting values only, never a constraint.
All output is separate from earlier studies. No fits are launched here.
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
R=P/'runs/20261031_shared_sr7_floatmu_fitlimits'
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

def _dcb_acceptance(w,proc,era,region,pdf):
    lo,hi,_=REGIONS[region]
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

def build(mass,degree=7,family='envelope',output=None,freeze_means=False,seed=None):
    if freeze_means:
        raise ValueError('This campaign requires all eight Upsilon means floating')
    mass=f'{float(mass):.1f}'
    output=Path(output or R/'observed'/f'm{mass}'/'joint/workspace_seeded.root')
    output.parent.mkdir(parents=True,exist_ok=True)
    src=source_path(mass);sf=ROOT.TFile.Open(str(src));old=sf.Get('w')
    assert old and family in ['envelope','cheb','bern']
    w=ROOT.RooWorkspace('w');keep=[]
    cr_result=json.loads(CR_REFERENCE.read_text())
    cr_seed={p['name']:p['value'] for p in cr_result['selected']['parameters']}
    means={p['name']:p['value'] for p in cr_result['means']}
    # Reuse prior CR-only central values as numerical starting values only.
    # Means float in the joint fit; no uncertainty/constraint is imported from
    # that earlier fit and no new CR pre-fit is performed by this builder.
    for reg,(lo,hi,nb) in REGIONS.items():
        x=ROOT.RooRealVar('mass_'+reg,'mass [GeV]',(lo+hi)/2,lo,hi)
        x.setBins(nb);x.setRange('global_SR',6,12);_import(w,x)
    cat=ROOT.RooCategory('CMS_channel','CMS_channel')
    for reg in REGIONS:
        for i in range(1,5):cat.defineType(f'{reg}_ch{i}')
    _import(w,cat);cat=w.cat('CMS_channel')
    old.var('r').setVal(1)
    nuis=list(old.obj('ModelConfig').GetNuisanceParameters())
    globs=list(old.obj('ModelConfig').GetGlobalObservables())
    constraints=list(old.obj('ModelConfig').GetPdf().getAllConstraints(old.obj('ModelConfig').GetObservables(),old.obj('ModelConfig').GetNuisanceParameters(),False))
    assert len(constraints)==10
    for con in constraints:_import(w,con)
    counts={f'{reg}_ch{i}':0. for reg in REGIONS for i in range(1,5)}
    original=old.data('data_obs')
    rows=[]
    for j in range(original.numEntries()):
        row=original.get(j);x=row.getRealValue('mass');ch=row.find('CMS_channel').getLabel();weight=original.weight()
        reg=next(k for k,(lo,hi,_) in REGIONS.items() if lo<=x<hi)
        counts[reg+'_'+ch]+=weight;rows.append((reg,ch,x,weight))
    positivity=[];integrals={};era_audits=[]
    for i,era in enumerate(ERAS,1):
        ch=f'ch{i}';seed_c=[float(old.var(f't{k-1}_{era}_envelope').getVal()) if k<=7 else 0. for k in range(1,degree+1)]
        shrink=1.
        while _minimum_cheb([1]+seed_c)<=1e-5:
            seed_c=[v*.98 for v in seed_c];shrink*=.98
        coeff=[]
        for j,v in enumerate(seed_c,1):coeff.append(_import(w,ROOT.RooRealVar(f'local_c{j}_{era}','',v,-2 if j==1 else -1,2 if j==1 else 1)))
        positivity.append(_positive_constraint(w,'positive_continuum_'+era,coeff))
        bs=[]
        if family in ['envelope','bern']:
            power=Chebyshev([1]+seed_c).convert(kind=Polynomial)(Polynomial([-1,2])).coef
            raw=np.array([sum(power[k]*math.comb(j,k)/math.comb(degree,k) for k in range(j+1)) for j in range(degree+1)])
            raw=np.maximum(raw,1e-5);raw/=raw[0]
            bs=[_import(w,ROOT.RooConstVar(f'local_b0_{era}','',1.))]
            bs.extend(_import(w,ROOT.RooRealVar(f'local_b{j}_{era}','',min(999.,max(1e-6,float(raw[j]))),1e-6,1000.)) for j in range(1,degree+1))
        if family=='envelope':
            idx=ROOT.RooCategory('pdfindex_'+era,'');_import(w,idx);idx=w.cat(idx.GetName())
        I={}
        for reg in ['SRL','SRH']:
            lo,hi,_=REGIONS[reg]
            ci=_poly_integral_node(w,f'integral_cheb_{era}_{reg}',coeff,_poly_integrals(degree,'cheb',lo,hi),True)
            if family in ['bern','envelope']:
                bi=_poly_integral_node(w,f'integral_bern_{era}_{reg}',bs,_poly_integrals(degree,'bern',lo,hi))
            if family=='envelope':I[reg]=_import(w,ROOT.RooFormulaVar(f'integral_selected_{era}_{reg}','@0==0?@1:@2',ROOT.RooArgList(idx,ci,bi)))
            else:I[reg]=ci if family=='cheb' else bi
        integrals[era]=I
        nSR=counts['SRL_'+ch]+counts['SRH_'+ch]
        yn=_import(w,ROOT.RooRealVar('local_yield_'+era,'continuum yield in both SRs',nSR,0,max(10,5*nSR)))
        crcoeff=[]
        for j in range(1,5):crcoeff.append(_import(w,ROOT.RooRealVar(f'local_c{j}_{era}_CR','',cr_seed[f'local_c{j}_{era}'],-2 if j==1 else -1,2 if j==1 else 1)))
        positivity.append(_positive_constraint(w,'positive_continuum_'+era+'_CR',crcoeff))
        nyCR=_import(w,ROOT.RooRealVar('local_yield_'+era+'_CR','CR continuum yield',cr_seed['local_yield_'+era],0,5*counts['CR_'+ch]))
        # Import original rates; only source mass>11 has explicitly extrapolated signal rate.
        fullrates={}
        for proc in ['Zd','upsilon1s','upsilon2s']:
            rate=old.function(f'extrapolated_rate_Zd_{era}') if proc=='Zd' and float(mass)>11 else old.function(f'n_exp_bin{ch}_proc_{proc}')
            fullrates[proc]=_import(w,rate)
            if proc!='Zd':
                scale=w.var(f'CMS_EXO25020_bkgScale_{proc}_{era}');scale.setVal(cr_seed[scale.GetName()]);scale.setConstant(False)
                original_pdf=old.pdf(f'shapeBkg_{proc}_{ch}')
                mn=ROOT.RooRealVar(f'Ymean_{proc}_{era}','',means[f'Ymean_{proc}_{era}'],original_pdf.x0().getVal()-.2,original_pdf.x0().getVal()+.2);mn.setConstant(freeze_means);_import(w,mn)
                _import(w,ROOT.RooConstVar(f'Ysigma_{proc}_{era}','',original_pdf.sigmaL().getVal()))
                for j,arg in enumerate([original_pdf.alphaL(),original_pdf.nL(),original_pdf.alphaR(),original_pdf.nR()]):_import(w,ROOT.RooConstVar(f'Ytail_{proc}_{era}_{j}','',arg.getVal()))
        era_record={'era':era,'seed_cheb_shrink':shrink,'processes':{}}
        for reg,(lo,hi,nb) in REGIONS.items():
            x=w.var('mass_'+reg)
            if reg=='CR':
                cont=ROOT.RooChebychev(f'local_continuum_{era}_{reg}','',x,ROOT.RooArgList(*crcoeff));cont=_import(w,cont);dy_rate=nyCR
            else:
                cheb=ROOT.RooChebychev(f'local_chebyshev{degree}_{era}_{reg}','',x,ROOT.RooArgList(*coeff));cheb.selectNormalizationRange('global_SR',True);cheb=_import(w,cheb)
                if family in ['bern','envelope']:
                    bern=ROOT.RooBernstein(f'local_bernstein{degree}_{era}_{reg}','',x,ROOT.RooArgList(*bs));bern.selectNormalizationRange('global_SR',True);bern=_import(w,bern)
                if family=='envelope':
                    env_idx=idx if reg=='SRL' else ROOT.RooCategory(idx.GetName(),'')
                    keep.append(env_idx)
                    cont=ROOT.RooMultiPdf(f'local_envelope_{era}_{reg}','',env_idx,ROOT.RooArgList(cheb,bern));cont.setCorrectionFactor(.25);cont=_import(w,cont)
                else:cont=cheb if family=='cheb' else bern
                dy_rate=_import(w,ROOT.RooFormulaVar(f'local_rate_dy_{era}_{reg}',('@0*@1/(@1+@2)' if reg=='SRL' else '@0*@2/(@1+@2)'),ROOT.RooArgList(yn,I['SRL'],I['SRH'])))
            dy_rate.setStringAttribute('combine.channel',reg+'_'+ch);dy_rate.setStringAttribute('combine.process','dy')
            pdfs={};rates={}
            for proc in ['Zd','upsilon1s','upsilon2s']:
                if proc=='Zd':
                    op=old.pdf('shapeSig_Zd_'+ch);cu=ROOT.RooCustomizer(op,reg);cu.replaceArg(old.var('mass'),x);pdf=cu.build();pdf.SetName(f'shape_{proc}_{era}_{reg}');keep.extend([cu,pdf]);pdf=_import(w,pdf)
                else:
                    sg=w.function(f'Ysigma_{proc}_{era}');tails=[w.function(f'Ytail_{proc}_{era}_{j}') for j in range(4)]
                    pdf=_import(w,ROOT.RooCrystalBall(f'shape_{proc}_{era}_{reg}','',x,w.var(f'Ymean_{proc}_{era}'),sg,sg,*tails))
                acceptance=_dcb_acceptance(w,proc,era,reg,pdf)
                rate=_import(w,ROOT.RooFormulaVar(f'local_rate_{proc}_{era}_{reg}','@0*@1',ROOT.RooArgList(fullrates[proc],acceptance)))
                rate.setAttribute('combine.signal',proc=='Zd');rate.setStringAttribute('combine.channel',reg+'_'+ch);rate.setStringAttribute('combine.process',proc)
                pdfs[proc]=pdf;rates[proc]=rate
                era_record['processes'].setdefault(proc,{})[reg]=dict(acceptance=acceptance.getVal(),rate_r1=rate.getVal())
            for stage in ['s','b']:
                pp=[cont,pdfs['upsilon1s'],pdfs['upsilon2s']];rr=[dy_rate,rates['upsilon1s'],rates['upsilon2s']]
                if stage=='s':pp=[pdfs['Zd']]+pp;rr=[rates['Zd']]+rr
                channel=ROOT.RooAddPdf(f'pdf_{stage}_{reg}_{ch}','',ROOT.RooArgList(*pp),ROOT.RooArgList(*rr));_import(w,channel)
        ySR=sum(era_record['processes'][proc][reg]['rate_r1'] for proc in ['upsilon1s','upsilon2s'] for reg in ['SRL','SRH'])
        yn.setVal(max(1,nSR-ySR))
        for ix in ([0,1] if family=='envelope' else [0]):
            if family=='envelope':idx.setIndex(ix)
            assert abs(sum(w.function(f'local_rate_dy_{era}_{rr}').getVal() for rr in ['SRL','SRH'])/yn.getVal()-1)<1e-12
        if family=='envelope':idx.setIndex(0)
        for proc in ['Zd','upsilon1s','upsilon2s']:
            assert abs(sum(v['acceptance'] for v in era_record['processes'][proc].values())-1)<1e-10
        era_audits.append(era_record)
    for stage,mcname in [('s','ModelConfig'),('b','ModelConfig_bonly')]:
        sim=ROOT.RooSimultaneousOpt('joint_'+mcname,'',cat)
        for reg in REGIONS:
            for i in range(1,5):sim.addPdf(w.pdf(f'pdf_{stage}_{reg}_ch{i}'),f'{reg}_ch{i}')
        for con in constraints:sim.addExtraConstraint(w.pdf(con.GetName()))
        for con in positivity:sim.addExtraConstraint(w.pdf(con.GetName()))
        _import(w,sim)
        mc=ROOT.RooStats.ModelConfig(mcname,w);mc.SetPdf(w.pdf(sim.GetName()))
        mc.SetObservables(ROOT.RooArgSet(*[w.var('mass_'+reg) for reg in REGIONS],cat));mc.SetParametersOfInterest(ROOT.RooArgSet(w.var('r')))
        mc.SetNuisanceParameters(ROOT.RooArgSet(*[w.var(v.GetName()) for v in nuis]));mc.SetGlobalObservables(ROOT.RooArgSet(*[w.var(v.GetName()) for v in globs]))
        getattr(w,'import')(mc);keep.append(mc)
    weight=ROOT.RooRealVar('evt_weight','',0,1e9)
    obs=ROOT.RooArgSet(*[w.var('mass_'+reg) for reg in REGIONS],cat,weight)
    data=ROOT.RooDataSet('data_obs','',obs,ROOT.RooFit.WeightVar(weight))
    for reg,ch,x,n in rows:
        cat.setLabel(reg+'_'+ch);w.var('mass_'+reg).setVal(x);weight.setVal(n);data.add(obs,n)
    getattr(w,'import')(data)
    for reg,(lo,hi,_) in REGIONS.items():w.var('mass_'+reg).setVal((lo+hi)/2)
    if family=='envelope':
        discretes=ROOT.RooArgSet(*[w.cat('pdfindex_'+era) for era in ERAS]);discretes.setName('discreteParams');getattr(w,'import')(discretes,'discreteParams')
    if seed:
        if isinstance(seed,(str,Path)):
            fseed=ROOT.TFile.Open(str(seed));sw=fseed.Get('w')
            for v in sw.allVars():
                if w.var(v.GetName()) and v.GetName().startswith(('local_','Ymean_','CMS_EXO25020_bkgScale_')):w.var(v.GetName()).setVal(v.getVal())
            fseed.Close()
        else:
            for k,v in seed.items():
                if w.var(k):w.var(k).setVal(v)
    w.var('r').setVal(0);w.var('r').setRange(-1000,1000);w.var('r').setConstant(False)
    assert data.numEntries()==1800 and data.sumEntries()==original.sumEntries()
    w.writeToFile(str(output))
    audit=dict(mass=float(mass),degree=degree,family=family,source=str(src),source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),
               regions=REGIONS,counts=counts,entries=data.numEntries(),events=data.sumEntries(),
               CR_initial_values_source=str(CR_REFERENCE),CR_initial_values_source_sha256=hashlib.sha256(CR_REFERENCE.read_bytes()).hexdigest(),
               CR_initial_values_only=True,mean_constraints_from_CR=False,
               means_fixed=freeze_means,means=means,eras=era_audits,
               shared_SR_continuum_intensity=True,polynomial_reference=[6,12],
               discrete_categories=4 if family=='envelope' else 0,
               correction_factor_per_SR_pdf=.25 if family=='envelope' else 0,
               nuisance_constraints=[a.GetName() for a in constraints],original_signal_nuisances=True,
               simultaneous_floating_means=[f'Ymean_{proc}_{era}' for era in ERAS for proc in ['upsilon1s','upsilon2s']],
               signal_rate_extrapolated=float(mass)>11)
    output.with_suffix('.audit.json').write_text(json.dumps(audit,indent=2))
    sf.Close()
    return output

def make_cr_workspace(source,output,means_floating=True):
    """Select CR channels without changing their likelihood or counted events."""
    source=Path(source);output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    f=ROOT.TFile.Open(str(source));src=f.Get('w');w=ROOT.RooWorkspace('w')
    cat=ROOT.RooCategory('CMS_channel','CMS_channel')
    for i in range(1,5):cat.defineType('CR_ch'+str(i))
    keep=[]
    for mcname in ['ModelConfig','ModelConfig_bonly']:
        sm=src.obj(mcname);sim=ROOT.RooSimultaneousOpt('control_'+mcname,'',cat)
        for i in range(1,5):sim.addPdf(sm.GetPdf().getPdf('CR_ch'+str(i)),'CR_ch'+str(i))
        for con in sm.GetPdf()._extraConstraints:
            if not con.GetName().startswith('positive_continuum_') or con.GetName().endswith('_CR'):sim.addExtraConstraint(con)
        _import(w,sim);mc=ROOT.RooStats.ModelConfig(mcname,w);mc.SetPdf(w.pdf(sim.GetName()))
        mc.SetObservables(ROOT.RooArgSet(w.var('mass_CR'),w.cat('CMS_channel')))
        for getter,setter in [('GetParametersOfInterest','SetParametersOfInterest'),('GetNuisanceParameters','SetNuisanceParameters'),('GetGlobalObservables','SetGlobalObservables')]:
            getattr(mc,setter)(ROOT.RooArgSet(*[w.arg(a.GetName()) for a in getattr(sm,getter)() if w.arg(a.GetName())]))
        getattr(w,'import')(mc);keep.append(mc)
    weight=ROOT.RooRealVar('evt_weight','',0,1e9);obs=ROOT.RooArgSet(w.var('mass_CR'),w.cat('CMS_channel'),weight)
    data=ROOT.RooDataSet('data_obs','',obs,ROOT.RooFit.WeightVar(weight));d=src.data('data_obs')
    for j in range(d.numEntries()):
        row=d.get(j);ch=row.find('CMS_channel').getLabel()
        if not ch.startswith('CR_'):continue
        w.cat('CMS_channel').setLabel(ch);w.var('mass_CR').setVal(row.getRealValue('mass_CR'));weight.setVal(d.weight());data.add(obs,d.weight())
    getattr(w,'import')(data)
    for v in w.allVars():
        if v.GetName().startswith('Ymean_'):v.setConstant(not means_floating)
    w.var('mass_CR').setVal(9.8);w.var('r').setVal(0);w.writeToFile(str(output));f.Close();return output

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mass');parser.add_argument('--degree',type=int,default=7);parser.add_argument('--family',choices=['envelope','cheb','bern'],default='envelope');parser.add_argument('--output');parser.add_argument('--seed');parser.add_argument('--float-means',action='store_true',help='Compatibility flag: means always float in this campaign');parser.add_argument('--cr-output')
    args=parser.parse_args();p=build(args.mass,args.degree,args.family,args.output,False,args.seed)
    if args.cr_output:make_cr_workspace(p,args.cr_output)
    print(p,flush=True)
