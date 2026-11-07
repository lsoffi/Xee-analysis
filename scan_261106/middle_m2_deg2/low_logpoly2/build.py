"""Native Combine binned workspace: exact bin integrals, frozen CR shapes, C1."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
from array import array
import numpy as np
from numpy.polynomial import Chebyshev,Polynomial
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
BASE=Path(__file__).resolve().parents[2]
OUT=Path(os.environ.get('XEE_BUILD_OUT',BASE/'runs/20261105_sliding95_cheb3_c1_combine'))
ERAS=['2022','2022EE','2023','2023BPix']
RECORDS=json.loads((BASE/'records.json').read_text())
DATA=list(csv.DictReader((BASE/'bins.csv').open()))
DIRS={'omega_phi':'runs/20261105_three_cr_fixed_tails','upsilon':'runs/20261105_three_cr_fixed_tails','charmonium':'runs/20261105_charmonium_decreasing_exp_alpha50_n50'}
NAMES=['CMS_eff_e_id','CMS_eff_e_reco','CMS_eff_trigger_e','lumi_1','lumi_2','CMS_scale_e']+['CMS_EXO25020_signalModelSigmaNuisance_'+e for e in ERAS]


def minimum(a):
    roots=Chebyshev(a).deriv().roots();ts=[-1.,1.]+[float(x.real) for x in roots if abs(x.imag)<1e-9 and -1<x.real<1]
    return min(float(Chebyshev(a)(t)) for t in ts)


def build(mass):
    degree=int(os.environ.get("XEE_DEGREE","3"));size=degree+1
    assert degree in (2,3,4,5,6)
    if mass >= 9.093333333333334:raise ValueError('Above the upper scan boundary')
    rec=next(r for r in RECORDS if r['eligible'] and abs(r['mass']-mass)<1e-8)
    target=OUT/f'm{mass:.1f}';target.mkdir(parents=True,exist_ok=True)
    w=ROOT.RooWorkspace('w');keep=[];frozen={};hashes={};constraints=[];global_names=set();channels=[];era_meta=[]
    def imp(obj):
        
        if obj.InheritsFrom('RooStats::ModelConfig'):getattr(w,'import')(obj)
        else:getattr(w,'import')(obj,ROOT.RooFit.RecycleConflictNodes(),ROOT.RooFit.Silence())
        return w.arg(obj.GetName()) or w.obj(obj.GetName())
    def const(name,value):
        return w.arg(name) or imp(ROOT.RooConstVar(name,'',float(value)))
    def var(name,value,low,high,error=None):
        x=ROOT.RooRealVar(name,'',float(value),float(low),float(high));x.setError(float(error if error is not None else max(.01,abs(value)*.05)));x.setAttribute('flatParam');return imp(x)
    def form(name,expr,args):return w.arg(name) or imp(ROOT.RooFormulaVar(name,expr,ROOT.RooArgList(*args)))
    def add(name,args):return imp(ROOT.RooAddition(name,'',ROOT.RooArgList(*args)))
    def sha(path):hashes[str(path.relative_to(BASE))]=hashlib.sha256(path.read_bytes()).hexdigest()
    def primitive(prefix,edge,p):
        name=prefix+'_F_'+f'{edge:.12f}'.replace('-','m').replace('.','p')
        t=f'(({edge:.17g}-@0)/abs(@1))';al='abs(@2)';ar='abs(@4)';nl='@3';nr='@5';bl=f'({nl}/{al}-{al})';br=f'({nr}/{ar}-{ar})'
        left=f'pow({nl}/{al},{nl})*exp(-0.5*{al}*{al})*(pow({bl}-{t},1-{nl})-pow({bl}+{al},1-{nl}))/({nl}-1)'
        right=f'pow({nr}/{ar},{nr})*exp(-0.5*{ar}*{ar})*(pow({br}+{t},1-{nr})-pow({br}+{ar},1-{nr}))/(1-{nr})'
        expr=f'({t}<-{al})?(sqrt(TMath::Pi()/2)*TMath::Erf(-{al}/sqrt(2))+({left})):(({t}>{ar})?(sqrt(TMath::Pi()/2)*TMath::Erf({ar}/sqrt(2))+({right})):(sqrt(TMath::Pi()/2)*TMath::Erf({t}/sqrt(2))))'
        return form(name,expr,p)
    def integral(prefix,a,b,p):return form(prefix+'_I_'+f'{a:.12f}_{b:.12f}'.replace('-','m').replace('.','p'),'@1-@0',[primitive(prefix,a,p),primitive(prefix,b,p)])
    r=var('r',0.,-30,30,.1)
    # Channel definitions are fixed before introducing the common integer-bin observable.
    raw={};max_bins=0
    for era in ERAS:
        groups={}
        for group in rec['selected_CRs']:
            path=BASE/"cr"/group/era/'display_fit.json';sha(path);groups[group]=json.loads(path.read_text())
        srrows=[x for x in DATA if x['era']==era and any(float(x['bin_low_GeV'])>=a-1e-9 and float(x['bin_high_GeV'])<=b+1e-9 for a,b in rec['sr_segments'])]
        raw[era]=(groups,srrows)
        for group,q in groups.items():max_bins=max(max_bins,len(q['counts']))
        for a,b in rec['sr_segments']:max_bins=max(max_bins,sum(float(x['bin_low_GeV'])>=a-1e-9 and float(x['bin_high_GeV'])<=b+1e-9 for x in srrows))
    axis=ROOT.RooRealVar('CMS_th1x','bin index',0.,float(max_bins));axis.setBins(max_bins);axis=imp(axis)
    cat=ROOT.RooCategory('CMS_channel','CMS_channel');cat=imp(cat)
    zero=const('zero',0.)
    def channel(label,era,region,edges,counts,processes):
        cat.defineType(label);nb=len(counts);template=ROOT.TH1D('template_'+label,'',max_bins,0,max_bins)
        for i,n in enumerate(counts,1):template.SetBinContent(i,n);template.SetBinError(i,math.sqrt(n))
        plist=[];clist=[];bplist=[];bclist=[];meta=[]
        for proc,is_signal,binpars,coef in processes:
            assert len(binpars)==nb
            pars=ROOT.RooArgList(*binpars,*([zero]*(max_bins-nb)))
            pdf=imp(ROOT.RooParametricHist(('shapeSig_' if is_signal else 'shapeBkg_')+proc+'_'+label,'',axis,pars,template))
            rate=form('n_exp_bin'+label+'_proc_'+proc,'@0',[coef]);rate.setAttribute('combine.signal',is_signal);rate.setStringAttribute('combine.channel',label);rate.setStringAttribute('combine.process',proc)
            plist.append(pdf);clist.append(rate)
            if not is_signal:bplist.append(pdf);bclist.append(rate)
            meta.append(dict(name=proc,signal=is_signal,bin_nodes=[x.GetName() for x in binpars],coefficient=rate.GetName(),shape=pdf.GetName(),shape_sum=sum(x.getVal() for x in binpars)))
        imp(ROOT.RooAddPdf('pdf_s_'+label,'',ROOT.RooArgList(*plist),ROOT.RooArgList(*clist)))
        imp(ROOT.RooAddPdf('pdf_b_'+label,'',ROOT.RooArgList(*bplist),ROOT.RooArgList(*bclist)))
        channels.append(dict(label=label,era=era,region=region,edges=edges,counts=counts,processes=meta,padded_bins=max_bins-nb));keep.append(template)
    for era in ERAS:
        groups,srrows=raw[era];a=min(float(x['bin_low_GeV']) for x in srrows);b=max(float(x['bin_high_GeV']) for x in srrows);width=b-a
        norms={};resonances=[];source_meta=[]
        for group,q in groups.items():
            lo,hi=q['edges'][0],q['edges'][-1];n=sum(q['counts']);bg=var('N_cont_'+group+'_'+era,math.exp(q['fit']['values'][0]),1e-5,5*n)
            off=7+q['tail_count'];slope=q['fit']['values'][off]/(hi-lo) if q['degree']==-1 else 0.
            unit_norm=math.expm1(slope*(hi-lo))/slope if abs(slope)>1e-12 else hi-lo
            norms[group]=dict(node=bg,lo=lo,hi=hi,slope=slope,integral=unit_norm)
            frozen['continuum_slope_'+group+'_'+era]=slope;const('continuum_slope_'+group+'_'+era,slope)
            for item in q['physical_parameters']:
                name=item['name'];p=item['fitted'];mid=sum(v['fitted']['mean'] for v in q['physical_parameters'])/2
                left=(item==q['physical_parameters'][0]);mu=var('mu_'+name+'_'+era,p['mean'],lo if left else mid,mid if left else hi,p['sigma']*.05)
                tail=[]
                for key in ['sigma','alphaL','nL','alphaR','nR']:
                    node=const(name+'_'+key+'_'+era,p[key]);tail.append(node);frozen[node.GetName()]=p[key]
                pp=[mu]+tail;den=integral(name+'_'+era,lo,hi,pp)
                yn=var('N_'+name+'_'+era,item['events_in_CR'],0,5*n)
                resonances.append(dict(group=group,name=name,params=pp,den=den,norm=yn))
            source_meta.append(dict(group=group,range=[lo,hi],continuum_norm=bg.GetName(),slope=slope))
        # C0/C1 are imposed only at physically shared SR/selected-CR edges.
        contacts=[]
        for group,d in norms.items():
            for edge in [d['lo'],d['hi']]:
                if any(abs(edge-lo)<1e-9 or abs(edge-hi)<1e-9 for lo,hi in rec['sr_segments']):contacts.append((group,edge))
        assert len(contacts)<=1
        rank=len(contacts); targets=[]; C=np.empty((0,3))
        anchor=0.
        if contacts:
            g,edge=contacts[0];d=norms[g]
            anchor=2*(edge-a)/width-1
            density=math.exp(d['slope']*(edge-d['lo']))/d['integral']
            targets=[(g,density)]
            amplitude=form('logpoly_anchor_'+era,f'{density:.17g}*@0',[d['node']])
        else:
            scale=sum(float(z['events']) for z in srrows)/width
            amplitude=var('logpoly_density_'+era,scale,scale*1e-6,scale*100)
        shape=[var(f'logpoly_{era}_c{j}',0,-8,8,.02) for j in (1,2)]
        coeff=[amplitude]+shape
        free=shape if contacts else coeff
        nodes,weights=np.polynomial.legendre.leggauss(8)
        # Original signal functions and their original ten Gaussian constraints.
        region=0 if mass<2.4 else 1 if mass<7 else 2
        if os.environ.get('XEE_SIGNAL_REGION'):region=int(os.environ['XEE_SIGNAL_REGION'])
        src=BASE/f'signal/cards_region{region}_data_envelope_allCorrections_binned/ee/{mass:.1f}/Xee_ee_4_{era}.root';sha(src)
        sf=ROOT.TFile.Open(str(src),'READ');sw=sf.Get('w');op=sw.pdf('shapeSig_Zd_Xee_ee_4_'+era)
        sp=[imp(v) for v in [{v.GetName().split('_')[0]:v for v in op.servers()}[k] for k in ('mean','sigma','alphaL','nL','alphaR','nR')]]
        rate=imp(sw.function('n_exp_binXee_ee_4_'+era+'_proc_Zd'))
        original=[sw.var('mass').getMin(),sw.var('mass').getMax()];sigden=integral('signal_'+era,*original,sp)
        for name in NAMES:
            gaussian=sw.pdf(name+'_Pdf')
            if gaussian:
                imp(gaussian)
                if gaussian.GetName() not in constraints:constraints.append(gaussian.GetName())
        for v in sw.obj('ModelConfig').GetGlobalObservables():global_names.add(v.GetName())
        # Build original-mass bins into a common integer index: exact Poisson means.
        definitions=[]
        for group,q in groups.items():definitions.append(('CR_'+group,[(lo,hi,n) for lo,hi,n in zip(q['edges'],q['edges'][1:],q['counts'])],group))
        for j,(lo,hi) in enumerate(rec['sr_segments']):definitions.append(('SR'+str(j),[(float(z['bin_low_GeV']),float(z['bin_high_GeV']),float(z['events'])) for z in srrows if float(z['bin_low_GeV'])>=lo-1e-9 and float(z['bin_high_GeV'])<=hi+1e-9],None))
        for reg,rows,group in definitions:
            label=reg+'_'+era;edges=[rows[0][0]]+[z[1] for z in rows];counts=[z[2] for z in rows];processes=[]
            if group:
                d=norms[group];parts=[]
                for k,(lo,hi,n) in enumerate(rows):
                    unit=((math.exp(d['slope']*(lo-d['lo']))*math.expm1(d['slope']*(hi-lo))/d['slope']) if abs(d['slope'])>1e-12 else hi-lo)/d['integral']
                    parts.append(const('contfrac_'+label+'_'+str(k),unit))
                processes.append(('continuum',False,parts,d['node']))
            else:
                parts=[]
                for k,(lo,hi,n) in enumerate(rows):
                    terms=[]
                    for u,wgt in zip(nodes,weights):
                        t=2*((lo+hi)/2+(hi-lo)*u/2-a)/width-1
                        exponent='+'.join(f'({t**j-anchor**j:.17g})*@{j}' for j in (1,2))
                        terms.append(f'({wgt*(hi-lo)/2:.17g})*exp({exponent})')
                    parts.append(form('contbin_'+label+'_'+str(k),'@0*('+'+'.join(terms)+')',coeff))
                processes.append(('continuum',False,parts,add('Ncont_'+label,parts)))
            for res in resonances:
                if group and res['group']!=group:continue
                fractions=[form('frac_'+res['name']+'_'+label+'_'+str(k),'@0/@1',[integral(res['name']+'_'+era,lo,hi,res['params']),res['den']]) for k,(lo,hi,n) in enumerate(rows)]
                acceptance=add('accept_'+res['name']+'_'+label,fractions);coef=form('N_'+res['name']+'_'+label,'@0*@1',[res['norm'],acceptance]);processes.append((res['name'],False,fractions,coef))
            fractions=[integral('signal_'+era,lo,hi,sp) for lo,hi,n in rows];totalI=add('signalI_'+label,fractions);coef=form('Nsignal_'+label,'@0*@1/@2',[rate,totalI,sigden]);processes.append(('Zd',True,fractions,coef))
            channel(label,era,reg,edges,counts,processes)
        era_meta.append(dict(era=era,SR_reference=[a,b],continuum_coefficients=[x.GetName() for x in coeff],C1_contacts=[dict(group=g,edge=x) for g,x in contacts],constraint_matrix=C.tolist(),constraint_targets=[dict(group=g,unit_density_or_derivative=v,norm=norms[g]['node'].GetName()) for g,v in targets],rank=rank,free_SR_parameters=len(free),CRs=source_meta,signal_source=str(src.relative_to(BASE)),signal_reference=original))
        keep.extend([sf,sw])
    for name in global_names:w.var(name).setConstant(True)
    for name in NAMES:
        if w.var(name):w.var(name).setVal(0.);w.var(name).setConstant(False);w.var(name).setError(1.)
    r.setRange(-30,30);r.setVal(0.);r.setConstant(False)
    nuisance=ROOT.RooArgSet(*[v for v in w.allVars() if not v.isConstant() and v.GetName() not in ['r','CMS_th1x']])
    assert all(w.function(name) and w.function(name).isConstant() for name in frozen)
    for stage,mcname in [('s','ModelConfig'),('b','ModelConfig_bonly')]:
        sim=ROOT.RooSimultaneousOpt('model_'+stage,'',cat)
        for ch in channels:sim.addPdf(w.pdf('pdf_'+stage+'_'+ch['label']),ch['label'])
        for name in constraints:sim.addExtraConstraint(w.pdf(name))
        sim=imp(sim);mc=ROOT.RooStats.ModelConfig(mcname,w);mc.SetPdf(sim);mc.SetObservables(ROOT.RooArgSet(axis,cat));mc.SetParametersOfInterest(ROOT.RooArgSet(r));mc.SetNuisanceParameters(nuisance);mc.SetGlobalObservables(ROOT.RooArgSet(*[w.var(n) for n in global_names]));imp(mc)
    histmap=ROOT.std.map('string','TH1*')()
    for ch in channels:
        h=ROOT.TH1D('data_'+ch['label'],'',max_bins,0,max_bins)
        for i,n in enumerate(ch['counts'],1):h.SetBinContent(i,n);h.SetBinError(i,math.sqrt(n))
        histmap[ch['label']]=h;keep.append(h)
    data=ROOT.RooDataHist('data_obs','original binned counts',ROOT.RooArgList(axis),cat,histmap);imp(data)
    expected=sum(sum(ch['counts']) for ch in channels);assert abs(data.sumEntries()-expected)<1e-8
    path=target/'workspace.root';w.writeToFile(str(path))
    for name,h in hashes.items():assert hashlib.sha256((BASE/name).read_bytes()).hexdigest()==h
    metadata=dict(degree=degree,mass=mass,geometry=rec,channels=channels,eras=era_meta,frozen_shape_values=frozen,original_signal_nuisances=NAMES,Gaussian_auxiliary_constraints=10,
                  constraints=constraints,free_parameters=[v.GetName() for v in nuisance],events=expected,integer_axis_bins=max_bins,
                  SR_model="exp(quadratic in scaled mass)",SR_quadrature_order=8,exact_bin_integration=False,signal_in_SR_and_all_selected_CRs=True,C1_only_actual_shared_boundaries=not bool(os.environ.get("XEE_C0")),physical_continuum_density=True,
                  source_hashes=hashes,source_hashes_unchanged=True,workspace_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    (target/'model.json').write_text(json.dumps(metadata,indent=2))
    print('BUILT',mass,len(channels),'channels',len(nuisance)+1,'parameters',[(q['rank'],q['free_SR_parameters']) for q in era_meta],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--mass',type=float,required=True);args=p.parse_args();build(args.mass)
