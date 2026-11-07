"""Read-only checks of the new two-region models; no likelihood fits."""
import argparse,json,struct,hashlib
import ROOT
import numpy as np
from model import R,P,ERAS,sha,regions_for

def audit(side,mass,deep=False):
    mass=f'{float(mass):.1f}';path=R/side/f'm{mass}/workspace.root';before=sha(path)
    meta=json.loads(path.with_suffix('.audit.json').read_text());regions=regions_for(side)
    handle=ROOT.TFile.Open(str(path));w=handle.Get('w');sf=ROOT.TFile.Open(meta['source']);old=sf.Get('w')
    assert sha(meta['source'])==meta['source_sha256'] and before==meta['workspace_sha256']
    assert meta['degree']==(6 if side=='low' else 4) and meta['CR_degree']==4
    means=[v for v in w.allVars() if v.GetName().startswith('Ymean_')]
    assert len(means)==8 and all(not v.isConstant() for v in means)
    assert not w.genobj('discreteParams') and not any(p.InheritsFrom('RooMultiPdf') for p in w.allPdfs())
    assert w.cat('CMS_channel').numTypes()==8
    mc=w.obj('ModelConfig');assert mc.GetGlobalObservables().getSize()==10 and mc.GetNuisanceParameters().getSize()==10
    constraints=list(mc.GetPdf()._extraConstraints);assert len(constraints)==18
    assert len(set(p.GetName() for p in constraints))==18
    assert {p.GetName() for p in constraints if not p.GetName().startswith('positive_continuum_')}==set(meta['nuisance_constraints'])
    assert w.data('data_obs').numEntries()==sum(4*r[2] for r in regions.values())
    basis_errors=[];acceptance_errors=[];dynamic_errors=[];rate_errors=[];c0_errors=[];c0_components=[]
    old.var('r').setVal(1);w.var('r').setVal(1)
    for i,era in enumerate(ERAS,1):
        for reg,(lo,hi,nb) in regions.items():
            x=w.var('mass_'+reg);assert abs(x.getMin()-lo)<1e-12 and abs(x.getMax()-hi)<1e-12 and x.getBins()==nb
            deg=meta['degree'] if reg=='SR' else 4;suffix='' if reg=='SR' else '_CR'
            c=[1.]+[w.var(f'local_c{k}_{era}{suffix}').getVal() for k in range(1,deg+1)]
            native=w.pdf(f'local_continuum_{era}_{reg}');norm=ROOT.RooArgSet(x)
            assert not native.dependsOn(w.var('mass_'+('CR' if reg=='SR' else 'SR')))
            integral=np.polynomial.Chebyshev(c).integ();total=(hi-lo)/2*(integral(1)-integral(-1))
            for frac in [.13,.53,.89]:
                xx=lo+frac*(hi-lo);x.setVal(xx);expected=np.polynomial.chebyshev.chebval(2*frac-1,c)/total
                error=abs(native.getVal(norm)/expected-1);basis_errors.append(error);assert error<1e-9
        for proc in ['Zd','upsilon1s','upsilon2s']:
            original=old.pdf(('shapeSig_' if proc=='Zd' else 'shapeBkg_')+proc+f'_ch{i}')
            if proc!='Zd':
                mu=w.var(f'Ymean_{proc}_{era}');assert abs(mu.getMin()-(original.x0().getVal()-.2))<1e-12 and abs(mu.getMax()-(original.x0().getVal()+.2))<1e-12
                sg=w.function(f'Ysigma_{proc}_{era}');assert sg.isConstant() and sg.getVal()==original.sigmaL().getVal()
                tails=[w.function(f'Ytail_{proc}_{era}_{k}') for k in range(4)]
                for a,b in zip(tails,[original.alphaL(),original.nL(),original.alphaR(),original.nR()]):assert a.isConstant() and a.getVal()==b.getVal()
                assert not w.var(f'CMS_EXO25020_bkgScale_{proc}_{era}').isConstant()
            pdf=w.pdf(f'shape_{proc}_{era}_SR');x=ROOT.RooRealVar(f'audit_mass_{proc}_{era}','',9.,6.,12.)
            full=ROOT.RooCrystalBall(f'audit_full_{proc}_{era}','',x,pdf.x0(),pdf.sigmaL(),pdf.sigmaR(),pdf.alphaL(),pdf.nL(),pdf.alphaR(),pdf.nR())
            variations=[None]
            if deep:variations += [-.03,.03] if proc!='Zd' else []
            central=pdf.x0().getVal()
            for delta in variations:
                if delta is not None:mu.setVal(min(mu.getMax()-.001,max(mu.getMin()+.001,central+delta)))
                for reg,(lo,hi,_) in regions.items():
                    x.setRange(reg,lo,hi);native=full.createIntegral(ROOT.RooArgSet(x),ROOT.RooFit.NormSet(ROOT.RooArgSet(x)),ROOT.RooFit.Range(reg)).getVal()
                    analytic=w.function(f'local_acceptance_{proc}_{era}_{reg}').getVal();error=abs(native-analytic)
                    assert error<1e-8;acceptance_errors.append(error)
                    if delta is not None:dynamic_errors.append(error)
                    assert w.pdf(f'shape_{proc}_{era}_{reg}').dependsOn(pdf.x0())
                    if delta is None:
                        record=next(q for q in meta['eras'] if q['era']==era)['processes'][proc][reg]
                        actual=w.function(f'local_rate_{proc}_{era}_{reg}').getVal();expected=record['full_rate_r1']*analytic
                        err=abs(actual-expected)/max(1.,abs(expected));rate_errors.append(err);assert err<1e-10
                        if proc=='Zd':
                            source_rate=old.function(f'extrapolated_rate_Zd_{era}') if float(mass)>11 else old.function(f'n_exp_binch{i}_proc_Zd')
                            assert abs(source_rate.getVal()-record['full_rate_r1'])<1e-8
            if proc!='Zd':mu.setVal(central)
    def floating_parameters(ww):
        return {v.GetName() for v in ww.obj('ModelConfig').GetPdf().getParameters(ww.obj('ModelConfig').GetObservables()) if v.InheritsFrom('RooRealVar') and not v.isConstant()}
    new_params=floating_parameters(w)
    expected_params={'r'}|{v.GetName() for v in mc.GetNuisanceParameters() if not v.isConstant()}
    for era in ERAS:
        expected_params.add('local_yield_'+era+'_CR')
        expected_params.update(f'local_c{k}_{era}' for k in range(1,meta['degree']+1))
        expected_params.update(f'local_c{k}_{era}_CR' for k in range(1,5))
        for proc in ['upsilon1s','upsilon2s']:
            expected_params.update([f'Ymean_{proc}_{era}',f'CMS_EXO25020_bkgScale_{proc}_{era}'])
    assert new_params==expected_params,dict(unexpected=sorted(new_params-expected_params),missing=sorted(expected_params-new_params))
    # The portable workflow does not build the historical model without C0.
    # Keep that extra regression check when available, and report its absence.
    previous_path=P/f'runs/20261031_separate_sr_cheb_limits/{side}/m{mass}/workspace.root'
    predecessor_checked=False
    if previous_path.is_file():
        previous_handle=ROOT.TFile.Open(str(previous_path));previous=previous_handle.Get('w')
        old_params=floating_parameters(previous)
        assert old_params-new_params=={'local_yield_'+era for era in ERAS}
        assert not new_params-old_params
        previous_handle.Close();predecessor_checked=True
    # Check actual normalized RooFit densities, not just the algebraic ratio.
    saved={v.GetName():v.getVal() for v in w.allVars()}
    for variant in (['nominal','normalizations','coefficients_means_nuisances'] if deep else ['nominal']):
        for name,value in saved.items():w.var(name).setVal(value)
        if variant=='normalizations':
            w.var('r').setVal(2.3)
            for era in ERAS:
                v=w.var('local_yield_'+era+'_CR');v.setVal(v.getVal()*1.13)
                for proc in ['upsilon1s','upsilon2s']:
                    v=w.var(f'CMS_EXO25020_bkgScale_{proc}_{era}');v.setVal(min(v.getMax()-.001,v.getVal()*1.05+.001))
        elif variant=='coefficients_means_nuisances':
            for v in w.allVars():
                name=v.GetName()
                if name.startswith('local_c'):v.setVal(v.getVal()+.002)
                elif name.startswith('Ymean_'):v.setVal(min(v.getMax()-.001,v.getVal()+.007))
            for v in w.obj('ModelConfig').GetNuisanceParameters():
                if not v.isConstant():v.setVal(min(v.getMax()-.001,max(v.getMin()+.001,v.getVal()+.15)))
        for era in ERAS:
            assert not w.var('local_yield_'+era) and w.function('local_yield_'+era).InheritsFrom('RooFormulaVar')
            assert not w.var('local_yield_'+era+'_CR').isConstant()
            values={}
            for reg,(lo,hi,_) in regions.items():
                x=w.var('mass_'+reg);x.setVal(meta['continuity_boundary']);norm=ROOT.RooArgSet(x)
                name='local_yield_'+era+('' if reg=='SR' else '_CR');rate=w.function(name) or w.var(name)
                values[reg]={'dy':rate.getVal()*w.pdf(f'local_continuum_{era}_{reg}').getVal(norm)}
                for proc in ['Zd','upsilon1s','upsilon2s']:
                    values[reg][proc]=w.function(f'local_rate_{proc}_{era}_{reg}').getVal()*w.pdf(f'shape_{proc}_{era}_{reg}').getVal(norm)
                values[reg]['total_b']=sum(v for k,v in values[reg].items() if k!='Zd')
                values[reg]['total_s']=values[reg]['total_b']+values[reg]['Zd']
                coeff=[1.]+[w.var(f'local_c{k}_{era}'+('' if reg=='SR' else '_CR')).getVal() for k in range(1,(meta['degree'] if reg=='SR' else 4)+1)]
                integral=np.polynomial.Chebyshev(coeff).integ();expected=(hi-lo)/2*(integral(1)-integral(-1))
                assert abs(w.function(f'c0_integral_{era}_{reg}').getVal()-expected)<1e-10
            for component in values['SR']:
                left=values['SR'][component];right=values['CR'][component];error=abs(left-right)/max(1.,abs(left),abs(right))
                c0_errors.append(error);assert error<1e-8,(side,mass,era,variant,component,left,right,error)
                if deep:c0_components.append(dict(era=era,variant=variant,component=component,SR_density=left,CR_density=right,error=error))
    for name,value in saved.items():w.var(name).setVal(value)
    endpoint_guard_checked=False
    if deep:
        era=ERAS[0];edge=1 if side=='low' else -1;degree=meta['degree']
        others=1+sum(w.var(f'local_c{k}_{era}').getVal()*edge**k for k in range(2,degree+1))
        w.var(f'local_c1_{era}').setVal(-others/edge)
        assert abs(w.function(f'c0_endpoint_{era}_SR').getVal())<1e-9
        assert w.pdf('positive_continuum_'+era).getVal()==0.
        assert w.function(f'c0_yield_ratio_{era}').getVal()==0.
        assert np.isfinite(w.function('local_yield_'+era).getVal())
        endpoint_guard_checked=True
        for name,value in saved.items():w.var(name).setVal(value)
    def digest(data,new):
        entries=[]
        for j in range(data.numEntries()):
            row=data.get(j);ch=row.find('CMS_channel').getLabel()
            if new:
                reg,ch=ch.split('_',1);xx=row.getRealValue('mass_'+reg)
            else:
                xx=row.getRealValue('mass')
                if not any(lo<=xx<hi for lo,hi,_ in regions.values()):continue
            entries.append((ch,xx,data.weight()))
        h=hashlib.sha256()
        for ch,xx,weight in sorted(entries):h.update(ch.encode()+struct.pack('dd',xx,weight))
        return h.hexdigest()
    assert digest(w.data('data_obs'),True)==digest(old.data('data_obs'),False)
    handle.Close();sf.Close();assert sha(path)==before
    return dict(side=side,mass=float(mass),passed=True,workspace_sha256=before,source_unchanged=True,data_exact=True,means_floating=8,widths_tails_fixed=True,
                bins=meta['entries'],events=meta['events'],discrete_categories=0,envelope_penalty=0,nuisance_constraints_once=10,
                local_basis_max_relative_error=max(basis_errors),DCB_acceptance_max_error=max(acceptance_errors),dynamic_mu_checks=len(dynamic_errors),
                dynamic_mu_max_error=max(dynamic_errors,default=0.),rate_max_relative_error=max(rate_errors),
                continuity='C0 physical event density',continuity_max_scaled_error=max(c0_errors),continuity_checks=len(c0_errors),
                floating_parameters_checked=len(new_params),historical_predecessor_comparison=predecessor_checked,
                continuum_yield_parameter_reduction=len(ERAS),continuity_components=c0_components,
                forbidden_zero_endpoint_guard_checked=endpoint_guard_checked)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--all',action='store_true');args=p.parse_args()
    points=[(q['side'],q['mass']) for q in json.loads((R/'workspace_manifest.json').read_text())['workspaces']] if args.all else [('low',7.5),('high',11.0)]
    rows=[audit(side,mass,deep=(side,mass) in [('low',7.5),('high',11.0)]) for side,mass in points]
    result=dict(passed=True,workspaces=len(rows),fits_launched=0,rows=rows,builder_sha256=sha(P/'analysis/separate_sr_c0_limits/model.py'))
    output=R/('all_mass_model_audit.json' if args.all else 'pilot_model_audit.json');output.write_text(json.dumps(result,indent=2)+'\n');print(output,flush=True)
