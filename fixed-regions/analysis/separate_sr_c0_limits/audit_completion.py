"""Independent read-only completion/numerical audit. Never fit or alter results."""
from pathlib import Path
import collections,datetime,hashlib,json,math,re
import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
P=Path(__file__).resolve().parents[2]
R=P/'runs/20261031_separate_sr_c0_limits'
ERAS=['2022','2022EE','2023','2023BPix']
QUANTILES={'observed':-1.,'expected_2p5':.025,'expected_16':.16,'expected_50':.5,'expected_84':.84,'expected_97p5':.975}

def read(path):return json.loads(Path(path).read_text())
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def close(a,b,rtol=1e-10,atol=1e-10):return math.isfinite(a) and math.isfinite(b) and abs(a-b)<=atol+rtol*max(abs(a),abs(b))
def native_file(path):
    f=ROOT.TFile.Open(str(path));assert f and not f.IsZombie(),str(path);return f

def inspect_significance(rows):
    """Audit newly requested native significances without modifying any result."""
    out=[];missing=[];receipts=[]
    for point in rows:
        side=point['side'];mass=point['mass'];D=R/side/f'm{mass:.1f}';J=D/'significance'
        path=D/'significance_results.json'
        if not path.exists():missing.append(dict(side=side,mass=mass));continue
        result=read(path);receipt=read(J/'combine.receipt.json');assert result['receipt']==receipt
        receipts.append(dict(side=side,mass=mass,method='Significance',**receipt))
        command=read(J/'combine.command.json');protocol=read(J/'protocol.json')
        assert command[command.index('-M')+1]=='Significance'
        assert command[command.index('-d')+1]==str(D/'workspace.root')
        assert '--freezeParameters' not in command and '--cminRunAllDiscreteCombinations' not in command
        assert protocol['parameter_freezes']==[] and protocol['attempts']==1
        assert protocol['workspace_sha256']==sha(D/'workspace.root')
        values=[];native=result.get('native_output')
        if native and Path(native).exists():
            f=native_file(native);tree=f.Get('limit')
            if tree:values=[float(q.limit) for q in tree]
            f.Close()
        if result['Z'] is not None:
            assert len(values)==1 and values[0]==result['Z']
            assert close(result['p0'],.5*math.erfc(values[0]/math.sqrt(2.)),rtol=1e-12,atol=0.)
            assert math.isfinite(result['log10_p0'])
            assert result['p0_underflow']==(result['p0']==0.)
            if result['p0_underflow']:assert 'p0_floating_point_underflow' in result['flags']
            else:assert close(result['log10_p0'],math.log10(result['p0']),rtol=1e-12,atol=1e-12)
        b=point['fit']['b'];s=point['fit']['s'];ref_valid=b['point_valid'] and s['point_valid']
        refZ=math.sqrt(max(0.,b['exact_deviance']-s['exact_deviance'])) if s['r']>0 else 0.
        assert result['reference_FD']['point_valid']==ref_valid
        assert close(refZ,result['reference_FD']['Z'],atol=1e-7)
        if not ref_valid:assert 'FD_reference_point_invalid' in result['flags']
        if b['covQual']!=3 or s['covQual']!=3:assert 'FD_covariance_warning' in result['flags']
        if ref_valid and result['Z'] is not None:
            delta=abs(result['Z']-refZ);agree=delta<=.1
            assert agree==result['FD_agreement']['agrees']
            assert close(delta,result['FD_agreement']['absolute_Z_difference'],atol=1e-7)
            if not agree:assert 'native_FD_significance_mismatch' in result['flags']
        out.append(dict(side=side,mass=mass,Z=result['Z'],p0=result['p0'],log10_p0=result['log10_p0'],p0_underflow=result['p0_underflow'],status=result['status'],flags=result['flags'],
                        numerical_ok=result['numerical_ok'],covariance_ok=result['covariance_ok'],FD_agreement=result['FD_agreement'],receipt_exit=receipt['exit_code']))
    counts=collections.Counter(flag for q in out for flag in q['flags'])
    result=dict(target=76,outcomes=len(out),missing=missing,production_closed=len(out)==76,
        native_complete=sum(q['status']=='complete' for q in out),native_receipts_exit0=sum(q['exit_code']==0 for q in receipts),
        numerical_ok=sum(q['numerical_ok'] for q in out),covariance_ok=sum(q['covariance_ok'] for q in out),
        flags=dict(counts),native_values_preserved=True,FD_used_only_as_check=True,local_one_sided_p0=True,rows=out)
    (R/'significance_completion.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return result,receipts

def inspect_selected(w,fd_path,stage,extracted,side,grouped):
    f=native_file(fd_path);fit=f.Get('fit_'+stage);assert fit
    for v in list(fit.floatParsFinal())+list(fit.constPars()):
        if v.InheritsFrom('RooRealVar') and w.var(v.GetName()):w.var(v.GetName()).setVal(v.getVal())
    if stage=='b':w.var('r').setVal(0.)
    native_parameters={v.GetName():v.getVal() for v in fit.floatParsFinal()}
    assert native_parameters.keys()=={q['name'] for q in extracted['parameters']}
    assert all(close(native_parameters[q['name']],q['value']) for q in extracted['parameters'])
    assert fit.status()==extracted['status'] and fit.covQual()==extracted['covQual'] and close(fit.edm(),extracted['edm'])
    means=[v for v in fit.floatParsFinal() if v.GetName().startswith('Ymean_')];assert len(means)==8
    degree=6 if side=='low' else 4;minima=[];continuity=[];D=0.;minimum_bin=math.inf
    model=w.obj('ModelConfig' if stage=='s' else 'ModelConfig_bonly').GetPdf()
    edge=w.var('mass_SR').getMax() if side=='low' else w.var('mass_SR').getMin()
    for i,era in enumerate(ERAS,1):
        density={}
        for region,suffix,deg in [('SR','',degree),('CR','_CR',4)]:
            coeff=[1.]+[w.var(f'local_c{k}_{era}{suffix}').getVal() for k in range(1,deg+1)]
            polynomial=np.polynomial.Chebyshev(coeff);roots=polynomial.deriv().roots()
            positions=[-1.,1.]+[q.real for q in roots if abs(q.imag)<1e-9 and -1<q.real<1]
            minima.append(float(min(polynomial(t) for t in positions)))
            x=w.var('mass_'+region);norm=ROOT.RooArgSet(x);cp=model.getPdf(f'{region}_ch{i}')
            assert cp.InheritsFrom('RooAddPdf')
            x.setVal(edge);components={}
            for pdf,rate in zip(cp.pdfList(),cp.coefList()):
                name=pdf.GetName();key='dy' if name.startswith('local_continuum_') else next(k for k in ['Zd','upsilon1s','upsilon2s'] if name.startswith('shape_'+k+'_'))
                components[key]=pdf.getVal(norm)*rate.getVal()
            components['total']=sum(components.values());density[region]=components
            expected=cp.expectedEvents(norm);width=(x.getMax()-x.getMin())/x.getBins();count=0.;term=0.;valid=True
            for xx,n in grouped[f'{region}_ch{i}']:
                x.setVal(xx);mu=cp.getVal(norm)*expected*width;minimum_bin=min(minimum_bin,mu);count+=n
                if mu<=0 or not math.isfinite(mu):valid=False
                elif n>0:term+=n*math.log(n/mu)
            D+=2*(expected-count+term) if valid else math.inf
        for key in density['SR']:
            left=density['SR'][key];right=density['CR'][key]
            continuity.append(abs(left-right)/max(abs(left),abs(right),1e-12))
    globals=list(w.obj('ModelConfig').GetGlobalObservables());assert len(globals)==10
    D+=sum((w.var(v.GetName()[:-3]).getVal()-v.getVal())**2 for v in globals)
    assert close(D,extracted['exact_deviance'],rtol=1e-10,atol=1e-7),(side,fd_path,stage,D,extracted['exact_deviance'])
    assert close(min(minima),extracted['minimum'],atol=1e-10)
    assert close(minimum_bin,extracted['minimum_bin'],atol=1e-9)
    assert max(continuity)<1e-8 and extracted['continuity_passed']
    rv=fit.floatParsFinal().find('r');r_boundary=bool(rv and min(abs(rv.getVal()-rv.getMin()),abs(rv.getVal()-rv.getMax()))<1e-4*max(1.,rv.getMax()-rv.getMin()))
    valid=fit.status()==0 and math.isfinite(fit.edm()) and 0<=fit.edm()<.01 and min(minima)>0 and minimum_bin>0 and math.isfinite(D) and not r_boundary
    assert valid==extracted['point_valid'] and (valid and fit.covQual()==3)==extracted['numerical_ok']
    out=dict(status=fit.status(),covQual=fit.covQual(),EDM=fit.edm(),exact_deviance=D,point_valid=valid,numerical_ok=valid and fit.covQual()==3,
             r=rv.getVal() if rv else 0.,r_boundary=r_boundary,continuity_max_relative_difference=max(continuity),
             minimum_continuum=min(minima),minimum_bin=minimum_bin,means_floating=8,means_at_boundary=extracted['mean_boundary'],
             anchor=fd_path.parent.name,source=str(fd_path),extraction_matches_native=True)
    f.Close();return out

def main():
    points=read(R/'manifest.json')['points'];target={('low',k/10) for k in range(60,106)}|{('high',k/10) for k in range(91,121)}
    assert len(points)==76 and {(q['region'],q['mass']) for q in points}==target
    wm=read(R/'workspace_manifest.json');workspaces={(q['side'],q['mass']):q for q in wm['workspaces']};assert set(workspaces)==target
    rows=[];source_hashes={};receipts=[];integrity=[]
    for point in points:
        side=point['region'];mass=point['mass'];D=R/side/f'm{mass:.1f}';meta=workspaces[(side,mass)]
        assert sha(D/'workspace.root')==meta['workspace_sha256']
        source_hashes.setdefault(meta['source'],sha(meta['source']));assert source_hashes[meta['source']]==meta['source_sha256']
        result=read(D/'results.json');fit_result=read(D/'fit_results.json');completion=read(D/'completion.json');fd_completion=read(D/'fits_completion.json')
        lr=read(D/'limits/combine.receipt.json');assert lr['exit_code']==0;receipts.append(dict(side=side,mass=mass,method='AsymptoticLimits',**lr))
        files=list((D/'limits').glob('higgsCombine*.root'));assert len(files)==1
        f=native_file(files[0]);tree=f.Get('limit');assert tree and tree.GetEntries()==6
        raw=[(float(v.quantileExpected),float(v.limit)) for v in tree];f.Close();limits={}
        for name,q in QUANTILES.items():
            hits=[v for z,v in raw if abs(z-q)<2e-6];assert len(hits)==1 and math.isfinite(hits[0]) and hits[0]>0
            limits[name]=hits[0];assert limits[name]==result['limits'][name]
        expected=[limits[k] for k in QUANTILES if k!='observed'];assert all(a<=b for a,b in zip(expected,expected[1:]))
        text=(D/'limits/combine.log').read_text(errors='replace')
        pairs=[(float(a),float(b)) for a,b in re.findall(r'At r =\s*([\d.eE+-]+):[^\n]*CLs\s*=\s*([\d.eE+-]+)',text)]
        crossing=dict(verified=False)
        if pairs:
            rr,cls=min(pairs,key=lambda v:abs(v[0]-limits['observed']));relative=abs(rr-limits['observed'])/max(abs(limits['observed']),1e-12)
            crossing=dict(r=rr,CLs=cls,r_relative_distance=relative,verified=abs(cls-.05)<.002,r_proximity_verified=relative<.001)
        assert crossing['verified']==result['observed_limit_check']['verified']
        trials=[]
        for directory in sorted(D.glob('fitdiagnostics_anchor_*')):
            receipt=read(directory/'combine.receipt.json');assert receipt['exit_code']==0
            command=read(directory/'combine.command.json');fd=native_file(directory/'fitDiagnostics_c0.root')
            required=['s']+([] if '--skipBOnlyFit' in command else ['b'])
            assert all(fd.Get('fit_'+stage) for stage in required)
            trials.append(dict(anchor=directory.name,required_stages=required,stage_statuses={stage:fd.Get('fit_'+stage).status() for stage in required},**receipt))
            receipts.append(dict(side=side,mass=mass,method='FitDiagnostics',anchor=directory.name,**receipt));fd.Close()
        assert {q['anchor'] for q in trials}=={'fitdiagnostics_anchor_zero','fitdiagnostics_anchor_expected'}
        wf=native_file(D/'workspace.root');w=wf.Get('w');data=w.data('data_obs');grouped=collections.defaultdict(list)
        for k in range(data.numEntries()):
            row=data.get(k);ch=row.find('CMS_channel').getLabel();region=ch.split('_')[0];grouped[ch].append((row.getRealValue('mass_'+region),data.weight()))
        selected={}
        for stage in ['b','s']:
            q=fit_result['fit'][stage];path=Path(q['source']);assert path.is_relative_to(D)
            selected[stage]=inspect_selected(w,path,stage,q,side,grouped)
        wf.Close()
        assert result['status']==completion['status']=='complete' and fit_result['status']==fd_completion['status']=='complete'
        assert all(selected[s]['point_valid'] for s in ['b','s'])==fit_result['point_valid']
        assert all(selected[s]['numerical_ok'] for s in ['b','s'])==fit_result['numerical_ok']
        delta=selected['b']['exact_deviance']-selected['s']['exact_deviance']
        assert close(delta,fit_result['delta_deviance_B_minus_SB'],atol=2e-7)
        rows.append(dict(side=side,mass=mass,complete=True,limits=limits,crossing=crossing,limit_flags=result['flags'],expected_bands_validated=False,
            fit=selected,fit_flags=fit_result['flags'],fit_trials=trials,delta_deviance_B_minus_SB=delta,workspace_source_hashes_unchanged=True))
        print('AUDITED',side,f'{mass:.1f}',flush=True)
    summaries={}
    for side in ['low','high']:
        q=[a for a in rows if a['side']==side];b=[a['fit']['b']['exact_deviance'] for a in q]
        summaries[side]=dict(points=len(q),B_only_deviance_min=min(b),B_only_deviance_max=max(b),B_only_deviance_span=max(b)-min(b),B_only_consistent=max(b)-min(b)<1e-6)
    assert all(q['B_only_consistent'] for q in summaries.values())
    invalid=[dict(side=q['side'],mass=q['mass'],stage=s,**{k:q['fit'][s][k] for k in ['status','covQual','EDM','r_boundary']}) for q in rows for s in ['b','s'] if not q['fit'][s]['point_valid']]
    covariance_only=[dict(side=q['side'],mass=q['mass'],stage=s,covQual=q['fit'][s]['covQual']) for q in rows for s in ['b','s'] if q['fit'][s]['point_valid'] and q['fit'][s]['covQual']!=3]
    unverified=[dict(side=q['side'],mass=q['mass'],**q['crossing']) for q in rows if not q['crossing']['verified']]
    counts=dict(target_points=76,completed_limits=len(rows),completed_fit_points=len(rows),native_receipts=len(receipts),native_receipts_exit0=len(receipts),
        native_limit_quantiles=6*len(rows),selected_fit_stages=2*len(rows),selected_fit_points_valid=sum(all(q['fit'][s]['point_valid'] for s in ['b','s']) for q in rows),
        selected_fit_points_cov3_valid=sum(all(q['fit'][s]['numerical_ok'] for s in ['b','s']) for q in rows),observed_crossings_verified=sum(q['crossing']['verified'] for q in rows),
        expected_band_bracket_flags=sum('expected_bands_bracket_hit' in q['limit_flags'] for q in rows),nonzero_asimov_bestfit_warnings=sum('nonzero_asimov_bestfit_warning' in q['limit_flags'] for q in rows))
    manifest=read(R/'manifest.json');significance_requested=bool(manifest.get('significance_requested') or 'Significance' in manifest.get('requested_methods',[]))
    significance=None
    if significance_requested:
        significance,new_receipts=inspect_significance(rows);receipts.extend(new_receipts)
        counts.update(significance_outcomes=significance['outcomes'],significance_native_complete=significance['native_complete'],
                      significance_numerical_ok=significance['numerical_ok'],native_receipts=len(receipts),native_receipts_exit0=sum(q['exit_code']==0 for q in receipts))
    output=dict(audit_passed=True,production_complete=not significance_requested or significance['production_closed'],all_results_numerically_valid=False,
        timestamp=datetime.datetime.now().astimezone().isoformat(),read_only=True,fits_launched=0,counts=counts,B_only_consistency=summaries,
        max_continuity_relative_difference=max(q['fit'][s]['continuity_max_relative_difference'] for q in rows for s in ['b','s']),
        point_invalid=invalid,covariance_only_flags=covariance_only,unverified_observed_crossings=unverified,source_workspaces_unchanged=True,
        limits_global_minimum_verified=False,expected_bands_validated=False,significance_requested=significance_requested,significance=significance,
        limitations=['A verified CLs crossing does not establish the global likelihood minimum used by native AsymptoticLimits.',
                     'Native expected bands retained unchanged; bracket and Asimov warnings are preserved.',
                     'A converged point with covariance quality 2 is distinct from a failed point; its parameter errors remain flagged.',
                     'Low+CR and high+CR reuse the same CR data and are correlated separate analyses, not a statistical combination.',
                     'Significance, when requested, is native and checked against FD; no new bias or goodness-of-fit test is included.'],rows=rows,receipts=receipts)
    (R/'completion_audit.json').write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    lines=['# Audit di completamento C0','',f'Campagna completa: {len(rows)}/76 limiti e {len(rows)}/76 punti FitDiagnostics. Audit in sola lettura; nessun nuovo fit.',
           '',f"Controllati {len(receipts)} ricevute native con exit 0, {6*len(rows)} quantili ROOT e {2*len(rows)} stati B/S+B selezionati. Valori estratti, parametri, devianze e stato numerico concordano con i file nativi. Input e workspace iniziali conservano gli hash del manifest.",
           '',f"Punti con B e S+B convergenti: {counts['selected_fit_points_valid']}/76. Con covarianze entrambe di qualità 3: {counts['selected_fit_points_cov3_valid']}/76. I {len(covariance_only)} flag di covarianza non sono classificati come fallimenti della minimizzazione.",
           '',f"Crossing osservati verificati: {counts['observed_crossings_verified']}/76. Bande attese native non validate: {counts['expected_band_bracket_flags']} flag di bracket e {counts['nonzero_asimov_bestfit_warnings']} avvisi sul best fit Asimov. Un crossing corretto non prova il minimo globale del denominatore.",
           '',f"Continuità C0 verificata nei 152 stati selezionati; massimo errore relativo {output['max_continuity_relative_difference']:.3g}.",'', '## Punti non convergenti','']
    lines+=['| Analisi | Massa | Fit | Status | Covarianza | EDM |','|---|---:|---|---:|---:|---:|']
    lines += [f"| {q['side']} | {q['mass']:.1f} | {q['stage']} | {q['status']} | {q['covQual']} | {q['EDM']:.6g} |" for q in invalid]
    lines+=['','## Coerenza B-only','']
    for side,q in summaries.items():lines.append(f"- {side}: D = {q['B_only_deviance_min']:.9f}–{q['B_only_deviance_max']:.9f}; dispersione {q['B_only_deviance_span']:.3g} fra le masse.")
    lines+=['','Crossing osservati non verificati: '+', '.join(f"{q['side']} {q['mass']:.1f}" for q in unverified)+'.','',
            'Le due analisi condividono la CR e non sono combinate. Completezza dei prodotti e validità numerica sono valutazioni distinte. Nessun risultato o flag precedente è stato corretto o rimosso.']
    if significance:
        lines+=['','## Significanze native aggiunte su richiesta','',
                f"Esiti {significance['outcomes']}/76; valori nativi completi {significance['native_complete']}/76; confronto FD numericamente accettabile {significance['numerical_ok']}/76. I p-value sono locali, a una coda. La qualità della covarianza FD resta un flag distinto dalla convergenza e dall'accordo della significanza.",
                '', 'Flag: '+json.dumps(significance['flags'],ensure_ascii=False)+'.',
                '', 'Nessuna significanza FD sostituisce il valore nativo e non è stato avviato alcun tentativo di correzione.']
    (R/'completion_audit.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(counts),flush=True)

if __name__=='__main__':main()
