# -*- coding: utf-8 -*-
"""Read-only, incremental report for the common SR-low/SR-high continuum scan.

Input files are never modified. Missing points remain missing; fit failures are
retained, and the primary bias estimate is explicitly conditional on valid fits.
"""
from pathlib import Path
import argparse
import csv
from collections import Counter
import datetime
import html
import hashlib
import json
import math
import os
import re
import time
from zoneinfo import ZoneInfo

os.environ.setdefault('MPLCONFIGDIR', str(Path.cwd() / 'work/mpl-cache'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import beta, norm
import uproot

P = Path.cwd()
R = P / 'runs/20261031_shared_sr_cr_scan'
MASSES = [k / 10 for k in range(60, 121)]
ERAS = ['2022', '2022EE', '2023', '2023BPix']
REGIONS = [('SRL', 'SR bassa'), ('CR', 'CR'), ('SRH', 'SR alta')]
INJECTIONS = ['r0', 'rExpected']
COLORS = ['#0072B2', '#D55E00']


def read(path, default=None):
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return default


def finite(x):
    return isinstance(x, (float, int)) and math.isfinite(x)


def write_json(path, value):
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, indent=2, allow_nan=False))
    temp.replace(path)


def csv_write(path, rows):
    if not rows:
        return
    keys = list(dict.fromkeys(k for q in rows for k in q))
    with path.open('w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def fmt(x, digits=4):
    return f'{x:.{digits}g}' if finite(x) else '—'


def fit_notes(q):
    notes=[]
    flags=q.get('flags', [])
    if any(f.startswith('fit_') or 'minimum' in f or 'Significance' in f for f in flags):
        notes.append('FLAG fit/significanza: vedere JSON')
    if any(not finite(stage.get('edm')) or stage['edm']>=.01 for stage in q.get('fit', {}).values()):
        notes.append('EDM fit non valido/non disponibile')
    if q.get('limits') and not q.get('expected_bands_verified', q.get('expected_validated', False)):
        notes.append('bande attese provvisorie')
    if q.get('status')!='complete':
        notes.append('scan parziale')
    return '; '.join(notes)


def observed_fit_ok(q):
    return all(q.get('fit', {}).get(stage, {}).get('numerical_ok', False)
               and finite(q['fit'][stage].get('edm')) and q['fit'][stage]['edm']<.01
               for stage in ['b', 's'])


def nominal_qc(diagnostics):
    """Summarize an audit outcome without retaining its parameter snapshot."""
    if not isinstance(diagnostics, dict) or not (diagnostics.get('outcome_complete') is True or diagnostics.get('available') is True):
        return dict(state='unchecked', checked=False, valid=False, status=None,
                    edm=None, flags=[], outcome_complete=False)
    d=diagnostics
    edm=d.get('edm')
    required=['minimizer_ok', 'fitresult_state_consistent', 'snapshot_complete',
              'continuum_positive', 'total_bins_positive', 'score_state_consistent']
    ok=(d.get('valid') is True and d.get('available') is True
        and finite(edm) and 0<=edm<.01
        and all(d.get(key) is True for key in required)
        and d.get('previous_score_parity') is not False
        and d.get('source_unchanged') is not False)
    return dict(state='checked_valid' if ok else 'checked_failed', checked=True,
                valid=ok, status=d.get('status'), edm=edm, flags=d.get('flags', []),
                outcome_complete=True)


def cp95(k, n):
    if not n:
        return [None, None]
    return [0. if k == 0 else float(beta.ppf(.025, k, n-k+1)),
            1. if k == n else float(beta.ppf(.975, k+1, n-k))]


def stats(values):
    a = np.array([x for x in values if finite(x)])
    if not len(a):
        return dict(n=0, mean=None, sd=None, se=None)
    sd = float(np.std(a, ddof=1)) if len(a) > 1 else None
    return dict(n=len(a), mean=float(a.mean()), sd=sd,
                se=sd / math.sqrt(len(a)) if sd is not None else None)


def bias_ratio(values, seed):
    a = np.array([x for x in values if finite(x)])
    if len(a) < 3 or np.std(a, ddof=1) <= 0:
        return dict(value=None, se=None, ci95=[None, None])
    sd = np.std(a, ddof=1)
    rng = np.random.default_rng(seed)
    samples = rng.choice(a, size=(2000, len(a)), replace=True)
    sample_sd = samples.std(axis=1, ddof=1)
    ok = sample_sd > max(1e-12, sd * 1e-10)
    ratios = samples.mean(axis=1)[ok] / sample_sd[ok]
    return dict(value=float(a.mean()/sd), se=float(ratios.std(ddof=1)),
                ci95=[float(x) for x in np.quantile(ratios, [.025, .975])])


def collect():
    observed = []
    for path in sorted((R/'observed').glob('m*/results.json')):
        q = read(path)
        if q is not None:
            observed.append(q)
    bias = []
    cr_info_cache={}
    for path in sorted((R/'bias/fits').glob('m*/*/toy*/result.json')):
        q = read(path)
        if q is None:
            continue
        mass = q.get('mass', float(path.parents[2].name[1:]))
        label = q.get('label', path.parents[1].name)
        toy = q.get('toy', int(path.parent.name[3:]))
        if toy < 1 or toy > 50 or label not in INJECTIONS:
            continue
        r, truth = q.get('r', q.get('rhat')), q.get('rtrue')
        flags = q.get('failure_flags', q.get('flags', []))
        edm=q.get('edm'); edm_ok=finite(edm) and edm<.01
        cr_edm=q.get('CR_edm')
        cache_directory=q.get('CR_cache_directory')
        cr_info={}
        if cache_directory:
            if cache_directory not in cr_info_cache:
                cr_result=read(Path(cache_directory)/'cr_result.json', {})
                selected=cr_result.get('selected', {})
                cr_info_cache[cache_directory]=dict(
                    edm=selected.get('edm'), valid=cr_result.get('valid'),
                    means={item['name']:item.get('value') for item in cr_result.get('means', [])},
                    boundary_parameters=selected.get('boundary_parameters', []))
            cr_info=cr_info_cache[cache_directory]
            if cr_edm is None:cr_edm=cr_info.get('edm')
        cr_edm_ok=finite(cr_edm) and cr_edm<.01
        bias.append(dict(mass=mass, injection=label, toy=toy, r=r, rtrue=truth,
                         bias=r-truth if finite(r) and finite(truth) else None,
                         valid=bool(q.get('fit_valid', q.get('valid', False))) and edm_ok and cr_edm_ok,
                         edm=edm, edm_valid=edm_ok, CR_edm=cr_edm, CR_edm_valid=cr_edm_ok,
                         numerical_validity_version=q.get('numerical_validity_version'),
                         CR_valid=q.get('CR_valid'), r_boundary=bool(q.get('r_boundary', False)),
                         CR_mu_boundary=any(str(name).startswith('Ymean_') for name in q.get('CR_boundary_parameters', [])),
                         CR_boundary_parameters=q.get('CR_boundary_parameters', []),
                         CR_means=cr_info.get('means', {}), CR_cache_directory=cache_directory,
                         CR_result_valid=cr_info.get('valid'), truth_source=q.get('truth_source'),
                         status=q.get('status'), covQual=q.get('covQual'),
                         flags=flags, result=str(path.relative_to(R))))
    gof = []
    observed_by_mass = {q['mass']: q for q in observed}
    for mass in MASSES:
        records = []
        for path in sorted((R/'gof').glob(f'toy*/m{mass:.1f}/result.json')):
            q = read(path)
            if q is not None:
                qc=nominal_qc(q.get('nominal_diagnostics'))
                records.append(dict(score=q.get('score'), qc=qc,
                                    CR_valid=q.get('CR_valid'), CR_edm=q.get('CR_edm'),
                                    r_boundary=bool(q.get('r_boundary', False)),
                                    source_valid=q.get('valid'), flags=q.get('flags', [])))
        qobs = observed_by_mass.get(mass, {})
        score = qobs.get('gof_observed')
        obs_diagnostics=qobs.get('gof_diagnostics')
        if not finite(score):
            standalone = read(R/f'gof/observed/m{mass:.1f}/result.json', {})
            score = standalone.get('score')
            if obs_diagnostics is None:
                obs_diagnostics=standalone.get('nominal_diagnostics')
        obs_qc=nominal_qc(obs_diagnostics)
        if not records and not finite(score) and not obs_qc['checked']:
            continue
        scores = [q for q in records if finite(q.get('score'))]
        n = len(scores); attempted=len(records)
        checked_valid=[q for q in records if q['qc']['valid']]
        # A failed CR calibration is also unresolved for the full two-stage protocol.
        resolved=[q for q in checked_valid if finite(q['score']) and q['CR_valid'] is True and not q['r_boundary']
                  and finite(q['CR_edm']) and 0<=q['CR_edm']<.01]
        unresolved=attempted-len(resolved)
        k = sum(q['score'] > score for q in scores) if finite(score) else None
        kr = sum(q['score'] > score for q in resolved) if finite(score) else None
        ci = cp95(k, n) if k is not None else [None, None]
        checked=sum(q['qc']['checked'] for q in records)
        failed=sum(q['qc']['state']=='checked_failed' for q in records)
        p_low=kr/attempted if kr is not None and attempted else None
        p_high=(kr+unresolved)/attempted if kr is not None and attempted else None
        gof.append(dict(mass=mass, observed=score, observed_qc=obs_qc['state'],
                        observed_nominal_status=obs_qc['status'], observed_nominal_edm=obs_qc['edm'],
                        observed_qc_flags=obs_qc['flags'], observed_diagnostic_completed=obs_qc['checked'],
                        observed_nominal_valid=obs_qc['valid'],
                        n_attempted=attempted, n_finite=n, k=k,
                        p=k/n if k is not None and n else None,
                        ci95_low=ci[0], ci95_high=ci[1],
                        n_nominal_checked=checked, n_nominal_valid=len(checked_valid),
                        n_nominal_failed=failed, n_nominal_unchecked=attempted-checked,
                        nominal_failure_flags=dict(Counter(flag for q in records if q['qc']['state']=='checked_failed' for flag in q['qc']['flags'])),
                        n_nominal_EDM_failed=sum(q['qc']['checked'] and (not finite(q['qc']['edm']) or not 0<=q['qc']['edm']<.01) for q in records),
                        n_resolved=len(resolved), n_unresolved=unresolved,
                        k_resolved=kr, p_qc_low=p_low, p_qc_high=p_high,
                        qc_bound_denominator=attempted,
                        qc_bound_interpretation='Conditional on the stored observed score; invalid, unchecked, nonfinite, CR-invalid or r-boundary completed toy outcomes are unresolved. Not a binomial confidence interval.',
                        calibrated_conclusion_available=obs_qc['valid'] and unresolved==0 and attempted==100,
                        n_nonfinite=attempted-n,
                        p_missing_score_low=k/attempted if k is not None and attempted else None,
                        p_missing_score_high=(k+attempted-n)/attempted if k is not None and attempted else None,
                        n_CR_invalid=sum(q['CR_valid'] is False for q in records),
                        n_CR_EDM_high=sum(finite(q['CR_edm']) and not 0<=q['CR_edm']<.01 for q in records),
                        n_CR_EDM_missing=sum(not finite(q['CR_edm']) for q in records),
                        n_r_boundary=sum(q['r_boundary'] for q in records)))
    profiles=[]
    for path in sorted((R/'bias/profiles').glob('m*/*/toy001/result.json')):
        q=read(path)
        if q is None or q.get('outcome_complete') is False: continue
        mass=q.get('mass', float(path.parents[2].name[1:]))
        label=q.get('injection', path.parents[1].name)
        if mass not in MASSES or label not in INJECTIONS: continue
        point_retained=bool(q.get('profile_failed_point_retained', False))
        edm=q.get('edm') if not point_retained else None
        edm_ok=finite(edm) and edm<.01
        exception=str(q.get('exception') or '')
        source_outcome=q.get('profile_outcome')
        skipped=str(source_outcome or '').startswith('skipped_') or bool(q.get('skipped_invalid_anchor')) or any(text in exception.lower() for text in ['no numerically valid immutable anchor', 'no readable point estimate', 'no valid anchor'])
        valid=bool(q.get('profile_valid')) and edm_ok
        profiles.append(dict(mass=mass, injection=label, toy=1, outcome_completed=True,
                             profile_valid=valid, source_outcome=source_outcome, outcome='valid_interval' if valid else 'skipped_invalid_or_missing_anchor' if skipped else 'failed_or_unavailable_interval',
                             edm=edm, edm_valid=edm_ok,
                             profile_found_lower_minimum=bool(q.get('profile_found_lower_minimum', False)),
                             point_r=q.get('point_r'), profile_r=q.get('r') if not point_retained else None,
                             point_retained_after_profile_failure=point_retained,
                             interval_low=q.get('interval_low'), interval_high=q.get('interval_high'),
                             primary_point_estimate_unchanged=q.get('primary_point_estimate_unchanged', False),
                             exception=exception, result=str(path.relative_to(R))))
    return sorted(observed, key=lambda q:q['mass']), bias, gof, profiles


def summarize_bias(records):
    summaries = []
    for mass in MASSES:
        for li, label in enumerate(INJECTIONS):
            rows = [q for q in records if q['mass'] == mass and q['injection'] == label]
            if not rows:
                continue
            valid = [q for q in rows if q['valid']]
            summaries.append(dict(mass=mass, injection=label, n=len(rows), n_valid=len(valid),
                                  n_invalid=len(rows)-len(valid),
                                  n_CR_invalid=sum(q['CR_valid'] is False for q in rows),
                                  n_CR_mu_boundary=sum(q['CR_mu_boundary'] for q in rows),
                                  n_EDM_failure=sum(not q['edm_valid'] for q in rows),
                                  n_EDM_high=sum(finite(q['edm']) and q['edm']>=.01 for q in rows),
                                  n_EDM_missing=sum(not finite(q['edm']) for q in rows),
                                  n_CR_EDM_failure=sum(not q['CR_edm_valid'] for q in rows),
                                  n_CR_EDM_high=sum(finite(q['CR_edm']) and q['CR_edm']>=.01 for q in rows),
                                  n_CR_EDM_missing=sum(not finite(q['CR_edm']) for q in rows),
                                  n_r_boundary=sum(q['r_boundary'] for q in rows),
                                  conditional=stats(q['bias'] for q in valid),
                                  all_finite=stats(q['bias'] for q in rows),
                                  normalized_bias=bias_ratio([q['bias'] for q in valid],
                                                            310000+round(mass*10)*10+li)))
    return summaries


def summarize_cr_mean_shifts(records):
    """Use recorded CR fits and generator truth; never rerun or alter fits."""
    generator=R/'bias/generator/results.json'
    gen=read(generator, {})
    gen_hash=hashlib.sha256(generator.read_bytes()).hexdigest() if generator.exists() else None
    truth_means={item['name']:item.get('value') for item in gen.get('means', [])}
    names=[f'Ymean_{res}_{era}' for res in ['upsilon1s','upsilon2s'] for era in ERAS]
    truth_cache={}; rows=[]; issues=[]
    for q in records:
        source=q.get('truth_source')
        if source not in truth_cache:
            path=Path(source) if source else None
            truth=read(path,{}) if path else {}
            truth_cache[source]=dict(matches=truth.get('generator_fit_sha256')==gen_hash,
                                     sha256=hashlib.sha256(path.read_bytes()).hexdigest() if path and path.exists() else None,
                                     recorded_generator_fit_sha256=truth.get('generator_fit_sha256'))
        metadata=truth_cache[source]
        if not metadata['matches']:
            if source not in issues:issues.append(source)
            continue
        cr_valid=(q.get('CR_valid') is True and q.get('CR_result_valid') is True
                  and finite(q.get('CR_edm')) and 0<=q['CR_edm']<.01)
        for name in names:
            truth=truth_means.get(name); fitted=q['CR_means'].get(name)
            if not finite(truth) or not finite(fitted):continue
            rows.append(dict(mass=q['mass'],injection=q['injection'],toy=q['toy'],parameter=name,
                             generator_mean_GeV=truth,fitted_mean_GeV=fitted,shift_MeV=1000*(fitted-truth),
                             CR_valid=cr_valid,point_valid=q['valid'],
                             parameter_at_boundary=name in q['CR_boundary_parameters'],
                             any_CR_mean_at_boundary=q['CR_mu_boundary'],
                             CR_result=str(Path(q['CR_cache_directory'])/'cr_result.json'),
                             point_result=q['result'],truth_source=source))
    summaries=[]
    for mass in MASSES:
        for label in INJECTIONS:
            cell=[q for q in rows if q['mass']==mass and q['injection']==label]
            for name in names:
                rr=[q for q in cell if q['parameter']==name]
                if not rr:continue
                valid=[q for q in rr if q['CR_valid']]
                point_valid=[q for q in rr if q['point_valid']]
                sv=stats(q['shift_MeV'] for q in valid); sa=stats(q['shift_MeV'] for q in rr)
                sp=stats(q['shift_MeV'] for q in point_valid)
                summaries.append(dict(mass=mass,injection=label,parameter=name,
                    generator_mean_GeV=rr[0]['generator_mean_GeV'],n_finite=len(rr),n_CR_valid=len(valid),
                    n_point_valid=len(point_valid),n_boundary=sum(q['parameter_at_boundary'] for q in rr),
                    n_boundary_CR_valid=sum(q['parameter_at_boundary'] for q in valid),
                    CR_valid_mean_shift_MeV=sv['mean'],CR_valid_SE_MeV=sv['se'],CR_valid_SD_MeV=sv['sd'],
                    all_finite_mean_shift_MeV=sa['mean'],all_finite_SE_MeV=sa['se'],all_finite_SD_MeV=sa['sd'],
                    point_valid_mean_shift_MeV=sp['mean'],point_valid_SE_MeV=sp['se']))
    provenance=dict(generator_results=str(generator),generator_results_sha256=gen_hash,
                    means_GeV=truth_means,truth_metadata=truth_cache,unmatched_truth_sources=issues,
                    truth_source='Existing generator/results.json means; generation metadata generator_fit_sha256 verified per cell.',
                    generation_code=str(P/'scripts/shared_sr_cr_scan/generator_toys.py'),
                    interpretation='Mean shifts from B-only CR calibration relative to the common generator truth. r=0 datasets/CR fits are reused across masses; expected-injection samples are independent by mass. These shifts do not prove the cause of any observed excess.',
                    validity='Primary subset: CR validity and 0<=EDM<0.01; all finite and final point-valid subsets retained separately. Boundary parameters are flagged, not silently removed.')
    return summaries,rows,provenance


def cr_mean_shift_plots(summaries):
    if not summaries:return
    signature=hashlib.sha256(json.dumps(summaries,sort_keys=True).encode()).hexdigest()
    cache=R/'cr_mean_shift_plot_cache.json'
    outputs=['cr_mean_shift','cr_mean_shift_all_finite']
    if read(cache,{}).get('signature')==signature and all((R/f'{name}.png').exists() for name in outputs):return
    names=[f'Ymean_{res}_{era}' for res in ['upsilon1s','upsilon2s'] for era in ERAS]
    labels=[f'{"Υ(1S)" if "upsilon1s" in name else "Υ(2S)"} · '+name.split('_',2)[2] for name in names]
    all_values=[q[key] for q in summaries for key in ['CR_valid_mean_shift_MeV','all_finite_mean_shift_MeV'] if finite(q[key])]
    vmax=max([1.]+[abs(v) for v in all_values])
    edges=np.arange(5.95,12.051,.1)
    cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#ededed')
    for filename,key,selection,nkey in [
        ('cr_mean_shift','CR_valid_mean_shift_MeV','CR numericamente valide','n_CR_valid'),
        ('cr_mean_shift_all_finite','all_finite_mean_shift_MeV','tutti i valori finiti','n_finite')]:
        fig,axes=plt.subplots(1,2,figsize=(16,5.6),sharey=True)
        for label,ax in zip(INJECTIONS,axes):
            lookup={(q['parameter'],q['mass']):q for q in summaries if q['injection']==label}
            values=np.array([[lookup.get((name,m),{}).get(key,np.nan) for m in MASSES] for name in names],dtype=float)
            im=ax.pcolormesh(edges,np.arange(9),np.ma.masked_invalid(values),cmap=cmap,vmin=-vmax,vmax=vmax,shading='flat')
            for yi,name in enumerate(names):
                for m in MASSES:
                    q=lookup.get((name,m),{})
                    boundary=q.get('n_boundary_CR_valid',0) if nkey=='n_CR_valid' else q.get('n_boundary',0)
                    if boundary and finite(q.get(key)):ax.plot(m,yi+.5,'x',color='black',ms=4,mew=.8)
            ax.axvline(9.093333333333334,color='black',ls='--',lw=.8)
            ax.axvline(10.506666666666668,color='black',ls='--',lw=.8)
            counts=[q[nkey] for q in lookup.values()]
            count_label=f'{min(counts)}–{max(counts)}' if counts else '0'
            ax.set_title(('r vero = 0' if label=='r0' else 'r vero = limite atteso mediano')+f' | n/cella: {count_label}',fontsize=11)
            ax.set_yticks(np.arange(8)+.5,labels=labels)
            ax.set_xlabel('Ipotesi di massa del segnale [GeV]');ax.set_xlim(5.95,12.05)
        axes[0].invert_yaxis()
        fig.suptitle('Spostamento medio delle otto medie Υ nella calibrazione B-only della CR\n'+selection+'; Δμ = μ stimata − μ generatore',fontsize=13)
        fig.colorbar(im,ax=axes,label='Spostamento medio [MeV]',fraction=.025,pad=.025)
        fig.text(.5,.035,'Grigio: dati ancora assenti. ×: almeno una stima della stessa media al bordo. Tratteggi: confini della CR.\nr=0: stessi 50 campioni/fit CR; a completamento il pannello è piatto per costruzione. Nessuna prova della causa degli eccessi.',ha='center',fontsize=9)
        fig.subplots_adjust(top=.81,bottom=.2,wspace=.12,right=.87);save(fig,filename)
    write_json(cache,dict(signature=signature))


def shade(ax):
    ax.axvspan(9.093333333333334, 10.506666666666668, color='#b0b0b0', alpha=.15,
               label='Ipotesi di massa nella CR')
    ax.axvspan(11, 12, color='#f6d17f', alpha=.13, label='Resa segnale estrapolata')
    ax.set_xlim(5.97, 12.03)
    ax.set_xlabel('Massa ipotizzata [GeV]')
    ax.grid(alpha=.18)


def save(fig, name):
    fig.savefig(R/f'{name}.png', dpi=160, bbox_inches='tight')
    fig.savefig(R/f'{name}.svg', bbox_inches='tight')
    plt.close(fig)


def observed_plots(observed):
    lookup = {q['mass']: q for q in observed}
    x = np.array(MASSES)
    def values(key):
        return np.array([lookup.get(m, {}).get('limits', {}).get(key, np.nan) for m in MASSES], dtype=float)
    any_limits = any(q.get('limits') for q in observed)
    if any_limits:
        fig, ax = plt.subplots(figsize=(12, 5.5))
        ax.fill_between(x, values('expected_2p5'), values('expected_97p5'), color='#ffcf40', label='Atteso ±2σ')
        ax.fill_between(x, values('expected_16'), values('expected_84'), color='#67bd63', label='Atteso ±1σ')
        ax.plot(x, values('expected_50'), '.--', ms=3, c='black', label='Atteso mediano')
        present=np.isfinite(values('expected_50'))
        isolated=present & ~np.r_[False,present[:-1]] & ~np.r_[present[1:],False]
        # Isolated completed hypotheses have no filled-band area yet. Show their
        # actual quantiles as vertical bars without bridging pending masses.
        for lo,hi,color,width in [('expected_2p5','expected_97p5','#ffcf40',7),('expected_16','expected_84','#67bd63',4)]:
            good=isolated & np.isfinite(values(lo)) & np.isfinite(values(hi))
            ax.vlines(x[good],values(lo)[good],values(hi)[good],color=color,lw=width,zorder=1)
        ax.plot(x, values('observed'), 'o-', ms=3, c='#16345d', label='Osservato')
        for q in observed:
            value = q.get('limits', {}).get('observed')
            if finite(value) and not (observed_fit_ok(q) and q.get('observed_limit_check', {}).get('verified', False)):
                ax.plot(q['mass'], value, 'x', c='red', ms=7)
        ax.set_yscale('log')
        ax.set_ylabel('Limite superiore 95% CLs su r')
        title = 'Limiti: continuo comune SR bassa/alta, envelope Cheb6/Bern6 + CR4'
        provisional = any(q.get('limits') and not q.get('expected_bands_verified', q.get('expected_validated', False)) for q in observed)
        ax.set_title(title + ('\nBande native provvisorie: crossing dei quantili non verificati' if provisional else ''))
        shade(ax); ax.legend(fontsize=8, ncol=2); save(fig, 'limits')
    if any(finite(q.get('p0')) for q in observed):
        fig, ax = plt.subplots(figsize=(12, 5.5))
        y = [lookup.get(m, {}).get('p0', np.nan) for m in MASSES]
        ax.plot(x, y, 'o-', color='#0072B2', ms=3)
        for sigma in [1, 2, 3, 4, 5]:
            ax.axhline(norm.sf(sigma), color='gray', lw=.6, ls='--')
            ax.text(12.02, norm.sf(sigma), f'{sigma}σ', va='bottom', ha='right', fontsize=8)
        for q in observed:
            if finite(q.get('p0')) and not (observed_fit_ok(q) and q.get('significance_consistent', False)):
                ax.plot(q['mass'], q['p0'], 'x', c='red', ms=7)
        ax.plot([], [], 'x', c='red', ms=7, label='Flag numerici del fit o della significanza')
        ax.set_yscale('log'); ax.set_ylim(bottom=min(1e-7, min([y for y in y if finite(y) and y > 0], default=1e-7)*.5), top=1)
        ax.set_ylabel('p-value locale, eccesso')
        ax.set_title('Significanza locale: unico r, segnale DCB in SR bassa, CR e SR alta')
        shade(ax); ax.legend(fontsize=8); save(fig, 'pvalues')


def old_vs_shared_limits(observed):
    """Overlay native stored limits only; no resampling or missing-mass bridges."""
    old_run=P/'runs/20261030_sr6_cr4_scan'
    old={}
    for path in old_run.glob('m*/results.json'):
        q=read(path)
        if q and 6<=q.get('mass',-1)<=9:old[q['mass']]=q
    if not old:return
    new={q['mass']:q for q in observed if 6<=q['mass']<=9}
    masses=[k/10 for k in range(60,91)]
    rows=[]
    for mass in masses:
        a=old.get(mass,{});b=new.get(mass,{})
        la=a.get('limits') or {};lb=b.get('limits') or {}
        old_ready=finite(la.get('observed')) and finite(la.get('expected_50'))
        new_ready=finite(lb.get('observed')) and finite(lb.get('expected_50'))
        rows.append(dict(mass=mass,
            old_observed=la.get('observed'),old_expected_50=la.get('expected_50'),
            new_observed=lb.get('observed'),new_expected_50=lb.get('expected_50'),
            old_ready=old_ready,new_ready=new_ready,
            old_numerical_flag=old_ready and not a.get('numerical_ok',False),
            new_numerical_flag=new_ready and not (observed_fit_ok(b) and b.get('observed_limit_check',{}).get('verified',False)),
            new_fit_numerical_ok=observed_fit_ok(b),
            new_observed_crossing_verified=b.get('observed_limit_check',{}).get('verified',False),
            old_expected_bands_provisional=old_ready and not a.get('expected_bands_verified',False),
            new_expected_bands_provisional=new_ready and not b.get('expected_bands_verified',b.get('expected_validated',False)),
            new_source_flags=b.get('flags',[]),
            old_source=str(old_run/f'm{mass:.1f}/results.json'),
            new_source=str(R/f'observed/m{mass:.1f}/results.json')))
    metadata=dict(mass_range=[6,9],mass_step=.1,unit='signal strength r',
        old_configuration='SR 6–9.093333 + CR 9.093333–10.506667; eight Upsilon means free in the simultaneous fit; SR-local Cheb6/Bern6 continuum.',
        new_configuration='SR low+high and CR; eight Upsilon means estimated B-only in CR then fixed; common global Cheb6/Bern6 intensity across the two SRs.',
        interpretation='Simultaneous changes in data regions, Upsilon-mean treatment and continuum constraints. The difference cannot be attributed uniquely to adding the high SR.',
        expected_bands='Native provisional bands omitted for readability; native observed and median expected values shown unchanged.',
        missing_mass_policy='Absent native values remain null; plotted curves break across missing hypotheses. No interpolation.',
        numerical_marker_policy='Preserve the existing reports: new markers indicate FD or observed crossing flags; historical markers use its saved numerical_ok flag. No new selection or quality gate.',
        rows=rows)
    signature=hashlib.sha256(json.dumps(metadata,sort_keys=True).encode()).hexdigest()
    cache=R/'old_vs_shared_limits_plot_cache.json'
    if read(cache,{}).get('signature')==signature and all((R/f'old_vs_shared_limits.{ext}').exists() for ext in ['png','pdf','json','csv']):return
    write_json(R/'old_vs_shared_limits.json',metadata)
    csv_write(R/'old_vs_shared_limits.csv',rows)
    fig,ax=plt.subplots(figsize=(12,6.7))
    for prefix,color,label in [('old','#0072B2','Storico'),('new','#7B3294','Nuovo')]:
        for key,ls,marker,meaning in [('observed','-','o','osservato'),('expected_50','--',None,'atteso mediano')]:
            values=[q[f'{prefix}_{key}'] if finite(q[f'{prefix}_{key}']) else np.nan for q in rows]
            ax.plot(masses,values,ls=ls,marker=marker,ms=3,lw=1.7,color=color,label=f'{label}: {meaning}')
        bad=[q for q in rows if q[f'{prefix}_numerical_flag']]
        for key in ['observed','expected_50']:
            xx=[q['mass'] for q in bad if finite(q[f'{prefix}_{key}'])]
            yy=[q[f'{prefix}_{key}'] for q in bad if finite(q[f'{prefix}_{key}'])]
            ax.plot(xx,yy,ls='none',marker='x' if prefix=='new' else 'D',ms=7,
                    color='#C00000' if prefix=='new' else '#C28E00',markerfacecolor='none',mew=1.5,
                    label=('Nuovo: flag fit/crossing osservato' if prefix=='new' else 'Storico: flag numerici') if key=='observed' and bad else None)
    nold=sum(q['old_ready'] for q in rows);nnew=sum(q['new_ready'] for q in rows)
    fig.suptitle(f'Limiti 95% CLs: confronto di due configurazioni | storico {nold}/31, nuovo {nnew}/31\nStorico: SR bassa + CR, μΥ libere | Nuovo: SR bassa + alta + CR, μΥ fissate dalla sola CR',fontsize=12,y=.98)
    ax.set_xlim(5.97,9.03);ax.set_yscale('log');ax.set_xlabel('Massa ipotizzata [GeV]')
    ax.set_ylabel('Limite superiore 95% CLs su r');ax.grid(alpha=.2,which='both')
    ax.legend(fontsize=9,ncol=2,loc='upper left')
    fig.text(.5,.065,'Il nuovo fit usa un continuo globale comune alle due SR: dati inclusi, trattamento delle μΥ e vincoli cambiano insieme.\nLa variazione non isola una sola causa. Bande native provvisorie omesse per leggibilità; punti mancanti non interpolati.',ha='center',fontsize=9)
    fig.subplots_adjust(top=.83,bottom=.2,left=.09,right=.98)
    for ext in ['png','pdf']:
        target=R/f'old_vs_shared_limits.{ext}';temp=R/f'.old_vs_shared_limits_{os.getpid()}.{ext}'
        fig.savefig(temp,dpi=160);temp.replace(target)
    plt.close(fig);write_json(cache,dict(signature=signature))


def regional_deviance(observed):
    """Decompose existing fitted data deviances; no new fits or regional tests."""
    rows=[]
    for q in observed:
        mass=q['mass']; path=R/f'observed/m{mass:.1f}/bins.json'
        bins=read(path, {})
        if not all(stage in bins for stage in ['b','s']): continue
        row=dict(mass=mass,r=q.get('r'),Z=q.get('Z'),fit_significance_numerical_ok=observed_fit_ok(q) and q.get('significance_consistent',False))
        for region,_ in REGIONS:
            for stage in ['b','s']:
                values=[bins[stage].get(f'{region}_ch{i}',{}).get('deviance') for i in range(1,5)]
                row[f'D_{stage}_{region}']=sum(values) if all(finite(v) for v in values) else None
            a,b=row[f'D_b_{region}'],row[f'D_s_{region}']
            row[f'delta_{region}']=a-b if finite(a) and finite(b) else None
        parts=[row[f'delta_{region}'] for region,_ in REGIONS]
        row['delta_data_sum']=sum(parts) if all(finite(v) for v in parts) else None
        a,b=[q.get('fit',{}).get(stage,{}).get('aux_penalty') for stage in ['b','s']]
        row['delta_aux']=a-b if finite(a) and finite(b) else None
        row['delta_total']=q.get('delta_deviance_B_minus_SB')
        row['closure_residual']=(row['delta_data_sum']+row['delta_aux']-row['delta_total']) if all(finite(row[k]) for k in ['delta_data_sum','delta_aux','delta_total']) else None
        row['decomposition_consistent']=finite(row['closure_residual']) and abs(row['closure_residual'])<1e-5
        row['bins_source']=str(path.relative_to(R))
        rows.append(row)
    return rows


def regional_deviance_plot(rows):
    if not rows:return
    fig,ax=plt.subplots(figsize=(12,6.5))
    lookup={q['mass']:q for q in rows}
    for key,label,color in [('delta_SRL','SR bassa','#0072B2'),('delta_CR','CR','#7B3294'),('delta_SRH','SR alta','#D55E00')]:
        y=[lookup.get(m,{}).get(key,np.nan) if lookup.get(m,{}).get('decomposition_consistent') else np.nan for m in MASSES]
        ax.plot(MASSES,y,'o-',ms=3,color=color,label=label)
    total=[lookup.get(m,{}).get('delta_total',np.nan) if lookup.get(m,{}).get('decomposition_consistent') else np.nan for m in MASSES]
    ax.plot(MASSES,total,'--',color='black',lw=1.3,label='Totale, inclusi vincoli ausiliari')
    bad=[q for q in rows if not q['fit_significance_numerical_ok'] and finite(q['delta_total']) and q['decomposition_consistent']]
    ax.scatter([q['mass'] for q in bad],[q['delta_total'] for q in bad],marker='x',color='red',s=45,label='Flag numerici del fit o della significanza')
    ax.axhline(0,color='gray',lw=.8)
    ax.set_ylabel(r'$\Delta D = D_{B}-D_{S+B}$')
    fig.suptitle('Dove cambia la descrizione dei dati: contributi alla devianza\nSomma delle quattro ere; valori positivi = miglioramento del fit S+B nella regione',fontsize=12,y=.99)
    shade(ax)
    handles,labels=ax.get_legend_handles_labels()
    fig.legend(handles,labels,fontsize=8,ncol=3,loc='upper center',bbox_to_anchor=(.5,.91),frameon=False)
    fig.text(.5,.005,'Sono contributi di devianza del fit simultaneo, non p-value regionali calibrati. Per r̂ < 0 il miglioramento riguarda un deficit.',ha='center',fontsize=8)
    fig.subplots_adjust(top=.78,bottom=.14);save(fig,'regional_delta_deviance')


def bias_plots(summaries):
    if not summaries:
        return
    for name, title, ylabel in [
        ('bias_normalized', 'Bias medio / dispersione empirica — condizionato ai fit numericamente validi', r'$(\overline{\hat r}-r_{true})/SD(\hat r)$'),
        ('bias', 'Bias medio — fit numericamente validi e tutti i risultati finiti', r'$\overline{\hat r}-r_{true}$'),
        ('failure_fraction', 'Frazione di risultati numericamente non validi — nessun toy sostituito', 'Risultati non validi / tentativi'),
        ('cr_mean_boundary_fraction', 'Medie delle Υ al bordo nel fit iniziale della CR', 'Toy con almeno una μ Υ al bordo / tentativi'),
        ('edm_failure_fraction', 'Controllo numerico EDM — soglia finita < 0,01', 'EDM non valido o non disponibile / tentativi')]:
        if name=='bias_normalized' and not any(finite(s['normalized_bias']['value']) for s in summaries):
            for ext in ['png', 'svg']:
                (R/f'{name}.{ext}').unlink(missing_ok=True)
            continue
        fig, axes = plt.subplots(1, 2, figsize=(14, 4.5), sharex=True)
        for li, (label, ax) in enumerate(zip(INJECTIONS, axes)):
            ss = [s for s in summaries if s['injection'] == label]
            xx, yy, ee = [], [], []
            for s in ss:
                if name == 'bias_normalized':
                    val, err = s['normalized_bias']['value'], s['normalized_bias']['se']
                elif name == 'bias':
                    val, err = s['conditional']['mean'], s['conditional']['se']
                elif name == 'failure_fraction':
                    val, err = s['n_invalid']/s['n'], None
                elif name == 'cr_mean_boundary_fraction':
                    val, err = s['n_CR_mu_boundary']/s['n'], None
                else:
                    val, err = s['n_EDM_failure']/s['n'], None
                if finite(val): xx.append(s['mass']); yy.append(val); ee.append(err)
            if xx:
                if name in ['failure_fraction', 'cr_mean_boundary_fraction', 'edm_failure_fraction']:
                    ax.plot(xx, yy, 'o', color=COLORS[li], ms=3)
                else:
                    ax.errorbar(xx, yy, yerr=[e if finite(e) else 0 for e in ee], fmt='o', color=COLORS[li], ms=3, capsize=2)
                points=dict(zip(xx,yy))
                ax.plot(MASSES,[points.get(m,np.nan) for m in MASSES],'-',color=COLORS[li],lw=1)
            if name == 'bias':
                ax.scatter([s['mass'] for s in ss], [s['all_finite']['mean'] for s in ss], marker='x', color='black', s=30, label='Tutti i valori finiti')
                ax.legend(fontsize=8)
            if name == 'edm_failure_fraction':
                ax.plot([s['mass'] for s in ss], [s['n_CR_EDM_failure']/s['n'] for s in ss], 'x', color='black', ms=5, label='Fit iniziale CR')
                ax.plot([], [], 'o', color=COLORS[li], label='Fit simultaneo')
                ax.legend(fontsize=8)
            if name not in ['failure_fraction', 'cr_mean_boundary_fraction', 'edm_failure_fraction']: ax.axhline(0, color='black', lw=.7, ls='--')
            else: ax.set_ylim(-.02, 1.02)
            nkey='n' if name in ['failure_fraction', 'cr_mean_boundary_fraction', 'edm_failure_fraction'] else 'n_valid'
            sample_counts=[s[nkey] for s in ss]
            sample_label=f'{min(sample_counts)}–{max(sample_counts)}' if sample_counts else '0'
            nlabel='tentativi' if nkey=='n' else 'stime valide'
            ax.set_title(('r vero = 0' if label == 'r0' else 'r vero = limite atteso mediano')+f' | {nlabel}: {sample_label}')
            ax.set_ylabel(ylabel); shade(ax)
        progress=sum(s['n'] for s in summaries)
        note='; barre = errore standard bootstrap' if name=='bias_normalized' else ''
        fig.suptitle(title+f'\nGeneratore Chebyshev7; {progress}/6100 toy, target 50 per massa e iniezione'+note, fontsize=12)
        fig.tight_layout(); save(fig, name)


def gof_plot(rows):
    rows = [q for q in rows if finite(q['p'])]
    if not rows:
        return
    fig, axes = plt.subplots(2,1,figsize=(12,9),sharex=True,gridspec_kw={'height_ratios':[1.4,1]})
    ax,bx=axes
    for state,color,marker,label in [
        ('checked_valid','#009E73','o','Osservato: minimo verificato'),
        ('checked_failed','#D55E00','X','Osservato: controllo fallito'),
        ('unchecked','#9B7400','s','Osservato: non ancora verificato')]:
        rr=[q for q in rows if q['observed_qc']==state]
        if not rr: continue
        x=np.array([q['mass'] for q in rr]); pp=np.array([q['p'] for q in rr])
        low=np.array([q['ci95_low'] for q in rr]); high=np.array([q['ci95_high'] for q in rr])
        ax.errorbar(x,pp,yerr=[pp-low,high-pp],fmt=marker,ms=4,capsize=2,color=color,label=label)
    points={q['mass']:q['p'] for q in rows}
    ax.plot(MASSES,[points.get(m,np.nan) for m in MASSES],'-',color='gray',lw=.8,alpha=.6)
    ax.axhline(.05, ls='--', color='red', lw=.8)
    ax.set_ylim(-.02,1.02); ax.set_ylabel('p grezzo; intervallo Monte Carlo 95%')
    nmin=min(q['n_finite'] for q in rows); nmax=max(q['n_finite'] for q in rows)
    fig.suptitle(f'GoF combinato saturated, r libero | {nmin}–{nmax} score finiti/massa, target 100\nTutti gli score finiti conservati nel p-value; campioni condivisi fra masse',fontsize=12)
    shade(ax); ax.set_xlabel(''); ax.legend(fontsize=8,ncol=3,loc='lower right')
    for q in rows:
        if finite(q['p_qc_low']):
            bx.plot([q['mass'],q['mass']],[q['p_qc_low'],q['p_qc_high']],color='#CC8800',lw=2,alpha=.65)
            bx.plot(q['mass'],q['p'],'o',color='black',ms=2)
    bx.plot([],[],color='#CC8800',lw=2,label='Limiti conservativi: toy non risolti sotto/sopra la soglia osservata')
    bx.plot([],[],'o',color='black',ms=3,label='p grezzo da tutti gli score finiti')
    bx.set_ylim(-.02,1.02); bx.set_ylabel('Intervallo dei toy non risolti')
    shade(bx); bx.legend(fontsize=8,loc='lower right')
    fig.text(.5,.015,'Limiti riferiti ai soli esiti completati e alla soglia osservata salvata; non sono intervalli MC.\nUn osservato con controllo fallito/non eseguito non consente una conclusione GoF validata.',ha='center',fontsize=9)
    fig.subplots_adjust(top=.89,bottom=.13,hspace=.18); save(fig,'gof')


def source_fd(mass):
    base = R/f'observed/m{mass:.1f}'
    candidates = list(base.glob('joint/fitdiagnostics/fitDiagnostics*.root')) + list(base.glob('fitdiagnostics/fitDiagnostics*.root'))
    return candidates[0] if candidates else None


def plot_full_spectrum(mass, q, curves, source, target):
    """Display the joint model; never bridge the SR polynomial through the CR."""
    image = target/'fit_full_log.png'
    if not all(f'{reg}_ch{i}' in curves['s'] for reg, _ in REGIONS for i in range(1, 5)):
        return None
    result_path = R/f'observed/m{mass:.1f}/results.json'
    result_time = result_path.stat().st_mtime if result_path.exists() else 0
    if image.exists() and image.stat().st_mtime > max(source.stat().st_mtime, result_time):
        return str(image.relative_to(R))
    fig = plt.figure(figsize=(14, 10))
    outer = fig.add_gridspec(2, 2, hspace=.32, wspace=.2)
    negative_signal = False
    for i, (era, cell) in enumerate(zip(ERAS, outer), 1):
        sub = cell.subgridspec(2, 1, height_ratios=[3, 1], hspace=.06)
        ax, rx = fig.add_subplot(sub[0]), fig.add_subplot(sub[1])
        ymax = 1.
        for region, _ in REGIONS:
            channel = f'{region}_ch{i}'
            item = curves['s'][channel]
            x = np.asarray(item['mass'], dtype=float)
            data = np.asarray(item['data'], dtype=float)
            pred = np.asarray(item['total'], dtype=float)
            width = item['bin_width']; edges = np.r_[x-width/2, x[-1]+width/2]
            ymax = max(ymax, float(np.max(data+np.sqrt(np.maximum(data, 0)))))
            positive_data = data > 0
            ax.errorbar(x[positive_data], data[positive_data], yerr=np.sqrt(data[positive_data]), fmt='.', color='black', ms=1.8, lw=.4, label='Dati' if region=='SRL' else None)
            ax.stairs(np.where(pred>0, pred, np.nan), edges, baseline=None, color='#009E73', lw=1.6, label='Totale S+B' if region=='SRL' else None)
            bonly = curves.get('b', {}).get(channel)
            if bonly:
                btotal = np.asarray(bonly['total'], dtype=float)
                ax.stairs(np.where(btotal>0, btotal, np.nan), edges, baseline=None, color='gray', ls='--', lw=1, label='Solo fondo' if region=='SRL' else None)
            for name, raw in item['components'].items():
                values = np.asarray(raw, dtype=float)
                if 'upsilon1s' in name: color, title = '#E69F00', 'Υ(1S)'
                elif 'upsilon2s' in name: color, title = '#CC79A7', 'Υ(2S)'
                elif 'Zd' in name:
                    color, title = '#D55E00', 'Segnale DCB (parte positiva)'
                    negative_signal |= bool(np.any(values < 0))
                else:
                    color, title = ('#7B3294', 'Continuo CR4') if region=='CR' else ('#0072B2', 'Continuo comune SR')
                show_label = region=='SRL' or title=='Continuo CR4'
                ax.stairs(np.where(values>0, values, np.nan), edges, baseline=None, color=color, lw=1.1, label=title if show_label else None)
            residual = np.full_like(pred, np.nan); ok = np.isfinite(pred) & (pred > 0)
            residual[ok] = (data[ok]-pred[ok])/np.sqrt(pred[ok])
            rx.plot(x, residual, '.', color='#009E73', ms=2)
        for axis in [ax, rx]:
            axis.axvspan(9.093333333333334, 10.506666666666668, color='#999999', alpha=.13)
            axis.axvline(9.093333333333334, color='#999999', lw=.6)
            axis.axvline(10.506666666666668, color='#999999', lw=.6)
            axis.set_xlim(6, 12); axis.grid(alpha=.12)
        ax.set_title(era); ax.set_yscale('log'); ax.set_ylim(.15, ymax*1.9)
        ax.set_ylabel('Eventi / 13,33 MeV'); ax.tick_params(labelbottom=False)
        rx.axhline(0, color='gray', lw=.7); rx.set_ylabel('(D−M)/√M'); rx.set_xlabel('Massa [GeV]')
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, ncol=4, fontsize=8.2, loc='upper center', bbox_to_anchor=(.5, .925), frameon=False)
    fig.suptitle(f'mX={mass:.1f} GeV | Fit simultaneo dell’intero spettro\nSR bassa e alta: stessa intensità Cheb6/Bern6; fascia grigia: CR con continuo Cheb4 indipendente', fontsize=13)
    fig.subplots_adjust(top=.85, bottom=.12)
    fit = q.get('fit', {})
    line = f'r̂={fmt(q.get("r"))}; Z locale={fmt(q.get("Z"))}; FD B/S status={fit.get("b", {}).get("status", "—")}/{fit.get("s", {}).get("status", "—")}, cov={fit.get("b", {}).get("covQual", "—")}/{fit.get("s", {}).get("covQual", "—")}'
    if fit_notes(q): line += ' | '+fit_notes(q)
    fig.text(.5, .04, line, ha='center', fontsize=9)
    log_note = 'Scala log: i bin a zero restano nei residui; componenti non positive non visualizzate.'
    if negative_signal: log_note += ' La DCB stimata è negativa: vedere i plot regionali lineari.'
    fig.text(.5, .018, log_note, ha='center', fontsize=8)
    fig.savefig(image, dpi=145); fig.savefig(image.with_suffix('.svg')); plt.close(fig)
    return str(image.relative_to(R))


def plot_bins(mass, q, errors):
    source = R/f'observed/m{mass:.1f}/bins.json'
    curves = read(source)
    if not curves or not curves.get('s'):
        return None
    target = R/f'fits/m{mass:.1f}'
    target.mkdir(parents=True, exist_ok=True)
    links = []
    result_path = R/f'observed/m{mass:.1f}/results.json'
    result_time = result_path.stat().st_mtime if result_path.exists() else 0
    for region, label in REGIONS:
        channels = [f'{region}_ch{i}' for i in range(1, 5)]
        if not all(ch in curves['s'] for ch in channels):
            continue
        image = target/f'fit_{region}.png'
        links.append(str(image.relative_to(R)))
        if image.exists() and image.stat().st_mtime > max(source.stat().st_mtime, result_time):
            continue
        fig = plt.figure(figsize=(13, 9))
        outer = fig.add_gridspec(2, 2, hspace=.32, wspace=.2)
        for era, channel, cell in zip(ERAS, channels, outer):
            item = curves['s'][channel]
            x = np.asarray(item['mass'], dtype=float)
            data = np.asarray(item['data'], dtype=float)
            pred = np.asarray(item['total'], dtype=float)
            width = item['bin_width']
            edges = np.r_[x-width/2, x[-1]+width/2]
            sub = cell.subgridspec(2, 1, height_ratios=[3, 1], hspace=.06)
            ax, rx = fig.add_subplot(sub[0]), fig.add_subplot(sub[1])
            ax.errorbar(x, data, yerr=np.sqrt(np.maximum(data, 0)), fmt='.', color='black', ms=2, lw=.5, label='Dati')
            ax.stairs(pred, edges, baseline=None, color='#009E73', lw=1.5, label='Totale S+B')
            bonly = curves.get('b', {}).get(channel)
            if bonly:
                ax.stairs(bonly['total'], edges, baseline=None, color='gray', ls='--', label='Solo fondo')
            for name, values in item['components'].items():
                if 'upsilon1s' in name: color, title = '#E69F00', 'Υ(1S)'
                elif 'upsilon2s' in name: color, title = '#CC79A7', 'Υ(2S)'
                elif 'Zd' in name: color, title = '#D55E00', 'Segnale DCB'
                else: color, title = '#0072B2', 'Continuo CR4' if region == 'CR' else 'Continuo comune SR'
                ax.stairs(values, edges, baseline=None, color=color, label=title)
            ax.set_title(era); ax.set_ylabel(f'Eventi / {width*1000:.2f} MeV'); ax.tick_params(labelbottom=False)
            ax.grid(alpha=.15); ax.legend(fontsize=7, ncol=2)
            residual = np.full_like(pred, np.nan)
            ok = np.isfinite(pred) & (pred > 0)
            residual[ok] = (data[ok]-pred[ok])/np.sqrt(pred[ok])
            rx.axhline(0, color='gray', lw=.7); rx.plot(x, residual, '.', color='#009E73', ms=2.5)
            rx.set_ylabel('(D−M)/√M'); rx.set_xlabel('Massa [GeV]'); rx.grid(alpha=.15)
            ax.set_xlim(edges[0], edges[-1]); rx.set_xlim(edges[0], edges[-1])
        fig.suptitle(f'mX={mass:.1f} GeV | {label} | Envelope comune Cheb6/Bern6, CR4\nΥ μ dalla sola CR, poi fisse; DCB in tutte e tre le regioni', fontsize=13)
        fig.subplots_adjust(top=.89, bottom=.11)
        fit = q.get('fit', {})
        annotation = f'r̂={fmt(q.get("r"))}; Z locale={fmt(q.get("Z"))}'
        if fit:
            annotation += f'; FD B/S status={fit.get("b", {}).get("status", "—")}/{fit.get("s", {}).get("status", "—")}, cov={fit.get("b", {}).get("covQual", "—")}/{fit.get("s", {}).get("covQual", "—")}'
        if fit_notes(q): annotation += ' | '+fit_notes(q)
        fig.text(.5, .025, annotation, ha='center', fontsize=9)
        fig.savefig(image, dpi=140); fig.savefig(image.with_suffix('.svg')); plt.close(fig)
    full = plot_full_spectrum(mass, q, curves, source, target)
    if full:
        links.insert(0, full)
    return links


def plot_fit(mass, q, errors):
    exact_curves = plot_bins(mass, q, errors)
    if exact_curves is not None:
        return exact_curves
    source = source_fd(mass)
    if source is None:
        return []
    target = R/f'fits/m{mass:.1f}'
    target.mkdir(parents=True, exist_ok=True)
    links = []
    result_time = (R/f'observed/m{mass:.1f}/results.json').stat().st_mtime if q else 0
    try:
        with uproot.open(source) as root:
            available = set(root['shapes_fit_s'].keys(cycle=False))
            for region, label in REGIONS:
                image = target/f'fit_{region}.png'
                channel_names = [f'{region}_ch{i}' for i in range(1, 5)]
                if not all(c in available for c in channel_names):
                    continue
                links.append(str(image.relative_to(R)))
                if image.exists() and image.stat().st_mtime > max(source.stat().st_mtime, result_time):
                    continue
                fig = plt.figure(figsize=(13, 9))
                outer = fig.add_gridspec(2, 2, hspace=.32, wspace=.2)
                for era, channel, cell in zip(ERAS, channel_names, outer):
                    d = root[f'shapes_fit_s/{channel}']
                    pred, edges = d['total'].to_numpy(); bw = np.diff(edges); x = (edges[1:]+edges[:-1])/2
                    pred = pred*bw; data = d['data'].values('y')*bw
                    sub = cell.subgridspec(2, 1, height_ratios=[3, 1], hspace=.06)
                    ax, rx = fig.add_subplot(sub[0]), fig.add_subplot(sub[1])
                    ax.errorbar(x, data, yerr=np.sqrt(np.maximum(data, 0)), fmt='.', color='black', ms=2, lw=.5, label='Dati')
                    ax.stairs(pred, edges, baseline=None, color='#009E73', lw=1.5, label='Totale S+B')
                    if f'shapes_fit_b/{channel}/total' in root:
                        ax.stairs(root[f'shapes_fit_b/{channel}/total'].values()*bw, edges, baseline=None, color='gray', ls='--', label='Solo fondo')
                    for process, color, pname in [('dy', '#0072B2', 'Continuo comune SR' if region != 'CR' else 'Continuo CR4'), ('upsilon1s', '#E69F00', 'Υ(1S)'), ('upsilon2s', '#CC79A7', 'Υ(2S)'), ('Zd', '#D55E00', 'Segnale DCB')]:
                        if process in d:
                            ax.stairs(d[process].values()*bw, edges, baseline=None, color=color, label=pname)
                    ax.set_title(era); ax.set_ylabel('Eventi / bin'); ax.tick_params(labelbottom=False)
                    ax.grid(alpha=.15); ax.legend(fontsize=7, ncol=2)
                    rx.axhline(0, color='gray', lw=.7)
                    rx.plot(x, (data-pred)/np.sqrt(np.maximum(pred, 1e-10)), '.', color='#009E73', ms=2.5)
                    rx.set_ylabel('(D−M)/√M'); rx.set_xlabel('Massa [GeV]'); rx.grid(alpha=.15)
                    ax.set_xlim(edges[0], edges[-1]); rx.set_xlim(edges[0], edges[-1])
                fig.suptitle(f'mX={mass:.1f} GeV | {label} | Envelope comune Cheb6/Bern6, CR4\nΥ μ dalla sola CR, poi fisse; DCB in tutte e tre le regioni', fontsize=13)
                fig.subplots_adjust(top=.89, bottom=.11)
                annotation = f'r̂={fmt(q.get("r"))}; Z locale={fmt(q.get("Z"))}' if q else 'Fit disponibile; limiti/GoF in elaborazione'
                fig.text(.5, .025, annotation, ha='center', fontsize=9)
                fig.savefig(image, dpi=140); fig.savefig(image.with_suffix('.svg')); plt.close(fig)
    except (OSError, KeyError, ValueError) as exc:
        errors.append(dict(path=str(source.relative_to(R)), reason=str(exc)))
    return links


def report():
    R.mkdir(parents=True, exist_ok=True)
    observed, bias, gof, profiles = collect()
    bias_summaries = summarize_bias(bias)
    cr_shifts,cr_shift_toys,cr_shift_provenance=summarize_cr_mean_shifts(bias)
    cr_mean_shift_plots(cr_shifts)
    write_json(R/'cr_mean_truth_provenance.json',cr_shift_provenance)
    regional=regional_deviance(observed)
    observed_plots(observed); old_vs_shared_limits(observed); bias_plots(bias_summaries); gof_plot(gof)
    regional_deviance_plot(regional)
    errors, gallery = [], {}
    lookup = {q['mass']: q for q in observed}
    for mass in MASSES:
        links = plot_fit(mass, lookup.get(mass, {}), errors)
        if links: gallery[mass] = links
    counts = dict(observed=sum(q.get('status')=='complete' for q in observed), observed_partial=len(observed), observed_target=61, fit_points=len(gallery),
                  fitdiagnostics=sum(len(q.get('fit', {}))==2 for q in observed),
                  significance=sum(finite(q.get('Z')) for q in observed),
                  limits=sum(len(q.get('limits', {}))==6 for q in observed),
                  gof_observed=sum(finite(q.get('gof_observed')) for q in observed),
                  gof_observed_diagnostics=sum(q['observed_diagnostic_completed'] for q in gof),
                  gof_observed_diagnostics_target=61,
                  gof_observed_checked_valid=sum(q['observed_qc']=='checked_valid' for q in gof),
                  gof_observed_checked_failed=sum(q['observed_qc']=='checked_failed' for q in gof),
                  gof_observed_unchecked=sum(q['observed_qc']=='unchecked' for q in gof),
                  bias=len(bias), bias_target=6100, bias_valid=sum(q['valid'] for q in bias),
                  bias_EDM_high=sum(finite(q['edm']) and q['edm']>=.01 for q in bias),
                  bias_EDM_missing=sum(not finite(q['edm']) for q in bias),
                  bias_CR_EDM_high=sum(finite(q['CR_edm']) and q['CR_edm']>=.01 for q in bias),
                  bias_CR_EDM_missing=sum(not finite(q['CR_edm']) for q in bias),
                  control_profiles=len(profiles), control_profiles_target=122,
                  control_profiles_valid=sum(q['profile_valid'] for q in profiles),
                  control_profiles_skipped=sum(q['outcome']=='skipped_invalid_or_missing_anchor' for q in profiles),
                  control_profiles_failed=sum(q['outcome']=='failed_or_unavailable_interval' for q in profiles),
                  control_profiles_lower_minimum=sum(q['profile_found_lower_minimum'] for q in profiles),
                  gof_scores=sum(q['n_attempted'] for q in gof), gof_target=6100,
                  gof_finite_scores=sum(q['n_finite'] for q in gof),
                  gof_missing_scores=sum(q['n_nonfinite'] for q in gof),
                  gof_points_100_outcomes=sum(q['n_attempted']==100 for q in gof),
                  gof_calibrated_points=sum(q['calibrated_conclusion_available'] for q in gof),
                  gof_nominal_diagnostics=sum(q['n_nominal_checked'] for q in gof), gof_nominal_diagnostics_target=6100,
                  gof_nominal_valid=sum(q['n_nominal_valid'] for q in gof),
                  gof_nominal_failed=sum(q['n_nominal_failed'] for q in gof),
                  gof_nominal_unchecked=sum(q['n_nominal_unchecked'] for q in gof))
    complete = counts['observed']==61 and counts['bias']==6100 and counts['gof_scores']==6100 and counts['control_profiles']==122 and counts['gof_observed_diagnostics']==61 and counts['gof_nominal_diagnostics']==6100
    now = datetime.datetime.now(ZoneInfo('Europe/Rome')).isoformat(timespec='seconds')
    summary = dict(updated_at=now, completed=complete, counts=counts, observed=observed,
                   bias=bias_summaries, CR_mean_shifts=cr_shifts, gof=gof, control_profiles=profiles, regional_deviance=regional, rendering_errors=errors,
                   interpretation={'bias': 'Conditional on numerical validity; all finite results and failures retained.',
                                   'coverage': 'Not measured. Profiles only on predetermined control toys.',
                                   'point_validity': 'Finite EDM <0.01 in the CR and simultaneous point fits, in addition to the worker status, positivity and boundary checks. Hesse covariance quality is a separate interval diagnostic.',
                                   'control_profiles': '122 preselected outcomes are required for completion. Failed or skipped-invalid-anchor controls count as completed documented outcomes, not valid intervals. Primary point estimates remain immutable.',
                                   'CR_mean_plugin': 'Requested procedure: estimate the eight Upsilon means in the B-only fit to the same CR, then fix them in the joint fit. Their plug-in estimation uncertainty is not profiled in limits or significance. Recalibrating every toy repeats the preliminary stage but is not equivalent to profiling mean nuisance parameters or universal validation.',
                                   'local_significance': 'Asymptotic local significance; no look-elsewhere correction for scanning masses.',
                                   'CR_mu_boundary': 'Reported separately; a converged CR mean at its bound is not automatically excluded by the current point-validity policy.',
                                   'injected_CR_contamination': 'Injected-signal ensembles also measure absorption into the B-only CR calibration step; bias is not attributable solely to continuum degree.',
                                   'gof': 'Conditional parametric bootstrap from the fitted B model, not an exact composite-null p-value. r profiled; 100 B-only generated samples reused across 61 mass tests. Primary p values use all finite raw scores. The score-parity-validated diagnostic runtime captures the nominal state before saturated minimization. QC states: checked-valid, checked-failed, unchecked; no covariance3 requirement. Invalid or unchecked observed nominal minima preclude a validated GoF conclusion.',
                                   'gof_qc_bounds': 'Bounds assign unresolved completed toy outcomes entirely below/above the stored observed score. Resolved requires nominal QC, finite score, valid CR with 0<=EDM<0.01, and r away from its artificial fit bound. Bounds use completed outcomes as denominator, are conditional on the stored observed threshold, and are not MC confidence intervals.',
                                   'gof_completion': 'Requires 61 observed and 6100 toy nominal diagnostic outcomes, including documented failed audits, in addition to all scores, bias points and 122 control-profile outcomes.',
                                   'signal_above11': 'Original parametric DCB; previously fitted normalization extrapolation without an extrapolation nuisance.'})
    # Strip nonfinite numbers from source diagnostics before JSON export.
    def clean(value):
        if isinstance(value, float) and not math.isfinite(value): return None
        if isinstance(value, dict): return {k:clean(v) for k,v in value.items()}
        if isinstance(value, list): return [clean(v) for v in value]
        return value
    write_json(R/'report_summary.json', clean(summary))
    write_json(R/'progress.json', dict(updated_at=now, completed=complete, counts=counts, rendering_errors=errors))
    csv_write(R/'observed_results.csv', [dict(mass=q['mass'], status=q.get('status'), r=q.get('r'), r_error=q.get('r_error'), Z=q.get('Z'), p0=q.get('p0'), numerical_ok=bool(q.get('numerical_ok')) and observed_fit_ok(q), source_numerical_ok=q.get('numerical_ok'), EDM_B=q.get('fit', {}).get('b', {}).get('edm'), EDM_SB=q.get('fit', {}).get('s', {}).get('edm'), **q.get('limits', {})) for q in observed])
    csv_write(R/'bias_results.csv', [dict(mass=q['mass'], injection=q['injection'], n=q['n'], n_valid=q['n_valid'], n_invalid=q['n_invalid'], n_CR_invalid=q['n_CR_invalid'], n_CR_mu_boundary=q['n_CR_mu_boundary'], n_r_boundary=q['n_r_boundary'], n_EDM_high=q['n_EDM_high'], n_EDM_missing=q['n_EDM_missing'], n_CR_EDM_high=q['n_CR_EDM_high'], n_CR_EDM_missing=q['n_CR_EDM_missing'], mean_bias=q['conditional']['mean'], estimator_sd=q['conditional']['sd'], normalized_bias=q['normalized_bias']['value'], normalized_bias_se=q['normalized_bias']['se'], all_finite_mean_bias=q['all_finite']['mean']) for q in bias_summaries])
    csv_write(R/'bias_toys.csv', [{k:v for k,v in q.items() if k!='CR_means'} for q in bias])
    csv_write(R/'cr_mean_shifts.csv',cr_shifts)
    csv_write(R/'cr_mean_shift_toys.csv',cr_shift_toys)
    csv_write(R/'gof_results.csv', gof)
    csv_write(R/'control_profiles.csv', profiles)
    csv_write(R/'regional_delta_deviance.csv',regional)
    rows=[]
    for mass in MASSES:
        q=lookup.get(mass, {}); limits=q.get('limits', {})
        plot_labels={'fit_full_log.png':'Spettro completo (log)', 'fit_SRL.png':'SR bassa', 'fit_CR.png':'CR', 'fit_SRH.png':'SR alta'}
        links=' · '.join(f'<a href="{p}">{plot_labels.get(p.rsplit("/",1)[-1], "Plot")}</a>' for p in gallery.get(mass, [])) or 'In attesa'
        fit_ok=observed_fit_ok(q)
        if not q: quality='In attesa'
        elif q.get('numerical_ok') and fit_ok: quality='OK numerico'
        elif q.get('status')!='complete' and fit_ok and q.get('significance_consistent', True): quality='Fit OK; scan parziale'
        else: quality='Flag / da verificare'
        rows.append(f'<tr><td>{mass:.1f}</td><td>{links}</td><td>{fmt(q.get("r"))}</td><td>{fmt(q.get("Z"))}</td><td>{fmt(q.get("p0"))}</td><td>{fmt(limits.get("observed"))}</td><td>{fmt(limits.get("expected_50"))}</td><td>{quality}</td></tr>')
    def panel(name, label, explanation):
        extra=''
        if name=='old_vs_shared_limits':
            extra='<p><a href="old_vs_shared_limits.pdf">PDF</a> · <a href="old_vs_shared_limits.csv">Valori CSV</a> · <a href="old_vs_shared_limits.json">Valori, configurazioni e provenienza JSON</a> · <a href="comparison_with_old_sr6_cr4.md">Confronto metodologico dettagliato</a></p>'
        if name=='limits':
            if (R/'comparison_with_old_sr6_cr4.md').exists():
                extra+='<p><a href="comparison_with_old_sr6_cr4.md">Perché questi limiti differiscono dal precedente scan SR bassa + CR</a>: confronto di dati inclusi, trattamento delle medie Υ e vincoli del continuo.</p>'
            audit=read(R/'limits_unconditional_minimum_audit.json', {})
            if audit:
                n=audit.get('completed_limit_files',0)
                refs=audit.get('positive_FD_references',0)+audit.get('physical_boundary_references',0)
                missing=', '.join(fmt(v,4) for v in audit.get('no_valid_FD_reference', [])) or 'nessuna'
                links=' · '.join(f'<a href="{name}">{text}</a>' for name,text in [('limits_unconditional_minimum_audit.md','Audit del minimo dei limiti: spiegazione'),('limits_unconditional_minimum_audit.json','Valori e provenienza JSON')] if (R/name).exists())
                extra+=f'<p>{links}</p><p>Audit parziale su {n} limiti: r̂ compatibile con {refs} riferimenti FitDiagnostics fisicamente ammissibili; per r̂ FD negativo il riferimento di q̃ è r=0. Masse senza riferimento FD valido: {missing} GeV. <b>Questo confronto di r̂ non certifica lo stesso minimo completo o globale:</b> parametri, indici, EDM e NLL assoluta del minimo AsymptoticLimits non sono registrati negli output conservati. Il crossing CLs verifica il limite rispetto al minimo nativo, non la globalità di quel minimo. Nessun nuovo criterio di selezione o flag viene introdotto.</p>'
        if (R/f'{name}.png').exists(): return f'<section><h2>{label}</h2><p>{explanation}</p><a href="{name}.png"><img src="{name}.png" alt="{label}"></a>{extra}</section>'
        return f'<section><h2>{label}</h2><p class="pending">In attesa dei primi risultati. {explanation}</p>{extra}</section>'
    panels=''.join([panel('limits', 'Limiti', 'Limiti 95% CLs con bande attese a 1 e 2 sigma. I quantili non verificati restano esplicitamente provvisori.'), panel('old_vs_shared_limits', 'Limiti 6–9 GeV: storico e nuovo fit sovrapposti', 'Stessa unità r e stessi assi. Linee continue: osservato; tratteggiate: atteso mediano. Nello storico si fittavano SR bassa + CR con medie Υ libere; nel nuovo si fittano SR bassa + alta + CR con medie Υ fissate dal fit alla sola CR e continuo globale comune. Queste modifiche avvengono insieme: la differenza non misura isolatamente l’effetto della SR alta. Le bande native provvisorie sono omesse per leggibilità. Nessun valore è interpolato o corretto; i marker mantengono visibili i flag numerici.'), panel('pvalues', 'p-value locale', 'Significanza e p-value locali asintotici: non includono la correzione per la ricerca su più masse (look-elsewhere effect).'), panel('regional_delta_deviance', 'Origine regionale del miglioramento del fit', 'ΔD = D_B − D_S+B, sommata sulle quattro ere. Le tre curve mostrano i contributi di SR bassa, CR e SR alta; il totale comprende anche i vincoli ausiliari. Non sono p-value regionali calibrati. La tabella dettagliata è sotto i grafici.'), panel('bias_normalized', 'Bias normalizzato', 'Bias medio diviso per la dispersione empirica dei toy. Medie condizionate ai fit numericamente validi; barre di errore bootstrap. Non è una significanza locale. I campioni r=0 sono condivisi fra masse, quelli con segnale iniettato sono indipendenti per massa. Il bias con segnale comprende anche la sua contaminazione nel passaggio iniziale di calibrazione B-only della CR: non misura soltanto il grado del polinomio.'), panel('bias', 'Bias e sensibilità ai fallimenti', 'Le croci mostrano la media di tutti i risultati finiti, compresi quelli con flag. Nessun toy viene sostituito.'), panel('failure_fraction', 'Fallimenti numerici', 'Le frazioni si riferiscono ai tentativi completati, senza rimpiazzo dei toy falliti.'), panel('cr_mean_boundary_fraction', 'Medie Υ al bordo nella CR', 'Frequenza di almeno una media Υ al limite del suo intervallo nel primo fit alla sola CR. È mostrata separatamente dai fallimenti: la politica attuale non esclude automaticamente un fit convergente con una media al bordo.'), panel('cr_mean_shift', 'Spostamento delle medie Υ nella calibrazione CR', 'Δμ medio rispetto alla verità del generatore, in MeV, per le otto medie e le due iniezioni. La selezione principale richiede il solo fit CR numericamente valido; non dipende dall’esito del successivo fit simultaneo. Le croci conservano i casi al bordo. Le dimensioni dei campioni e gli errori standard sono nel CSV. Questo diagnostico può descrivere la risposta della calibrazione B-only al segnale iniettato, ma non dimostra la causa degli eccessi osservati.'), panel('cr_mean_shift_all_finite', 'Medie Υ: confronto con tutti i valori finiti', 'Stessa scala cromatica, includendo anche le calibrazioni CR con flag. I CSV conservano inoltre il sottoinsieme con fit simultaneo numericamente valido. Per r=0, gli stessi 50 dataset e fit CR sono riutilizzati a ogni massa: i pannelli CR-valido e tutti-finiti saranno piatti per costruzione a campagna completa, non 61 validazioni indipendenti. Eventuali differenze parziali riflettono campioni incompleti; il sottoinsieme selezionato per validità del fit simultaneo può variare con la massa. Le iniezioni di segnale usano campioni distinti per massa.'), panel('edm_failure_fraction', 'Controllo EDM', 'La validità numerica delle stime centrali richiede EDM finito inferiore a 0,01 sia nella CR sia nel fit simultaneo, oltre agli altri controlli di stato, positività e bordi. Le tabelle distinguono EDM sopra soglia da EDM non disponibile. La qualità Hesse degli intervalli è un controllo separato.'), panel('gof', 'Goodness of fit', 'Combine saturated con r libero; bootstrap parametrico con generazione solo fondo dal modello B stimato. Il risultato è condizionato a quel modello generatore e non è un p-value esatto per un’ipotesi nulla composita con parametri ignoti. Il p-value usa tutti gli score finiti, inclusi quelli con flag, e mostra intervalli binomiali Monte Carlo al 95%. I campioni sono condivisi fra masse. Il runtime diagnostico verificato salva stato, EDM e parametri prima della minimizzazione saturata, senza cambiare lo score. Si distinguono esiti controllati validi, controllati falliti e non ancora controllati. Gli osservati non validi sono evidenziati. Il pannello inferiore assegna tutti i toy non risolti sotto/sopra lo score osservato salvato, usando come denominatore gli esiti completati: questi limiti non sono intervalli Monte Carlo e non validano una soglia osservata non affidabile.')])
    refresh='' if complete else '<meta http-equiv="refresh" content="60">'
    regional_table=''
    if regional:
        rr=[]
        for q in regional:
            quality='OK numerico' if q['fit_significance_numerical_ok'] else 'Flag fit/significanza'
            if not q['decomposition_consistent']:quality+='; esportazione da verificare'
            rr.append(f'<tr><td>{q["mass"]:.1f}</td><td>{fmt(q["r"])}</td><td>{fmt(q["delta_SRL"])}</td><td>{fmt(q["delta_CR"])}</td><td>{fmt(q["delta_SRH"])}</td><td>{fmt(q["delta_aux"])}</td><td>{fmt(q["delta_total"])}</td><td>{quality}</td></tr>')
        regional_table='<section><h2>Contributi regionali alla variazione di devianza</h2><p>ΔD = D<sub>B</sub> − D<sub>S+B</sub>, sommata sulle quattro ere in ciascuna regione. Un contributo positivo indica una migliore descrizione locale nel fit S+B; uno negativo indica un peggioramento compensato altrove. I parametri sono stimati simultaneamente: <b>questi numeri non sono test indipendenti né p-value regionali calibrati</b>. La somma dei tre contributi e della variazione dei vincoli ausiliari riproduce la variazione totale. Un r negativo descrive un deficit, anche se ΔD è positivo.</p><p><a href="regional_delta_deviance.csv">Tabella completa CSV, inclusi D prima/dopo e verifica della somma</a></p><div style="overflow:auto"><table><thead><tr><th>Massa</th><th>r̂</th><th>ΔD SR bassa</th><th>ΔD CR</th><th>ΔD SR alta</th><th>ΔD vincoli</th><th>ΔD totale</th><th>Qualità</th></tr></thead><tbody>'+''.join(rr)+'</tbody></table></div></section>'
    cr_shift_table=''
    if cr_shifts:
        sr=[]
        for mass in MASSES:
            for label in INJECTIONS:
                cell=[q for q in cr_shifts if q['mass']==mass and q['injection']==label]
                if not cell:continue
                bounds=sum(q['n_boundary'] for q in cell)
                sr.append(f'<tr><td>{mass:.1f}</td><td>{label}</td><td>{min(q["n_finite"] for q in cell)}–{max(q["n_finite"] for q in cell)}</td><td>{min(q["n_CR_valid"] for q in cell)}–{max(q["n_CR_valid"] for q in cell)}</td><td>{min(q["n_point_valid"] for q in cell)}–{max(q["n_point_valid"] for q in cell)}</td><td>{bounds}</td></tr>')
        cr_shift_table='<section><h2>Campioni e provenienza degli spostamenti Υ</h2><p><a href="cr_mean_shifts.csv">Statistiche per massa, iniezione e media: media, SD, errore standard e bordi</a> · <a href="cr_mean_shift_toys.csv">Valori di ogni toy e riferimento al fit CR</a> · <a href="cr_mean_truth_provenance.json">Verità generata e verifica degli hash</a> · <a href="bias/generator/results.json">Risultato originale del generatore</a></p><details><summary>Conteggi per massa e iniezione</summary><p>I conteggi indicano il minimo–massimo fra le otto medie; i bordi sono il numero di coppie toy–parametro al bordo, non il numero di toy distinti. Il sottoinsieme CR valido è distinto da quello del fit simultaneo valido. I valori al bordo restano nella media se numericamente validi.</p><table><thead><tr><th>Massa</th><th>Iniezione</th><th>Valori finiti / media</th><th>CR valide / media</th><th>Fit simultanei validi / media</th><th>Coppie toy–media al bordo</th></tr></thead><tbody>'+''.join(sr)+'</tbody></table></details></section>'
    profile_rows=[]
    outcome_labels={'valid_interval':'Intervallo valido', 'skipped_invalid_or_missing_anchor':'Saltato: stima/CR assente o non valida', 'failed_or_unavailable_interval':'Intervallo fallito/non disponibile'}
    for q in profiles:
        profile_rows.append(f'<tr><td>{q["mass"]:.1f}</td><td>{q["injection"]}</td><td>{outcome_labels[q["outcome"]]}</td><td>{fmt(q["point_r"])}</td><td>{fmt(q["profile_r"])}</td><td>[{fmt(q["interval_low"])}, {fmt(q["interval_high"])}]</td><td>{fmt(q["edm"])}</td><td>{"Sì" if q["profile_found_lower_minimum"] else "No"}</td><td><a href="{q["result"]}">JSON</a></td></tr>')
    profile_table='<section><h2>Profili sui controlli prestabiliti</h2><p>Un controllo sul toy 1 per ogni massa e iniezione: <b>'+str(counts['control_profiles'])+'/122 esiti documentati</b>, '+str(counts['control_profiles_valid'])+' intervalli validi, '+str(counts['control_profiles_failed'])+' falliti/non disponibili e '+str(counts['control_profiles_skipped'])+' saltati per stima/CR assente o non valida. Tutti gli esiti documentati contribuiscono al completamento della coda; solo quelli validi sono intervalli utilizzabili. Questi controlli non misurano la copertura dell’intero insieme. Il valore centrale primario resta quello dei due tentativi fissi; i minimi più bassi trovati durante il profilo vengono segnalati.</p>'
    if profile_rows:
        profile_table+='<div style="overflow:auto"><table><thead><tr><th>Massa</th><th>Iniezione</th><th>Esito</th><th>r punto</th><th>r profilo</th><th>Intervallo</th><th>EDM profilo</th><th>Minimo più basso</th><th>Dettagli</th></tr></thead><tbody>'+''.join(profile_rows)+'</tbody></table></div>'
    else: profile_table+='<p class="pending">Coda dei controlli ancora in attesa. Lo scan non sarà dichiarato completo prima di tutti i 122 esiti.</p>'
    profile_table+='</section>'
    qc_labels={'checked_valid':'Verificato valido','checked_failed':'CONTROLLO FALLITO','unchecked':'Non verificato'}
    gof_table=''
    if gof:
        gr=[]
        for g in gof:
            ci=f'[{fmt(g["ci95_low"],3)}, {fmt(g["ci95_high"],3)}]'
            unresolved=f'[{fmt(g["p_qc_low"],3)}, {fmt(g["p_qc_high"],3)}]'
            toy_flags=html.escape('; '.join(f'{key}: {value}' for key,value in g['nominal_failure_flags'].items())) or '—'
            rowstyle=' style="background:#ffe9e4"' if g['observed_qc']=='checked_failed' else ''
            gr.append(f'<tr{rowstyle}><td>{g["mass"]:.1f}</td><td>{qc_labels[g["observed_qc"]]}<br>stato {fmt(g["observed_nominal_status"])}; EDM {fmt(g["observed_nominal_edm"],3)}</td><td>{g["n_attempted"]}/100</td><td>{g["n_finite"]}</td><td>{g["k"] if g["k"] is not None else "—"}</td><td>{fmt(g["p"],3)}</td><td>{ci}</td><td>{g["n_nominal_valid"]} / {g["n_nominal_failed"]} / {g["n_nominal_unchecked"]}</td><td>{g["n_unresolved"]}</td><td>{unresolved}</td><td>{toy_flags}</td></tr>')
        gof_table='<section><h2>Dettaglio del GoF e controlli del minimo nominale</h2><p>Il p-value grezzo resta la frazione di <b>tutti gli score finiti</b> sopra quello osservato, senza selezionare i toy convergenti. Il controllo richiede minimizzazione riuscita, 0 ≤ EDM &lt; 0,01, uguaglianza dello stato salvato, snapshot completo, positività del continuo e dei bin totali e accordo dello score ricostruito. Non richiede covQual=3. Un audit concluso con errore conta come esito controllato fallito; la sua assenza conta come non verificato.</p><p class="notice"><b>Se il controllo osservato fallisce o manca, il p-value riportato non supporta una conclusione GoF validata.</b> I limiti conservativi assegnano tutti i toy non risolti prima sotto, poi sopra la soglia osservata salvata. Sono condizionati a tale soglia e ai tentativi completati, non sono intervalli Monte Carlo. Un toy è risolto solo se ha score finito, controllo nominale valido, calibrazione CR valida con 0 ≤ EDM &lt; 0,01 e r lontano dal bordo artificiale di fit. I toy mancanti, non controllati e con controllo fallito restano nel denominatore di questi limiti.</p><p><a href="../../work/shared_gof_diagnostic_runtime/v2/validation/parity.json">Verifica di identità degli score</a> · <a href="../../work/shared_gof_diagnostic_runtime/v2/validation/parser_validation.json">Verifica del parser diagnostico</a> · <a href="../../work/shared_gof_diagnostic_runtime/v2/build_manifest.json">Provenienza della build</a> · <a href="../../work/shared_gof_diagnostic_runtime/v2/IMMUTABLE.json">Immutabilità del runtime v2</a></p><div style="overflow:auto"><table><thead><tr><th>Massa</th><th>GoF osservato</th><th>Esiti</th><th>Score finiti</th><th>Sopra obs.</th><th>p grezzo</th><th>MC 95%</th><th>Toy QC validi / falliti / non controllati</th><th>Non risolti</th><th>Limiti conservativi</th><th>Flag nominali dei toy</th></tr></thead><tbody>'+''.join(gr)+'</tbody></table></div></section>'
    protocol_audit_links=' · '.join(f'<a href="{name}">{label}</a>' for name,label in [('gof_protocol_independent_audit.json','Audit indipendente del protocollo GoF (JSON)'),('gof_protocol_independent_audit.md','Audit del protocollo GoF: note'),('gof/backfill_parity_audit.json','Identità degli score originali dopo il recupero diagnostico')] if (R/name).exists())
    if protocol_audit_links:
        gof_table+='<section><h2>Verifica indipendente del protocollo GoF</h2><p>'+protocol_audit_links+'</p><p>Controllo in sola lettura di generazione, calibrazione preliminare, limiti dei parametri, score e rendicontazione; nessun nuovo fit o modifica della statistica.</p></section>'
    gof_alert=''
    obs_failed=[q for q in gof if q['observed_qc']=='checked_failed']
    if obs_failed:
        masses=', '.join(f'{q["mass"]:.1f}' for q in obs_failed)
        gof_alert='<p class="notice"><b>GoF osservato: controllo del minimo fallito a '+masses+' GeV.</b> Gli score finiti sono conservati, ma i p-value associati non costituiscono risultati GoF validati. Consultare stato, EDM e diagnostica nella tabella GoF.</p>'
    page=f'''<!doctype html>
<html lang="it">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">{refresh}<title>Xee — due SR con continuo comune</title>
<style>body{{font:16px/1.5 system-ui,sans-serif;color:#192639;max-width:1250px;margin:auto;padding:28px;background:#f4f6fa}}h1{{font-size:30px;margin-bottom:8px}}h2{{font-size:22px}}section,.card{{background:white;padding:22px;border-radius:12px;margin:18px 0;border:1px solid #d8e1ec}}.cards{{display:grid;grid-template-columns:repeat(4,1fr);gap:15px}}.cards .card{{margin:10px 0}}.count{{font-size:28px;font-weight:700}}.pending{{color:#68798d}}.notice{{border-left:5px solid #d18a18;background:#fff6e7;padding:14px}}table{{border-collapse:collapse;width:100%;font-size:14px}}th,td{{border-bottom:1px solid #dde4ec;padding:8px;text-align:left}}th{{background:#edf2f8;position:sticky;top:0}}img{{max-width:100%;height:auto}}a{{color:#0758a2}}code{{word-break:break-all}}@media(max-width:700px){{.cards{{grid-template-columns:1fr}}body{{padding:12px}}table{{font-size:11px}}}}</style>
</head>
<body>
<h1>Due SR con un unico continuo, CR simultanea</h1>
<p>{'Esecuzioni completate; verifiche e flag riportati di seguito.' if complete else 'Campagna in corso — risultati pubblicati progressivamente, nessun punto mancante è interpolato.'}<br>Aggiornato: {now}</p>
<div class="cards">
<div class="card">
<div class="count">{counts['observed']} / 61</div>Punti con tutti i risultati osservati<br>Fit: {counts['fitdiagnostics']}/61 · p-value: {counts['significance']}/61<br>Limiti: {counts['limits']}/61 · GoF osservati: {counts['gof_observed']}/61<br>{counts['fit_points']} punti con plot disponibili<br>QC GoF osservati: {counts['gof_observed_diagnostics']}/61; {counts['gof_observed_checked_failed']} falliti</div>
<div class="card">
<div class="count">{counts['bias']} / 6.100</div>Fit dei toy di bias<br>{counts['bias_valid']} stime numericamente valide</div>
<div class="card">
<div class="count">{counts['gof_scores']} / 6.100</div>Esiti GoF raccolti<br>{counts['gof_finite_scores']} score finiti · {counts['gof_missing_scores']} mancanti<br>{counts['gof_points_100_outcomes']} punti con 100 esiti<br>QC nominali: {counts['gof_nominal_diagnostics']}/6.100<br>{counts['gof_nominal_valid']} validi · {counts['gof_nominal_failed']} falliti · {counts['gof_nominal_unchecked']} non controllati</div>
<div class="card">
<div class="count">{counts['control_profiles']} / 122</div>Esiti dei profili sui controlli<br>{counts['control_profiles_valid']} intervalli validi<br>{counts['control_profiles_failed']} falliti · {counts['control_profiles_skipped']} saltati</div>
</div>
<section>
<h2>Configurazione</h2>
<table>
<tr>
<th>Regione</th>
<th>Range [GeV]</th>
<th>Continuo</th>
</tr>
<tr>
<td>SR bassa</td>
<td>6–9,093333</td>
<td rowspan="2">Envelope Chebyshev6/Bernstein6: coefficienti, indice di famiglia e intensità comuni fra le due SR, indipendenti per era</td>
</tr>
<tr>
<td>SR alta</td>
<td>10,506667–12</td>
</tr>
<tr>
<td>CR</td>
<td>9,093333–10,506667</td>
<td>Chebyshev4 indipendente</td>
</tr>
</table>
<p>Quattro ere, binning originale. Υ(1S) e Υ(2S) esplicite in tutte le regioni: medie determinate nella sola CR e poi fissate; sigma e code MC fisse, normalizzazioni condivise e libere. Segnale DCB in tutte le regioni con unico r e nuisance originali. La CR non è esclusa dalla likelihood. La penalità discreta dell’envelope è contata una sola volta per era. Chebyshev e Bernstein dello stesso grado descrivono lo stesso spazio polinomiale con domini dei coefficienti differenti: non sono due famiglie fisiche indipendenti.</p>
<p><b>Procedura richiesta per le medie Υ:</b> le otto μ sono stimate con un fit di solo fondo nella stessa CR e poi fissate nel fit simultaneo. Sono stime inserite come valori fissi (plug-in): <b>la loro incertezza di stima non viene profilata nei limiti e nelle significanze</b>. Ricalibrarle per ogni toy riproduce il passaggio preliminare della procedura; non equivale a trattarle come nuisance profilate e non costituisce una validazione universale.</p>
<p>Bias: generatore alternativo Chebyshev7 globale nelle SR, CR4; 50 toy per massa e iniezione r=0 / limite atteso mediano del nuovo modello. Le stime centrali richiedono EDM finito &lt;0,01 nel fit CR e nel simultaneo. I 122 profili dei toy di controllo prestabiliti sono elaborati in una coda separata, mantenendo immutati i valori centrali; nessuna misura complessiva della copertura.</p>
<p class="notice">Sopra 11 GeV: shape DCB parametriche originali e normalizzazioni estrapolate nello studio precedente, senza nuova incertezza di estrapolazione. Convergenza numerica e completamento delle esecuzioni non costituiscono validazione scientifica del modello.</p>
</section>{gof_alert}{panels}{regional_table}{gof_table}{cr_shift_table}{profile_table}<section>
<h2>Tutti i fit</h2>
<p>Ogni link mostra le quattro ere e le singole componenti. I residui sono diagnostici, non significanze locali.</p>
<div style="overflow:auto">
<table>
<thead>
<tr>
<th>Massa</th>
<th>Plot</th>
<th>r̂</th>
<th>Z locale</th>
<th>p₀</th>
<th>Limite obs.</th>
<th>Atteso med.</th>
<th>Qualità</th>
</tr>
</thead>
<tbody>{''.join(rows)}</tbody>
</table>
</div>
</section>
<section>
<h2>Dati del report</h2>
<p>
<a href="progress.json">Avanzamento JSON</a> · <a href="report_summary.json">Riepilogo JSON</a> · <a href="configuration.json">Configurazione</a> · <a href="pilot_independent_model_audit.json">Verifica indipendente del continuo</a> · <a href="pilot_independent_DCB_audit.json">Verifica DCB e dati</a>{' · <a href="completion_audit.json">Audit indipendente di completezza e provenienza</a>' if (R/'completion_audit.json').exists() else ''}{' · <a href="observed_results.csv">Osservati CSV</a>' if observed else ''}{' · <a href="bias_results.csv">Bias CSV</a> · <a href="bias_toys.csv">Tutti i toy CSV</a>' if bias else ''}{' · <a href="gof_results.csv">GoF CSV</a>' if gof else ''}{' · <a href="control_profiles.csv">Controlli profilo CSV</a>' if profiles else ''}</p>
<p>Errori di rendering temporanei: {len(errors)}. I risultati sorgente non vengono modificati dal report.</p>
</section>
</body>
</html>'''
    temp=R/'index.html.tmp'; temp.write_text(page); temp.replace(R/'index.html')
    print(json.dumps(dict(updated_at=now, completed=complete, counts=counts, rendering_errors=len(errors))), flush=True)
    return complete


if __name__ == '__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--watch', action='store_true'); parser.add_argument('--interval', type=int, default=60)
    args=parser.parse_args()
    while True:
        done=report()
        if done or not args.watch: break
        time.sleep(args.interval)
