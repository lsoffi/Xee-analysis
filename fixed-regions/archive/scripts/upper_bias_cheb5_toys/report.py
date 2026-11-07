"""Progress and final reports for Chebyshev5 -> cubic-envelope toy bias.

Read-only regarding fit/generation artifacts. Scientific summaries explicitly retain
failure counts and distinguish valid-fit summaries from all finite point estimates.
"""
from pathlib import Path
import argparse, collections, csv, datetime, html, json, math, os, time
from zoneinfo import ZoneInfo
os.environ.setdefault('MPLCONFIGDIR', str(Path.cwd()/'work/mpl-cache'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import beta, chi2

P = Path.cwd()
R = P/'runs/20261030_upper_bias_cheb5_toys'
MASSES = [k/10 for k in range(105, 121)]
LABELS = ['r0', 'rExpected']
NTOYS = 50
TARGET = len(MASSES) * len(LABELS) * NTOYS
COLORS = ['#0072B2', '#D55E00']
TITLES = [r'$r_{true}=0$', r'$r_{true}=L_{expected,50\%}$']


def finite(x):
    return isinstance(x, (int,float)) and math.isfinite(x)


def read(path, default=None):
    try: return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError): return default


def write_json(path, obj):
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(obj, indent=2, allow_nan=False))
    temporary.replace(path)


def interval95(k, n):
    if not n: return [None,None]
    return [0.0 if k==0 else float(beta.ppf(.025,k,n-k+1)),
            1.0 if k==n else float(beta.ppf(.975,k+1,n-k))]


def stats(values):
    a=np.asarray([x for x in values if finite(x)], dtype=float)
    n=len(a)
    if not n: return dict(n=0,mean=None,sd=None,se=None,median=None)
    sd=float(np.std(a,ddof=1)) if n>1 else None
    return dict(n=n,mean=float(np.mean(a)),sd=sd,se=sd/math.sqrt(n) if sd is not None else None,
                median=float(np.median(a)))


def normalized_bias(values, seed):
    a=np.asarray([x for x in values if finite(x)],dtype=float)
    if len(a)<3 or np.std(a,ddof=1)<=0: return dict(value=None,se=None,ci95=[None,None])
    value=float(np.mean(a)/np.std(a,ddof=1))
    rng=np.random.default_rng(seed)
    sample=rng.choice(a,size=(2000,len(a)),replace=True)
    sd=sample.std(axis=1,ddof=1)
    nondegenerate=sd>max(1e-12,float(np.std(a,ddof=1))*1e-10)
    bs=sample.mean(axis=1)[nondegenerate]/sd[nondegenerate]
    return dict(value=value,se=float(np.std(bs,ddof=1)),ci95=[float(x) for x in np.quantile(bs,[.025,.975])])


def normalize(q, path):
    # The aliases allow partial pilot outputs while the worker settles its schema.
    try: mass=float(path.parents[2].name[1:])
    except ValueError: mass=q.get('mass')
    label=path.parents[1].name
    r=q.get('r',q.get('rhat'))
    rtrue=q.get('rtrue')
    err=q.get('error',q.get('r_error'))
    bias=r-rtrue if finite(r) and finite(rtrue) else None
    lo=q.get('interval_low',q.get('r_lo',q.get('r_lower')))
    hi=q.get('interval_high',q.get('r_hi',q.get('r_upper')))
    interval_kind=q.get('interval_kind',q.get('interval_method',q.get('error_method')))
    if (lo is None or hi is None) and interval_kind != 'profile68':
        if finite(r) and finite(err) and err>0:
            lo,hi=r-err,r+err
            interval_kind='Hesse conditional'
    good_interval=finite(lo) and finite(hi) and finite(r) and lo<r<hi
    protocol=q.get('protocol','legacy_v1')
    protocol_final=protocol=='two_fixed_anchors_v2'
    strict_valid=bool(q.get('valid',False)) and protocol_final
    valid=bool(q.get('fit_valid',q.get('valid',False))) and protocol_final
    interval_valid=bool(q.get('interval_valid',good_interval)) and good_interval and valid
    pull=q.get('pull')
    if pull is None and finite(bias) and (interval_kind != 'profile68' or good_interval):
        # Use the interval side pointing towards truth, where an interval exists.
        sigma=(r-lo if bias>=0 else hi-r) if good_interval else err
        pull=bias/sigma if finite(sigma) and sigma>0 else None
    flags=q.get('flags',q.get('failure_flags',q.get('failure_reasons',[])))
    if isinstance(flags,str): flags=[flags]
    if not protocol_final: flags=list(flags)+['pending_protocol_v2_upgrade']
    elif not valid and not flags: flags=['numerical_invalid']
    row=dict(mass=mass,injection=label,toy=q.get('toy',int(path.parent.name.replace('toy',''))),
             rtrue=rtrue,r=r,error=err,bias=bias,pull=pull,valid=valid,strict_valid=strict_valid,
             interval_low=lo,interval_high=hi,interval_kind=interval_kind,
             interval_valid=interval_valid,covered=bool(lo<=rtrue<=hi) if interval_valid and finite(rtrue) else None,
             status=q.get('status'),covQual=q.get('covQual'),
             CR_valid=q.get('CR_valid'),CR_status=q.get('CR_status'),CR_covQual=q.get('CR_covQual'),
             r_boundary=q.get('r_boundary',False),flags=flags,
             protocol=protocol,protocol_final=protocol_final,hesse_valid=bool(q.get('hesse_valid',q.get('covQual')==3 and finite(err) and err>0)) and valid,hesse_pull=bias/err if finite(bias) and finite(err) and err>0 else None,result=str(path.relative_to(R)))
    return row


def load_rows():
    rows=[]
    for p in sorted((R/'fits').glob('m*/*/toy*/result.json')):
        q=read(p)
        if q is not None:
            normalized=normalize(q,p)
            if 1 <= normalized['toy'] <= NTOYS: rows.append(normalized)
    return rows


def summarize(rows):
    summaries=[]
    for mass in MASSES:
        for li,label in enumerate(LABELS):
            rr=[r for r in rows if r['mass']==mass and r['injection']==label]
            vv=[r for r in rr if r['valid']]
            iv=[r for r in vv if r['hesse_valid']]
            covered=0
            n=len(rr);nv=len(vv);ni=len(iv)
            sb=stats([r['bias'] for r in vv]); sp=stats([r['hesse_pull'] for r in iv]); allb=stats([r['bias'] for r in rr])
            truth=read(R/f'generation/m{mass:.1f}/{label}/truth.json',{})
            injection=next((r['rtrue'] for r in rr if finite(r['rtrue'])),truth.get('rtrue'))
            sums=dict(mass=mass,injection=label,rtrue=injection,n=n,target=NTOYS,n_valid=nv,n_invalid=sum(r['protocol_final'] and not r['valid'] for r in rr),n_pending_upgrade=sum(not r['protocol_final'] for r in rr),n_strict_valid=sum(r['strict_valid'] for r in rr),
                      n_finite=sum(finite(r['r']) for r in rr),n_hesse_valid=ni,n_interval=0,n_covered=0,
                      bias=sb,bias_all_finite=allb,pull=sp,
                      normalized_bias=normalized_bias([r['bias'] for r in vv],105000+round(mass*10)*10+li),
                      coverage=None,coverage_ci95=[None,None],coverage_all_attempted_low=None,coverage_all_attempted_high=None,coverage_not_computed='Profiles are requested only for predetermined toy1 per cell; mixed intervals are not aggregated.',
                      interval_kinds=dict(collections.Counter(r['interval_kind'] or 'unavailable' for r in rr)),
                      flags=dict(collections.Counter(f for r in rr for f in r['flags'])),
                      n_CR_invalid=sum(r['CR_valid'] is False for r in rr),
                      n_r_boundary=sum(r['r_boundary'] for r in rr))
            summaries.append(sums)
    return summaries


def save_plot(fig, name):
    fig.savefig(R/(name+'.png'),dpi=160,bbox_inches='tight')
    fig.savefig(R/(name+'.svg'),bbox_inches='tight')
    plt.close(fig)


def panel_plot(summaries, name, title, ylabel, accessor, err_accessor=None, reference=None, ylim=None, all_accessor=None):
    fig,axes=plt.subplots(1,2,figsize=(13,4.3),sharex=True)
    for li,(label,ax) in enumerate(zip(LABELS,axes)):
        ss=[s for s in summaries if s['injection']==label]
        sample_key='pull' if name.startswith('pull_') else 'bias'
        use=[s for s in ss if finite(accessor(s)) and s[sample_key]['n']>=3]
        if use:
            x=[s['mass'] for s in use]; y=[accessor(s) for s in use]
            errors=[err_accessor(s) or 0 for s in use] if err_accessor else None
            ax.errorbar(x,y,yerr=errors,color=COLORS[li],marker='o',capsize=3,label='Fit numericamente validi')
        if all_accessor:
            aa=[s for s in ss if finite(all_accessor(s)) and s['bias_all_finite']['n']>=3]
            if aa: ax.plot([s['mass'] for s in aa],[all_accessor(s) for s in aa],ls='none',marker='x',color='#555555',label='Tutti i valori finiti (anche flag)')
        if reference is not None: ax.axhline(reference,color='black',lw=1,ls='--')
        nlabel=f" · n={min(s[sample_key]['n'] for s in use)}–{max(s[sample_key]['n'] for s in use)}" if use else ''
        ax.set_title(TITLES[li]+nlabel);ax.set_xlabel('Massa ipotizzata [GeV]');ax.set_ylabel(ylabel)
        ax.grid(alpha=.22);ax.set_xlim(10.45,12.05)
        if ylim:ax.set_ylim(*ylim)
        if ax.get_legend_handles_labels()[0]:ax.legend(fontsize=8)
        if not use:ax.text(.5,.5,'In attesa di almeno 3 toy validi per punto',transform=ax.transAxes,ha='center',color='#777777')
    fig.suptitle(title+'\nChebyshev5 → envelope Chebyshev3/Bernstein3; punti mostrati da 3 toy validi, target 50',fontsize=12)
    fig.tight_layout();save_plot(fig,name)


def plots(summaries,rows):
    panel_plot(summaries,'bias','Bias medio CONDIZIONATO alla riuscita del fit (barre: errore MC)',r'$\langle\hat r-r_{true}\rangle$',
               lambda s:s['bias']['mean'],lambda s:s['bias']['se'],reference=0,
               all_accessor=lambda s:s['bias_all_finite']['mean'])
    panel_plot(summaries,'bias_normalized','Bias medio / dispersione empirica — solo stime valide (condizionato)',r'$\langle\hat r-r_{true}\rangle / SD(\hat r)$',
               lambda s:s['normalized_bias']['value'],lambda s:s['normalized_bias']['se'],reference=0)
    panel_plot(summaries,'pull_mean','Diagnostica secondaria: media del pull Hesse (solo covQual=3)',r'$\langle(\hat r-r_{true})/\sigma_{Hesse}\rangle$',
               lambda s:s['pull']['mean'],lambda s:s['pull']['se'],reference=0)
    panel_plot(summaries,'pull_width','Diagnostica secondaria: larghezza del pull Hesse (solo covQual=3)',r'$SD[(\hat r-r_{true})/\sigma_{Hesse}]$',
               lambda s:s['pull']['sd'],reference=1)
    panel_plot(summaries,'estimator_sd','Dispersione empirica della stima di r',r'$SD(\hat r)$',lambda s:s['bias']['sd'])
    fig,axes=plt.subplots(1,2,figsize=(13,4.2),sharey=True)
    for li,(label,ax) in enumerate(zip(LABELS,axes)):
        ss=[s for s in summaries if s['injection']==label]
        ax.bar([s['mass'] for s in ss],[s['n_valid'] for s in ss],width=.075,color=COLORS[li],label='Stima valida')
        ax.bar([s['mass'] for s in ss],[s['n_invalid'] for s in ss],bottom=[s['n_valid'] for s in ss],width=.075,color='#d9a3a3',label='Flag/falliti')
        ax.bar([s['mass'] for s in ss],[s['n_pending_upgrade'] for s in ss],bottom=[s['n_valid']+s['n_invalid'] for s in ss],width=.075,color='#b7bdc5',label='Da aggiornare a due inizializzazioni')
        ax.axhline(NTOYS,color='black',ls='--',lw=1);ax.set_ylim(0,NTOYS*1.08);ax.set_title(TITLES[li]);ax.set_xlabel('Massa [GeV]');ax.set_ylabel(f'Toy completati / {NTOYS}');ax.legend(fontsize=8)
    fig.suptitle(f'Avanzamento: {len(rows)}/{TARGET} toy — tutti i fallimenti restano contabilizzati');fig.tight_layout();save_plot(fig,'progress')
    fig,axes=plt.subplots(1,2,figsize=(13,4.3),sharey=True)
    largest=max([s['n_invalid']/max(1,s['n']-s['n_pending_upgrade']) for s in summaries]+[0.])
    for label,ax,title in zip(LABELS,axes,TITLES):
        ss=[s for s in summaries if s['injection']==label]
        nn=[s['n']-s['n_pending_upgrade'] for s in ss]
        ff=[s['n_invalid']/n if n else 0. for s,n in zip(ss,nn)]
        ax.bar([s['mass'] for s in ss],ff,width=.075,color='#B44B4B')
        for x,y,s,n in zip([s['mass'] for s in ss],ff,ss,nn):
            ax.text(x,y+.008,f"{s['n_invalid']}/{n}",ha='center',va='bottom',fontsize=7,rotation=90)
        ax.set_title(title);ax.set_xlabel('Massa [GeV]');ax.set_ylabel('Stime fallite / toy completati con protocollo finale')
        ax.set_ylim(0,min(1.,max(.25,largest+.1)));ax.grid(axis='y',alpha=.2)
    fig.suptitle('Frazione di fit numericamente non validi — nessun toy sostituito o rimosso\nConteggi sulle esecuzioni finali già completate; i risultati ai bordi restano segnalati',fontsize=12)
    fig.tight_layout();save_plot(fig,'failure_fraction')
    gallery=[]
    directory=R/'distributions';directory.mkdir(exist_ok=True)
    for s in summaries:
        rr=[r for r in rows if r['mass']==s['mass'] and r['injection']==s['injection']]
        if not rr:continue
        fig,axes=plt.subplots(1,2,figsize=(10,3.6))
        for ax,key,xlabel in [(axes[0],'bias',r'$\hat r-r_{true}$'),(axes[1],'hesse_pull','Pull Hesse (diagnostico)')]:
            vv=[r[key] for r in rr if (r['hesse_valid'] if key=='hesse_pull' else r['valid']) and finite(r[key])];flag=[r[key] for r in rr if not (r['hesse_valid'] if key=='hesse_pull' else r['valid']) and finite(r[key])]
            aa=vv+flag
            if aa:
                # Bound bin count BEFORE allocation: pathological flagged Hesse pulls
                # can make automatic Freedman-Diaconis edges require huge arrays.
                nbins=min(30,max(8,math.ceil(math.sqrt(len(aa)))))
                lower,upper=min(aa),max(aa)
                if lower==upper:lower-=.5;upper+=.5
                bins=np.linspace(lower,upper,nbins+1)
                ax.hist([vv,flag],bins=bins,stacked=True,color=['#0072B2','#C76D6D'],label=['Validi','Flag, finiti'])
                ax.legend(fontsize=8)
            ax.axvline(0,color='black',ls='--',lw=1);ax.set_xlabel(xlabel);ax.set_ylabel('Toy');ax.grid(alpha=.2)
        fig.suptitle(f"m={s['mass']:.1f} GeV, {s['injection']}: {s['n']}/{NTOYS}, {s['n_valid']} validi")
        fig.tight_layout();name=f"distributions/m{s['mass']:.1f}_{s['injection']}";save_plot(fig,name);gallery.append((s,name+'.png'))
    return gallery




def cr_mean_diagnostic(rows):
    values=collections.defaultdict(list)
    eras=['2022','2022EE','2023','2023BPix']
    parameters=[]
    for row in rows:
        if not row['protocol_final'] or row['CR_valid'] is not True:continue
        truth=read(R/f"generation/m{row['mass']:.1f}/{row['injection']}/truth.json",{})
        reference=truth.get('means',{})
        cr=read((R/row['result']).parent/'cr_result.json',{})
        for mean in cr.get('means',[]):
            name=mean['name']
            if name not in reference:continue
            values[(row['mass'],row['injection'],name)].append(1000*(mean['value']-reference[name]))
            if name not in parameters:parameters.append(name)
    parameters.sort(key=lambda name: (0 if '1s' in name.lower() else 1,next((i for i,e in enumerate(eras) if name.endswith('_'+e)),99),name))
    # Keep a deterministic order from truth even before completed records exist.
    if not parameters:
        truth=read(R/'generation/m10.8/r0/truth.json',{})
        parameters=sorted(truth.get('means',{}))
    summary=[]
    for mass in MASSES:
        for label in LABELS:
            for name in parameters:
                st=stats(values[(mass,label,name)])
                summary.append(dict(mass=mass,injection=label,parameter=name,n=st['n'],mean_shift_MeV=st['mean'],se_MeV=st['se'],sd_MeV=st['sd']))
    if summary:
        with (R/'cr_mean_shifts.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(summary[0]));w.writeheader();w.writerows(summary)
    maximum=max([abs(q['mean_shift_MeV']) for q in summary if finite(q['mean_shift_MeV'])]+[1.])
    fig,axes=plt.subplots(1,2,figsize=(14,5),sharey=True)
    im=None
    for label,ax,title in zip(LABELS,axes,TITLES):
        grid=np.full((len(parameters),len(MASSES)),np.nan)
        counts=[]
        for j,mass in enumerate(MASSES):
            subset=[q for q in summary if q['mass']==mass and q['injection']==label]
            for i,name in enumerate(parameters):
                q=next((q for q in subset if q['parameter']==name),None)
                if q and finite(q['mean_shift_MeV']):grid[i,j]=q['mean_shift_MeV']
            counts.append(max([q['n'] for q in subset]+[0]))
        im=ax.imshow(grid,origin='upper',aspect='auto',cmap='RdBu_r',vmin=-maximum,vmax=maximum)
        ax.set_xticks(range(len(MASSES)),[f'{m:.1f}\n(n={n})' for m,n in zip(MASSES,counts)],rotation=90)
        ax.set_yticks(range(len(parameters)),[n.replace('Ymean_','') for n in parameters]);ax.set_title(title);ax.set_xlabel('Massa [GeV] e numero di fit CR validi')
    fig.suptitle('Stima delle medie Υ nella sola CR: scostamento medio dalla verità generata\nSolo CR valide e toy aggiornati al protocollo finale; incertezza MC e conteggi nel CSV',fontsize=12)
    fig.subplots_adjust(left=.15,right=.85,bottom=.25,top=.84,wspace=.15)
    color_axis=fig.add_axes([.88,.25,.018,.59])
    fig.colorbar(im,cax=color_axis,label='Media(μ stimata − μ generata) [MeV]')
    save_plot(fig,'cr_mean_bias')
    return '<h2>Diagnostica delle medie Υ stimate nella CR</h2><p>La mappa mostra lo scostamento medio rispetto alle medie usate per generare i toy. Sono inclusi i fit CR validi dei toy con protocollo finale. Il segnale iniettato è presente anche nella CR, ma il primo stadio fitta solo fondo: un eventuale spostamento delle medie fa quindi parte del bias dell’intera procedura. A r=0 possono contribuire anche fluttuazioni, non linearità e modello del continuo. Non attribuire automaticamente ogni scostamento alla contaminazione. <a href="cr_mean_shifts.csv">Medie, errori MC e numerosità per parametro</a>.</p><img src="cr_mean_bias.png">'


def generator_plots():
    generator=read(R/'generator/results.json')
    if not generator:return ''
    for region in ['SR','CR']:
        fig,axes=plt.subplots(2,2,figsize=(12,8),sharex=True)
        for ax,ch in zip(axes.ravel(),[q for q in generator['channels'] if q['region']==region]):
            points=ch['points'];x=np.array([p['mass'] for p in points]);n=np.array([p['observed'] for p in points])
            low=np.where(n>0,.5*chi2.ppf(.158655254,n*2),0.);high=.5*chi2.ppf(.841344746,(n+1)*2)
            ax.errorbar(x,n,yerr=[n-low,high-n],fmt='.',color='black',ms=3,elinewidth=.7,label='Dati osservati')
            ax.plot(x,[p['expected'] for p in points],color='#0072B2',lw=1.7,label='Generatore B totale')
            for key,color,label in [('dy','#D55E00','Continuo'),('upsilon1s','#009E73','Υ(1S)'),('upsilon2s','#CC79A7','Υ(2S)')]:
                ax.plot(x,[p['components'].get(key,0.) for p in points],color=color,lw=1.2,ls='--',label=label)
            ax.set_title(ch['era']+f" · D={ch['deviance']:.2f}");ax.set_xlabel('Massa [GeV]');ax.set_ylabel('Eventi / bin originale (13,33 MeV)')
            ax.set_xlim(*ch['range']);ax.set_ylim(bottom=0);ax.grid(alpha=.2);ax.legend(fontsize=8)
        degree=5 if region=='SR' else 4
        fig.suptitle(f'Generatore dei toy — {region}: Chebyshev{degree} + Υ(1S) + Υ(2S)\nFit B-only sui dati; medie Υ stimate nella sola CR e poi fisse',fontsize=13)
        fig.tight_layout();save_plot(fig,'generator_'+region)
    return f'<h2>Fondo usato per generare i toy</h2><p>Tre inizializzazioni confrontate con devianza senza offset. Fit selezionato: status {generator["status"]}, covQual {generator["covQual"]}, devianza {generator["deviance"]:.6f}. Massima differenza tra inizializzazioni valide: {generator["max_valid_deviance_difference"]:.3g}. La positività è verificata agli estremi e a tutti i punti stazionari del polinomio. <a href="generator/results.json">Parametri e audit del generatore</a>.</p><details><summary>Generatore nella SR, quattro ere</summary><img loading="lazy" src="generator_SR.png"></details><details><summary>Generatore nella CR, quattro ere</summary><img loading="lazy" src="generator_CR.png"></details>'


def fmt(x, digits=3):
    return '—' if not finite(x) else f'{x:.{digits}g}'


def export_csv(rows,summaries):
    if rows:
        out=[]
        for r in rows:
            q=dict(r);q['flags']=';'.join(q['flags']);out.append(q)
        with (R/'toys.csv').open('w',newline='') as h:
            w=csv.DictWriter(h,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
    ss=[]
    for s in summaries:
        q={k:v for k,v in s.items() if not isinstance(v,(dict,list))}
        for name in ['bias','bias_all_finite','pull','normalized_bias']:
            q.update({name+'_'+k:v for k,v in s[name].items() if not isinstance(v,list)})
        q['coverage_ci95_low'],q['coverage_ci95_high']=s['coverage_ci95'];q['flags']=json.dumps(s['flags']);q['interval_kinds']=json.dumps(s['interval_kinds']);ss.append(q)
    with (R/'summary.csv').open('w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=list(ss[0]));w.writeheader();w.writerows(ss)


def render(active=False):
    R.mkdir(parents=True,exist_ok=True)
    rows=load_rows();summaries=summarize(rows);gallery=plots(summaries,rows);export_csv(rows,summaries);generator_section=generator_plots();cr_section=cr_mean_diagnostic(rows)
    audits=sorted((R/'pilot_profile').glob('m*/*/toy*/profile_endpoint_audit/audit.json'))
    profile_rows=[r for r in rows if r['toy']==1 and r['protocol_final']]
    profile_section='<h2>Controlli profilati: solo toy 1 prestabilito per massa/iniezione</h2><p>Questi controlli verificano la forma locale della likelihood e le asimmetrie. Non sono un campione per stimare la copertura. I profili extra già calcolati restano nel CSV grezzo.</p><div class="scroll"><table><tr><th>Massa</th><th>Iniezione</th><th>r vero</th><th>r stimato</th><th>Intervallo profilato 68%</th><th>Esito</th><th>Risultato</th></tr>'+''.join('<tr>'+''.join('<td>'+v+'</td>' for v in [fmt(r['mass']),r['injection'],fmt(r['rtrue']),fmt(r['r']),('['+fmt(r['interval_low'])+', '+fmt(r['interval_high'])+']') if r['interval_kind']=='profile68' else 'Non disponibile',str(r['interval_valid']) if r['interval_kind']=='profile68' else 'Non profilato','<a href="'+r['result']+'">JSON</a>'])+'</tr>' for r in profile_rows)+'</table></div>'
    audit_section='<p>Controlli indipendenti dei bordi degli intervalli (profilo con tutte le combinazioni discrete): '+', '.join(f'<a href="{a.relative_to(R)}">{a.parents[3].name}, {a.parents[2].name}</a>' for a in audits)+'. Il livello nativo è CL=0,68, corrispondente a ΔNLL≈0,494473.</p>' if audits else ''
    timestamp=datetime.datetime.now(datetime.timezone.utc).isoformat()
    display_time=datetime.datetime.now(ZoneInfo('Europe/Rome')).strftime('%d/%m/%Y %H:%M:%S %Z')
    generated_cells=sum((R/f'generation/m{mass:.1f}/{label}/truth.json').exists() for mass in MASSES for label in LABELS)
    counts=dict(target=TARGET,generated_toys=100*generated_cells,generated_cells=generated_cells,selected_toy_ids=[1,NTOYS],completed=len(rows),valid=sum(r['valid'] for r in rows),invalid=sum(r['protocol_final'] and not r['valid'] for r in rows),pending_upgrade=sum(not r['protocol_final'] for r in rows),protocol_counts=dict(collections.Counter(r['protocol'] for r in rows)),strict_valid=sum(r['strict_valid'] for r in rows),interval_valid=sum(r['interval_valid'] for r in rows),profile_controls=sum(r['toy']==1 and r['interval_kind']=='profile68' and r['interval_valid'] for r in rows),profiles_total=sum(r['interval_kind']=='profile68' and r['interval_valid'] for r in rows),any_flags=sum(bool(r['flags']) for r in rows),
                complete_cells=sum(s['n']==NTOYS for s in summaries),target_cells=32)
    complete=len(rows)==TARGET and all(s['n']==NTOYS for s in summaries) and not counts['pending_upgrade']
    status='completed_with_flags' if complete and counts['any_flags'] else ('completed' if complete else 'running' if active else 'partial')
    write_json(R/'progress.json',dict(status=status,updated_at_utc=timestamp,**counts))
    write_json(R/'summary.json',dict(status=status,updated_at_utc=timestamp,counts=counts,summaries=summaries,
        interpretation='Alternative-background toy bias, Chebyshev5 generator to Chebyshev3/Bernstein3 envelope; conditional fit point estimates after per-toy CR-only mean estimation; profile controls only predeterminedtoy1, no ensemblecoverage',
        methodology=['Generated100 per cell before user reduction; only predetermined IDs1..50 analyzed, independent across mass/injection. IDs51..100 are not failures.',
                     'Only final two_fixed_anchors_v2 results enter valid-fit moments; legacy records are pending upgrades. All finite point estimates are also shown transparently.',
                     'No ensemblecoverage is computed: profile intervals are evaluated only for predeterminedtoy1 per cell. Extra previously computedprofiles remain inrawrecords.',
                     'Mean bias divided by ensemble standard deviation differs from the mean pull.',
                     'The expected-median injection is a fixed archived numerical benchmark, not an asserted validated expected limit.',
                     'With50 successful approximately unit-width pulls, the Monte Carlo error on mean pull is about0.14; tiny biases cannot be conclusively excluded.',
                     'Hesse pulls are secondary diagnostics restricted to covQual3; the primarybias/empiricalSD does not use individualfituncertainties.']))
    table=[]
    for s in summaries:
        table.append('<tr>'+''.join(f'<td>{x}</td>' for x in [f"{s['mass']:.1f}",s['injection'],fmt(s['rtrue']),f"{s['n']}/{NTOYS}",f"{s['n_valid']}",
                      fmt(s['bias']['mean'])+' ± '+fmt(s['bias']['se']),fmt(s['bias_all_finite']['mean']),
                      fmt(s['normalized_bias']['value']),fmt(s['pull']['mean'])+' ± '+fmt(s['pull']['se']),fmt(s['pull']['sd']),
                      str(s['n_hesse_valid']),html.escape(json.dumps(s['interval_kinds'],ensure_ascii=False)),
                      html.escape(json.dumps(s['flags'],ensure_ascii=False))])+'</tr>')
    dist=''.join(f'<details><summary>m={s["mass"]:.1f} GeV · {s["injection"]} · {s["n_valid"]}/{s["n"]} validi</summary><img loading="lazy" src="{name}"></details>' for s,name in gallery)
    figures=''.join(f'<figure><img src="{name}.png"><figcaption>{text}</figcaption></figure>' for name,text in [
        ('progress','Ogni combinazione massa/iniezione considera 50 toy (ID 1–50). I flag sono conservati.'),
        ('failure_fraction','Frazione di stime fallite tra i toy già completati con il protocollo finale; le etichette mostrano numeratore e denominatore. Un punto al bordo numerico non viene trasformato in una stima valida.'),
        ('bias','Bias medio sui fit validi; croci: media di tutti i risultati finiti, inclusi i flag.'),
        ('bias_normalized','Bias medio normalizzato alla dispersione empirica di r: misura distinta dalla media del pull. Errore bootstrap su 2000 ricampionamenti.'),
        ('estimator_sd','Dispersione empirica della stima di r nei toy con stima valida.'),
        ('pull_mean','Diagnostica secondaria: (r̂−r vero)/σHesse, solo con covQual=3. Non mescola errori Hesse e profilati.'),
        ('pull_width','Larghezza del pull Hesse, diagnostica secondaria. Non costituisce una misura di copertura.')])
    refresh='<meta http-equiv="refresh" content="30">' if active and not complete else ''
    document=f'''<!doctype html><html lang="it"><meta charset="utf-8">{refresh}<title>Bias toy Chebyshev5 → envelope3</title>
<style>body{{font:16px system-ui,sans-serif;max-width:1450px;margin:30px auto;padding:0 24px;color:#17212b}}h1{{font-size:29px}}.box{{background:#eef4f8;border-left:4px solid #0072b2;padding:14px}}.warn{{background:#fff4dd;padding:14px}}img{{width:100%;height:auto}}table{{border-collapse:collapse;font-size:12px;width:100%}}th,td{{padding:7px;border:1px solid #ddd;text-align:right}}th{{background:#eef4f8;position:sticky;top:0}}td:nth-child(2){{text-align:left}}figure{{margin:28px 0}}figcaption{{font-size:14px;color:#455}}details{{margin:14px 0}}.scroll{{overflow:auto}}</style>
<h1>Bias con toy: Chebyshev5 → envelope Chebyshev3/Bernstein3</h1>
<p class="box"><b>{counts['completed']}/{TARGET} toy completati</b> · {counts['valid']} stime valide · {counts['profile_controls']}/32 controlli profilati sul toy 1 · {counts['invalid']} stime fallite · {counts['any_flags']} risultati con almeno un flag · {counts['pending_upgrade']} in attesa del protocollo finale a due inizializzazioni · {counts['complete_cells']}/32 combinazioni complete.<br>Generati e conservati: {counts['generated_toys']}; selezionati per il fit: ID 1–{NTOYS} per combinazione.<br>Stato: {status}. Aggiornato: {display_time}.</p>
<p>16 masse, 10,5–12,0 GeV a passi di 0,1; 50 toy per massa e per iniezione: r=0 e r=limite atteso mediano del precedente envelope di grado 3. La generazione iniziale di 100 toy per combinazione è conservata; dopo la riduzione richiesta si analizzano solo gli ID 1–50, scelti prima dei fit. Gli ID 51–100 non sono fallimenti. Le quattro ere sono combinate. Il continuo del generatore è Chebyshev5 nella SR; nel fit è profilata la scelta Chebyshev3/Bernstein3 per era.</p>
<p>SR 10,506667–12 GeV; CR 9,093333–10,506667 GeV, Chebyshev4. In ogni toy si ripete prima il fit B-only della sola CR per stimare le medie Υ(1S)/Υ(2S), poi le medie sono fissate nel fit simultaneo; normalizzazioni Υ condivise e libere, σ e code MC fisse. Segnale DCB in entrambe le regioni, con unico r. Il segnale eventualmente iniettato nella CR è presente nei toy anche se il primo stadio è B-only.</p>
<p class="warn"><b>Il bias mostrato per le stime valide è CONDIZIONATO alla riuscita numerica, non una validazione incondizionata dell’ensemble.</b> I fit falliti o al bordo possono avere una distribuzione diversa dai fit validi: escluderli può cambiare molto la media. Le croci mostrano tutti i valori finiti, compresi quelli con flag; una stima al bordo numerico è censurata e non equivale al vero minimo. I limiti numerici di r non sono stati allargati per nascondere questo problema. Leggere insieme bias, frazione di fallimenti e distribuzioni. Questo studio non autorizza ad accettare il modello basandosi soltanto sulle curve dei toy validi.</p>
<p class="warn">Il protocollo finale esegue due inizializzazioni prestabilite (r=0 e r=expected) per ogni toy, in entrambe le iniezioni, e seleziona il minimo della devianza senza offset. Su richiesta dell’utente, il profilo al 68% viene calcolato solo sul toy 1 prestabilito per ogni combinazione; per tutti i 1600 toy si misura il bias della stima puntuale. I risultati del protocollo iniziale restano visibili come provvisori e sono esclusi dalle medie finali finché non aggiornati. Le statistiche dei soli fit validi sono condizionate alla riuscita numerica: non vanno interpretate come una validazione dell’intero ensemble se molti fit falliscono. Le croci e i conteggi rendono visibili i risultati problematici. Non viene misurata la copertura dell’ensemble: gli intervalli profilati di controllo e gli errori Hesse non sono mescolati per produrre una curva di copertura. I profili extra già calcolati sono conservati nei dati grezzi. Sopra 11 GeV restano la parametrizzazione DCB e l’estrapolazione delle rese già adottate, senza nuova incertezza di estrapolazione. I valori expected sono benchmark provvisori archiviati.</p>
<p><a href="summary.csv">Riassunto CSV</a> · <a href="summary.json">Riassunto JSON</a> · <a href="progress.json">Avanzamento JSON</a>{' · <a href="toys.csv">Tutti i toy CSV</a>' if rows else ''}</p>
{generator_section}{audit_section}{figures}{cr_section}{profile_section}<h2>Risultati per massa e iniezione</h2><div class="scroll"><table><thead><tr><th>Massa</th><th>Iniezione</th><th>r vero</th><th>Completati</th><th>Stime valide</th><th>Bias ± SE</th><th>Bias tutti finiti</th><th>Bias/SD</th><th>Pull Hesse ± SE</th><th>SD pull Hesse</th><th>Hesse validi</th><th>Tipo intervallo</th><th>Flag</th></tr></thead><tbody>{''.join(table)}</tbody></table></div>
<h2>Distribuzioni per combinazione</h2>{dist or '<p>Le distribuzioni appariranno con i primi risultati.</p>'}
<p>Con 50 toy validi e un pull largo circa 1, l’errore statistico sulla media del pull è circa 0,14: questo test non permette di escludere con precisione bias molto piccoli. Il bias di modello è studiato con un solo generatore alternativo specifico, non con tutte le possibili descrizioni del fondo. Nessun risultato precedente è modificato da questo report.</p></html>'''
    tmp=R/'index.html.tmp';tmp.write_text(document);tmp.replace(R/'index.html')
    print(json.dumps(dict(status=status,**counts)),flush=True)
    return complete


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--watch',action='store_true');parser.add_argument('--interval',type=int,default=30);parser.add_argument('--max-seconds',type=int,default=7200)
    args=parser.parse_args();start=time.monotonic()
    while True:
        complete=render(active=args.watch)
        if not args.watch or complete:break
        if time.monotonic()-start>=args.max_seconds:
            render(active=False)
            break
        time.sleep(args.interval)
