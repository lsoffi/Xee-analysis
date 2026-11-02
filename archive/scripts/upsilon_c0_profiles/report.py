"""Read-only plots of fixed-r conditional fits under the Upsilon peaks."""
from pathlib import Path
import argparse,csv,datetime,fcntl,hashlib,html,json,math,os
from zoneinfo import ZoneInfo

P=Path(__file__).resolve().parents[2]
R=P/'runs/20261102_upsilon_c0_profiles'
os.environ.setdefault('MPLCONFIGDIR',str(P/'work/mpl-cache'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SIDES=['low','high'];MASSES=[9.4,9.5,10.0,10.1]
LABELS={'low':'SR bassa Chebyshev6 + CR4','high':'SR alta Chebyshev4 + CR4'}
COLORS={9.4:'#0072B2',9.5:'#D55E00',10.0:'#009E73',10.1:'#CC79A7'}
PROCESSES=[('upsilon1s','Υ(1S)','#E69F00'),('upsilon2s','Υ(2S)','#CC79A7'),('Zd','Segnale DCB','#D55E00')]
VERSION=5

def read(path,default=None):
    try:return json.loads(Path(path).read_text())
    except (OSError,json.JSONDecodeError):return default

def finite(x):return isinstance(x,(float,int)) and math.isfinite(x)
def number(x):return f'{x:.5g}' if finite(x) else '—'
def profile_valid(point):return point.get('profile_point_valid',point.get('point_valid')) is True
def yield_value(point,key):
    value=(point.get('yields') or {}).get(key)
    return value if finite(value) else np.nan
def clean(x):
    if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [clean(v) for v in x]
    if isinstance(x,float) and not math.isfinite(x):return None
    return x

def dump(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+f'.tmp.{os.getpid()}')
    temp.write_text(json.dumps(clean(obj),indent=2,allow_nan=False)+'\n');temp.replace(path)

def csvwrite(path,rows):
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with Path(path).open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader();writer.writerows(rows)

def save(fig,name,directory=R):
    directory.mkdir(parents=True,exist_ok=True)
    for ext in ['png','pdf']:
        target=directory/f'{name}.{ext}';temp=target.with_name(f'.{name}_{os.getpid()}.{ext}')
        fig.savefig(temp,dpi=160);temp.replace(target)
    plt.close(fig)

def trigger_parameter(point,profile):
    """Read the selected fit's nuisance, without altering a native fit or scan."""
    native=Path(point.get('source') or '.')
    raw=read(native.parent/'result.json',{})
    parameters={p['name']:p.get('value') for p in raw.get('parameters',[])}
    value=parameters.get('CMS_eff_trigger_e')
    if finite(value):return value,str(native.parent/'result.json')
    original=read(Path(profile.get('source_workspace') or '.').parent/'fit_results.json',{})
    stage='b' if abs(point['r'])<1e-9 else 's'
    previous=original.get('fit',{}).get(stage,{})
    if previous.get('source')==str(native):
        parameters={p['name']:p.get('value') for p in previous.get('parameters',[])}
        value=parameters.get('CMS_eff_trigger_e')
        if finite(value):return value,str(Path(profile['source_workspace']).parent/'fit_results.json')
        if abs(point['r'])<1e-9 and abs(previous.get('aux_penalty',float('inf')))<1e-10:
            return 0.,'Audited previous B-only fit: zero auxiliary penalty and nominal global observables'
    return None,None

def collect():
    rows=[]
    for side in SIDES:
        for mass in MASSES:
            path=R/side/f'm{mass:.1f}'/'profile.json';q=read(path)
            if q:
                if q.get('side')!=side or abs(q.get('mass',-100)-mass)>1e-7:raise ValueError('Profile identity mismatch: '+str(path))
                q=dict(q,source=str(path.relative_to(R)))
                q['points']=sorted([dict(p,yields=p.get('yields') or {},flags=p.get('flags') or []) for p in q.get('points',[])],key=lambda x:x['r'])
                for p in q['points']:
                    p['CMS_eff_trigger_e'],p['trigger_parameter_source']=trigger_parameter(p,q)
                    p['signal_yield_per_r']=yield_value(p,'Zd')/p['r'] if p['r']>0 else None
                if len({p['r'] for p in q['points']})!=len(q['points']):raise ValueError('Duplicate scan coordinate: '+str(path))
                rows.append(q)
    return rows

def point_deviance(point,profile):
    d=point.get('exact_deviance');d0=profile.get('reference_min_deviance')
    return d-d0 if finite(d) and finite(d0) else np.nan

def arrays(profile):
    points=profile['points'];r=np.array([p['r'] for p in points],float)
    q=np.array([point_deviance(p,profile) for p in points],float)
    good=np.array([profile_valid(p) and finite(p.get('exact_deviance')) for p in points],bool)
    return r,q,good

def q_limits(profiles):
    values=[point_deviance(p,q) for q in profiles for p in q['points'] if profile_valid(p)]
    values=[x for x in values if math.isfinite(x)]
    if not values:return (-.15,5.)
    return (min(-.15,min(values)-.1),max(4.5,min(12.,max(values)*1.08+.05)))

def draw_profile(ax,profile,ylimits=None):
    r,q,good=arrays(profile);mass=profile['mass'];color=COLORS[mass]
    label=f'mX = {mass:.1f} GeV'+(' (parziale)' if profile.get('status')!='complete' else '')
    ax.plot(r,np.where(good,q,np.nan),'o-',color=color,ms=3.5,lw=1.6,label=label)
    lowcov=np.array([p.get('covQual')!=3 for p in profile['points']],bool)&good
    if lowcov.any():ax.plot(r[lowcov],q[lowcov],'o',color=color,mfc='white',ms=4,lw=.7)
    bad=~good
    if bad.any():
        y0,y1=ylimits or q_limits([profile]);positions=np.where(np.isfinite(q),np.clip(q,y0+.03*(y1-y0),y1-.03*(y1-y0)),y0+.03*(y1-y0))
        ax.plot(r[bad],positions[bad],'x',color='#B2182B',ms=6,mew=1.5,zorder=6)
    return bool(bad.any()),bool(lowcov.any())

def style_q(ax,profiles):
    limits=q_limits(profiles);ax.set_ylim(*limits)
    for value in [1.,3.84]:
        ax.axhline(value,color='#777777',ls=':',lw=.8)
        ax.text(.985,value,f'{value:g}',ha='right',va='bottom',color='#666666',fontsize=8,transform=ax.get_yaxis_transform())
    ax.set_xlabel('Intensità del segnale fissata, r');ax.set_ylabel('−2 Δ ln L (profilo)');ax.grid(alpha=.16)
    return limits

def profile_plot(profiles):
    fig,axes=plt.subplots(1,2,figsize=(14,6.5));flag=False;lowcov=False
    for side,ax in zip(SIDES,axes):
        selected=[q for q in profiles if q['side']==side];limits=style_q(ax,selected)
        for q in selected:
            b,c=draw_profile(ax,q,limits);flag|=b;lowcov|=c
        if selected:ax.legend(fontsize=9,loc='best')
        else:ax.text(.5,.5,'In attesa dei profili',ha='center',va='center',transform=ax.transAxes)
        ax.set_title(LABELS[side])
    fig.suptitle('Profili sotto le Υ — quattro ere combinate, continuità C0 SR/CR\nr fissato; normalizzazioni e otto medie Υ rifittate insieme agli altri parametri',fontsize=13)
    notes=['Livelli 1 e 3,84: riferimenti numerici, non intervalli calibrati in presenza di degenerazione.',
           'Le due analisi condividono i dati della CR e restano separate. Ogni curva ha il proprio riferimento di likelihood.']
    if flag:notes.append('× rosse: punto non valido, curva interrotta; valori fuori scala riportati al margine. Dati completi nei CSV.')
    if lowcov:notes.append('Cerchi vuoti: minimo convergente con covarianza segnalata; non sono barre di incertezza.')
    fig.text(.5,.018,'\n'.join(notes),ha='center',va='bottom',fontsize=8.4)
    fig.subplots_adjust(top=.79,bottom=.23,wspace=.23)
    save(fig,'profiles')

def draw_yields(ax,profile,delta=False):
    r,_,good=arrays(profile);points=profile['points']
    zeros=[p for p in points if abs(p['r'])<1e-9 and profile_valid(p)]
    baseline=zeros[0].get('yields',{}) if zeros else None
    processes=PROCESSES+([('dy','Continuo','#0072B2')] if delta else [])
    if delta and baseline is None:
        ax.text(.5,.5,'Manca un punto valido a r=0',ha='center',transform=ax.transAxes)
        return
    for proc,label,color in processes:
        y=np.array([yield_value(p,proc) for p in points],float)
        if delta:y-=baseline.get(proc,np.nan)
        ax.plot(r,np.where(good,y,np.nan),'o-',color=color,ms=2.3,lw=1.6,label=label)
    if not delta:
        total=np.array([sum(yield_value(p,proc) for proc,_,_ in PROCESSES) for p in points])
        ax.plot(r,np.where(good,total,np.nan),'o--',color='#333333',ms=1.6,lw=1,label='Υ(1S)+Υ(2S)+segnale')
    else:
        total=np.array([sum(yield_value(p,proc)-baseline.get(proc,np.nan) for proc,_,_ in processes) for p in points])
        ax.plot(r,np.where(good,total,np.nan),'o--',color='#333333',ms=1.6,lw=1,label='Somma di tutte le componenti')
        ax.axhline(0,color='#888888',lw=.7)
    ax.set_xlabel('Intensità del segnale fissata, r');ax.set_ylabel('Variazione degli eventi rispetto a r=0' if delta else 'Eventi attesi nelle regioni selezionate')
    ax.grid(alpha=.16);ax.legend(fontsize=7.5,loc='best')

def compensation_plot(profiles,side=None):
    selected={(q['side'],q['mass']):q for q in profiles}
    order=[(s,m) for s in SIDES for m in [9.5,10.0]] if side is None else [(side,m) for m in MASSES]
    fig,axes=plt.subplots(2,2,figsize=(13,9));anydata=False
    for pair,ax in zip(order,axes.ravel()):
        q=selected.get(pair)
        if q:draw_yields(ax,q);anydata=True
        else:ax.text(.5,.5,'In attesa del profilo',ha='center',transform=ax.transAxes)
        ax.set_title(f'{"SR bassa" if pair[0]=="low" else "SR alta"} + CR | mX = {pair[1]:.1f} GeV')
    fig.suptitle('Come il fit redistribuisce gli eventi fra segnale e Υ\nRese condizionali integrate su SR+CR e sommate sulle quattro ere',fontsize=13)
    fig.text(.5,.025,'Rese condizionali, non intervalli. Anche le nuisance del segnale variano: Nsegnale non è necessariamente proporzionale a r.\nIl continuo resta libero; i punti non validi interrompono le curve. Le due analisi SR+CR non sono combinate.',ha='center',fontsize=9)
    fig.subplots_adjust(top=.87,bottom=.14,hspace=.34,wspace=.24)
    save(fig,'compensation'+('_'+side if side else ''))

def nuisance_plot(profiles):
    lookup={(q['side'],q['mass']):q for q in profiles}
    fig,axes=plt.subplots(2,2,figsize=(13,9))
    for pair,ax in zip([(s,m) for s in SIDES for m in [9.5,10.0]],axes.ravel()):
        profile=lookup.get(pair);right=ax.twinx()
        ax.set_title(f'{"SR bassa" if pair[0]=="low" else "SR alta"} + CR | mX={pair[1]:.1f} GeV')
        if profile:
            r,_,good=arrays(profile)
            theta=np.asarray([p.get('CMS_eff_trigger_e') for p in profile['points']],dtype=float)
            response=np.asarray([p.get('signal_yield_per_r') for p in profile['points']],dtype=float)
            ax.plot(r,np.where(good,theta,np.nan),'o-',ms=2.5,lw=1.5,color='#0072B2',label='θ trigger elettronico')
            right.plot(r,np.where(good,response,np.nan),'o-',ms=2.3,lw=1.6,color='#D55E00',label='Nsegnale / r dopo il fit')
            h1,l1=ax.get_legend_handles_labels();h2,l2=right.get_legend_handles_labels();ax.legend(h1+h2,l1+l2,fontsize=7.5,loc='best')
        else:ax.text(.5,.5,'In attesa del profilo',ha='center',transform=ax.transAxes)
        ax.axhline(0,color='#999999',lw=.7);ax.grid(alpha=.12)
        ax.set_xlabel('Intensità del segnale fissata, r');ax.set_ylabel('θ trigger [unità del vincolo gaussiano]',color='#0072B2');ax.tick_params(axis='y',colors='#0072B2')
        right.set_ylabel('Eventi di segnale per unità r',color='#D55E00');right.tick_params(axis='y',colors='#D55E00');right.set_ylim(bottom=0)
    fig.suptitle('La profilazione modifica anche la risposta del segnale\nNuisance del trigger e resa effettiva per unità r nei fit selezionati',fontsize=13)
    fig.text(.5,.022,'Nsegnale/r è mostrato soltanto per r>0. Tutte le altre nuisance restano profilate.\nInput originale a 10 GeV (2022): −1σ del trigger moltiplica la resa per 0,2504. La nuisance era già nella datacard.\nLa coda del profilo include sia lo scambio segnale–Υ sia la riduzione della resa tramite le incertezze del segnale.',ha='center',fontsize=8.7)
    fig.subplots_adjust(top=.87,bottom=.14,hspace=.34,wspace=.43);save(fig,'signal_nuisance')

def point_rows(profiles):
    rows=[]
    for q in profiles:
        for p in q['points']:
            row=dict(side=q['side'],mass=q['mass'],r=p['r'],profile_q=point_deviance(p,q),
                exact_deviance=p.get('exact_deviance'),reference_min_deviance=q.get('reference_min_deviance'),
                point_valid=p.get('point_valid'),profile_point_valid=profile_valid(p),minimum_check_passed=p.get('minimum_check_passed'),
                lower_feasible_candidate_gap=p.get('lower_feasible_candidate_gap'),status=p.get('status'),edm=p.get('edm'),covQual=p.get('covQual'),
                anchor=p.get('anchor'),flags='; '.join(p.get('flags',[])),source=p.get('source'),
                data_deviance=p.get('data_deviance'),aux_penalty=p.get('aux_penalty'),CMS_eff_trigger_e=p.get('CMS_eff_trigger_e'),
                signal_yield_per_r=p.get('signal_yield_per_r'),trigger_parameter_source=p.get('trigger_parameter_source'))
            row.update({'yield_'+key:p.get('yields',{}).get(key) for key in ['dy','upsilon1s','upsilon2s','Zd']})
            rows.append(clean(row))
    return rows

def compensation_examples(profiles):
    out=[]
    for mass in [9.5,10.0]:
        profile=next((q for q in profiles if q['side']=='low' and q['mass']==mass),None)
        if profile is None:continue
        zeros=[p for p in profile['points'] if abs(p['r'])<1e-9 and profile_valid(p)]
        candidates=[p for p in profile['points'] if profile_valid(p) and finite(p.get('exact_deviance'))]
        if not zeros or not candidates:continue
        zero=zeros[0];best=min(candidates,key=lambda p:p['exact_deviance'])
        out.append(dict(side='low',mass=mass,r=best['r'],baseline_yields=zero['yields'],selected_yields=best['yields'],
            changes={k:yield_value(best,k)-yield_value(zero,k) for k in ['dy','upsilon1s','upsilon2s','Zd']},
            selection='Lowest deviance among usable scanned points; not a new global minimum claim.'))
    return out

def single_profile(profile):
    directory=R/profile['side']/f'm{profile["mass"]:.1f}';stamp=hashlib.sha256(json.dumps(clean([VERSION,profile]),sort_keys=True).encode()).hexdigest()
    if read(directory/'plot_cache.json',{}).get('signature')==stamp and all((directory/f'{name}.{ext}').exists() for name in ['profile_and_yields','yield_changes'] for ext in ['png','pdf']):return
    fig,axes=plt.subplots(1,2,figsize=(13,5.5));limits=style_q(axes[0],[profile]);draw_profile(axes[0],profile,limits);draw_yields(axes[1],profile)
    axes[0].legend(fontsize=9);fig.suptitle(f'{LABELS[profile["side"]]} | mX={profile["mass"]:.1f} GeV\nProfilo con r fissato e compensazione delle rese Υ/segnale',fontsize=13)
    fig.text(.5,.025,'Quattro ere; C0; medie e normalizzazioni Υ variabili. Livelli 1/3,84 descrittivi, non intervalli calibrati.\n×: punto non valido; cerchio vuoto: covarianza segnalata. Le rese includono soltanto le due regioni di questa analisi.',ha='center',fontsize=8.5)
    fig.subplots_adjust(top=.78,bottom=.22,wspace=.28);save(fig,'profile_and_yields',directory)
    fig,ax=plt.subplots(figsize=(9,5.5));draw_yields(ax,profile,delta=True)
    fig.suptitle(f'{LABELS[profile["side"]]} | mX={profile["mass"]:.1f} GeV\nVariazione di ciascuna resa rispetto al fit r=0',fontsize=12)
    fig.text(.5,.02,'Rese integrate nelle regioni selezionate e sommate sulle quattro ere. Nessuna banda di incertezza.',ha='center',fontsize=8.5)
    fig.subplots_adjust(top=.81,bottom=.17);save(fig,'yield_changes',directory)
    csvwrite(directory/'profile_points.csv',point_rows([profile]));dump(directory/'plot_cache.json',dict(signature=stamp))

def render():
    R.mkdir(parents=True,exist_ok=True);profiles=collect();expected=[(s,m) for s in SIDES for m in MASSES]
    signature=hashlib.sha256(json.dumps(clean([VERSION,profiles]),sort_keys=True).encode()).hexdigest()
    if profiles and read(R/'plot_cache.json',{}).get('signature')!=signature:
        profile_plot(profiles);compensation_plot(profiles);nuisance_plot(profiles)
        for side in SIDES:compensation_plot(profiles,side)
        dump(R/'plot_cache.json',dict(signature=signature))
    for q in profiles:single_profile(q)
    rows=point_rows(profiles)
    if rows:csvwrite(R/'profile_points.csv',rows)
    counts=dict(expected_profiles=len(expected),available_profiles=len(profiles),completed_profiles=sum(q.get('status')=='complete' for q in profiles),
        points=len(rows),valid_points=sum(p['profile_point_valid'] is True for p in rows),flagged_points=sum(p['profile_point_valid'] is not True for p in rows),
        native_valid_points=sum(p['point_valid'] is True for p in rows),native_invalid_points=sum(p['point_valid'] is not True for p in rows),
        minimum_check_failed_points=sum(p['minimum_check_passed'] is False for p in rows))
    audit=read(R/'completion_audit.json',{})
    audited=bool(audit.get('passed',audit.get('audit_passed',False)))
    summary=dict(updated_at=datetime.datetime.now(ZoneInfo('Europe/Rome')).isoformat(timespec='seconds'),counts=counts,
        completed=counts['completed_profiles']==len(expected) and audited,production_complete=counts['completed_profiles']==len(expected),final_audit_passed=audited,
        profiles=[dict(side=q['side'],mass=q['mass'],status=q.get('status'),points=len(q['points']),
            valid_points=sum(profile_valid(p) for p in q['points']),native_valid_points=sum(p.get('point_valid') is True for p in q['points']),
            reference_min_deviance=q.get('reference_min_deviance'),reference_source=q.get('reference_source'),source=q['source']) for q in profiles],
        conventions=dict(statistic='Exact deviance difference from the per-profile reliable reference supplied by the fit runner. No negative values are clipped.',
            likelihood_levels='1 and 3.84 are numerical reference levels only; no coverage or nominal confidence interval is claimed.',
            yields='Analytical expected events integrated over the selected SR+CR, summed over all four eras.',
            covariance='A convergent conditional point with covQual != 3 remains plotted as an open marker. A known lower feasible candidate beyond the runner tolerance also excludes the point from the profile, independently of native convergence; profile_point_valid controls plots.',
            shared_data='Low/high analyses reuse the same CR observations and remain separate correlated fits.'))
    summary['compensation_examples']=compensation_examples(profiles)
    dump(R/'report_summary.json',summary)
    gallery=[]
    for side in SIDES:
        cards=[]
        for mass in MASSES:
            q=next((q for q in profiles if q['side']==side and q['mass']==mass),None)
            if q is None:cards.append(f'<section><h3>mX={mass:.1f} GeV</h3><p>In attesa del profilo.</p></section>');continue
            folder=f'{side}/m{mass:.1f}';cells=[]
            for p in q['points']:
                css=' class="bad"' if not profile_valid(p) else ''
                values=[p['r'],point_deviance(p,q),p.get('status'),p.get('covQual'),p.get('edm')]+[p.get('yields',{}).get(k) for k in ['upsilon1s','upsilon2s','Zd','dy']]
                cells.append(f'<tr{css}>'+''.join(f'<td>{number(v)}</td>' for v in values)+f'<td>{html.escape("; ".join(p.get("flags",[])))}</td></tr>')
            table='<div class="scroll"><table><thead><tr>'+''.join(f'<th>{h}</th>' for h in ['r','−2ΔlnL','status','covQual','EDM','Υ1','Υ2','segnale','continuo','flag'])+'</tr></thead><tbody>'+''.join(cells)+'</tbody></table></div>'
            cards.append(f'<section><h3>mX={mass:.1f} GeV — {html.escape(str(q.get("status","in corso")))}</h3><p>Riferimento devianza: {number(q.get("reference_min_deviance"))}. {sum(profile_valid(p) for p in q["points"])}/{len(q["points"])} punti utilizzabili nel profilo; {sum(p.get("point_valid") is True for p in q["points"])} punti con convergenza nativa.</p><a href="{folder}/profile_and_yields.png"><img loading="lazy" src="{folder}/profile_and_yields.png" alt="Profilo e rese {side} {mass:.1f}"></a><p><a href="{folder}/profile_and_yields.pdf">PDF</a> · <a href="{folder}/yield_changes.png">Variazioni di tutte le rese rispetto a r=0</a> · <a href="{folder}/yield_changes.pdf">PDF variazioni</a> · <a href="{folder}/profile_points.csv">Punti CSV</a> · <a href="{q["source"]}">Risultato originale JSON</a></p><details><summary>Stato numerico e rese per punto</summary>{table}</details></section>')
        gallery.append(f'<h2>{LABELS[side]}</h2>'+''.join(cards))
    overview=''
    if profiles:overview=''.join(f'<section><h2>{title}</h2><a href="{name}.png"><img src="{name}.png" alt="{title}"></a><p><a href="{name}.pdf">PDF</a></p></section>' for name,title in [('profiles','Profili di likelihood'),('compensation','Compensazione sotto i due picchi'),('signal_nuisance','Effetto delle incertezze sulla resa del segnale')])
    if profiles:overview='<section><p class="notice"><b>Il profilo include tutte le incertezze del modello.</b> A r fissato, le nuisance di efficienza e forma del segnale vengono rifittate. La resa effettiva del segnale può quindi stabilizzarsi anche mentre r aumenta. La coda larga del profilo non va attribuita soltanto allo scambio con la Υ: il grafico del trigger e di Nsegnale/r rende visibile anche questa compensazione. Non sono state congelate nuisance aggiuntive.</p></section>'+overview
    if profiles:overview='<section><p>La grande incertezza del trigger <b>è già presente negli input originali</b>. Per il segnale del 2022, la datacard a 10 GeV contiene <code>CMS_eff_trigger_e lnN 0.2504/1.7496</code>: una variazione di −1σ moltiplica la resa per 0,2504; a 9,5 GeV il fattore inferiore è 0,30417. Questo permette a r di aumentare mentre la resa effettiva cambia molto meno. Non è stata introdotta una nuova nuisance in questo studio. Fonti: <a href="../../data/mass_scan/261028/cards_region2_data_envelope_allCorrections_binned/ee/10.0/Xee_ee_4_allYears.txt">datacard originale a 10 GeV, riga 37</a> · <a href="../../data/mass_scan/261028/cards_region2_data_envelope_allCorrections_binned/ee/9.5/Xee_ee_4_allYears.txt">datacard originale a 9,5 GeV, riga 37</a>.</p></section>'+overview
    now=summary['updated_at'];all_links='<a href="report_summary.json">Configurazione e riepilogo</a>'
    if rows:all_links+=' · <a href="profile_points.csv">Tutti i punti CSV</a>'
    if (R/'completion_audit.json').exists():all_links+=' · <a href="completion_audit.json">Audit finale</a>'
    if profiles:overview='<section><p>'+('Verifica finale di integrità e completezza disponibile; i singoli flag numerici restano visibili.' if audited else 'Profili disponibili; verifica finale e controllo dei riferimenti di likelihood ancora in corso.')+'</p></section>'+overview
    if profiles:overview=f'<section><p>Convergenza nativa: {counts["native_valid_points"]}/{counts["points"]} punti. Utilizzabili nella curva dopo il controllo dei minimi: {counts["valid_points"]}/{counts["points"]}. Un punto convergente è comunque segnalato se esiste un candidato fisicamente ammissibile con devianza inferiore oltre la tolleranza 0,02. Il report conserva separatamente entrambe le valutazioni.</p></section>'+overview
    if summary['compensation_examples']:
        example_rows=[]
        for case in summary['compensation_examples']:
            d=case['changes'];mass=case['mass'];folder=f'low/m{mass:.1f}'
            example_rows.append(f'<tr><td>{mass:.1f}</td><td>{number(case["r"])}</td><td>{d["Zd"]:+.0f}</td><td>{d["upsilon1s"]:+.0f}</td><td>{d["upsilon2s"]:+.0f}</td><td>{d["dy"]:+.0f}</td><td><a href="{folder}/yield_changes.png">Grafico variazioni</a></td></tr>')
        overview='<section><h2>Quale componente compensa il segnale?</h2><p>Esempi SR bassa+CR: variazioni degli eventi attesi dal fit r=0 al punto scansionato di devianza minima fra quelli utilizzabili. A 9,5 GeV diminuiscono sia Υ(1S) sia continuo; a 10 GeV la resa Υ(2S) varia pochissimo e il segnale è compensato soprattutto dal continuo. Il modello reale non impone shape identiche a segnale e Υ.</p><div class="scroll"><table><thead><tr><th>Massa</th><th>r del punto</th><th>Δ segnale</th><th>Δ Υ1</th><th>Δ Υ2</th><th>Δ continuo</th><th>Figura</th></tr></thead><tbody>'+''.join(example_rows)+'</tbody></table></div></section>'+overview
    page=f'''<!doctype html><html lang="it"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Xee — profili sotto le Υ</title><style>body{{font:16px/1.5 system-ui,sans-serif;color:#172b40;background:#f4f6f9;max-width:1320px;margin:auto;padding:25px}}section{{background:white;border:1px solid #d8e0e9;border-radius:10px;padding:18px;margin:18px 0}}img{{max-width:100%}}.notice{{background:#fff5df;padding:14px;border-left:4px solid #bb861c}}.scroll{{overflow:auto}}table{{border-collapse:collapse;font-size:12px;width:100%}}td,th{{padding:7px;border-bottom:1px solid #ddd;text-align:left}}.bad{{background:#fff0ed}}a{{color:#0758a2}}summary{{cursor:pointer}}@media(max-width:700px){{body{{padding:12px}}}}</style></head><body><h1>Profili di likelihood sotto Υ(1S) e Υ(2S)</h1><p>Aggiornato: {now}. Profili disponibili: {counts['available_profiles']}/8; completati: {counts['completed_profiles']}/8. Punti validi: {counts['valid_points']}/{counts['points']}.</p><section><h2>Che cosa mostriamo</h2><p>Ipotesi di massa 9,4, 9,5, 10,0 e 10,1 GeV. Ogni punto fissa r e rifitta gli altri parametri, incluse le normalizzazioni delle due Υ e le otto medie. Il modello conserva le shape DCB originali e la continuità C0 del continuo, senza imporre una stessa pendenza. SR bassa 6–9,093333 GeV, SR alta 10,506667–12 GeV; CR 9,093333–10,506667 GeV. Le due analisi SR+CR sono separate e correlate tramite i dati della CR.</p><p>Le rese sono eventi attesi integrati nelle regioni effettivamente incluse e sommati sulle quattro ere. Le shape di segnale e Υ sono quelle del modello: non vengono rese artificialmente identiche. Il continuo resta rifittato; il grafico delle variazioni mostra anche il suo contributo.</p><p class="notice">Le linee a −2ΔlnL=1 e 3,84 sono soltanto riferimenti descrittivi: non si dichiarano intervalli di confidenza calibrati in presenza di degenerazione. Minimi non validi sono segnalati e non collegati alla curva. Il riferimento di ciascun profilo è conservato nel JSON originale. Nessun limite, p-value o bias precedente viene sostituito da questi grafici.</p></section>{overview}{''.join(gallery)}<section><h2>Dati</h2><p>{all_links}</p></section></body></html>'''
    temp=R/f'.index_{os.getpid()}.html';temp.write_text(page);temp.replace(R/'index.html')
    print(json.dumps(dict(counts=counts,completed=summary['completed'])),flush=True)
    return summary

if __name__=='__main__':
    argparse.ArgumentParser(description=__doc__).parse_args();R.mkdir(parents=True,exist_ok=True)
    with (R/'report.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX);render()
