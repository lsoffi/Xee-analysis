#!/usr/bin/env python3
"""Portable, explicit entry points for the current separate-SR C0 analysis.

No step runs implicitly. Statistical production is deliberately separate from
input import, model building, plots, and numerical audits.
"""
from pathlib import Path
import argparse, hashlib, importlib, json, os, shutil, subprocess, sys

P = Path(__file__).resolve().parents[1]
R = P/'runs/20261031_separate_sr_c0_limits'
INPUT = P/'data/mass_scan/261028/cards_region2_data_envelope_allCorrections_binned/ee'
S = P/'analysis/separate_sr_c0_limits'

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,q):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(q,indent=2)+'\n')
def run(command):
    print('RUN', ' '.join(map(str,command)),flush=True)
    subprocess.run(list(map(str,command)),cwd=P,check=True)
def points(args):
    pp=json.loads((P/'configs/current_scan.json').read_text())['points']
    if args.side: pp=[q for q in pp if q['region']==args.side]
    if args.masses: pp=[q for q in pp if round(q['mass'],1) in args.masses]
    if not pp: raise SystemExit('No configured mass/side combinations selected.')
    return pp

def doctor():
    import ROOT,numpy
    result=ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
    if result<0: raise RuntimeError('Cannot load libHiggsAnalysisCombinedLimit')
    print(json.dumps(dict(python=sys.version,ROOT=ROOT.gROOT.GetVersion(),numpy=numpy.__version__,
                         RooSimultaneousOpt=hasattr(ROOT,'RooSimultaneousOpt'),repository=str(P)),indent=2))
    run([P/'bin/combine','--help'])

def import_inputs(args):
    source=args.source.resolve()
    manifest=json.loads((P/'configs/input_manifest.json').read_text())['files']
    selected=[q for q in manifest if Path(q['relative_path']).parts[0] in
              {f'{m:.1f}' for m in args.masses}] if args.masses else manifest
    if not selected: raise SystemExit('No mass input paths selected.')
    missing=[]
    for q in selected:
        src=source/q['relative_path'];dst=INPUT/q['relative_path']
        if not src.is_file(): missing.append(str(src));continue
        if sha(src)!=q['sha256']: raise RuntimeError('Input checksum mismatch: '+str(src))
        if dst.exists():
            if sha(dst)!=q['sha256']: raise RuntimeError('Refusing to replace different input: '+str(dst))
            continue
        dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst)
    if missing: raise RuntimeError('Missing original inputs:\n'+'\n'.join(missing))
    print('Verified/imported',len(selected),'input files; original files not modified.')

def module(name,folder=S):
    sys.path.insert(0,str(folder));return importlib.import_module(name)

def build(args):
    model=module('model')
    R.mkdir(parents=True,exist_ok=True)
    manifest=json.loads((P/'configs/current_scan.json').read_text())
    manifest.update(requested_methods=['FitDiagnostics','AsymptoticLimits','Significance'],fits_requested=True,significance_requested=True,completed=False)
    path=R/'manifest.json'
    if not path.exists():dump(path,manifest)
    else:
        existing=json.loads(path.read_text());existing.update(fits_requested=True,significance_requested=True);dump(path,existing)
    for q in points(args):
        side,m=q['region'],q['mass'];out=R/side/f'm{m:.1f}/workspace.root'
        if out.exists():
            audit=json.loads(out.with_suffix('.audit.json').read_text())
            assert sha(out)==audit['workspace_sha256'];assert sha(audit['source'])==audit['source_sha256']
            print('Existing immutable workspace',out)
        else: model.build(side,m)
    rows=[]
    for ws in sorted(R.glob('*/m*/workspace.root')):
        audit=json.loads(ws.with_suffix('.audit.json').read_text())
        rows.append(dict(side=audit['side'],mass=audit['mass'],workspace=str(ws),
              workspace_sha256=sha(ws),source=audit['source'],source_sha256=audit['source_sha256']))
    dump(R/'workspace_manifest.json',dict(count=len(rows),workspaces=rows,builder=str(S/'model.py'),builder_sha256=sha(S/'model.py')))

def native(args):
    modules={'fit':'run_fits','limits':'run_limits','significance':'run_significance'}
    m=module(modules[args.step])
    for point in points(args):
        d=R/point['region']/f"m{point['mass']:.1f}"
        if not (d/'workspace.root').is_file():raise RuntimeError('Build workspace first: '+str(d))
        if args.step=='significance' and not (d/'fit_results.json').is_file():raise RuntimeError('FitDiagnostics required first: '+str(d))
        result=m.job(point)
        print(result,flush=True)
        if (isinstance(result,tuple) and result[-1]=='failed') or (isinstance(result,dict) and result.get('status')=='failed'):
            raise RuntimeError('Native method failed; preserve logs and inspect before retrying.')

def extend(args):
    out=P/'runs/20261030_upper_envelope_scan'
    stamp=out/'extension_manifest.json'
    if stamp.exists():
        for q in json.loads(stamp.read_text())['files']:
            assert sha(P/q['path'])==q['sha256'],'Changed signal extension: use a fresh clone for rebuilding.'
        print('Existing signal extension verified; no files rewritten.');return
    if list((out/'extrapolated_inputs').glob('m*/Xee_ee_4_allYears.root')):
        raise RuntimeError('Partial or untracked signal extension exists. Inspect it and use a fresh clone; refusing to rewrite ROOT inputs.')
    import ROOT
    ROOT.gSystem.Load('libHiggsAnalysisCombinedLimit')
    out.mkdir(parents=True,exist_ok=True)
    rates=[]
    for k in range(60,111):
        path=INPUT/f'{k/10:.1f}/Xee_ee_4_allYears.root'
        f=ROOT.TFile.Open(str(path));w=f.Get('w');w.var('r').setVal(1)
        rates.append(dict(mass=k/10,rates=[w.function(f'n_exp_binch{i}_proc_Zd').getVal() for i in range(1,5)]));f.Close()
    dump(out/'signal_rates_existing.json',rates)
    run([P/'bin/plot-python',P/'analysis/signal_extension/extrapolate_rates.py'])
    run([P/'bin/python',P/'analysis/signal_extension/extend_shapes.py'])
    run([P/'bin/python',P/'analysis/signal_extension/build_extrapolated_inputs.py'])
    files=sorted((out/'extrapolated_inputs').rglob('*'))+sorted((out/'extended_signal_shapes').rglob('*'))
    files += [out/'signal_rates_existing.json',out/'signal_rate_extrapolation.json']
    dump(stamp,dict(files=[dict(path=str(p.relative_to(P)),sha256=sha(p)) for p in files if p.is_file()]))

def profiles(args):
    folder=P/'analysis/upsilon_c0_profiles';sys.path.insert(0,str(S))
    m=module('run',folder);cases=[(q['region'],q['mass']) for q in points(args)]
    if any(mass not in [9.4,9.5,10.,10.1] for side,mass in cases):
        raise RuntimeError('Historical profiles were defined at 9.4,9.5,10.0,10.1. Select --masses explicitly.')
    m.R.mkdir(parents=True,exist_ok=True)
    manifest=m.R/'manifest.json'
    present={(q['side'],q['mass']) for q in json.loads(manifest.read_text())['cases']} if manifest.exists() else set()
    present.update(cases)
    dump(manifest,dict(cases=[dict(side=s,mass=mass) for s,mass in sorted(present)]))
    if args.profile_pass=='initial':
        results=[m.run_case(s,mass) for s,mass in cases]
        dump(m.R/'runner.json',dict(status='completed',results=results))
    elif args.profile_pass=='reverse':
        rev=module('reverse',folder)
        results=[rev.job(s,mass) for s,mass in cases]
        dump(m.R/'reverse_runner.json',dict(status='completed',results=results))
    elif args.profile_pass=='branches':
        issues=[]
        for side,mass in cases:
            pd=m.R/side/f'm{mass:.1f}';profile=json.loads((pd/'profile.json').read_text())
            for point in profile['points']:
                directory=pd/('r_'+format(point['r'],'.12g').replace('.','p'))
                trials=[json.loads(x.read_text()) for x in directory.glob('*/result.json')]
                feasible=[q for q in trials if q.get('exact_deviance') is not None and q.get('minimum_bin',0)>0 and min(q.get('continuum_minima',{'x':-1}).values())>0 and q.get('continuity_passed')]
                if not feasible:continue
                best=min(feasible,key=lambda q:q['exact_deviance'])
                if point.get('exact_deviance') is not None and point['exact_deviance']-best['exact_deviance']>.02:
                    issues.append(dict(side=side,mass=mass,r=point['r'],bestsource=best['source']))
        dump(m.R/'branch_check.json',dict(issues=issues,threshold=.02))
        run([P/'bin/python',folder/'branch_refit.py'])
    elif args.profile_pass=='report':
        run([P/'bin/plot-python',folder/'report.py'])
    elif args.profile_pass=='reference':
        for side,mass in cases:
            run([P/'bin/python',folder/'refine_reference.py','--run','--side',side,'--mass',str(mass),'--workers','1'])

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('step',choices=['doctor','import-inputs','extend-signal','build','fit','limits','significance','report','profiles','yields','audit'])
    ap.add_argument('--side',choices=['low','high']);ap.add_argument('--masses',nargs='+',type=float)
    ap.add_argument('--source',type=Path,help='Original ee directory, with mass subdirectories; required by import-inputs.')
    ap.add_argument('--profile-pass',choices=['initial','reverse','branches','reference','report'],default='initial')
    ap.add_argument('--dry-run',action='store_true',help='Print operations; import no ROOT and do not create files.')
    a=ap.parse_args();os.chdir(P)
    if a.masses:a.masses=[round(x,1) for x in a.masses]
    if a.dry_run:
        print(json.dumps(dict(step=a.step,side=a.side,masses=a.masses,profile_pass=a.profile_pass,
             points=points(a) if a.step in ['build','fit','limits','significance','profiles'] else None),indent=2));return
    if a.step=='doctor':doctor()
    elif a.step=='import-inputs':
        if a.source is None:ap.error('import-inputs requires --source /path/to/ee')
        import_inputs(a)
    elif a.step=='extend-signal':extend(a)
    elif a.step=='build':build(a)
    elif a.step in ['fit','limits','significance']:native(a)
    elif a.step=='report':run([P/'bin/plot-python',S/'report.py'])
    elif a.step=='profiles':profiles(a)
    elif a.step=='yields':
        run([P/'bin/python',P/'analysis/bonly_c0_yields/extract.py'])
        run([P/'bin/python',P/'analysis/bonly_c0_yields/compare_sb.py'])
        run([P/'bin/plot-python',P/'analysis/bonly_c0_yields/plot.py'])
    elif a.step=='audit':run([P/'bin/python',S/'audit_completion.py'])

if __name__=='__main__':main()
