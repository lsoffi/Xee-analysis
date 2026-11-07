"""Verify exported profile tables and report links against immutable scan results."""
from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import urlsplit,unquote
import csv,datetime,json,math
from PIL import Image

P=Path(__file__).resolve().parents[2];R=P/'runs/20261102_upsilon_c0_profiles'
def read(p):return json.loads(Path(p).read_text())
def same(a,b):
    if a is None:return b in ['',None]
    try:return math.isclose(float(a),float(b),rel_tol=1e-12,abs_tol=1e-8)
    except (TypeError,ValueError):return False

class Links(HTMLParser):
    def __init__(self):super().__init__();self.links=[]
    def handle_starttag(self,tag,attrs):self.links.extend(v for k,v in attrs if k in ['src','href'])

def audit():
    problems=[];sources={};profiles=[]
    for side in ['low','high']:
        for mass in [9.4,9.5,10.,10.1]:
            path=R/side/f'm{mass:.1f}/profile.json'
            if not path.exists():problems.append('missing_profile:'+str(path));continue
            q=read(path);profiles.append(q)
            for p in q['points']:
                key=(side,mass,p['r'])
                if key in sources:problems.append('duplicate_native_coordinate:'+str(key))
                sources[key]=(p,q)
    with (R/'profile_points.csv').open() as f:rows=list(csv.DictReader(f))
    seen=set()
    for row in rows:
        key=(row['side'],float(row['mass']),float(row['r']))
        if key in seen:problems.append('duplicate_csv_coordinate:'+str(key))
        seen.add(key)
        if key not in sources:problems.append('unexpected_csv_coordinate:'+str(key));continue
        p,q=sources[key]
        for name in ['r','exact_deviance','status','covQual','edm','data_deviance','aux_penalty']:
            if not same(p.get(name),row.get(name)):problems.append(f'csv_value_mismatch:{key}:{name}')
        d=p.get('exact_deviance');ref=q.get('reference_min_deviance');delta=d-ref if isinstance(d,(int,float)) and isinstance(ref,(int,float)) else None
        if not same(delta,row.get('profile_q')):problems.append('q_mismatch:'+str(key))
        if not same(ref,row.get('reference_min_deviance')):problems.append('reference_mismatch:'+str(key))
        if str(p.get('point_valid'))!=row['point_valid']:problems.append('validity_mismatch:'+str(key))
        profile_ok=p.get('profile_point_valid',p.get('point_valid')) is True
        if str(profile_ok)!=row.get('profile_point_valid'):problems.append('profile_validity_mismatch:'+str(key))
        if not same(p.get('lower_feasible_candidate_gap'),row.get('lower_feasible_candidate_gap')):problems.append('minimum_gap_mismatch:'+str(key))
        for key_y in ['dy','upsilon1s','upsilon2s','Zd']:
            if not same((p.get('yields') or {}).get(key_y),row.get('yield_'+key_y)):problems.append(f'yield_mismatch:{key}:{key_y}')
        y=(p.get('yields') or {}).get('Zd');response=y/p['r'] if p['r']>0 and isinstance(y,(int,float)) else None
        if not same(response,row.get('signal_yield_per_r')):problems.append('response_mismatch:'+str(key))
    if seen!=set(sources):problems.append('csv_coordinates_differ_from_sources')
    parser=Links();parser.feed((R/'index.html').read_text());local=[]
    for link in parser.links:
        u=urlsplit(link)
        if u.scheme or not u.path:continue
        path=(R/unquote(u.path)).resolve();local.append(str(path))
        if not path.exists():problems.append('broken_link:'+str(path))
    pngs=list(R.glob('*.png'))+list(R.glob('*/*/*.png'))
    pdfs=list(R.glob('*.pdf'))+list(R.glob('*/*/*.pdf'))
    for p in pngs:
        try:
            with Image.open(p) as im:im.verify()
        except Exception as error:problems.append('invalid_png:'+str(p)+':'+str(error))
    for p in pdfs:
        raw=p.read_bytes()
        if not raw.startswith(b'%PDF-') or b'%%EOF' not in raw[-200:]:problems.append('invalid_pdf_container:'+str(p))
    if len(pngs)!=21 or len(pdfs)!=21:problems.append('unexpected_plot_count')
    summary=read(R/'report_summary.json');native_audit=read(R/'completion_audit.json') if (R/'completion_audit.json').exists() else {}
    native_pass=bool(native_audit.get('passed',native_audit.get('audit_passed',False)))
    expected_completed=len(profiles)==8 and all(q.get('status')=='complete' for q in profiles) and native_pass
    if summary.get('completed')!=expected_completed:problems.append('report_completion_gate_mismatch')
    result=dict(passed=not problems,timestamp=datetime.datetime.now().astimezone().isoformat(),read_only_on_scientific_outputs=True,
        profiles=len(profiles),native_coordinates=len(sources),csv_rows=len(rows),valid_points=sum(p.get('profile_point_valid',p.get('point_valid')) is True for p,q in sources.values()),
        native_valid_points=sum(p.get('point_valid') is True for p,q in sources.values()),
        flagged_points=sum(p.get('profile_point_valid',p.get('point_valid')) is not True for p,q in sources.values()),png_files=len(pngs),pdf_files=len(pdfs),local_links_checked=len(local),
        native_audit_passed=native_pass,report_completed=summary.get('completed'),problems=problems,
        policy='Profiles and yields are read from saved native fits. Failed points are not connected; reference levels do not claim calibrated intervals. No fit or previous result is changed.')
    (R/'report_render_audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
    return result

if __name__=='__main__':
    out=audit();raise SystemExit(0 if out['passed'] else 1)
