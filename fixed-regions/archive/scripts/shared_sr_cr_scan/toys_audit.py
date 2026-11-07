"""Read-only production audit for the shared-SR toy campaign."""
from pathlib import Path
import argparse, hashlib, json, math, struct
import ROOT
from toys_helpers import atomic_json, apply_globals, channel_layout
from generator_toys import R, sha


def audit_workspace(path, toy, snapshot, cr):
    f=ROOT.TFile.Open(str(path));w=f.Get('w');data=w.data('toy_data');assert data
    assert data.numEntries()==toy.numEntries() and abs(data.sumEntries()-toy.sumEntries())<1e-8
    a=[];b=[]
    for sample,out in [(data,a),(toy,b)]:
        for k in range(sample.numEntries()):
            row=sample.get(k)
            out.append((row.find('CMS_channel').getLabel(),row.getRealValue('mass_SRL'),row.getRealValue('mass_CR'),row.getRealValue('mass_SRH'),sample.weight()))
    assert a==b,'Toy data changed while cloning into fit workspace'
    for glob in w.obj('ModelConfig').GetGlobalObservables():
        assert glob.getVal()==snapshot.find(glob.GetName()).getVal()
    for q in cr['means']:
        assert w.var(q['name']).isConstant() and w.var(q['name']).getVal()==q['value']
    assert len(channel_layout(w))==12
    assert w.genobj('discreteParams').getSize()==4
    f.Close()


def main(args):
    checks=[];issues=[]
    for path in sorted((R/'toy_templates').glob('m*/provenance.json')):
        a=json.loads(path.read_text());good=sha(a['immutable_copy'])==a['sha256']
        if not good:issues.append(dict(path=str(path),issue='immutable_fit_template_changed'))
    generated=[]
    # Common0 copied metadata is intentionally not counted as new independent ensembles.
    truths=list((R/'bias/generation').glob('m*/rExpected/truth.json'))
    truths += [p for p in [R/'bias/generation/common_r0/truth.json',R/'gof/generation/truth.json'] if p.exists()]
    for path in truths:
        a=json.loads(path.read_text());ok=sha(a['toy_file'])==a['toy_file_sha256'] and sha(a['source'])==a['source_sha256']
        if not ok:issues.append(dict(path=str(path),issue='generated_input_hash_changed'))
        generated.append(dict(path=str(path),seed=a['seed'],ntoys=a['ntoys'],hashes_verified=ok))
    paths=sorted((R/'bias/fits').glob('m*/*/toy*/result.json'))
    if args.all or len(paths)<=args.sample:selected=paths
    elif args.sample<=1:selected=paths[:args.sample]
    else:selected=[paths[round(k*(len(paths)-1)/(args.sample-1))] for k in range(args.sample)]
    for path in selected:
        a=json.loads(path.read_text())
        if a.get('exception') or not a.get('directory'):
            checks.append(dict(path=str(path),readable=False,exception=a.get('exception')));continue
        truth=json.loads(Path(a['truth_source']).read_text());f=ROOT.TFile.Open(truth['toy_file']);toy=f.Get(f"toys/toy_{a['toy']}");snap=f.Get(f"toys/toy_{a['toy']}_snapshot")
        cr=json.loads((Path(a['CR_cache_directory'])/'cr_result.json').read_text())
        try:
            audit_workspace(Path(a['directory'])/'workspace.root',toy,snap,cr)
            checks.append(dict(path=str(path),data_unchanged=True,aux_snapshot_exact=True,eight_means_from_CR_fixed=True,twelve_channels=True,four_shared_indices=True))
        except Exception as error:issues.append(dict(path=str(path),issue=repr(error)))
        f.Close()
    point_valid=0
    for path in paths:
        q=json.loads(path.read_text());point_valid+=bool(q.get('fit_valid'))
        if q.get('fit_valid') and not (q.get('edm') is not None and math.isfinite(q['edm']) and 0<=q['edm']<.01):issues.append(dict(path=str(path),issue='point_marked_valid_with_bad_EDM'))
    gof_paths=list((R/'gof').glob('toy*/m*/result.json'));gof_qc=0;gof_finite=0;gof_valid=0
    for path in gof_paths:
        q=json.loads(path.read_text());a=q.get('nominal_diagnostics',{});score=q.get('score')
        gof_qc+=bool(a.get('outcome_complete') or a.get('available'))
        gof_finite+=bool(score is not None and math.isfinite(score));gof_valid+=bool(q.get('valid'))
        if a.get('valid') and (not a.get('available') or not a.get('edm_valid') or not a.get('continuum_positive') or not a.get('score_state_consistent')):
            issues.append(dict(path=str(path),issue='GoF_nominal_marked_valid_without_required_checks'))
    profile_paths=list((R/'bias/profiles').glob('m*/*/toy001/result.json'))
    out=dict(bias_results=len(paths),bias_target=6100,point_valid=point_valid,point_invalid=len(paths)-point_valid,
             gof_scores=len(gof_paths),gof_target=6100,gof_finite_scores=gof_finite,gof_quality_valid=gof_valid,gof_QC_outcomes=gof_qc,
             immutable_templates=len(list((R/'toy_templates').glob('m*/provenance.json'))),
             independent_generation_ensembles=generated,independent_generated_datasets=sum(a['ntoys'] for a in generated),
             sampled_semantic_checks=checks,issues=issues,passed=not issues,
             coverage_evaluated=False,profile_outcomes=len(profile_paths),profile_target=122,
             points_and_gof_complete=len(paths)==6100 and len(gof_paths)==6100,
             complete=len(paths)==6100 and len(gof_paths)==6100 and gof_qc==6100 and len(profile_paths)==122)
    atomic_json(R/'toys_audit.json',out)
    print(json.dumps({k:out[k] for k in ['bias_results','point_valid','gof_scores','independent_generated_datasets','passed','complete','issues']},indent=2))
    return 0 if not issues else 1


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--sample',type=int,default=6);p.add_argument('--all',action='store_true');raise SystemExit(main(p.parse_args()))
