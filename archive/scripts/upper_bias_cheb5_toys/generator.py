"""Generate independent 100-toy ensembles from the fitted Chebyshev5 truth."""
from pathlib import Path
import argparse, json, subprocess, time, hashlib, math, struct, fcntl
import ROOT
from prepare_generator import P,R,OLD,ERAS,customize,sha

def generate(mass,label):
    D=R/'generation'/f'm{mass}'/label;D.mkdir(parents=True,exist_ok=True)
    with (D/'generation.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        return _generate(mass,label)

def _generate(mass,label):
    D=R/'generation'/f'm{mass}'/label;D.mkdir(parents=True,exist_ok=True)
    source=OLD/f'degree3/m{mass}/joint/workspace_seeded.root'
    fit_source=R/'generator/results.json'
    expected_source=OLD/f'degree3/m{mass}/results.json'
    if (D/'truth.json').exists():
        truth=json.loads((D/'truth.json').read_text())
        assert truth['source_sha256']==sha(source) and truth['generator_fit_sha256']==sha(fit_source)
        assert truth['expected_source_sha256']==sha(expected_source)
        assert truth['toy_file_sha256']==sha(truth['toy_file'])
        print('ALREADY GENERATED',mass,label,flush=True);return
    generator=json.loads(fit_source.read_text());expected=float(json.loads(expected_source.read_text())['limits']['expected_50'])
    assert math.isfinite(expected) and expected>0
    rtrue=0. if label=='r0' else expected
    seed=851000+round(float(mass)*10)*10+(0 if label=='r0' else 1)
    f=ROOT.TFile.Open(str(source));w=f.Get('w')
    w.var('r').setVal(1)
    original_signal={f'{era}_{reg}':w.function(f'local_rate_Zd_{era}{suffix}').getVal() for era in ERAS for reg,suffix in [('SR',''),('CR','_CR')]}
    keep=customize(w,fixed=True,parameters=generator['parameters'])
    w.var('r').setVal(1)
    signal={f'{era}_{reg}':w.function(f'local_rate_Zd_{era}{suffix}').getVal() for era in ERAS for reg,suffix in [('SR',''),('CR','_CR')]}
    assert max(abs(signal[k]-original_signal[k]) for k in signal)<1e-8
    means={v.GetName():v.getVal() for v in w.allVars() if v.GetName().startswith('Ymean_')}
    assert all(abs(means[a['name']]-a['value'])<1e-12 for a in generator['means'])
    w.var('r').setRange(-10000,10000);w.var('r').setVal(rtrue);w.var('mass').setVal(11.2);w.var('mass_CR').setVal(9.8)
    w.saveSnapshot('generation',w.allVars());w.writeToFile(str(D/'generator.root'));f.Close()
    cmd=[str(P/'bin/combine'),'-M','GenerateOnly',str(D/'generator.root'),'--snapshotName','generation','--saveToys','--toysFrequentist','--bypassFrequentistFit','--expectSignal',str(rtrue),'--genBinnedChannels',','.join(f'{reg}_ch{i}' for reg in ['SR','CR'] for i in range(1,5)),'-t','100','-s',str(seed),'-m',mass,'-n','_generated']
    (D/'command.json').write_text(json.dumps(cmd,indent=2));start=time.time()
    with (D/'generate.log').open('w') as log:rc=subprocess.run(cmd,cwd=D,stdout=log,stderr=subprocess.STDOUT,timeout=600).returncode
    (D/'receipt.json').write_text(json.dumps(dict(exit_code=rc,seconds=time.time()-start),indent=2));assert rc==0
    files=list(D.glob('higgsCombine*.root'));assert len(files)==1,files
    ff=ROOT.TFile.Open(str(files[0]));toy_audit=[]
    for i in range(1,101):
        toy=ff.Get(f'toys/toy_{i}');snap=ff.Get(f'toys/toy_{i}_snapshot');assert toy and snap
        counts={};digest=hashlib.sha256()
        for j in range(toy.numEntries()):
            row=toy.get(j);ch=row.find('CMS_channel').getLabel();n=toy.weight()
            assert math.isfinite(n) and n>=0 and abs(n-round(n))<1e-8
            counts[ch]=counts.get(ch,0)+n
            digest.update(ch.encode());digest.update(struct.pack('ddd',row.getRealValue('mass'),row.getRealValue('mass_CR'),n))
        assert len(counts)==8
        toy_audit.append(dict(toy=i,events=toy.sumEntries(),entries=toy.numEntries(),counts=counts,data_sha256=digest.hexdigest(),snapshot_parameters=snap.getSize()))
    ff.Close();assert len(set(a['data_sha256'] for a in toy_audit))==100
    truth=dict(mass=float(mass),label=label,injection=label,rtrue=rtrue,expected_median=expected,expected_source=str(expected_source),expected_source_sha256=sha(expected_source),expected_status='Previously computed native degree3 median: provisional numerical anchor, held fixed for this injection campaign',seed=seed,ntoys=100,generator_degree=5,fit_degree=3,CR_degree=4,source=str(source),source_sha256=sha(source),generator_fit_source=str(fit_source),generator_fit_sha256=sha(fit_source),generator_workspace_sha256=sha(D/'generator.root'),toy_file=str(files[0]),toy_file_sha256=sha(files[0]),signal_rates_r1=signal,means=means,generation='Frequentist binned Poisson toys including independent auxiliary measurements; fixed conditional truth from observed Cheb5 B-only fit',independent_seeds_across_masses_and_injections=True,signal_in_SR_and_CR=True,SR_range=[10.506666666666668,12],CR_range=[9.093333333333334,10.506666666666668],toy_audit=toy_audit)
    (D/'truth.json').write_text(json.dumps(truth,indent=2));print('GENERATED',mass,label,rtrue,'seed',seed,'100 toys',flush=True)

def audit():
    rows=[];seeds=[];allhash=[]
    for path in sorted((R/'generation').glob('m*/*/truth.json')):
        a=json.loads(path.read_text());assert a['ntoys']==100
        assert sha(a['source'])==a['source_sha256']
        assert sha(a['generator_fit_source'])==a['generator_fit_sha256']
        assert sha(a['expected_source'])==a['expected_source_sha256']
        assert sha(a['toy_file'])==a['toy_file_sha256']
        f=ROOT.TFile.Open(a['toy_file']);assert f and not f.IsZombie() and not f.TestBit(ROOT.TFile.kRecovered)
        for item in a['toy_audit']:
            i=item['toy'];toy=f.Get(f'toys/toy_{i}');snapshot=f.Get(f'toys/toy_{i}_snapshot');assert toy and snapshot
            assert toy.numEntries()==872 and toy.sumEntries()==item['events']
            digest=hashlib.sha256()
            for j in range(toy.numEntries()):
                row=toy.get(j);ch=row.find('CMS_channel').getLabel();n=toy.weight()
                digest.update(ch.encode());digest.update(struct.pack('ddd',row.getRealValue('mass'),row.getRealValue('mass_CR'),n))
            assert digest.hexdigest()==item['data_sha256'],(path,i)
            assert snapshot.getSize()==10
            allhash.append(digest.hexdigest())
        f.Close();seeds.append(a['seed']);rows.append(dict(mass=a['mass'],injection=a['injection'],rtrue=a['rtrue'],seed=a['seed'],toys=100,verified=True,toy_file=a['toy_file']))
    assert len(set(seeds))==len(seeds)
    assert len(set(allhash))==len(allhash)
    result=dict(cells=len(rows),target_cells=32,toys=len(allhash),target_generated_toys=3200,generated_toys_per_cell=100,selected_toys_per_cell=50,selected_toy_ids=list(range(1,51)),selected_toys=len(rows)*50,target_selected_toys=1600,selection='User reduced the fit campaign to50toys/cell after100/cell had already been generated; preserve all generated files and select exactlyIDs1..50 without fit-dependent filtering.',unique_seeds=len(set(seeds)),unique_datasets=len(set(allhash)),source_and_result_hashes_verified=True,all_872_bins_and_10_global_observables=True,rows=rows,priority_generation_note='The final12.0rExpected cell was generated twice with identical command and deterministic seed by overlapping main and priority launches. Final ROOT reopened without recovery; all100 data hashes and snapshots checked against saved audit, no malformed/changed dataset. Subsequent generator launches protected by per-cell filelock. Loaded12edgepilot dataset also matched finaltoy1bitwise: priority_semantic_audit.json.')
    (R/'generation/audit.json').write_text(json.dumps(result,indent=2));(R/'generation_summary.json').write_text(json.dumps(result,indent=2));print('AUDIT',len(rows),'cells',len(allhash),'generated unique toys VERIFIED; selected',len(rows)*50,flush=True)
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--mass',nargs='+',default=['10.8']+[f'{k/10:.1f}' for k in range(105,121) if k!=108]);parser.add_argument('--label',nargs='+',choices=['r0','rExpected'],default=['r0','rExpected']);parser.add_argument('--audit-only',action='store_true');args=parser.parse_args()
    if not args.audit_only:
        for mass in args.mass:
            for label in args.label:generate(mass,label)
    audit()
