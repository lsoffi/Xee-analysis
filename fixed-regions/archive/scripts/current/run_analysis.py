#!/usr/bin/env python3
"""Reproduce the era-specific B -> S+B workflow in a fresh output directory."""
import argparse,datetime,pathlib,subprocess,sys,json,os
root=pathlib.Path(__file__).resolve().parents[2]
p=argparse.ArgumentParser();p.add_argument('--dry-run',action='store_true');p.add_argument('--run-name');args=p.parse_args()
name=args.run_name or datetime.datetime.now().strftime('%Y%m%d_%H%M%S_%f');run=root/'runs'/name
if '/' in name or name in ['.','..']:raise SystemExit('Invalid run name')
if not args.dry_run:
 if run.exists():raise SystemExit('Run directory already exists; choose a new name')
 (run/'outputs').mkdir(parents=True);(run/'work').mkdir()
commands=[]
def task(label,cmd,cwd):
 commands.append({'label':label,'command':list(map(str,cmd)),'cwd':str(cwd)})
 print(label,flush=True)
 if args.dry_run:return
 with (run/'work'/(label+'.log')).open('w') as out:subprocess.run(list(map(str,cmd)),cwd=cwd,stdout=out,stderr=subprocess.STDOUT,check=True)
def python(label,script):task(label,[sys.executable,root/'scripts/current'/script],run)
def fit(label,era,workspace,extra=(),strategy=0,tolerance=None):
 cmd=[root/'bin/combine','-M','FitDiagnostics',workspace,*extra,'--saveNormalizations','--saveShapes','--keepFailures','--rMin','-30','--rMax','30','--cminDefaultMinimizerStrategy',str(strategy),'--robustFit','1','--X-rtd','MINIMIZER_freezeDisassociatedParams','--cminRunAllDiscreteCombinations','-n',label,'-v','3']
 if tolerance:cmd+=['--cminDefaultMinimizerTolerance',str(tolerance)]
 task(label,cmd,run/'outputs')
python('prepare_degree5_era_start','prepare_degree5_era_start.py')
for era in ['2022','2022EE','2023','2023BPix']:fit('_degree5_stageB_'+era,era,'Xee_degree5_eraStart_'+era+'.root',['--skipSBFit'])
python('seed_B_to_SB','seed_degree5_sb_from_b.py')
for era in ['2022','2022EE','2023','2023BPix']:fit('_degree5_Binitialized_'+era,era,'Xee_degree5_Binitialized_'+era+'.root',['--skipBOnlyFit','--preFitValue','0'])
fit('_degree5_Binitialized_combined','combined','Xee_degree5_Binitialized_combined.root',['--preFitValue','0'])
fit('_degree5_Binitialized_combined_check','combined','Xee_degree5_Binitialized_combined.root',['--preFitValue','0'],strategy=1,tolerance=.001)
task('significance',[root/'bin/combine','-M','Significance','-d','Xee_degree5_Binitialized_combined.root','--redefineSignalPOIs','r','--rMin','-50','--rMax','50','-n','.observed_significance_degree5_Binitialized','-v','1','--cminDefaultMinimizerStrategy','0','--X-rtd','MINIMIZER_freezeDisassociatedParams','--cminRunAllDiscreteCombinations'],run/'outputs')
for label,script in [('extract_report','extract_degree5_sequence.py'),('plot_fits','plot_degree5_sequence.py'),('plot_B_to_SB_functions','plot_degree5_sequence_functions.py')]:python(label,script)
if not args.dry_run:
 (run/'commands.json').write_text(json.dumps(commands,indent=2)+'\n');(root/'runs/LATEST.txt').write_text(str(run)+'\n')
print(('Dry-run verified: ' if args.dry_run else 'Results: ')+str(run))
