import pathlib,json,subprocess,time,concurrent.futures,os,signal,re
P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_gaussian_signal_scan'
def run(a):
 w=R/f"limit_repairs/{a['model']}/m{a['mass_label']}";w.mkdir(parents=True,exist_ok=True);t=json.loads((R/f"fits/{a['model']}/m{a['mass_label']}/limits/command.json").read_text());cmd=t['command'];cmd[cmd.index('--cminDefaultMinimizerStrategy')+1]='1';cmd+=['--cminDefaultMinimizerTolerance','0.001','--rMax',str(max(1,5*a['expected50']))];cmd[cmd.index('-n')+1]='_crossing';t.update(cwd=str(w),command=cmd);(w/'command.json').write_text(json.dumps(t,indent=2));start=time.time()
 with (w/'combine.log').open('w') as f:
  p=subprocess.Popen(cmd,cwd=w,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
  try:rc=p.wait(timeout=300)
  except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGTERM);p.wait();rc=124
 s=(w/'combine.log').read_text();matches=re.findall(r'At r =\s*([\d.eE+-]+):[^\n]*CLs\s*=\s*([\d.eE+-]+)',s);out={**t,'exit_code':rc,'elapsed_seconds':time.time()-start,'last_cls':float(matches[-1][1]) if matches else None};(w/'receipt.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:out[k] for k in ['model','mass','exit_code','elapsed_seconds','last_cls']}),flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=5) as pool:list(pool.map(run,json.loads((R/'limit_crossing_warnings.json').read_text())))
