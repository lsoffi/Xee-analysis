import pathlib,json,subprocess,concurrent.futures,time,os,signal
P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_alternative_background_scan';tasks=[]
for a in json.loads((R/'refinement_points.json').read_text()):
 if not a['reason']['covariance_or_status']:continue
 for role in ['limits','significance','fitdiagnostics']:
  src=R/f"fits/{a['family']}/m{a['mass_label']}/{role}/command.json";t=json.loads(src.read_text());cmd=t['command'];cmd[cmd.index('--cminDefaultMinimizerStrategy')+1]='1';cmd+=['--cminDefaultMinimizerTolerance','0.001'];cmd[cmd.index('-n')+1]='_accurate';wd=R/f"numerical_checks/{a['family']}/m{a['mass_label']}/{role}";wd.mkdir(parents=True,exist_ok=True);t.update(command=cmd,cwd=str(wd));tasks.append(t)
(R/'commands_numerical_checks.json').write_text(json.dumps(tasks,indent=2))
def work(t):
 wd=pathlib.Path(t['cwd']);(wd/'command.json').write_text(json.dumps(t,indent=2));start=time.time()
 with (wd/'combine.log').open('w') as f:
  proc=subprocess.Popen(t['command'],cwd=wd,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
  try:rc=proc.wait(timeout=300)
  except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGTERM);proc.wait();rc=124
 r={**t,'exit_code':rc,'elapsed_seconds':time.time()-start};(wd/'receipt.json').write_text(json.dumps(r,indent=2));return r
out=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
 for fu in concurrent.futures.as_completed([pool.submit(work,t) for t in tasks]):
  r=fu.result();out.append(r);progress=dict(completed=len(out),total=len(tasks),execution_failures=sum(a['exit_code']!=0 for a in out));(R/'progress_numerical.json').write_text(json.dumps(progress,indent=2));print(json.dumps(progress),flush=True)
