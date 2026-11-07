import pathlib,json,subprocess,concurrent.futures,time,os,signal
P=pathlib.Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');R=P/'runs/20261029_gaussian_signal_scan';tasks=[]
for a in json.loads((R/'fit_quality_detailed.json').read_text()):
 if a['fit_b_covQual']==3 and a['fit_s_covQual']==3 and a['fit_b_status']==0 and a['fit_s_status']==0:continue
 for role in ['fitdiagnostics','significance']:
  t=json.loads((R/f"fits/{a['model']}/m{a['mass_label']}/{role}/command.json").read_text());c=t['command'];c[c.index('--cminDefaultMinimizerStrategy')+1]='1';c+=['--cminDefaultMinimizerTolerance','0.001'];c[c.index('-n')+1]='_accurate';w=R/f"numerical_checks/{a['model']}/m{a['mass_label']}/{role}";w.mkdir(parents=True,exist_ok=True);t.update(command=c,cwd=str(w));tasks.append(t)
def run(t):
 w=pathlib.Path(t['cwd']);(w/'command.json').write_text(json.dumps(t,indent=2));start=time.time()
 with (w/'combine.log').open('w') as f:
  p=subprocess.Popen(t['command'],cwd=w,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
  try:rc=p.wait(timeout=180)
  except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGTERM);p.wait();rc=124
 t.update(exit_code=rc,elapsed_seconds=time.time()-start);(w/'receipt.json').write_text(json.dumps(t,indent=2));print(t['model'],t['mass'],t['role'],rc,flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:list(pool.map(run,tasks))
