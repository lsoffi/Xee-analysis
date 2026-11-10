"""Independent report workers, atomic generations, five-minute refresh, no fit submissions."""
from pathlib import Path
import subprocess,os,time,json,datetime,fcntl,shutil,socket
root=Path(__file__).resolve().parent
scan=Path('/afs/cern.ch/work/s/soffi/DIELE-WORK/Analisi-statistica-Xee/scan_261106')
web=Path('/eos/user/s/soffi/www/Xee/261106/scan_SR_CR/high_live')
lock=open('/tmp/xee_liveplots.lock','w')
try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
except BlockingIOError:raise SystemExit('Updater already active on this host')
web.mkdir(parents=True,exist_ok=True)
state={};active={};due={k:0 for k in ['bias','coverage','gof']}
def stamp():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def atomic(p,text):
 t=p.with_suffix('.tmp');t.write_text(text);t.replace(p)
while True:
 now=time.time()
 for name in due:
  if name in active:
   proc,start,stage,log=active[name]
   rc=proc.poll()
   if rc is None:
    if now-start>900:state[name]['warning']='Report running over 15 min; possible AFS delay. No duplicate worker launched.'
    continue
   log.close();state[name]['exit_code']=rc;state[name]['finished_utc']=stamp()
   if rc==0:
    dest=web/name;dest.mkdir(exist_ok=True)
    # Publish each file through rename; status goes last and identifies this generation.
    for f in sorted(stage.iterdir(),key=lambda p:p.name=='status.json'):
     tmp=dest/(f.name+'.tmp');shutil.copyfile(f,tmp);tmp.replace(dest/f.name)
    state[name]['published_utc']=stamp();state[name]['result']=json.loads((stage/'status.json').read_text());state[name]['warning']=None
   else:state[name]['warning']='Report failed; previous published figures preserved. See '+str(root/(name+'.log'))
   active.pop(name);due[name]=now+300
  elif now>=due[name]:
   stage=root/('stage_'+name);stage.mkdir(exist_ok=True)
   cmd=['python3','-u',str(root/(name+'.py'))]
   if name=='gof':cmd+=['--results',str(scan/'gof_high_sr/results'),'--observed',str(scan/'gof_high_sr/observed'),'--output',str(stage)]
   env=dict(os.environ,PLOT_OUTPUT=str(stage));log=open(root/(name+'.log'),'w')
   proc=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,env=env)
   state[name]={**state.get(name,{}),'started_utc':stamp(),'pid':proc.pid,'running':True};active[name]=(proc,now,stage,log)
 for name in state:state[name]['running']=name in active
 snapshot={'heartbeat_utc':stamp(),'host':socket.gethostname(),'pid':os.getpid(),'refresh_seconds':300,'reports':state}
 atomic(root/'watch_status.json',json.dumps(snapshot,indent=2));atomic(web/'watch_status.json',json.dumps(snapshot,indent=2))
 page='<html><head><meta http-equiv="refresh" content="300"></head><body><h1>High SR: live partial validation</h1><p>Report refresh every five minutes after completion; AFS delays can postpone publication. Last heartbeat UTC: '+snapshot['heartbeat_utc']+'</p><p>Incomplete samples are provisional and may depend on completion time. Invalid toys are not counted as successes. These figures do not establish final validation.</p>'
 for name,figs in [('bias',['bias_r0','bias_expected']),('coverage',['coverage_vs_mass']),('gof',['gof'])]:
  page+='<h2>'+name+'</h2><pre>'+json.dumps(state.get(name,{}),indent=2)+'</pre>'
  for f in figs:page+=f'<p><a href="{name}/{f}.png">PNG</a> | <a href="{name}/{f}.pdf">PDF</a></p><img style="max-width:100%" src="{name}/{f}.png">'
  page+=f'<p><a href="{name}/{name}.csv">CSV</a></p>'
 page+='</body></html>';atomic(web/'index.html',page)
 time.sleep(20)
