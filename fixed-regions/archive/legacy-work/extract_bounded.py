import subprocess,sys
for n in [1,2]:
 tag='rebin8_degree5_limit'+str(n)
 for script in ['work/extract_degree4.py','work/export_rebin8_full_bins.py']:
  with open('work/'+tag+'_'+script.split('/')[-1]+'.log','w') as log:subprocess.run([sys.executable,script,tag],stdout=log,stderr=subprocess.STDOUT,check=True)
