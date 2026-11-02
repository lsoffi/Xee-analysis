import subprocess,sys
for degree in [5,6]:
 for stem in ['rebin8','to11_rebin8']:
  tag=stem+'_degree'+str(degree)
  for script in ['work/extract_degree4.py','work/export_rebin8_full_bins.py']:
   with open('work/'+tag+'_'+script.split('/')[-1]+'.log','w') as log:subprocess.run([sys.executable,script,tag],stdout=log,stderr=subprocess.STDOUT,check=True)
