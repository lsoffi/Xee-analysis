from pathlib import Path
import json,shutil,os,re,datetime
P=Path('/Users/imb/Documents/Codex/Analisi-statistica-Xee');OLD=P/'runs/20261029_alternative_background_scan';ACTIVE=P/'runs/20261029_bernstein_residuals';archive=OLD/'archive/retired_models';manifest=[]
def move_with_link(src,dest):
 if not src.exists() or src.is_symlink():return
 n=sum(1 for p in src.rglob('*') if p.is_file()) if src.is_dir() else 1;dest.parent.mkdir(parents=True,exist_ok=True);assert not dest.exists();shutil.move(str(src),str(dest));src.symlink_to(os.path.relpath(dest,src.parent),target_is_directory=dest.is_dir());assert src.resolve()==dest.resolve();manifest.append(dict(original=str(src),archived=str(dest),compatibility_link=True,files=n))
for family in ['logquad','expquad','doubleexp']:
 for bucket in ['fits','workspaces','workspaces/seed','seeds','numerical_checks','limit_repairs','limit_diagnostics']:
  move_with_link(OLD/bucket/family,archive/family/bucket)
 for ext in ['png','pdf']:move_with_link(OLD/f'plots/limits_{family}.{ext}',archive/family/f'plots/limits_{family}.{ext}')
comp=OLD/'archive/comparison';comp.mkdir(parents=True,exist_ok=True)
if not (comp/'index.html').exists():
 raw=(OLD/'index.html').read_text();(comp/'index.html').write_text(raw.replace('<meta charset="utf-8">','<meta charset="utf-8"><base href="../../">',1))
if (OLD/'REPORT.md').exists() and not (OLD/'REPORT.md').is_symlink():
 text=(OLD/'REPORT.md').read_text();text=re.sub(r'(\]\()((?!https?://|/|#)[^)]+)(\))',lambda m:m[1]+str(OLD/m[2])+m[3],text);(comp/'REPORT.md').write_text('> ARCHIVIO: confronto precedente. Analisi attiva limitata al Bernstein su richiesta dell’utente.\n\n'+text);(OLD/'REPORT.md').unlink();(OLD/'REPORT.md').symlink_to('archive/comparison/REPORT.md')
for p in list((OLD/'plots').glob('*comparison*')):move_with_link(p,comp/'plots'/p.name)
(OLD/'index.html').write_text('''<!doctype html><meta charset="utf-8"><title>Studio precedente archiviato</title><h1>Studio precedente: confronto archiviato</h1><p>L’analisi attiva è limitata al Bernstein grado 7 identificabile. Le shape delle Υ sono assunte valide.</p><p><a href="../20261029_bernstein_residuals/index.html">Apri i residui del Bernstein</a></p><p><a href="archive/comparison/index.html">Confronto storico completo</a></p><p>I tre modelli alternativi sono conservati in archive/retired_models, con link di compatibilità ai vecchi percorsi.</p>''')
(OLD/'README.md').write_text('# Confronto precedente archiviato\n\nAnalisi attiva: [residui del solo Bernstein](../20261029_bernstein_residuals/README.md).\n\n[Confronto storico](archive/comparison/index.html) · [Report storico](archive/comparison/REPORT.md).\n\nlogquad, expquad e doubleexp sono archiviati in archive/retired_models; tutti i file e i log sono conservati con link di compatibilità. Il Bernstein rimane disponibile nel percorso fits/bernstein7. Nessun input originale è stato modificato.\n')
record=dict(retired_models=['logquad','expquad','doubleexp'],active_model='bernstein7',upsilon_shapes='assumed valid by explicit user instruction; unchanged',entries=manifest,archived_comparison=str(comp),active_directory=str(ACTIVE));(ACTIVE/'archive_manifest.json').write_text(json.dumps(record,indent=2));(OLD/'archive_manifest.json').write_text(json.dumps(record,indent=2))
link=ACTIVE/'fit_results'
if not link.exists():link.symlink_to(os.path.relpath(OLD/'fits/bernstein7',ACTIVE),target_is_directory=True)
print('Archived',sum(m['files'] for m in manifest),'files; preserved compatibility links')
