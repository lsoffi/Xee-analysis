# Archivio sorgenti originali

Questa directory conserva i file così come usati durante gli studi. `source_manifest.json` riporta percorso originale, SHA256, dimensione e dereferenziazione dei symlink. Non contiene cache, librerie ROOT/Combine, ROOT di dati o output di fit, grafici, log, credenziali o memoria della chat.

- `scripts/`: tutti gli script del progetto, incluse le vecchie versioni e i duplicati storici nominati.
- `legacy-work/`: primi script e macro del workspace precedente, compresi i file non collegati in `scripts/experiments`.
- `bin/`: wrapper originali, con percorsi storici.
- `config/`: configurazioni e provenienza originali.

I percorsi assoluti sono conservati per tracciabilità. **Non eseguire i wrapper originali come se fossero portabili.** La versione portabile corrente è in `analysis/`, con i comandi in `scripts/workflow.py` e la guida `docs/REPRODUCTION.md`. I vecchi test envelope/μ congelate non validano automaticamente l'attuale C0 con μ libere.
