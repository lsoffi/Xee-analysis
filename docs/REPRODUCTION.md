# Riprodurre l'analisi Xee

Questa guida riguarda la configurazione corrente **C0**, con analisi separate SR bassa+CR e SR alta+CR. L'archivio storico conserva il codice dei confronti precedenti: non significa che quelle validazioni si trasferiscano automaticamente a C0.

## 0. Ambiente e dati

Eseguire dalla radice del repository. Ambiente originario verificato: Combine v11.1.0, commit `6c543b6c0a29df51d390652c2b0ff559c8a211f1`, ROOT 6.34.10, Python 3.12.13, numpy 2.4.6. Il lock Pixi originale è incluso; i grafici usano un ambiente separato, numpy 2.0.2, matplotlib 3.9.4, scipy 1.13.1. Il bootstrap su un computer nuovo richiede rete e un compilatore C++ compatibile; non viene eseguito automaticamente dal clone.

```bash
bash environment/setup.sh
./bin/python scripts/workflow.py doctor
```

Se Combine è già installato e caricato, per esempio in CMSSW:

```bash
export XEE_USE_SYSTEM=1
./bin/python scripts/workflow.py doctor
```

Per usare interpreti/installazioni specifici, impostare `XEE_PYTHON`, `XEE_COMBINE`, `XEE_PLOT_PYTHON` ai rispettivi eseguibili. Le variabili sono percorsi di eseguibili, non stringhe di comando. Non fissare le directory del computer originale.

**I ROOT e le datacard originali non sono pubblicati in questo repository di codice.** Occorre la cartella `ee` originale con sottocartelle `6.0`, `6.1`, ..., `11.0`. Il manifest elenca hash e dimensioni attese. Importare una copia verificata:

```bash
./bin/python scripts/workflow.py import-inputs --source /percorso/alla/cartella/ee
```

È possibile importare soltanto un punto per la prova:

```bash
./bin/python scripts/workflow.py import-inputs --source /percorso/ee --masses 9.5
```

I file vengono copiati in `data/mass_scan/261028/cards_region2_data_envelope_allCorrections_binned/ee/`; un hash diverso causa un errore, senza sovrascrittura. Il file `configs/input_manifest.json` identifica esattamente la versione dei dati. Gli output sono esclusi da git.

## 1. Costruzione del modello corrente

```bash
./bin/python scripts/workflow.py build --masses 9.5
./bin/python tests/smoke.py
```

Questo costruisce entrambe le analisi a 9,5. Lo smoke test non minimizza: controlla dati, parametri liberi, C0, positività, assenza di envelope e convenzioni statistiche. Per il solo lato basso:

```bash
./bin/python scripts/workflow.py build --side low --masses 7.5
```

La costruzione usa `analysis/separate_sr_c0_limits/model.py`. Gli starting point numerici delle vecchie soluzioni sono già esportati, dopo le conversioni esatte Bernstein6→Chebyshev6 o grado3→grado4, in `configs/initial_values.json`. **Sono inizializzazioni, non prior né parametri congelati.** Non occorrono i vecchi ROOT dei fit per ricostruire il modello.

Le configurazioni sono:

- SR bassa 6–9,093333333 GeV, Chebyshev6; oppure SR alta 10,506666667–12 GeV, Chebyshev4.
- CR 9,093333333–10,506666667 GeV, Chebyshev4.
- Quattro ere, otto canali per analisi; nessun dato ripetuto entro una likelihood.
- Otto medie Υ floating; sigma e code MC fisse; rese globali condivise con accettanze regionali dinamiche.
- DCB di segnale in entrambe le regioni, unico `r`, nuisance originali.
- C0 del continuo per era: `N_SR=N_CR*I_SR*P_CR(edge)/(I_CR*P_SR(edge))`; C1 non imposto.

Ogni workspace ha un `.audit.json`. Il `workspace_manifest.json` cresce con i modelli costruiti. Non modificare i workspace già costruiti: usare un altro checkout/output per un nuovo modello.

## 2. Segnale sopra 11 GeV

Per costruire 11,1–12,0 servono prima tutti gli input originali fino a 11:

```bash
./bin/python scripts/workflow.py extend-signal
```

La DCB mantiene il grafo parametrico originale, sostituendo il valore di massa. Le rese nominali vengono rilette dagli input, fittate in 10–11 con `N(m)=N(11)*exp[b*(m-11)+c*(m-11)^2]` e poi estrapolate. Le incertezze relative di normalizzazione sono mantenute a quelle di 11 GeV. L'incertezza aggiuntiva di estrapolazione non è assegnata: questo è un punto aperto, non una misura MC sopra 11. I controlli sui punti lasciati fuori dal fit riguardano masse ≤11.

## 3. Limiti, fit e significanza a un punto

La sequenza di riproduzione usa i valori numerici archiviati in `configs/numerical_scales.json`, che fissano gli intervalli iniziali e la seconda partenza del fit. Non si tratta di vincoli statistici sulle rese.

```bash
./bin/python scripts/workflow.py limits --side low --masses 9.5
./bin/python scripts/workflow.py fit --side low --masses 9.5
./bin/python scripts/workflow.py significance --side low --masses 9.5
./bin/plot-python analysis/separate_sr_c0_limits/report.py
```

`limits` chiama AsymptoticLimits nativo, con r≥0, salva sei quantili e controlla il crossing osservato. `fit` esegue due FitDiagnostics, partenza zero e partenza expected, con la stessa prescrizione per tutti i punti; sceglie il minimo convergente usando la devianza senza offset. `significance` usa la snapshot S+B come inizializzazione delle nuisance, poi esegue Significance nativo e confronta il risultato con la differenza B/S+B di FitDiagnostics. Non sostituisce il risultato nativo quando i due non concordano.

Per leggere i comandi esatti, aprire `combine.command.json` nella sottocartella del metodo. Ogni comando salva log, ricevuta, tempo, codice di uscita; gli estrattori preservano i flag. Un `exit_code=0` non dimostra un minimo corretto.

La root dei risultati è `runs/20261031_separate_sr_c0_limits/`: il nome identifica la strategia storica; nel nuovo checkout contiene le nuove esecuzioni. Ogni punto è `low/m9.5` oppure `high/m9.5`.

## 4. Scan completo

```bash
./bin/python scripts/workflow.py build
./bin/python scripts/workflow.py limits
./bin/python scripts/workflow.py fit
./bin/python scripts/workflow.py significance
./bin/python scripts/workflow.py report
./bin/python scripts/workflow.py audit
```

Lo scan comprende 76 configurazioni: 46 masse low 6,0–10,5 e 30 high 9,1–12,0, passo 0,1. La root CLI esegue in sequenza per evitare di sovraccaricare il computer. È possibile selezionare gruppi disgiunti con `--side` e `--masses`, eseguendoli in processi separati; non eseguire due volte lo stesso punto contemporaneamente. I metodi nativi hanno lock per punto, ma gli aggiornamenti dei manifest non sono pensati per build concorrenti.

Usare `--dry-run` per vedere la selezione senza creare file o importare ROOT. La verifica finale `audit` richiede lo scan completo, non solo il pilot.

Nel report, distinguere sempre: completamento tecnico; convergenza; qualità covarianza; minimo e crossing; validazione statistica. Le bande native precedenti erano provvisorie e molti p-value vicino alle Υ erano incoerenti. La riproduzione non corregge automaticamente questi problemi.

## 5. Profili vicino alle Υ

Prima servono fit e limiti nelle due analisi alle quattro masse. La ricetta storica è:

```bash
./bin/python scripts/workflow.py profiles --masses 9.4 9.5 10.0 10.1 --profile-pass initial
./bin/python scripts/workflow.py profiles --masses 9.4 9.5 10.0 10.1 --profile-pass reverse
./bin/python scripts/workflow.py profiles --masses 9.4 9.5 10.0 10.1 --profile-pass branches
./bin/python scripts/workflow.py profiles --masses 9.4 9.5 10.0 10.1 --profile-pass reference
./bin/python scripts/workflow.py profiles --masses 9.4 9.5 10.0 10.1 --profile-pass report
```

`initial` usa 25 punti quadratici più il valore rhat e tre inizializzazioni; `reverse` percorre il profilo al contrario; `branches` identifica candidati fisicamente ammissibili più bassi di 0,02 nella devianza e applica una sola ripartenza dedicata. I punti non risolti restano marcati. Il risultato è un profilo osservato, non un intervallo calibrato. `reference` riproduce inoltre il raffinamento libero del riferimento, a partire dal miglior punto condizionato; appartiene alla sequenza storica finale dei 204/208 punti utilizzabili. Tutte queste operazioni sono esplicite e non si avviano durante lo scan dei limiti.

Nel metodo condizionale, `FitDiagnostics --customStartingPoint --skipSBFit` salva un oggetto chiamato `fit_b`, ma minimizza il modello completo al valore r fissato. L'estrattore verifica questa convenzione: non usare l'estrattore B-only standard, che forzerebbe r=0.

## 6. Rese fisiche e incertezze

Dopo i fit low/high a 9,5 e 10,0:

```bash
./bin/python scripts/workflow.py yields
```

Si estraggono i coefficienti estesi delle componenti, cioè numeri di eventi integrati nei range indicati. Le scale `bkgScale` da sole non sono rese. Le incertezze sono `sqrt(J V J^T)` con la covarianza completa del fit libero, includendo le medie floating e le accettanze. La covarianza di un punto del profilo a r fissato non viene usata per stimare l'incertezza marginale delle rese S+B. I confini dei parametri sono registrati; errori Hesse simmetrici non sostituiscono intervalli profilati.

## 7. NUOVA validazione C0: generatore candidato, bias e GoF

**Questa sezione fornisce codice eseguibile per il lavoro ancora da fare. Non è la certificazione che i test C0 siano già stati completati.** Nessuna campagna viene avviata senza `--execute`.

Un candidato alternativo può essere costruito e fittato con un grado SR diverso, lasciando CR4 e C0 invariati:

```bash
./bin/python scripts/generator.py --side low --mass 7.5 --degree 7 --tag cheb7_candidate
./bin/python scripts/generator.py --side low --mass 7.5 --degree 7 --tag cheb7_candidate --execute
```

Il runner richiede almeno due partenze B-only convergenti con covQual3, conserva tutte e tre le partenze e salva `runs/generators/cheb7_candidate/low/m7.5/workspace_postfit_b.root`. **Questo non basta per accettare il generatore**: vanno controllati la sua GoF, i residui, la positività e la stabilità. Un generatore della stessa famiglia non esaurisce le alternative possibili. I generatori Bernstein storici sono nell'archivio e non vengono trasformati tacitamente in generatori C0.

Bias r=0 e r=expected, 50 toy ciascuno:

```bash
./bin/python scripts/validation.py bias --side low --mass 7.5 --toys 50 --seed 913075 --injection 0 --tag bias0_cheb7 --generator runs/generators/cheb7_candidate/low/m7.5/workspace_postfit_b.root --execute
./bin/python scripts/validation.py bias --side low --mass 7.5 --toys 50 --seed 923075 --injection expected --tag biasE_cheb7 --generator runs/generators/cheb7_candidate/low/m7.5/workspace_postfit_b.root --execute
```

Il generatore deve contenere la DCB della massa iniettata e lo stesso layout SR/CR. Il programma controlla osservabili, range/bin e global observables; la scelta scientifica della verità generatrice resta esplicita. Il valore expected proviene dal limite mediano salvato e può essere provvisorio: serve come benchmark d'iniezione, non come validazione del limite.

Protocollo implementato:

1. Copia separata del workspace generatore; valori veri fissati, snapshot `generation`.
2. `GenerateOnly --saveToys --toysFrequentist --bypassFrequentistFit`, generazione binned degli otto canali, segnale presente anche nella CR, fluttuazioni Poisson e osservazioni ausiliarie.
3. Importazione di ciascun toy e del suo snapshot delle global observables in una copia del modello di fit; le otto medie Υ restano libere. Nessun passaggio CR-only per congelarle nella configurazione corrente.
4. Due FitDiagnostics con partenze 0 e expected, intervallo r simmetrico, stesso protocollo per tutte le iniezioni. Confronto mediante devianza senza offset, prior incluse una sola volta.
5. Salvataggio di tutti i fit, fallimenti e bordi. Media/SD sui validi e confronto con tutti i valori finiti; nessun toy sostituito.
6. Bias `mean(rhat-rtrue)`, SD campionaria con N−1, errore Monte Carlo SD/sqrt(N); bias/SD con 2.000 bootstrap. Pull Hesse separati con covQual3. **Nessun fit gaussiano alla distribuzione, nessuna misura completa di copertura.**

GoF corrente con r libero, 100 toy generati da B-only:

```bash
./bin/python scripts/validation.py gof --side low --mass 7.5 --toys 100 --seed 914075 --tag gof100 --execute
```

Il runner usa `GoodnessOfFit --algo saturated` sui dati e su ogni toy con r libero. Per ciascun toy confronta anche lo score nativo con il minimo FitDiagnostics e conserva i disaccordi. Il confronto dello score osservato con il FitDiagnostics osservato resta un controllo esplicito da completare prima di validare la nuova campagna. Il p-value grezzo è la frazione di score toy > score osservato; sono salvati anche `(k+1)/(N+1)`, numerosità e fallimenti. Se ci sono score mancanti il p-value è incompleto; non ignorarli. I controlli di convergenza di FitDiagnostics non certificano automaticamente la minimizzazione interna di GoodnessOfFit.

**Dipendenza dalla massa:** cambia la DCB S+B usata per fittare ciascun toy e i dati. Con r=0 fissato, range e modello di fondo invariati, il test sarebbe invece comune alle masse: piccole variazioni sarebbero solo Monte Carlo. I p-value di discovery e quelli di GoF misurano domande diverse.

Ripetere i comandi per le masse previste e poi:

```bash
./bin/plot-python scripts/validation_report.py --tag bias0_cheb7
./bin/plot-python scripts/validation_report.py --tag biasE_cheb7
./bin/plot-python scripts/validation_report.py --tag gof100
```

I tag distinti separano generatore, iniezione e protocollo; non mescolare campagne diverse. Usare seed distinti se si vogliono ensemble indipendenti. Toy appaiati fra masse richiedono un disegno esplicito: questo runner usa ensemble specifici per cella, mentre alcune campagne storiche condividevano i dataset.

## 8. Codice storico e percorso delle scelte

L'archivio contiene 628 file originali con hash verificati: script, configurazioni, wrapper e primi macro ROOT/C++. I symlink sono stati risolti in file reali. Il codice non è stato riscritto o ripulito cambiandone il significato.

| Studio | Cartella originale in `archive/scripts/` | Scopo |
|---|---|---|
| Input originali, scan e primi controlli | `mass_scan`, `current`, `experiments`, `validation` | Riprodurre il punto di partenza; primi bias/range/binning |
| Risoluzione e finestre | `common95_cheb2`, `upsilon_peak_windows` | Accettanze DCB e finestre comuni |
| CR e funzione sotto le Υ | `cr_optimization`, `cr_float_shapes`, `cr_float_mu`, `cr_degree4_mu`, `cr_windows` | Scelta range/modello, μ vs σ, stabilità e toy |
| Grado SR basso | `sr_degree_selection`, `sr_bernstein7_bias`, `sr_cheb7_bern_bias`, `sr_upper_mass_bias` | Generatori alternativi, grado minimo e segnale nella CR |
| Envelope basso | `sr6_cr4_scan`, `gof_profiled_combined`, `gof_thesis48` | Scan, bias Asimov, GoF storico |
| SR alta | `upper_sr_cr`, `upper_envelope_scan`, `upper_bias_cheb5_toys` | Gradi 1–3, segnale estrapolato e toy |
| SR unite, campagna interrotta | `shared_sr_cr_scan`, `shared_sr7_floatmu_fitlimits` | Polinomio globale e studi non completati |
| Analisi separate e C0 | `separate_sr_cheb_limits`, `separate_sr_c0_limits` | Ultima configurazione e confronto del raccordo |
| Profili e rese | `upsilon_c0_profiles`, `bonly_c0_yields` | Identificabilità e normalizzazioni |

**Limite di riproducibilità storica:** gli script originali possono dipendere da percorsi assoluti, runtime locali e output di campagne precedenti. L'archivio da solo non è una pipeline portabile di tutte le 50+ campagne. `archive/source_manifest.json` identifica il file originale; occorrono le dipendenze specifiche elencate nei rispettivi script e i risultati intermedi, oppure rigenerarli in ordine. Il percorso corrente portabile evita questa catena tramite snapshot numerici di inizializzazione con provenienza. Non lanciare uno script storico sul nuovo checkout presumendo che usi la configurazione corrente.

## 9. Controlli realmente eseguiti sul rilascio

`docs/release_smoke_test.json` registra i controlli: sintassi dei sorgenti attivi; import con hash; ricostruzione low/high a 9,5; C0 e positività; conteggi; parità delle predizioni iniziali con i workspace originali; convenzioni media/SD e conservazione dei fit falliti. Non è stata avviata una nuova produzione statistica o una nuova campagna toy durante il confezionamento. Il bootstrap completo su un computer diverso resta da verificare.
