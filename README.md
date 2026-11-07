# Xee: analisi simultanea SR + CR con CMS Combine

## Branch e struttura del repository

`main` è l’unico branch pubblicato in `lsoffi/Xee-analysis`. Contiene entrambi i workflow. La directory [sliding-window/](sliding-window/) conserva lo scan sviluppato sul branch storico `sliding-window` e integrato in `main`. Per eseguire entrambe le analisi basta il checkout di `main`.

| Workflow su `main` | Dove iniziare | Documentazione |
|---|---|---|
| **Sliding window**, 69 ipotesi da 0,4 a 9,0 GeV. La SR varia con la massa e il fit usa le CR selezionate per quel punto. | [Istruzioni per lxplus e Condor](sliding-window/README.md) · [Finestre e CR per massa](sliding-window/scan_261106/records.json) | [Nota di analisi in inglese (PDF)](sliding-window/docs/Xee_analysis_note.pdf). A pagina 6, il test sul bias di segnale e sull’efficienza del trigger a 1,4 GeV. |
| **Regioni fisse**, spettro 6–12 GeV. Due analisi separate condividono la CR delle Υ. | [Workflow portabile](scripts/workflow.py) · [Istruzioni di riproduzione](docs/REPRODUCTION.md) | [Guida alle regioni fisse (PDF)](docs/Xee_guida_SR_CR.pdf) · [Metodo](docs/METODO.md) |

Nello sliding window, [cr_261106/](sliding-window/cr_261106/) contiene il codice dei fit delle risonanze SM, [scan_261106/](sliding-window/scan_261106/) contiene il fit simultaneo, i job Condor e gli studi di bias, e [provenance/](sliding-window/provenance/) contiene i checksum degli input.

## Configurazione di riferimento a regioni fisse

Si eseguono due analisi separate, ciascuna simultanea sulle quattro ere (2022, 2022EE, 2023, 2023BPix):

| Analisi | Signal region | Control region | Continuo SR / CR |
|---|---|---|---|
| Bassa + CR | 6-9,093333 GeV | 9,093333-10,506667 GeV | Chebyshev 6 / 4 |
| Alta + CR | 10,506667-12 GeV | Stessa CR | Chebyshev 4 / 4 |

- In ogni canale: continuo + Upsilon(1S) + Upsilon(2S) + DCB del segnale.
- Le due Upsilon hanno medie e normalizzazioni libere, condivise fra SR e CR della stessa era; larghezze e code fissate ai modelli MC.
- Il segnale compare anche nella CR, con un unico parametro `r` e accettanze regionali calcolate dagli integrali delle shape.
- Il raccordo C0 collega esattamente le densita del continuo al confine. I coefficienti dei due polinomi restano distinti e liberi. Nessun raccordo della derivata e nessun envelope in questa configurazione.
- Le due analisi condividono i dati CR: sono correlate e non si combinano moltiplicandone le likelihood.

**Configurazione concordata non significa validazione conclusa.** Fit, limiti e p-value sono stati calcolati, con problemi numerici conservati e documentati. I precedenti bias e GoF con toy riguardano altri modelli; non validano automaticamente il C0. Il PDF termina ogni studio con le verifiche mancanti.

## Cartelle del workflow a regioni fisse

| Percorso | Uso |
|---|---|
| `scripts/workflow.py` | Entrata unica per importazione, build, fit, limiti, significanza, profili, rese e report |
| [bin/](bin/) e [tests/](tests/) | Comandi per usare l’ambiente locale e test del workflow |
| `analysis/` | Adattamento portabile dei moduli della configurazione C0 |
| `configs/` | Griglia, checksum degli input, inizializzazioni e impostazioni numeriche riproducibili |
| `environment/` | Versione Combine, ambiente Pixi e dipendenze dei grafici |
| `archive/` | Codice originale degli studi, preservato; include le prime macro C++ |
| `docs/` | Manuale, istruzioni operative e delimitazione dei risultati riproducibili |
| `evidence/` | Risultati compatti originali, con provenienza e checksum, usati nel manuale |

Gli input ROOT, gli ambienti installati, i toy e gli output pesanti non sono versionati. I file di configurazione riportano i checksum degli input attesi. La directory `archive/` conserva anche dipendenze e percorsi storici: leggere la guida prima di eseguire direttamente uno di quegli script.

## Ambiente

Ambiente di riferimento verificato: **Combine v11.1.0**, commit `6c543b6c0a29df51d390652c2b0ff559c8a211f1`, **ROOT 6.34.10**, Python **3.12.13**. Installare [Pixi](https://pixi.sh), quindi:

```bash
git clone https://github.com/lsoffi/Xee-analysis.git
cd Xee-analysis
bash environment/setup.sh
./bin/python scripts/workflow.py doctor
```

Il setup scarica Combine dal repository ufficiale e fissa il commit; non modifica una precedente installazione CMSSW. Per usare un ambiente gia disponibile consultare `docs/REPRODUCTION.md`.

## Primo punto riproducibile

Fornire la directory originale `ee`, contenente le sottodirectory delle masse. Non rinominare i file né sostituire un workspace con un altro che abbia lo stesso nome ma un checksum diverso.

```bash
./bin/python scripts/workflow.py import-inputs --source /percorso/ee --masses 7.5
./bin/python scripts/workflow.py build --side low --masses 7.5
./bin/python scripts/workflow.py limits --side low --masses 7.5
./bin/python scripts/workflow.py fit --side low --masses 7.5
./bin/python scripts/workflow.py significance --side low --masses 7.5
./bin/python scripts/workflow.py report
```

L'ordine limiti → fit è una convenzione operativa: le inizializzazioni usano il benchmark atteso archiviato in `configs/numerical_scales.json`, indipendente dal limite appena ricalcolato. Comandi, log, output e controlli restano sotto `runs/20261031_separate_sr_c0_limits/`. La data nel nome identifica la campagna originale, non l'orario della nuova esecuzione. I risultati esistenti sono preservati: per una campagna completamente nuova usare un clone/cartella di esecuzione distinto.

## Scan completo

```bash
./bin/python scripts/workflow.py import-inputs --source /percorso/ee
./bin/python scripts/workflow.py extend-signal
./bin/python scripts/workflow.py build
./bin/python scripts/workflow.py limits
./bin/python scripts/workflow.py fit
./bin/python scripts/workflow.py significance
./bin/python scripts/workflow.py report
./bin/python scripts/workflow.py audit
```

Sono **76 combinazioni**: 46 masse nella bassa+CR (6-10,5), 30 nell'alta+CR (9,1-12), passo 0,1 GeV. Le masse nella CR sono analizzate con entrambe le likelihood. L'estensione sopra 11 GeV usa le formule DCB originali e rese estrapolate: non include una nuova incertezza di estrapolazione.

Le istruzioni dettagliate per profili vicino alle Upsilon, rese regionali, bias e GoF sono in **[REPRODUCTION.md](docs/REPRODUCTION.md)**. Le ricette di validazione C0 sono distinte dalle campagne storiche completate. Nessun comando di questa sezione avvia automaticamente i toy.

## Come leggere i risultati

- Conservare gli output nativi, gli stati di convergenza, la qualita della covarianza, l'EDM, i bordi dei parametri e gli hash.
- Il crossing `CLs = 0.05` non dimostra da solo di avere trovato il minimo globale.
- Un buon GoF non dimostra assenza di bias sul segnale. Un bias medio piccolo non dimostra copertura corretta.
- Vicino alle Upsilon, shape quasi indistinguibili e nuisance del segnale possono produrre profili molto ampi. Non usare `resa / errore Hesse` come significanza.
- Le medie dei toy validi sono condizionate al successo numerico: riportare anche i fallimenti e non sostituirli silenziosamente.

## Fonti e ricostruzione

Gli algoritmi nativi appartengono a [CMS Combine](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/v11.1.0/). La parametrizzazione C0, i trasferimenti regionali, le campagne, gli estrattori, i controlli e i report sono documentati nel codice di questo repository. L'archivio e i suoi manifest permettono di identificare gli script originali; gli adattamenti portabili sono separati per evitare di alterare retroattivamente la provenienza dei risultati.

Per rigenerare il PDF installare `docs/requirements.txt`, quindi eseguire `python scripts/render_manual.py`. La verifica visiva del PDF resta necessaria dopo modifiche sostanziali al testo.
