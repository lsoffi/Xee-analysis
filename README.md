# Xee: analisi simultanea SR + CR con CMS Combine

Tabella delle sliding window e delle CR per massa · [PDF](sliding-window/docs/tabella_SR_CR_5MeV.pdf) · [CSV](sliding-window/docs/tabella_SR_CR_5MeV.csv)

## Branch e struttura del repository

`main` è l’unico branch pubblicato in `lsoffi/Xee-analysis`. Contiene entrambi i workflow. I workflow si trovano nelle directory [fixed-regions/](fixed-regions/README.md) e [sliding-window/](sliding-window/README.md). La directory [sliding-window/](sliding-window/) conserva lo scan sviluppato sul branch storico `sliding-window` e integrato in `main`. Per eseguire entrambe le analisi basta il checkout di `main`.

| Workflow su `main` | Dove iniziare | Documentazione |
|---|---|---|
| **Sliding window**, 69 ipotesi da 0,4 a 9,0 GeV. La SR varia con la massa e il fit usa le CR selezionate per quel punto. | [Istruzioni per lxplus e Condor](sliding-window/README.md) · [Finestre e CR per massa](sliding-window/scan_261106/records.json) | [Nota di analisi in inglese (PDF)](sliding-window/docs/Xee_analysis_note.pdf). Sezione 5: bias e trigger, copertura e GoF nella regione φ–J/ψ (GoF alle pagine 12–13). |
| **Regioni fisse**, spettro 6–12 GeV. Due analisi separate condividono la CR delle Υ. | [Workflow portabile](fixed-regions/scripts/workflow.py) · [Istruzioni di riproduzione](fixed-regions/docs/REPRODUCTION.md) | [Guida alle regioni fisse (PDF)](fixed-regions/docs/Xee_guida_SR_CR.pdf) · [Metodo](fixed-regions/docs/METODO.md) |

Nello sliding window, [cr_261106/](sliding-window/cr_261106/) contiene il codice dei fit delle risonanze SM, [scan_261106/](sliding-window/scan_261106/) contiene il fit simultaneo, i job Condor e gli studi di bias, e [provenance/](sliding-window/provenance/) contiene i checksum degli input.
