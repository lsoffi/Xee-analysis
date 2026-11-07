# Content uses helpers from build.py. All numerical tables are drawn from saved results.
for f in ['config/cr_definitive.json','runs/20261030_cr_optimization/REPORT.md','runs/20261030_cr_float_mu/REPORT.md','runs/20261030_cr_degree4_mu/REPORT.md','runs/20261030_sr_degree_selection/summary.json','runs/20261030_sr_bernstein7_bias/REPORT.md','runs/20261030_sr_cheb7_bern_bias/REPORT.md','runs/20261030_sr_upper_mass_bias/REPORT.md','runs/20261030_sr6_cr4_scan/summary.json','runs/20261030_sr6_cr4_scan/bias_summary.json','runs/20261030_sr6_cr4_scan/validation.json','runs/20261030_gof_profiled_combined/summary.json','runs/20261030_gof_thesis48/pooled.json','runs/20261029_signal_widths/common95.json','runs/20261030_gof_profiled_combined/REPORT.md','runs/20261030_sr6_cr4_scan/REPORT.md']:source(f)
source('runs/20261029_upsilon_peak_windows/summary.json')
source('/Users/imb/.codex/attachments/e3760e8e-42e3-4036-a278-5a903f7d6500/thesis_v322.pdf')
page('Fit simultaneo SR+CR: strategia e risultati','REPORT DI ANALISI', 'Ricerca X → ee | Quattro ere: 2022, 2022EE, 2023, 2023BPix | Stato al 30 ottobre 2026')
p('Report schematico e documentato della configurazione adottata: come sono stati scelti range e modelli, come vengono trasferite le code delle Υ e come entra il segnale in entrambe le regioni. Seguono risultati, controlli, limiti di validazione e atlante completo dei fit.',11)
table(['Configurazione adottata','Implementazione'],[
 ['Ipotesi di massa','31 punti: 6,0–9,0 GeV, passo 0,1 GeV; DCB originale del workspace'],
 ['SR unica','6,000000–9,093333 GeV; envelope Chebyshev6 / Bernstein6'],
 ['CR delle Υ','9,093333–10,506667 GeV; Chebyshev4; configurazione fissata'],
 ['Fit simultaneo','8 canali: SR e CR per ogni era. Segnale e Υ presenti in entrambe le regioni.'],
 ['Risultati principali','Massimo Z locale = 1,24 a 7,7 GeV; pGoF S+B = 0,47–0,68 su tutte le masse.']], [.25,.75],10)
box('La configurazione è definita; la validazione statistica non è conclusa. Il GoF è soddisfacente nel test combinato, ma alcuni bias superano 0,2σ; 10/124 fit di bias hanno flag numerici. Le bande attese dei limiti restano <b>provvisorie</b>, con estremi superiori troncati numericamente. Questi punti sono espliciti in tutto il report.','Come leggere le conclusioni',RED)
p('Nessun nuovo fit è stato eseguito per preparare questo documento. I numeri e i grafici provengono dalle campagne archiviate; le nuove figure descrittive sono derivate dai loro output. [S1–S14]',8.8,color=GRAY)

page('Mappa del report e stato delle conclusioni','GUIDA ALLA LETTURA')
table(['Parte','Contenuto e scopo'],[
 ['Strategia','Dati, regioni, likelihood, parametri condivisi e indipendenti.'],
 ['Scelte CR e SR','Percorso delle prove, evidenza numerica, scelta operativa e limiti delle conclusioni.'],
 ['Risonanze e bordi','Shape Υ, rese, medie, DCB di segnale, frazioni in CR, accettanza e contaminazione.'],
 ['Risultati','FitDiagnostics, bias, GoF, limiti CLs e p-value locale; tabelle per tutte le masse.'],
 ['Validazione e tracciabilità','Controlli effettuati, verifiche aperte, comandi e riferimenti agli artefatti.'],
 ['Appendici','Tabelle complete dei 124 bias e 62 tavole di fit: SR e CR a ogni massa, quattro ere per tavola.']], [.23,.77],9.6)
p('<b>Tre livelli da non confondere.</b> (1) <b>Scelta fissata:</b> range e forma della CR, SR unica ed envelope di grado 6. (2) <b>Controllo numerico:</b> convergenza, covarianza, positività, conteggio dei toy e integrità degli input. (3) <b>Validazione statistica:</b> robustezza rispetto al generatore, bias e copertura. Il superamento del livello 2 non implica il livello 3.')
p('<b>Tre probabilità diverse.</b> pGoF confronta la descrizione globale dei dati con i toy; p<sub>0</sub> quantifica un eccesso locale di segnale; CLs determina il limite su r. Il medesimo dataset può avere un buon GoF e una stima di segnale sensibile al modello di fondo.')
p('Nel lettore PDF usare i segnalibri per accedere direttamente a una sezione, a una massa e alle sue tavole SR/CR. Le pagine dei fit mantengono il formato vettoriale originale per consentire lo zoom.',9,color=GRAY)

page('Dati, binning e geometria delle regioni','CONFIGURAZIONE')
fig(TMP/'regions.png',150)
audit=json.loads((R/'m7.5/joint/audit.json').read_text())
table(['Era','Eventi SR','Eventi CR','Totale'],[[e,f'{audit["counts"][f"SR_ch{i}"]:,.0f}'.replace(',',' '),f'{audit["counts"][f"CR_ch{i}"]:,.0f}'.replace(',',' '),f'{audit["counts"][f"SR_ch{i}"]+audit["counts"][f"CR_ch{i}"]:,.0f}'.replace(',',' ')] for i,e in enumerate(eras,1)]+[['Totale','266 819','48 038','314 857']], [.25,.25,.25,.25],9,pad=4)
p('Si conserva il binning di partenza: 450 bin uniformi tra 6 e 12 GeV, ampiezza 13,333 MeV. La SR usa 232 bin e la CR 106: 338 bin per era, 1.352 osservazioni binned nel fit combinato. Nessun rebin nello scan finale. [S1,S9]',9.5)
p('I range nominali «6–9,1» e «9,1–10,5» sono allineati ai bordi esistenti: 9,093333 e 10,506667 GeV. Le regioni sono adiacenti e disgiunte: ogni bin entra una volta sola. Il tratto 10,506667–12 GeV non entra in questo fit. «Ipotesi di massa fino a 9» non significa che si usino dati soltanto fino a 9.',9.5)

page('Un’unica likelihood, non due fit indipendenti','MODELLO STATISTICO')
box('L(r,θ,j) = ∏<sub>era, regione, bin</sub> Pois(n | λ) × ∏<sub>k</sub> C<sub>k</sub>(θ<sub>k</sub>)<br/>λ = r · S<sup>ref</sup> · a<sub>sig,bin</sub> + B<sub>cont,regione</sub> · f<sub>cont,bin</sub> + Σ<sub>Υ1S,Υ2S</sub> κ<sub>Υ</sub> · N<sub>Υ</sub><sup>ref</sup> · a<sub>Υ,bin</sub>', 'Struttura del fit esteso')
p('θ raccoglie i parametri continui di fondo e le nuisance; j indica le scelte discrete del modello SR. I fattori C vincolano soltanto i parametri con incertezza assegnata.',9)
table(['Elemento','Legame tra SR e CR / tra ere'],[
 ['r, intensità del segnale','Un unico r per tutte le ere e per entrambe le regioni; il segnale cambia con l’ipotesi mX.'],
 ['Υ(1S), Υ(2S)','Per ciascuna era, stessa media μ e stessa scala di normalizzazione κ in SR e CR. Ere indipendenti.'],
 ['Continuo','Coefficienti e normalizzazioni indipendenti tra SR/CR e tra ere. Nessun trasferimento del continuo dalla CR.'],
 ['Nuisance del segnale','Vincoli originali del workspace conservati e profilati una sola volta nella likelihood combinata.'],
 ['Scelta discreta SR','Un indice Cheb6/Bern6 per era; 2⁴ = 16 combinazioni possibili profilate.']], [.23,.77],9.5)
p('Il fit alla CR usato in preparazione fornisce <b>valori iniziali</b>. Non si fissa il suo risultato come valore esatto nella SR e non si aggiunge una prior ricavata dagli stessi eventi. Nel fit finale le due regioni aggiornano insieme i parametri condivisi, propagandone le correlazioni. [S1,S9]',9.5)

page('Il percorso seguito: dalle anomalie al modello attuale','PERCORSO DELLE SCELTE')
table(['Passaggio','Domanda','Indicazione ottenuta'],[
 ['Fit globali iniziali 6–12','Perché il segnale modifica il continuo anche lontano dalla massa testata?','Polinomi globali e componenti Υ possono compensarsi: una variazione locale del segnale può cambiare coefficienti e forma al bordo.'],
 ['Range, binning, gradi e shape','Il comportamento dipende dal modello o dalla risoluzione?','Provati rebin, gradi minori, gaussiana/DCB, fit locali e omissione delle Υ. Risultati storici diagnostici, non baseline finale.'],
 ['Finestre al 95% della DCB','Come scegliere un range legato alla shape?','Intervallo minimo al 95% per era, unione tra ere e allineamento ai bin. Non è un criterio sufficiente per la flessibilità del fondo.'],
 ['Separazione SR+CR','Come usare i picchi Υ senza ignorarne le code?','CR intorno a entrambe le risonanze; parametri Υ condivisi con la SR.'],
 ['Ottimizzazione CR','Quanto restringere e quale continuo usare?','Scelta fissata: 9,0933–10,5067 GeV, Cheb4; μ e rese liberi, σ/code MC fisse.'],
 ['Ottimizzazione SR','Qual è il grado minimo con bias contenuto?','Confronti incrociati Bernstein7↔Cheb: grado6 minimo nei test condizionali richiesti.'],
 ['Scan finale','Come si comporta l’envelope con nuisance originali profilate?','31 masse, 124 bias Asimov, fit/limiti/p0 e GoF con r libero; risultati e criticità separati.']], [.18,.30,.52],9,pad=4)
p('La scelta non viene giustificata dalla scomparsa di un eccesso osservato: la motivazione deve poggiare su descrizione dei dati, stabilità e bias. Non tutti i test storici hanno la stessa likelihood o le stesse nuisance; non si sovrappongono le loro conclusioni senza specificarne le condizioni. [S2–S10]',9)

page('CR: cercare il range più stretto senza perdere informazione','SCELTA DEL RANGE CR')
p('La CR deve determinare le due Υ abbastanza bene da descriverne le code nella SR. Ridurla guadagna spazio di ricerca, ma può rendere poco distinguibili continuo, Υ(1S) e Υ(2S). Un buon GoF sotto i picchi non garantisce una buona estrapolazione delle rese. [S2,S3]')
table(['Studio','Disegno del test','Esito rilevante'],[
 ['Griglia di finestre','16 range ottenuti combinando bordi inferiori ≈9,093/9,2/9,293/9,4 e superiori ≈10,507/10,4/10,307/10,2 GeV.','Gradi1–3; shape Υ dapprima fisse, poi variazioni di μ e σ.'],
 ['Confronto a quantità fissa','Previsione delle code Υ negli intervalli fissi 6–9,0933 e 10,5067–12.','Non si cambia la quantità estrapolata quando si stringe la CR. Non si confrontano direttamente NLL di dataset diversi.'],
 ['Toy e generatore alternativo','200 toy per configurazione selezionata, GoF/LR e test di resa con Bernstein4 alternativo; 2 seed.','GoF anche buono ma bias di resa e problemi di copertura Hesse; nessun range minimo universalmente validato.'],
 ['Finestre dei soli nuclei','Valutate anche CR attorno ai singoli picchi e contenimenti DCB.','Le code rendono le finestre di contenimento larghe e sovrapposte. Una finestra stretta perde informazioni sulla separazione delle componenti.']], [.2,.43,.37],9.3)
box('Si mantiene la CR 9,093333–10,506667 GeV. È il compromesso operativo fissato dall’utente; non è presentato come minimo matematico né come soluzione che abbia già superato ogni test di estrapolazione.','Decisione sul range')

page('CR: perché solo μ liberi e perché Chebyshev4','SCELTA DEL MODELLO CR')
p('La prova con μ e σ delle due Υ liberi mostrava problemi di identificabilità. Congelando σ e code al MC e lasciando libere solo le medie, i fit di era numericamente validi sono passati da 128/192 a 171/192; 34/48 combinazioni hanno tutte le ere valide. Rimanevano parametri ai bordi in alcuni test. [S3]')
table(['CR [GeV]','Grado','Ere valide','Ere con μ al bordo','Devianza CR','Toy GoF4 validi'],[
 ['9,093–10,507','3','4/4','2','381,05','-'],['9,093–10,507','4','4/4','0','369,25','50/50'],
 ['9,200–10,400','3','4/4','0','316,18','-'],['9,200–10,400','4','4/4','0','305,51','36/50'],
 ['9,293–10,307','3','4/4','1','252,27','-'],['9,293–10,307','4','4/4','1','247,12','13/50']], [.22,.10,.14,.20,.14,.20],9,pad=4)
p('Nel range mantenuto, grado4 dà pGoF = 36/50 = 0,72 nel pilot e migliora la devianza di 11,803 rispetto al grado3. Il test LR aveva 45/50 toy validi, zero superamenti e 5 fallimenti: <b>non</b> significa p = 0, né consente da solo una decisione conclusiva al 5%. [S4]',9.5)
p('La scelta Cheb4 evita le instabilità più evidenti delle CR strette e migliora la separazione delle componenti. È stata poi fissata come configurazione definitiva. Coefficienti, μ e normalizzazioni restano liberi nel fit finale; «congelare la configurazione» non significa congelare i parametri.',9.5)

page('CR: cosa è migliorato e cosa resta da verificare','LIMITI DELLA VALIDAZIONE CR')
table(['Test nel range CR adottato','Cheb3','Cheb4','Lettura corretta'],[
 ['Toy Bernstein4: bias resa Υ sotto CR / σ','3,18','0,21','Miglioramento; il caso Cheb4 è vicino a una chiusura nella stessa famiglia polinomiale.'],
 ['Toy Bernstein5: bias resa Υ sotto CR / σ','3,40','0,73','Migliora, ma non raggiunge l’obiettivo esplorativo0,2.'],
 ['Toy Bernstein5: bias resa Υ sopra CR / σ','2,31','2,28','L’estrapolazione superiore non è validata da quella inferiore.'],
 ['Asimov Bernstein4: bias sotto/sopra, fit 4','-','≈0 / ≈0','La chiusura deterministica non elimina effetti di bordo e fluttuazioni nei toy.'],
 ['Asimov Bernstein5: bias sotto/sopra, fit 4','-','0,219 / 1,400','Rimane una componente di mismatch senza fluttuazioni.']], [.37,.09,.15,.39],9.2)
p('I bias di resa qui riguardano <b>CR-only</b>, con errori Hesse diagnostici, e in alcuni casi sono calcolati sui soli fit numericamente validi. Non sono direttamente il bias di r o dei limiti del fit simultaneo: la SR aggiunge informazione. Le percentuali di fallimento sono state conservate nei report originali. [S4]')
box('Restano separati dalla scelta fissata: profili di likelihood di μ e rese, robustezza con più seed/toy, identificabilità soprattutto della Υ(2S), generatori alternativi e copertura degli intervalli. La CR non viene riaperta in questo report.','Verifiche ulteriori, non cambi di configurazione',RED)

page('SR: perché una regione unica e come si definiscono i bordi','SCELTA DEL RANGE SR')
p('Le prime finestre mobili erano costruite dalle DCB nominali: intervallo più corto contenente il 95% del segnale, normalizzato su 6–12 GeV, per ogni era; poi minimo bordo inferiore e massimo bordo superiore tra ere. I fit arrotondavano verso l’esterno ai bin originali. Questo collegava la finestra alla risoluzione e alle code, senza assumere una gaussiana simmetrica. [S14]')
table(['Strategia storica','Configurazione attuale'],[
 ['Una finestra diversa per massa, con coefficienti adattati localmente.','Un’unica SR6–9,093333 GeV, stessi dati e grado del continuo per tutte le masse.'],
 ['Disgiunzione della finestra95% dalla CR: limitava la griglia storica a circa8,5GeV.','Disgiunzione dei <i>dati</i> SR/CR, non obbligo di tenere il segnale fuori dalla CR.'],
 ['Gradi bassi possono bastare localmente; confronti di eccessi sensibili a range e modello.','La SR larga richiede una scelta del grado mediante GoF e bias dedicati.'],
 ['Il 95% è un criterio geometrico, non dimostra ottimalità statistica del range.','La DCB entra anche in CR con lo stesso r; lo scan può arrivare a 9 GeV, con sensibilità degradata vicino al bordo.']], [.5,.5],9.6)
p('Il bordo inferiore6GeV deriva dal dominio degli input di questa regione; il superiore è il bordo inferiore della CR fissata. Le ipotesi6–9GeV sono la griglia richiesta ed effettivamente testata: <b>non</b> una dimostrazione che9GeV sia il massimo punto fisicamente analizzabile. Nessun risultato di questo scan copre la SR oltre10,5067GeV.',9.6)

page('SR: GoF e test di grado non bastano da soli','SELEZIONE DEL CONTINUO SR')
p('Con CR fissata a Cheb4, sono stati confrontati Chebyshev3–7 nella SR unica. Si cercano gradi bassi che descrivano i dati e non producano un segnale spurio nei test alternativi. Il grado8 come modello di fit è stato escluso dall’utente. [S5]')
table(['Grado SR','Devianza totale','Toy GoF validi','pGoF sui validi','Osservazione'],[
 ['3','1.415,97','100/100','0,020','Descrizione insufficiente in questo test.'],
 ['4','1.290,64','99/100','0,606','GoF buono; non garantisce bias contenuto.'],
 ['5','1.282,10','99/100','0,657','Un ulteriore grado cambia poco il GoF.'],
 ['6','1.279,88','97/100','0,567','Scelta da completare con test di segnale spurio.'],
 ['7','1.277,58','97/100','0,588','Maggiore flessibilità non è automaticamente migliore.']], [.13,.20,.20,.20,.27],9.3)
p('Il salto3→4 aveva ΔD≈125,32;0 superamenti su98 toy LR validi e2 fallimenti. Per4→5:ΔD≈8,55 e9/99 superamenti. Sono test preliminari con fallimenti espliciti; i p condizionati ai successi non vanno promossi a probabilità calibrate senza riserve.',9.5)
box('Il grado minimo suggerito dal GoF non coincide necessariamente con quello richiesto dal bias. La selezione successiva usa generatori alternativi e più ipotesi di massa, senza scegliere il grado perché riduce una significanza osservata.','Criterio adottato')

page('SR: scelta incrociata del grado minimo','BIAS CONDIZIONALE DI SELEZIONE')
p('Generatori fittati B-only ai dati: Bernstein7 → fit Chebyshev3–7 e, in senso inverso, Chebyshev7 → fit Bernstein3–7. CR invariata.21 masse6,5–8,5GeV; Asimov con r=0 e r pari al limite mediano atteso di uno scan <b>storico</b>, usato solo come ampiezza di riferimento. Nuisance del segnale tenute nominali in questa fase. [S6,S7]',9.7)
table(['Grado','Bern7→Cheb: massimo |bias|/σ','Cheb7→Bern: massimo |bias|/σ','Esito nei test condizionali'],[
 ['3','7,086 (42/42)','7,102 (36/42)','Non passa'],['4','1,736 (42/42)','1,550 (24/42)','Non passa'],['5','0,568 (42/42)','0,552 (28/42)','Non passa'],['6','0,135 (42/42)','0,113 (42/42)','Primo grado sotto0,2 in entrambe le direzioni'],['7','≈0 (42/42)','0,031 (42/42)','Controllo di chiusura più flessibile']], [.10,.27,.27,.36],9.3)
p('Il grado6 è quindi il minimo che supera il criterio Asimov nelle condizioni richieste. Questa conclusione non annulla altri scenari: con un precedente generatore Cheb8, il massimo bias di Cheb6 era circa0,36σ. Il generatore scelto conta. [S5–S7]',9.6)
p('Un pilot con30 dataset fluttuati, riusati sulle masse e sui gradi, dava perCheb6 un massimo bias assoluto≈0,45±0,20σ; il confronto appaiatoCheb6−Cheb7≈0,154±0,117σ. Le curve sono correlate e il campione non certifica |bias|&lt;0,2σ ovunque né la copertura.',9.6)

page('SR: evidenza grafica della selezione di grado','RISULTATO CONDIZIONALE, NON VALIDAZIONE FINALE')
fig(TMP/'selection_degrees.png',335,'Massimo bias normalizzato sulle21masse e sulle due ampiezze di iniezione. I punti di grado inferiore con fallimenti non costituiscono un test completo superato. Le due direzioni individuano grado6 nel test condizionale. [S6,S7]')
box('Il successivo scan Combine usa l’envelope e profila le nuisance originali del segnale. Non è statisticamente identico a questa selezione condizionale: i nuovi bias devono essere letti autonomamente.','Passaggio al modello finale')

page('Estensione della ricerca da 8,5 a 9 GeV','SEGNALI VICINI ALLA CR')
p('L’estensione ha mantenuto la DCB in SR e CR, con r comune.6 masse8,5–9,0GeV,2 generatori×5 gradi×2 iniezioni =120 fit Asimov. Le iniezioni sono r=0 e r=1 <b>fisso</b>, non il nuovo limite atteso. Nuisance del segnale nominali nel test storico. [S8]')
table(['Grado','Bern7→Cheb: max |bias|/σ','Cheb7→Bern: max |bias|/σ'],[['3','6,966','6,951'],['4','2,255','2,247'],['5','0,629','0,609'],['6','0,191','0,168'],['7','≈0','0,037']], [.2,.4,.4],9.5)
p('Tutti120 fit erano numericamente validi. Grado6 rimaneva sotto0,2σ nei due test condizionali. Contemporaneamente l’errore Hesse su r cresceva da≈0,17 a8,5GeV a≈1,07 a 9 GeV: bias ridotto non significa sensibilità uniforme.',10)
box('La massa di9GeV dista soltanto93,3MeV dal bordo della SR. Una parte consistente del segnale entra in CR; viene fittata, non eliminata. Il limite finale aumenta per l’insieme di accettanze, rese nominali e correlazioni con continuo/Υ. Non attribuiamo automaticamente tutto il degrado a una sola componente.','Bordo e sensibilità')

page('Envelope SR: cosa varia e cosa viene profilato','IMPLEMENTAZIONE FINALE')
table(['Elemento','Scelta e significato'],[
 ['Chebyshev6','T₀=1; sei coefficienti liberi per era. c₁∈[−2,2], c₂…c₆∈[−1,1]. Normalizzazione del continuo separata.'],
 ['Bernstein6','b₀=1; sei rapporti b₁…b₆∈[10⁻⁶,1000], tutti positivi. Nessuna ridondanza di scala con il yield separato.'],
 ['Positività','Cheb: controllo su257punti durante il fit; verifica esatta agli estremi e punti stazionari dopo i fit dati. Bernstein positivo per costruzione.'],
 ['Profiling discreto','RooMultiPdf per era, indice indipendente. Tutte16combinazioni profilate. I parametri della funzione inattiva sono disassociati, non nuove sorgenti di informazione.'],
 ['Penalità','Correzione RooMultiPdf uguale per le due opzioni:3,5 nella NLL per era. Costante totale28 in−2lnL; stesso trattamento in dati e toy.'],
 ['Continuo CR','Chebyshev4 e yield propri di ciascuna era, indipendenti dai coefficienti SR.']], [.22,.78],9.5)
box('Chebyshev e Bernstein allo stesso grado descrivono lo stesso spazio polinomiale, con domini diversi; i coefficienti Bernstein positivi lo restringono. L’envelope tratta la scelta di parametrizzazione/vincoli, ma non rappresenta due famiglie funzionali indipendenti e non esaurisce l’incertezza di forma.','Limite del confronto',RED)


page('Parametri attivi e nuisance originali','STRUTTURA DELLA LIKELIHOOD')
table(['Gruppo','Parametri nella configurazione attiva','Vincoli e correlazioni'],[
 ['Continuo SR','6 coefficienti + 1 yield per era','Liberi; nessun vincolo derivato dalla CR; scelta discreta indipendente per era.'],
 ['Continuo CR','4 coefficienti + 1 yield per era','Liberi e indipendenti dal continuo SR.'],
 ['Risonanze Υ','2 medie + 2 scale di resa per era','Condivise SR/CR, indipendenti fra ere; intervalli tecnici, non prior di misura.'],
 ['Segnale','Un unico r per tutte le ere e regioni','Intensità rispetto al benchmark originale; nominale DCB da MC.'],
 ['Efficienze e luminosità','CMS_eff_e_id, CMS_eff_e_reco, CMS_eff_trigger_e; lumi_1, lumi_2','Nodi e correlazioni originali del workspace conservati.'],
 ['Deformazioni del segnale','CMS_scale_e e 4 CMS_EXO25020_signalModelSigmaNuisance_ERA','Scala condivisa secondo il modello originale; risoluzione specifica per era; vincoli originali profilati.']], [.21,.40,.39],9.2)
p('Per una scelta attiva dell’envelope:16 parametri di fondo per era (64 nelle quattro ere), più r e le10 nuisance originali. Le coordinate della forma SR inattiva non sono una seconda componente sommata al continuo. Gli indici discreti selezionano una funzione per era, non una frazione di miscela.',9.5)
p('Gli otto parametri μ e le otto normalizzazioni Υ sono distinti dalle nuisance gaussiane originali del segnale. La likelihood finale contiene una sola copia di ogni vincolo originale. Non si aggiungono vincoli50%,20% o10% sulle Υ: la loro informazione arriva dai dati SR+CR. [S1,S9]',9.5)

page('Υ(1S) e Υ(2S): shape, normalizzazioni e trasferimento','FONDO RISONANTE')
p('Le Υ sono <b>fondi risonanti</b>; il polinomio rappresenta il fondo continuo. Entrambe le risonanze sono incluse esplicitamente in entrambe le regioni, comprese le code della2S sotto la1S e nella SR. Le shape originali sono considerate valide come richiesto. [S1,S9]')
table(['Parametro','Trattamento nel fit simultaneo'],[
 ['μ di ciascuna Υ','Libero per risonanza ed era; stesso parametro in SR e CR. Intervallo tecnico μMC±0,2GeV; non è una prior gaussiana.'],
 ['σ del nucleo e parametri di coda','Costanti ai valori MC originali, distinti per risonanza ed era. Non vengono rifittati nello scan.'],
 ['Scala κ della resa','Libera in[0,10] rispetto alla resa di riferimento6–12, per risonanza ed era; condivisa SR/CR. Nessuna prior50/20/10% adottata.'],
 ['Resa nella regione R','NΥ,R = κΥ × NΥ,6–12ref × AΥ,R(μΥ; σMC, codeMC). L’accettanza cambia con μ.'],
 ['Informazione della CR','I picchi aiutano a determinare μ e rese; la SR contribuisce nello stesso fit. L’incertezza e le correlazioni sono propagate profilando i parametri.']], [.28,.72],9.4)
p('Il fit preparatorio CR-only non diventa una misura esterna aggiuntiva sugli stessi dati. Non si duplica il suo vincolo nella SR. Analogamente, mantenere shape MC fisse è un’ipotesi del modello, non una dimostrazione che la loro incertezza sia nulla.',9.5)

page('Υ: esempio numerico nel fit S+B a 7,5 GeV','PARAMETRI POST-FIT, ERRORI HESSE')
pars=list(csv.DictReader((R/'m7.5/joint/parameters.csv').open()));pd={x['name']:x for x in pars if x['stage']=='s'}
yt=[]
for e in eras:
 for pr,lab in [('upsilon1s','Υ(1S)'),('upsilon2s','Υ(2S)')]:
  mu=pd[f'Ymean_{pr}_{e}'];k=pd[f'CMS_EXO25020_bkgScale_{pr}_{e}'];yt.append([e,lab,f'{float(mu["value"]):.4f} ± {float(mu["error"]):.4f}',f'{float(k["value"]):.3f} ± {float(k["error"]):.3f}'])
table(['Era','Componente','μ [GeV]','κ rispetto alla resa nominale6–12'],yt,[.18,.17,.30,.35],9.8)
p('Questi errori sono Hesse locali, non intervalli profilati. La Υ(2S) è meno precisamente determinata in alcune ere. Non si deve fissarne la resa nella SR ignorando questa incertezza. I valori della tabella non sono imposti agli altri punti di massa: il fit viene ripetuto a ogni ipotesi.',10)
p('Controllo aggiunto ai risultati salvati: nelle31soluzioni S+B finali nessuna delle16variabili μ/κ delle Υ è entro10⁻⁴ dell’ampiezza del proprio intervallo dal bordo. Questo è un controllo tecnico; non sostituisce profili, identificabilità e copertura. Tabelle complete dei parametri B-only/S+B nelle directory di massa. [S9]',9.4)

page('Segnale DCB e contaminazione della CR','TRATTAMENTO DEL SEGNALE')
box('Nsig,R(r,θ,mX) = r × Nsig,6–12ref(θ,mX) × A<sub>R</sub>(θ,mX)<br/>A<sub>R</sub> = ∫<sub>R</sub> DCB(m; mX,θ)dm / ∫<sub>6–12</sub> DCB(m; mX,θ)dm','Accettanza e normalizzazione')
p('La DCB nominale è quella originale del workspace alla massa mX; non una gaussiana e non una shape appresa dai dati di ricerca. I parametri nominali MC e le nuisance originali di scala/risoluzione sono conservati. «Shape da MC» non significa congelare tutte le nuisance del segnale nel fit finale. [S9]')
p('Le rese di segnale in SR e CR sono derivate dalla stessa normalizzazione e dallo stesso r. Le accettanze sono ricalcolate quando le nuisance che deformano la DCB cambiano. Primitive analitiche e integrali ROOT sono stati confrontati nominalmente e con variazioni di nuisance.')
p('<b>Contaminazione della CR</b> significa che un segnale di ricerca può produrre eventi anche nei bin usati per vincolare le Υ. Non si presume contaminazione nulla: la DCB è inclusa esplicitamente in CR. Ometterla permetterebbe al segnale di modificare rese/medie Υ e continuo senza essere attribuito coerentemente a r.')
p('Le frazioni riportate sono riferite al segnale <b>nel dominio originale6–12GeV</b>, non a tutta la DCB sull’asse reale. Le rese di riferimento del workspace vengono preservate; non si inventano correzioni per eventi fuori dal suo dominio. Gli eventi oltre10,506667 non sono utilizzati, e la loro quota non è rinormalizzata nella SR.',9.5)

page('Frazioni di segnale: cosa accade verso il bordo alto','CONTAMINAZIONE E ACCETTANZA')
fig(TMP/'acceptance.png',285,'Frazioni nominali per era derivate dagli audit degli input, normalizzate su6–12GeV. La CR contiene una quota crescente della DCB; questa quota entra nella likelihood con lo stesso r. [S9]')
rows=[]
for m in [6.0,7.5,8.5,8.8,8.9,9.0]:
 a=[v for v in S if v['mass']==m];rr=[]
 for k in ['SR_fraction','CR_fraction','unused_fraction']:rr.append(f'{100*min(x[k] for x in a):.2f}–{100*max(x[k] for x in a):.2f}%')
 rows.append([f'{m:.1f}',*rr])
table(['mX [GeV]','In SR','In CR','Oltre10,5067, non usato'],rows,[.15,.28,.28,.29],8.9,pad=3)

page('Segnali ai bordi: protezioni e limiti reali','INTERPRETAZIONE DELLA COPERTURA GEOMETRICA')
table(['Problema','Trattamento adottato','Limite residuo'],[
 ['Taglio della DCB al bordo SR','Rese trasferite con integrali; segnale incluso anche in CR.','Un buon contenimento SR+CR non garantisce che r sia separabile dal fondo.'],
 ['Degenerazione con Υ','μ e rese condivisi; r comune; σ/code Υ MC fisse.','Le correlazioni devono essere profilate; non si risolvono congelando a posteriori le normalizzazioni.'],
 ['Massa 9 GeV vicina alla CR','In CR finisce15,87–17,17% del segnale nominale6–12; in SR82,57–83,82%.','I limiti diventano meno stringenti: r95 osservato≈7,91 contro≈1,11 a8,5GeV.'],
 ['Bordo inferiore6GeV','Si usa la definizione di segnale e resa del workspace nel dominio6–12.','I nuovi bias completi a6,0 e6,1GeV superano0,2σ anche con rtrue=0.'],
 ['Ricerca sopra la CR','Non inclusa nella campagna6–9GeV presentata qui.','Il comportamento dell’estrapolazione superiore CR-only non autorizza uno scan della SR alta.']], [.22,.41,.37],9.5)
box('Avere signal leakage modellato in CR evita un’omissione della likelihood. Non elimina la perdita di sensibilità, l’incertezza di forma o il rischio di segnale spurio. Il criterio operativo è controllare bias, profili e copertura nel modello effettivamente usato.','Cosa significa “trattare la contaminazione”')

page('Esecuzione dello scan: cosa è stato profilato','PROCEDURA COMPUTAZIONALE')
table(['Metodo','Configurazione numerica finale','Controlli effettuati'],[
 ['FitDiagnostics','31fit B-only e31S+B. r può essere negativo; range adattato per massa. Rifit strategia2, tolleranza0,001, precisione10⁻¹² ove promossi.','31/31 punti con status 0/covQual 3 in B e S; r non al bordo; continuo positivo; componenti e normalizzazioni salvate.'],
 ['AsymptoticLimits','95%CLs; strategia2, tolleranza10⁻⁵; rRelAcc0,001, rAbsAcc10⁻⁵.','31crossing osservati vicini aCLs=0,05. Bande attese superiori non validate.'],
 ['Significance','r POI, strategia1, tolleranza0,001; test locale di eccesso.','Confronto con variazioneNLL di FitDiagnostics; tutti coerenti entro il criterio impostato.'],
 ['Bias Asimov','2generatori×31masse×2iniezioni; nuisance originali profilate nel fit.','114/124 status 0/covQual 3;10flag residui espliciti.'],
 ['GoF finale','Combinato quattro ere; saturated, strategia0, tolleranza default, r libero[-100,100];100toy/massa.','3.100 score finiti,100 ID/massa,nessunflaglog o rtoy al bordo. Non equivale a audit globale dei minimi.']], [.18,.43,.39],9.1)
p('In tutti i metodi del modello finale: profiling delle16combinazioni dell’envelope e congelamento dei parametri disassociati dalla funzione inattiva. Le impostazioni del primo comando storico (p.es. Υ1S congelata in FitDiagnostics) non descrivono questa configurazione finale: qui entrambe le normalizzazioni Υ sono libere. [S9,S11]',9.3)

page('Come leggere le tavole dei fit','GUIDA ALL’ATLANTE')
p('L’atlante finale contiene <b>62 pagine vettoriali</b>: una SR e una CR per ciascuna delle31masse. Ogni pagina mostra separatamente2022,2022EE,2023,2023BPix. La legenda è conservata insieme al binning originale e al pannello dei residui.')
table(['Elemento del grafico','Significato'],[
 ['Punti neri','Dati binned; barre visuali √N.'],
 ['Curva totale S+B','Somma post-fit del continuo, Υ(1S), Υ(2S) e segnale alla massa testata.'],
 ['Curva B-only tratteggiata','Risultato del fit con r=0, per confronto visivo.'],
 ['Componenti','Continuo SR dell’envelope o Cheb4 in CR, ciascunaΥ, DCB del segnale. Il segnale può avere resa negativa se r fitted&lt;0.'],
 ['Pannello inferiore','(D−M)/√M con M totale post-fit. È un residuo visuale, non il GoF saturated e non un pull che includa tutta la covarianza del fit.'],
 ['Scala e normalizzazioni','Eventi per13,33MeV; ordinate lineari. Le componenti possono compensarsi senza produrre grosse differenze nel totale.']], [.25,.75],10)
box('Per vedere come cambia il fit al bordo, confrontare le tavole7,5;7,7;8,5;9,0GeV, siaSR siaCR. La CR va guardata insieme alla SR perché il segnale e le rese Υ sono condivisi. I valori r, Z e pGoF sono richiamati nel piè di pagina di ogni tavola.','Percorso di lettura suggerito')

page('Bias finale: definizione e disegno del test','SCAN COMBINE CON NUISANCE PROFILATE')
box('bias assoluto = r̂ − rtrue<br/>bias normalizzato = (r̂ − rtrue) / σ<sub>Hesse</sub>(r)<br/>Il criterio |bias|/σ&lt;0,2 è un obiettivo diagnostico scelto per il confronto, non una legge statistica.','Quantità riportate')
p('Sono stati generati Asimov dai fit B-only Bernstein7 e Chebyshev7, con CR4. Ogni generatore è fittato con l’envelope Cheb6/Bern6. Si coprono31masse6–9GeV e due ampiezze: <b>rtrue=0 e rtrue=1</b>. Totale124fit. Il segnale è presente anche in CR. [S9,S10]')
p('r=1 è il benchmark del workspace. <b>Non</b> è il limite atteso specifico per massa e non è10volte quel limite. Le iniezioni storiche basate sul limite atteso appartengono ai test condizionali di selezione, non a questa campagna. Le precedenti campagne interrotte non vengono usate come risultati finali.')
p('La generazione Asimov non include fluttuazioni statistiche; il fit profila le nuisance originali. La differenza rispetto agli studi con nuisance nominali può cambiare sia il minimo sia l’errore e quindi il bias normalizzato. Non si può trasferire automaticamente la validazione condizionale del grado6 a questo modello completo.')
p('Si riporta bias/σ come quantitativo affidabile solo per status 0/covQual 3 e r non al bordo. I flag residui sono mostrati come croci, non eliminati dal resoconto. Una media dei soli successi non certifica la campagna intera; nessuna prova di copertura è inclusa qui.',9.5)

page('Bias finale: tutte le masse e le due iniezioni','RISULTATI ASIMOV')
fig(R/'bias.png',370,'Righe:rtrue=0 e1. Colonne:bias normalizzato e assoluto. Blu:generatoreBern7; arancio:Cheb7. Banda grigia:±0,2σ. Croci:verifica numerica aperta. I valori r sono intensità rispetto al riferimento, non eventi. [S10]')
p('Il bias con rtrue=0 supera0,2σ ai bordi6,0;6,1;9,0GeV per entrambi i generatori. Con rtrue=1 vi sono superamenti anche all’interno. Un GoF buono non cancella questi risultati.',9.5,color=RED)

page('Bias finale: riepilogo dei superamenti e dei flag','ESITI E LIMITAZIONI')
summary=[]
for fam,label in [('bern7','Bernstein7'),('cheb7','Chebyshev7')]:
 for rt in [0,1]:
  q=[v for v in B if v['family']==fam and v['rtrue']==rt];valid=[v for v in q if v['valid']];worst=max(valid,key=lambda v:abs(v['bias_sigma']));outside=[f'{v["mass"]:.1f}' for v in valid if abs(v['bias_sigma'])>.2];summary.append([label,str(rt),f'{len(valid)}/31',f'{abs(worst["bias_sigma"]):.3f} a {worst["mass"]:.1f}',', '.join(outside)])
table(['Generatore','rtrue','Validi','Max |bias|/σ nei validi','Masse oltre0,2σ [GeV]'],summary,[.18,.09,.1,.22,.41],9.2)
flags=[v for v in B if not v['valid']]
table(['Generatore','Iniezione','Masse con flag numerico'],[[fam,'r=1',', '.join(f'{v["mass"]:.1f}' for v in flags if v['family']==fam)] for fam in ['bern7','cheb7']], [.2,.15,.65],9.3)
p('Tutti i10flag residui appartengono all’iniezione r=1. I124valori, errori, status e covQual sono riportati nelle quattro tabelle analitiche in appendice. Ricalcoli numerici e seed migliorati sono conservati; non si sono cambiati range o gradi per rimuovere un bias sfavorevole. [S10]',9.6)
box('La scelta del grado6 resta quella adottata. La campagna completa non supporta l’affermazione “bias sotto0,2σ su tutto6–9 per tutte le ampiezze”. Occorrono test ulteriori se questo è il requisito per un risultato fisico finale.','Conclusione sul bias',RED)

page('Goodness of fit: procedura finale, come nel comando fornito','GoF SATURATED CON r LIBERO')
p('Per ciascuna massa si esegue il fit GoF sui dati e su100toy frequentist generati dal modello B-only fittato ai dati. Durante la valutazione del GoF si usa il modello <b>S+B con r libero</b>, su dati e toy. È la distinzione chiarita dal comando dell’utente: manca --fixedSignalStrength 0. [S11]')
box('t = −2 ln[L(θ̂,r̂) / L<sub>sat</sub>]<br/>pGoF(mX) = N(ttoy &gt; tobs) / 100<br/>Il modello saturated assegna a ciascun bin la sua aspettativa indipendente.','Definizione operativa')
p('La massa entra attraverso la DCB nel fit S+B: cambia dove il modello può assorbire una fluttuazione. Generare toy B-only non significa fissare r=0 quando li si rifitta. La distribuzione dello score viene quindi calibrata a ogni massa; non si aggregano i toy di masse diverse per ottenere un unico pGoF.')
p('Combine viene eseguito con strategia0 e tolleranza default, r[-100,100] ereditato dai workspace; tutte16combinazioni discrete profilate.100toy per massa,31masse,3.100 score. Gli intervalli mostrati sono Clopper–Pearson95% per l’incertezza Monte Carlo del p-value, non bande di sensibilità del segnale.')
p('Il salvataggio dei toy aggiunge nel TTree righe di generazione con quantileExpected=−2. Sono state escluse: si usano soltanto score GoF con quantileExpected=−1 e100 ID distinti. La penalità costante dell’envelope è trattata allo stesso modo in dati e toy.',9.4)

page('GoF combinato: risultato su tutte le masse','RISULTATO FINALE')
fig(G/'gof_vs_mass.png',350,'Combinato quattro ere, SR+CR:100toy per massa, fit S+B con r libero. Tutti i pGoF sono compresi tra0,47 e0,68 e superano0,05. Barre:precisione Monte Carlo95%. [S11]')
p('Valori guida:6,1GeV→0,47;7,5→0,60;7,7→0,50;7,8→0,68;9,0→0,59. Nessun errore rilevato nei log finali, tutti 3.100 score finiti e nessun rtoy al bordo. Non è stato eseguito un multistart sistematico di ogni toy; il test non dimostra assenza di bias o copertura dei limiti.',9.4)


page('GoF: le distribuzioni dei toy','DUE ESEMPI DEL TEST FINALE')
twofig(G/'m7.5/gof.png',G/'m9.0/gof.png',235,'A sinistra 7,5 GeV:60score su100 superano l’osservato; a destra 9 GeV:59su100. I toy sono generati B-only e rifittati con il segnale libero alla rispettiva massa. Non si usa una distribuzione asintotica χ² per attribuire questi p-value. [S11]')
p('Gli istogrammi mostrano la distribuzione nulla empirica della statistica saturated; la linea verticale è il risultato sui dati. La porzione oltre la linea determina pGoF. La risoluzione della frequenza empirica è1/100 e gli intervalli binomiali quantificano l’incertezza del campione finito.')
p('Le fluttuazioni della curva pGoF fra masse includono sia la dipendenza del fit S+B da mX sia la statistica Monte Carlo. Con100toy per punto non è corretto attribuire significato fisico a differenze di pochi centesimi. Le distribuzioni per tutte 31masse sono salvate nella campagna S11.',9.7)

page('GoF: confronto con i controlli precedenti','PROBABILITÀ DA NON CONFONDERE')
table(['Controllo','r nel fit','Risultato','Interpretazione'],[
 ['Combinato B-only storico','Fissato a0','p=0,47,100toy;MC95%0,369–0,572','Un solo test comune alle masse: modello e dati non dipendono da mX.'],
 ['GoF finale combinato','Libero','p(mX)=0,47–0,68;100toy/massa','Procedura del comando fornito; dipendenza reale dalla massa attraverso il segnale nel fit.'],
 ['B-only per era','Fissato a0','2022:0,869;2022EE:0,444;2023:0,868;2023BPix:0,0468','31campioni da100toy dello stesso test per era, aggregati in3.100toy/era.'],
 ['p0 di significanza','Test di r=0 contro eccesso','Minimo p0≈0,107 a7,7GeV','Probabilità locale di una fluttuazione almeno altrettanto compatibile con un eccesso; non misura la qualità globale.']], [.24,.12,.30,.34],9.2)
box('2023BPix è vicino alla soglia nel controllo B-only separato: p≈0,0468, intervallo MC95%[0,0396;0,0548]. Il buon GoF combinato non autorizza a dire che ogni era passi chiaramente ogni controllo. È un indizio da approfondire con residui e stabilità; non identifica una massa specifica problematica.','Informazione che il combinato può rendere meno evidente',RED)
p('Il riferimento alla figura4.8 della tesi riguarda la procedura di GoF, non l’adozione dei suoi range o dei suoi numeri. La figura mostrava2022; l’utente ha poi richiesto esplicitamente il combinato delle quattro ere e100toy/massa. [S11,S12,S15]',9.4)

page('Limiti al 95% CLs: risultato e stato delle bande','LIMITE SULL’INTENSITÀ r')
fig(R/'limits.png',350,'Limite osservato e mediano atteso in unità di r del workspace. Bande1/2σ originali mantenute come richiesto, chiaramente PROVVISORIE. [S9]')
p('<b>Problema noto:</b> la ricerca numerica dei quantili superiori attesi ha restituito estremi troncati dal bracket senza verifica del crossing. Il ricalcolo è stato interrotto per esplicita richiesta dell’utente. Non usare le bande superiori per conclusioni quantitative. I31limiti osservati hanno un controllo separato del crossingCLs≈0,05. Il mediano atteso resta quello nominale dello scan.',9.4,color=RED)

page('Come interpretare i limiti e il degrado al bordo','CLs, NORMALIZZAZIONE E ACCETTANZA')
p('Il limite è calcolato con AsymptoticLimits al 95%CLs. La curva è un limite su <b>r</b>, intensità rispetto alla resa di riferimento del workspace a ogni massa, con nuisance e indice dell’envelope profilati. r non è un numero di eventi; senza la definizione del benchmark non è una sezione d’urto assoluta.')
table(['mX [GeV]','r95 osservato','Mediano atteso nominale','Segnale in CR, frazione6–12'],[[f'{m:.1f}',fmt(next(v for v in A if v['mass']==m)['limits']['observed'],5),fmt(next(v for v in A if v['mass']==m)['limits']['expected_50'],5),f'{100*min(v["CR_fraction"] for v in S if v["mass"]==m):.2f}–{100*max(v["CR_fraction"] for v in S if v["mass"]==m):.2f}%'] for m in [6.0,7.0,7.5,7.7,8.5,8.8,8.9,9.0]],[.16,.26,.29,.29],9.4,pad=4)
p('L’aumento del limite verso9GeV riflette la sensibilità del modello completo: r viene correlato con continuo e risonanze, e le rese di riferimento e le accettanze cambiano con la massa. Non è una perdita automatica di tutti gli eventi entrati nella CR: quei bin restano nel fit.',9.8)
p('Non è stata misurata la copertura frequentista dei limiti in questa campagna. Il controllo del crossing osservato verifica l’esito numerico della ricercaCLs, non la correttezza del modello di generazione né dell’approssimazione asintotica.',9.5)

page('p-value locale e significanza del segnale','SCAN DELLE IPOTESI DI MASSA')
fig(R/'pvalues.png',350,'p0 locale da Combine Significance con lo stesso modello simultaneo. Non è un p-value globale corretto per lo scan in massa. [S9]')
p('Il massimo è Z=1,241 a7,7GeV, p0=0,1073. A7,8GeV Z=1,102, p0=0,1352; a 7,5 GeV Z=0 e p0=0,5. Il test di eccesso è unilaterale: una stima r negativa non produce una significanza positiva. Non emerge un eccesso locale rilevante nello scan corrente; ciò non costituisce da solo validazione del fondo.',9.6)

rows=[]
for a in A:
 g=gb[a['mass']];rows.append([f'{a["mass"]:.1f}',fmt(a['r'],4),fmt(a['Z'],4),fmt(a['p0'],4),fmt(a['limits']['observed'],5),fmt(a['limits']['expected_50'],5),f'{g["p"]:.2f}',f'{g["CP95"][0]:.3f}–{g["CP95"][1]:.3f}'])
for part,rr in enumerate([rows[:16],rows[16:]],1):
 page(f'Risultati numerici completi: parte {part}/2','FIT, LIMITE, SIGNIFICANZA E GoF')
 table(['mX GeV','r̂ FD','Z locale','p0','r95 obs','r95 exp50*','pGoF','MC95% pGoF'],rr,[.09,.12,.10,.11,.14,.14,.09,.21],9.0,pad=4)
 p('*Mediano nominale; bande attese superiori provvisorie. r̂ è quello di FitDiagnostics con i controlli numerici finali; il GoF usa il proprio fit con strategia 0 e può dare piccoli scostamenti numerici di r̂. [S9,S11]',9.0,color=GRAY)

page('Cosa risulta verificato, cosa non è ancora dimostrato','BILANCIO DELLA VALIDAZIONE')
table(['Aspetto','Stato','Conseguenza'],[
 ['Geometria e dati','Verificati','Binning originale, regioni disgiunte, stessi dati alle31masse, input preservati via SHA256.'],
 ['PDF e trasferimenti','Verificati numericamente','Integrali DCB confrontati con ROOT; parametri Υ condivisi; σ/code MC fisse; positività controllata.'],
 ['Fit dati e p0','31/31 controlli superati','B/S status 0,covQual 3; r non al bordo; consistenzaNLL/Significance. Nessun μ/κΥ vicino ai limiti tecnici nel controllo salvato.'],
 ['GoF finale combinato','Soddisfacente nel test eseguito','pGoF 0,47–0,68.100toy/massa danno una precisione limitata; non multistart globale dei toy.'],
 ['Bias completo','Non supera tutti i criteri','Superamenti0,2σ e10flag numerici. Non confondere con selezione Asimov condizionale passata.'],
 ['Bande attese','Provvisorie','Estremi superiori troncati numericamente; non usare per valutazioni quantitative.'],
 ['Copertura/profili','Non completati','Servono toy dedicati per limiti e intervalli, soprattutto ai bordi e per Υ2S.'],
 ['Robustezza della forma','Condizionata ai generatori provati','Envelope a pari grado non copre tutte le forme plausibili; non prova universale.']], [.23,.24,.53],9.2,pad=4)

page('Conclusioni operative e verifiche ulteriori','CHIUSURA DEL REPORT')
box('La configurazione adottata è: SR unica6–9,093333GeV con envelope Cheb6/Bern6; CR9,093333–10,506667GeV con Cheb4; μ e rese Υliberi e condivisi; σ/codeΥMCfisse; DCB originale con nuisance e r comune inSR/CR. La presenza del segnale in CR è parte del modello finale.','Strategia fissata')
p('<b>Risultato sui dati.</b> I fit superano i controlli numerici eseguiti e il GoF combinato è buono. Lo scan non mostra un eccesso locale significativo. I limiti osservati sono disponibili a tutte 31masse; le bande attese vanno considerate provvisorie.')
p('<b>Da verificare senza riaprire automaticamente la configurazione:</b> risolvere o comprendere i10flag di bias; verificare il bias completo con ampiezze motivate e più generatori; profili di r, μ e rese Υ; copertura dei limiti; robustezza dei minimi con seed alternativi; residui e comportamento di2023BPix; validazione delle code e delle incertezze MC se richiesta dall’analisi fisica.')
p('<b>Non inferire:</b> “GoF alto ⇒ nessun bias”, “grado6 passa ogni generatore”, “CR definitiva ⇒ tutti i suoi parametri noti”, oppure “segnale in CR ⇒ perso”. Queste equivalenze non sono supportate dai test.')
p('Le tabelle seguenti riportano tutti i124fit di bias, inclusi quelli con flag. L’atlante finale mostra ogni fit SR/CR. Gli artefatti e i comandi originali restano disponibili nelle directory elencate nei riferimenti.',10)

page('Riproducibilità: comandi e dettagli essenziali','APPENDICE TECNICA')
p('Percorso del progetto: <font face="Mono" size="8.3">/Users/imb/Documents/Codex/Analisi-statistica-Xee</font>. Input originali: <font face="Mono" size="8">data/mass_scan/261028/cards_region2_data_envelope_allCorrections_binned/ee/MASS/Xee_ee_4_allYears.root</font>. Runtime usato: ROOT6.34.10, Combine11.1 tramite wrapper locali. [S9,S11]',9.5)
table(['Operazione','Schema del comando; opzioni complete in command.json per massa'],[
 ['FitDiagnostics','combine -M FitDiagnostics -d workspace.root --saveShapes --saveNormalizations --keepFailures --robustFit 0 --minos none --cminDefaultMinimizerStrategy 2 --cminDefaultMinimizerTolerance 0.001; range r specifico.'],
 ['Limiti','combine -M AsymptoticLimits -d workspace.root --cminDefaultMinimizerStrategy 2 --cminDefaultMinimizerTolerance 0.00001 --rRelAcc 0.001 --rAbsAcc 0.00001; rMax specifico.'],
 ['Significanza','combine -M Significance -d workspace.root --redefineSignalPOIs r --cminDefaultMinimizerStrategy 1 --cminDefaultMinimizerTolerance 0.001; range r specifico.'],
 ['GoF osservato','combine -M GoodnessOfFit workspace.root --algo saturated --cminDefaultMinimizerStrategy 0; senza fixedSignalStrength.'],
 ['GoF toy','Stesso comando più -t 100 --toysFrequentist --saveToys -s SEED --genBinnedChannels SR_ch1,…,CR_ch4 --trackParameters r.']], [.19,.81],9.2)
p('Gli schemi non sostituiscono gli array esatti dei file command.json, che includono seed, nome di output, range r e opzioni specifiche per massa. Tutti i metodi includono --X-rtd MINIMIZER_freezeDisassociatedParams e --cminRunAllDiscreteCombinations. I risultati non sono stati ricalcolati per il PDF.',9.1,color=GRAY)

page('Riferimenti agli studi e agli artefatti','TRACCIABILITÀ 1/2')
refs=[
 ('S1','Configurazione CR definitiva','config/cr_definitive.json'),
 ('S2','Ottimizzazione range CR: gradi1–3, shapeΥfisse','runs/20261030_cr_optimization/REPORT.md'),
 ('S3','CR con soli μ liberi; confronto con σ floating','runs/20261030_cr_float_mu/REPORT.md'),
 ('S4','CR Cheb4:range,stabilità e bias delle rese','runs/20261030_cr_degree4_mu/REPORT.md'),
 ('S5','Prima selezione SR:GoF,LR e generatoriCheb7/8','runs/20261030_sr_degree_selection/summary.json'),
 ('S6','Bernstein7→Cheb3–7,Asimov e pilot30toy','runs/20261030_sr_bernstein7_bias/REPORT.md'),
 ('S7','Cheb7→Bernstein3–7,Asimov','runs/20261030_sr_cheb7_bern_bias/REPORT.md'),
 ('S8','Estensione condizionale8,5–9GeV','runs/20261030_sr_upper_mass_bias/REPORT.md')]
table(['ID','Contenuto','Percorso relativo al progetto'],refs,[.07,.35,.58],9,pad=6)
p('Le conclusioni storiche sono riportate con le condizioni del rispettivo studio. I report precedenti possono descrivere uno stato “candidato” prima della scelta esplicita dell’utente; per la configurazione attiva fa fede S1, integrata dalla sceltaSRenvelope dello scan S9.',9.5)

page('Riferimenti, dati tabellari e atlante','TRACCIABILITÀ 2/2')
refs=[('S9','Scan finale:fit/limiti/p0,componenti,qualità','runs/20261030_sr6_cr4_scan/REPORT.md; summary.json; validation.json; mMASS/joint/parameters.csv'),('S10','124bias Asimov e flag completi','runs/20261030_sr6_cr4_scan/bias_summary.json; bias_summary.csv; bias/'),('S11','GoF combinato r libero,100toy/massa','runs/20261030_gof_profiled_combined/REPORT.md; summary.json; mMASS/observed/ e toys/'),('S12','Controlli B-only perera','runs/20261030_gof_thesis48/REPORT.md; pooled.json'),('S13','Geometria finestreΥ e signal leakage','runs/20261029_upsilon_peak_windows/summary.json'),('S14','Finestre storiche 95%DCB','runs/20261029_signal_widths/common95.json'),('S15','Tesi allegata,fig4.8 e sez4.2.3','thesis_v322.pdf, pagina stampata109; comandoGoF chiarito direttamente dall’utente nella conversazione')]
table(['ID','Contenuto','Percorso / riferimento'],refs,[.07,.35,.58],8.1,pad=3)
p('Tutti i fit dell’atlante provengono dai PDF vettoriali mMASS/fit_SR.pdf e fit_CR.pdf dello scan S9, rigenerati dai FitDiagnostics finali durante la campagna. Le tavole mantengono legenda e componenti originali; la fascia aggiunta riporta r̂, Z locale e pGoF finale.',9.4)
p('Il manifest associato al PDF registra sorgenti e SHA256, indice delle pagine e corrispondenza massa/regione. I file di input e i risultati originali non sono modificati dalla produzione del report.',9.3,color=GRAY)

for fam,label in [('bern7','Bernstein7'),('cheb7','Chebyshev7')]:
 for rt in [0,1]:
  page(f'Bias analitici: generatore {label}, rtrue={rt}','APPENDICE: TUTTI I124FIT')
  q=sorted([v for v in B if v['family']==fam and v['rtrue']==rt],key=lambda v:v['mass']);rows=[]
  for v in q:
   rows.append([f'{v["mass"]:.1f}',fmt(v['rhat'],5),fmt(v['error'],5),fmt(v['bias'],5),fmt(v['bias_sigma'],4),f'{v["status"]}/{v["covQual"]}','OK' if v['valid'] else 'FLAG'])
  table(['mX GeV','r̂','σHesse(r)','r̂−rtrue','bias/σ','status/covQual','Numerica'],rows,[.10,.17,.17,.17,.15,.15,.09],7.7,pad=1.6)
  p('OK:controlli numerici del fit superati; non significa |bias|/σ&lt;0,2. FLAG:errori e bias normalizzati solo diagnostici, non validati. Nessun toy fluttuato né test di copertura in questa tabella. [S10]',8.4,color=GRAY)
