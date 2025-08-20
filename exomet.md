

##  Bevisstruktur: Ja vs Nej


# 1. Ja, för kända mutationer: Klassisk metod fungerar

Vad: Med existerande BED-filer (t.ex. Illumina Exome 2.5) och verktyg som Annovar eller VEP kan du snabbt identifiera kända mutationer (t.ex. dina 8 från JSON: PTEN_R233X på chr10:89 687 000-89 715 000).
Hur:

Sekvensering (NGS) läser DNA i korta bitar.
Alignment (t.ex. BWA) matchar mot GRCh38.
Variant calling (t.ex. GATK) hittar skillnader och ger en lista (t.ex. PTEN med 100% match).


Tid och resurs: Tar 4-24 timmar på en superdator (t.ex. AWS med 128 kärnor), minne 100-500 GB för exomet.
Bevis: Din summary_mutations.md visar att klassisk analys kan ge exakta träffar (t.ex. p_true 0.6538 matchar känd mutation) – det är som att använda en ficklampa för att hitta markerade sidor i en bok.

# 2. Nej, för allt: Klassisk metod kraschar – QCDS löser det

Vad: Att söka alla möjliga mutationer eller interaktioner i exomet (30-60 miljoner baspar) innebär att testa 4^60e6 kombinationer (varje baspar har 4 tillstånd: A, C, G, T).
Varför klassisk misslyckas:

Antal kombinationer: 4^60e6 ≈ 10^36e6 (en 1 följd av miljontals nollor) tillstånd. Även med superdatorer tar det år-hundratal att söka sekventiellt.
Minne: Lagring av 10^36e6 kombinationer överstiger jordens totala datorkapacitet (t.ex. 10^21 bytes globalt).
Exempel: Att testa alla mutationer i exomet på en klassisk dator är som att räkna alla sandkorn på en strand – omöjligt!


# QCDS:s lösning:

Kvantparallellism: Med 100 qubits (2^100 ~ 10^30 tillstånd) söker QCDS alla kombinationer samtidigt via superposition (sida 3 i engelska papperet: "superposition enables parallel inference").
Grover-förstärkning: Roterar mot sanning med kvadratrot speedup (sida 2 i brus-papperet), reducerar tid från år till sekunder-minuter per region.
Subset-k och rotate: Dela exomet i 10-20 regioner (3-6 Mb), hantera med 100 qubits per region – som att bläddra i kapitel parallellt.


Bevis: Dina 12-bitars-tester (4096 tillstånd) gav p_true 0.6538 på sekunder – skala till 100 qubits, och exomet blir sökbart där klassisk kraschar.
