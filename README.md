# QCDS
Quantum Condition-Driven Synthesis

QCDS är byggt för att låta många små bevisbitar “rösta” samtidigt och låta våginterferens (kvant-grejen) förstärka det som är konsekvent och dämpa det som är motsägelsefullt. Om/när det fungerar som tänkt blir det här överlägset därför att du får både bättre signal-till-brus och färre falska larm – utan att tappa spårbarhet. Här är varför, på ren svenska:

1) Alla hypoteser testas samtidigt

Klassiskt: du provar en regel i taget (trösklar, if-satser) och riskerar ordnings- och interaktions-effekter.
QCDS: lägger många möjliga evidenskombinationer i superposition och låter systemet “pröva allt på en gång”. Det gör att mönster som verkligen hänger ihop poppar upp, medan slumpmässiga mönster tar ut varandra.

2) “Sanning förstärks, motsägelser släcks”

Vågor som är i fas förstärker varandra (konsistent evidens), ur fas släcker varandra (motsägelsefullt). Vår 4-lagers pipeline (L1→L4) är i praktiken en domino av sanity-checks: går signalen igenom alla steg får du en robust “låst” topp-stat; om inte, försvinner den.

3) Tre oberoende mätningar per par (TG/AC/CE)

Truth Gradient (TG): “tänder/släcker lampor ihop?” (samvariation i Z-bas).

Amplitude Confirmation (AC): kikar på fas → amplitud, ett helt annat perspektiv.

Conditional Enhancement (CE): riktad förstärkning av logiskt villkor.
Tre olika vinklar på samma fråga ger intern korsvalidering redan på 2 qubits.

4) Entropi + bias = inbyggd kvalitetsmätare

Vi loggar entropi (H) och bias för varje körning:

H≈2.0 (2 qubits) eller H≈4.0 (4 qubits) betyder “inget överdrivet mönsterläckage”.

Låg bias ⇒ ingen qubit dominerar.
Det här är din mätare för pålitlighet i realtid, inte bara ett resultat “på tro”.

5) Evidensmotor som är konservativ men transparent

De 4 bitarna (PVS1_LOF, PM2_Rare, PP3_DamagingPred, ClinVar_P/LP) + extra flaggor (BS1, BP4, ClinVar_B/LB) gör att:

Varje ‘1’ har tydlig innebörd som går att läsa i JSON och lekmannatexten.

Inga hemliga trösklar; du ser exakt varför en poäng tippade åt ett håll.

Konservativ tolkning ⇒ hellre VUS än felaktigt “patogen”.

6) Färre falska positiva genom konsekvenskrav

För att ett 4-bitars mönster ska vinna måste det överleva fyra lager av villkor och interferens. Brus och spretiga prediktorer tenderar att självdämpa; bara robusta kombinationer förstärks.

7) Skalar från laptop → hårdvara (utan kodändring)

Du kan mass-köra AER över hela BRCA2 (tätt rutnät eller alla VCF-varianter) för fart och debugging, och sedan slå på IBM för verklighetskontroll. Samma fil, samma API, med flaggan som växlar backend. Det gör att du kan:

Emulera snabbt, iterera regler, fintrimma.

Verifiera på riktig kvantdator när du vill.

8) Robust mot “tröskelkampsport”

I klassiska pipelines bråkar ofta AF-gräns, prediktorskalning och ClinVar-taggar med varandra. QCDS låter helheten avgöra via våglogik, i stället för att du måste “mikro-tuna” varje if-sats. Resultatet blir stabilare mellan dataset.

9) Spårbarhet för både expert och lekman

Allt sparas i två JSON-rapporter + lättförståelig text:

Onkolog får siffror (⟨Z⊗Z⟩, H, bias, topp-tillstånd).

Lekman får vardagsmetaforer (“två lampor som tänds ihop”, “vi nyper systemet och ser hur det svarar”).
Det gör resultaten kommunicerbara utan att tappa djupet.

10) Utbyggbart utan att riva allt

Vill du lägga till fler regler (t.ex. spliceprediktorer, domain-hits, co-segregation)? Lägg till fler bitar/lagervillkor. QCDS-mönstret (superposition → interferens → förstärkning) återanvänds. Du skalar komplexitet horisontellt, inte med en växande djungel av ad-hoc-if.

Vad betyder “vi hittar inget evidens nu?”

Det betyder inte att QCDS “misslyckas” – det betyder att den konservativa evidensmotorn inte fann en robust, konsekvent kombination (på de regler vi valt) som överlever 4 lager. Det är exakt så vi vill att det beter sig när materialet inte stödjer ett starkt påstående: hellre “inget” än ett falskt “något”.

Vill du öka träffchansen (utan att offra konservatismen) kan vi:

Lägga till fler tillåtna evidenstyper (fler bitar) – t.ex. spliceprediktorer, domänspecifik kunskap (BRCA2 BRC-repeats), som extra “röster”.

Göra region-scan tätare i AER, men samma setup ska också kunna köras mot IBM (färre punkter, men samma logik).

Justera PM2/BS1-trösklar och vilka INFO-fält vi accepterar (VEP/ANN/CADD/SIFT/PolyPhen/gnomAD/ClinVar).

“QCDS kombinerar flera oberoende kvantprover och en konservativ evidensmotor; endast mönster som är konsekventa över flera lager förstärks. Det ger hög motståndskraft mot brus och lägre risk för falska positiva jämfört med sekventiella tröskelpipelines, samtidigt som resultat och osäkerhet kan förklaras begripligt för både experter och lekmän.”
