# QCDS

Quantum Condition-Driven Synthesis

QCDS is built to let many small pieces of evidence “vote” at the same time and let wave interference (the quantum thing) amplify what is consistent and dampen what is contradictory. If/when it works as intended, this will be superior because you get both better signal-to-noise and fewer false alarms – without losing traceability. Here’s why, in plain English:

All hypotheses are tested simultaneously

Classical: you try one rule at a time (thresholds, if-statements) and risk order and interaction effects. QCDS: superposes many possible evidence combinations and lets the system “try everything at once”. This makes patterns that are truly connected pop up, while random patterns cancel each other out.

“Truth is amplified, contradictions are extinguished”

Waves that are in phase amplify each other (consistent evidence), out of phase extinguish each other (contradictory). Our 4-layer pipeline (L1→L4) is in effect a domino of sanity checks: if the signal passes all the steps, you get a robust “locked” peak state; if not, it disappears.

Three independent measurements per pair (TG/AC/CE)

Truth Gradient (TG): “do lights turn on/off together?” (covariation in Z-basis).

Amplitude Confirmation (AC): looks at phase → amplitude, a completely different perspective.

Conditional Enhancement (CE): targeted reinforcement of logical conditions. Three different angles on the same question provide internal cross-validation already on 2 qubits.

Entropy + bias = built-in quality meter

We log entropy (H) and bias for each run:

H≈2.0 (2 qubits) or H≈4.0 (4 qubits) means “no excessive pattern leakage”.

Low bias ⇒ no qubit dominates. This is your real-time reliability gauge, not just a “faith” result.

Evidence engine that is conservative but transparent

The 4 bits (PVS1_LOF, PM2_Rare, PP3_DamagingPred, ClinVar_P/LP) + extra flags (BS1, BP4, ClinVar_B/LB) mean that:

Every ‘1’ has a clear meaning that can be read in JSON and layman’s text.

No secret thresholds; you see exactly why a score tipped in one direction.

Conservative interpretation ⇒ rather VUS than falsely “pathogenic”.

Fewer false positives through consistency requirements

For a 4-bit pattern to win, it must survive four layers of conditions and interference. Noise and sprawling predictors tend to self-attenuate; only robust combinations are amplified.

Scales from laptop → hardware (without code change)

You can mass-run AER over the entire BRCA2 (dense grid or all VCF variants) for speed and debugging, and then turn on IBM for reality check. Same file, same API, with the flag that switches the backend. This allows you to:

Emulate quickly, iterate rules, fine-tune.

Verify on a real quantum computer whenever you want.

Robust against “threshold martial arts”

In classic pipelines, AF limit, predictor scaling and ClinVar tags often fight each other. QCDS lets the whole decide via wave logic, instead of you having to “micro-tune” each if statement. The result is more stable between datasets.

Traceability for both expert and layperson

Everything is saved in two JSON reports + easy-to-understand text:

Oncologist gets numbers (⟨Z⊗Z⟩, H, bias, peak state).

Layman gets everyday metaphors (“two lights that turn on together”, “we pinch the system and see how it responds”). It makes the results communicable without losing depth.

Extensible without tearing everything down

Want to add more rules (e.g. splice predictors, domain hits, co-segregation)? Add more bits/layer conditions. The QCDS pattern (superposition → interference → amplification) is reused. You scale complexity horizontally, not with a growing jungle of ad-hoc-ifs.

What does “we find no evidence now” mean?

It doesn’t mean that QCDS “fails” – it means that the conservative evidence engine didn’t find a robust, consistent combination (of the rules we chose) that survives 4 layers. That’s exactly how we want it to behave when the material doesn’t support a strong claim: “nothing” is better than a false “something”.

If you want to increase the hit rate (without sacrificing conservatism), we can:

Add more allowed evidence types (more bits) – e.g. splice predictors, domain-specific knowledge (BRCA2 BRC repeats), as extra “votes”.

Make region-scan denser in AER, but the same setup should also be able to run against IBM (fewer points, but same logic).

Adjust PM2/BS1 thresholds and which INFO fields we accept (VEP/ANN/CADD/SIFT/PolyPhen/gnomAD/ClinVar).

“QCDS combines multiple independent quantum samples and a conservative evidence engine; only patterns that are consistent across multiple layers are amplified. This provides high resilience to noise and lower risk of false positives compared to sequential threshold pipelines, while results and uncertainty can be explained in a way that is understandable to both experts and laypeople.”



#############################################

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
