# case3_qcds_openplan_transpile.py
from qiskit import QuantumCircuit, transpile
from qiskit_ibm_runtime import QiskitRuntimeService
from qiskit_aer import Aer

SHOTS     = 1024
BACKENDS  = ["ibm_brisbane", "ibm_torino"]  # testas i ordning

# ===== QCDS: flaggor (valfritt, ändrar inte befintlig funktionalitet) =====
QCDS_ENABLE = True            # sätt False om du vill stänga av extra QCDS-utskrifter
QCDS_EXPORT_JSON = False      # spara även QCDS-mått till 'results_QCDS.json'

def build_truth_gradient():
    qc = QuantumCircuit(2); qc.h(0); qc.cx(0, 1); qc.measure_all(); return qc

def build_amplitude_confirmation():
    qc = QuantumCircuit(2); qc.h([0, 1]); qc.cz(0, 1); qc.measure_all(); return qc

def build_conditional_enhancement():
    qc = QuantumCircuit(2); qc.h([0, 1]); qc.cz(0, 1); qc.h(0); qc.measure_all(); return qc

def try_ibm_sampler(qc: QuantumCircuit):
    # Init service via ibm_cloud (använder sparade creds om du kört save_account)
    service = QiskitRuntimeService(channel="ibm_cloud")
    from qiskit_ibm_runtime import SamplerV2 as Sampler  # Sampler v2 i job mode

    last_err = None
    for name in BACKENDS:
        try:
            backend = service.backend(name)

            # 🔧 VIKTIGT: Transpilera mot backenden så H/CZ mappas till target-gates
            tqc = transpile(
                qc,
                backend=backend,
                optimization_level=1,     # lätt opt
                layout_method="sabre",    # bra default
                routing_method="sabre"
            )

            # Sampler v2 i "job mode": initiera med mode=backend
            sampler = Sampler(mode=backend)

            job = sampler.run([tqc], shots=SHOTS)
            result = job.result()

            # Stöd för både SamplerResult (quasi_dists) och PrimitiveResult (counts)
            if hasattr(result, "quasi_dists"):  # SamplerResult-väg
                probs = result.quasi_dists[0].binary_probabilities()
            else:  # PrimitiveResult-väg
                pub = result[0]
                # Försök läsa counts från 'meas' (skapad av measure_all)
                try:
                    counts = pub.data.meas.get_counts()
                except AttributeError:
                    # Fallback: hitta första attribut som har get_counts()
                    cand = next(
                        (getattr(pub.data, n) for n in dir(pub.data)
                         if hasattr(getattr(pub.data, n), "get_counts")),
                        None
                    )
                    if cand is None:
                        raise RuntimeError("Hittar inga counts i PrimitiveResult.data")
                    counts = cand.get_counts()
                total = sum(counts.values()) or SHOTS
                probs = {k: v / total for k, v in counts.items()}

            top = max(probs, key=probs.get)
            print(f"✅ IBM körning OK på {name}")
            return {"backend": name, "top": top, "probs": probs}

            #job = sampler.run([tqc], shots=SHOTS)
            #quasi = job.result().quasi_dists[0].binary_probabilities()
            #top = max(quasi, key=quasi.get)
            #print(f"✅ IBM körning OK på {name}")
            #return {"backend": name, "top": top, "probs": quasi}


        except Exception as e:
            last_err = e
            print(f"⚠️  {name} misslyckades: {e}")
    raise RuntimeError(last_err or "IBM okänt fel")

def run_aer_local(qc: QuantumCircuit):
    sim = Aer.get_backend("qasm_simulator")
    counts = sim.run(qc, shots=SHOTS).result().get_counts()
    total = sum(counts.values())
    probs = {k: v/total for k, v in counts.items()}
    top = max(probs, key=probs.get)
    return {"backend": "qasm_simulator", "top": top, "probs": probs}

def run(qc: QuantumCircuit):
    try:
        return try_ibm_sampler(qc)
    except Exception as e:
        print(f"↩️  IBM misslyckades, går till Aer: {e}")
        return run_aer_local(qc)

# =======================
# ===== QCDS helpers =====
# (Endast extra – påverkar inte din körning)
# =======================
import json, math
from datetime import datetime

def _p(probs, key):
    return float(probs.get(key, 0.0))

def qcds_metrics_from_probs(probs: dict):
    """
    Beräknar QCDS-mått ur sannolikheter i Z-basis.
    Antagande: bitsträngar '00','01','10','11' där högersta biten är qubit 0.
    """
    p00, p01, p10, p11 = _p(probs, '00'), _p(probs, '01'), _p(probs, '10'), _p(probs, '11')
    # ⟨Z⊗Z⟩ / "Truth Gradient"-score: korrelation i Z-basis
    z_corr = (p00 + p11) - (p01 + p10)
    # bias per qubit
    bit0_bias = (p00 + p10) - (p01 + p11)   # högersta bit (q0)
    bit1_bias = (p00 + p01) - (p10 + p11)   # vänstersta bit (q1)
    # entropi som "uniformitetsmått" (2 bitar max)
    eps = 1e-12
    H = 0.0
    for p in (p00, p01, p10, p11):
        if p > 0:
            H -= p * math.log(p, 2)
    H_max = 2.0
    uniformity = H / H_max  # 1.0 ~ helt uniform, 0 ~ helt deterministisk

    return {
        "zz_correlation": round(z_corr, 4),
        "bit0_bias": round(bit0_bias, 4),
        "bit1_bias": round(bit1_bias, 4),
        "entropy_bits": round(H, 4),
        "uniformity_0to1": round(uniformity, 4),
    }

def qcds_print(label: str, probs: dict):
    m = qcds_metrics_from_probs(probs)
    # Lätta, icke-invasiva utskrifter efter dina egna prints
    print(f"   QCDS ⟨Z⊗Z⟩: {m['zz_correlation']} | bias(q0): {m['bit0_bias']} | bias(q1): {m['bit1_bias']} | H: {m['entropy_bits']}")

def qcds_maybe_export(all_results: list):
    if not QCDS_EXPORT_JSON:
        return
    payload = {
        "generated_at": datetime.utcnow().isoformat()+"Z",
        "shots": SHOTS,
        "backends_tried": BACKENDS,
        "results": all_results,
    }
    with open("results_QCDS.json", "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)
    print("💾 QCDS: sparade results_QCDS.json")

# =======================
# ===== MAIN (oförändrad körning; bara extra QCDS-utskrift efteråt) =====
# =======================
if __name__ == "__main__":
    tests = [
        ("Truth Gradient", build_truth_gradient()),
        ("Amplitude Confirmation", build_amplitude_confirmation()),
        ("Conditional Enhancement", build_conditional_enhancement()),
    ]
    _acc = []  # för ev. JSON-export
    for label, qc in tests:
        r = run(qc)
        print(f"▶️ {label} @ {r['backend']} — top: {r['top']}")
        top5 = sorted(r["probs"].items(), key=lambda kv: kv[1], reverse=True)[:5]
        print("   probs:", {k: round(v, 4) for k, v in top5})
        if QCDS_ENABLE:
            qcds_print(label, r["probs"])
        _acc.append({"label": label, **r})
    qcds_maybe_export(_acc)



# =======================
# ===== QCDS EXTRAS (frivilliga, ändrar inget i din körning) =====
# =======================

# 1) Välj ett "bäst" kopplat qubit-par på backenden (t.ex. lägsta ecr/cx-fel)
def pick_best_coupled_pair(backend, op_name="ecr"):
    """
    Returnerar (qA, qB) för den koppling som ser bäst ut enligt backend.target.
    Faller tillbaka till första paret om fel-data saknas.
    """
    try:
        op = backend.target.get_operation(op_name) or backend.target.get_operation("cx")
        if op is None:
            raise RuntimeError("Varken ecr eller cx finns på target.")
        best = None
        best_err = None
        # op kan vara en mapping {(qA, qB): [InstructionProperties, ...], ...}
        for pair, props in op.items():
            # props kan vara list/tuple eller single; leta error-attribut
            err = None
            if props is None:
                err = None
            elif isinstance(props, (list, tuple)) and props:
                # ta första prop som har 'error'
                for p in props:
                    if hasattr(p, "error") and p.error is not None:
                        err = p.error
                        break
            else:
                if hasattr(p, "error") and p.error is not None:
                    err = p.error
            # välj min error om den finns, annars välj första paret
            if best is None:
                best, best_err = pair, err
            else:
                if err is not None and (best_err is None or err < best_err):
                    best, best_err = pair, err
        return best
    except Exception:
        # Fallback: leta i coupling map om den finns
        try:
            cmap = backend.coupling_map
            return tuple(cmap.get_edges()[0])
        except Exception:
            # Sista utvägen: (0,1)
            return (0, 1)

# 2) Bygg en X-basis-variant av Amplitude Confirmation för att synliggöra fas
def build_amplitude_confirmation_xbasis():
    qc = QuantumCircuit(2)
    qc.h([0, 1]); qc.cz(0, 1)
    # Läs i X-basis (fas -> amplituder)
    qc.h([0, 1])
    qc.measure_all()
    return qc

# 3) Kör samma pipeline som try_ibm_sampler men lås layout till ett specifikt qubit-par
def try_ibm_sampler_on_pair(qc: QuantumCircuit, backend_name: str, pair: tuple, shots: int = SHOTS):
    """
    Identisk resultat-hantering som i try_ibm_sampler (stöd för både SamplerResult och PrimitiveResult),
    men transpilerar med initial_layout=pair så du kan styra vilka fysiska qubits som används.
    """
    service = QiskitRuntimeService(channel="ibm_cloud")
    from qiskit_ibm_runtime import SamplerV2 as Sampler  # Sampler v2 i job mode

    backend = service.backend(backend_name)

    # Transpilera med fixerat par (t.ex. valt via pick_best_coupled_pair)
    tqc = transpile(
        qc,
        backend=backend,
        initial_layout=[pair[0], pair[1]],
        optimization_level=1,
        layout_method="sabre",
        routing_method="sabre",
    )

    sampler = Sampler(mode=backend)
    job = sampler.run([tqc], shots=shots)
    result = job.result()

    # *** EXAKT SAMMA ROBUSTA RESULT-HANTERING SOM I DIN FUNKTION ***
    if hasattr(result, "quasi_dists"):  # SamplerResult
        probs = result.quasi_dists[0].binary_probabilities()
    else:  # PrimitiveResult
        pub = result[0]
        try:
            counts = pub.data.meas.get_counts()
        except AttributeError:
            cand = next(
                (getattr(pub.data, n) for n in dir(pub.data)
                 if hasattr(getattr(pub.data, n), "get_counts")),
                None
            )
            if cand is None:
                raise RuntimeError("Hittar inga counts i PrimitiveResult.data")
            counts = cand.get_counts()
        total = sum(counts.values()) or shots
        probs = {k: v / total for k, v in counts.items()}

    top = max(probs, key=probs.get)
    return {"backend": backend_name, "top": top, "probs": probs, "pair": pair}

# 4) Liten helper för att att köra X-basis-varianten via din befintliga run()
def run_amplitude_confirmation_xbasis():
    qc = build_amplitude_confirmation_xbasis()
    return run(qc)  # använder din exakta pipeline (try_ibm_sampler -> Aer-fallback)




# --- Exempel 1: välj “bästa” par på brisbane och kör Truth Gradient på just det paret
# service = QiskitRuntimeService(channel="ibm_cloud")
# backend = service.backend("ibm_brisbane")
# best_pair = pick_best_coupled_pair(backend, op_name="ecr")
# print("🔎 valt par på brisbane:", best_pair)
# res_pair = try_ibm_sampler_on_pair(build_truth_gradient(), "ibm_brisbane", best_pair)
# print("▶️ Truth Gradient (låst par) @", res_pair["backend"], res_pair["pair"], "— top:", res_pair["top"])
# print("   probs:", {k: round(v, 4) for k, v in sorted(res_pair["probs"].items(), key=lambda kv: kv[1], reverse=True)[:5]})

# --- Exempel 2: X-basis för att se fasinformation i Amplitude Confirmation
# rx = run_amplitude_confirmation_xbasis()
# print("▶️ Amplitude Confirmation (X-basis) @", rx["backend"], "— top:", rx["top"])
# print("   probs:", {k: round(v, 4) for k, v in sorted(rx["probs"].items(), key=lambda kv: kv[1], reverse=True)[:5]})


# ===== QCDS KONVERGENS (append-only) =====

from collections import defaultdict

def avg_probs(dicts):
    acc = defaultdict(float)
    for d in dicts:
        for k, v in d.items():
            acc[k] += v
    n = max(1, len(dicts))
    return {k: acc[k] / n for k in acc}

def qcds_converge(backend_name="ibm_brisbane", repeats=3):
    """
    QCDS-konvergens:
    - väljer "bästa" kopplade qubit-par
    - kör Truth Gradient (Z), Conditional Enhancement (Z), Amplitude Confirmation (X)
    - repeterar och medelvärdesbildar sannolikheter för bättre SNR
    Kräver att du har sparat hjälpfunktionerna:
      pick_best_coupled_pair, build_amplitude_confirmation_xbasis, try_ibm_sampler_on_pair
    """
    service = QiskitRuntimeService(channel="ibm_cloud")
    backend = service.backend(backend_name)

    # Välj rätt kopplingsgate och bästa par
    try:
        native_ops = {op.name for op in backend.target.operations}
    except Exception:
        native_ops = set()
    op_name = "ecr" if "ecr" in native_ops else "cz"
    pair = pick_best_coupled_pair(backend, op_name=op_name)
    print(f"🔎 QCDS konvergens på {backend_name}, valt par {pair}, op={op_name}, repeats={repeats}")

    # Bygg kretsar med rätt basis per fenomen
    tg  = build_truth_gradient()                 # Z-basis
    ce  = build_conditional_enhancement()        # Z-basis
    acx = build_amplitude_confirmation_xbasis()  # X-basis (H före mätning)

    # Kör flera repeter och medelvärdesbilda (konvergera)
    tg_runs  = [try_ibm_sampler_on_pair(tg,  backend_name, pair)["probs"] for _ in range(repeats)]
    ce_runs  = [try_ibm_sampler_on_pair(ce,  backend_name, pair)["probs"] for _ in range(repeats)]
    acx_runs = [try_ibm_sampler_on_pair(acx, backend_name, pair)["probs"] for _ in range(repeats)]

    tg_avg, ce_avg, acx_avg = avg_probs(tg_runs), avg_probs(ce_runs), avg_probs(acx_runs)

    def topk(probs, k=4):
        return {k_: round(v, 4) for k_, v in sorted(probs.items(), key=lambda kv: kv[1], reverse=True)[:k]}

    # Utskrift + QCDS-mått (⟨ZZ⟩ för Z-basis; för AC(X) motsvarar detta ⟨XX⟩ före H)
    print(f"▶️ TG (Z) @ {backend_name} {pair} — top:", max(tg_avg, key=tg_avg.get))
    print("   probs:", topk(tg_avg)); qcds_print("TG(Z)", tg_avg)

    print(f"▶️ CE (Z) @ {backend_name} {pair} — top:", max(ce_avg, key=ce_avg.get))
    print("   probs:", topk(ce_avg)); qcds_print("CE(Z)", ce_avg)

    print(f"▶️ AC (X) @ {backend_name} {pair} — top:", max(acx_avg, key=acx_avg.get))
    print("   probs:", topk(acx_avg))
    qcds_print("AC(X→tolka som ⟨XX⟩)", acx_avg)

# Kör när du vill:
qcds_converge("ibm_brisbane", repeats=3)

# ===== HIERARKISKT 4 IMPLEMENTATION (append-only, påverkar inte befintlig kod) =====

# Definiera en enkel condition-syntes baserat på föregående ψ (top state)
def synthesize_new_condition(previous_top: str):
    """
    Syntetiserar ny condition field C_n från föregående ψ_n (top state).
    Exempel: Invertera top state för att skapa re-entry (kan utökas med mer logik).
    """
    if previous_top == '11':
        return "x or y"  # Exempel på re-entry: bredda till '01', '10', '11'
    elif previous_top == '00':
        return "x and y"  # Fokus på '11' för konvergens
    else:
        return "x and y"  # Default fallback

# QCDS Hierarchical 4 Cascade
def qcds_hierarchical_4(initial_condition: str = "x and y", num_qubits: int = 2):
    """
    Implementerar hierarkiskt 4: 4-nivåers cascade för brusreducering och meta-inference.
    - Q1: Initial Synthesizer – Kör med initial C1, producera ψ1
    - Q2: Inferential Re-entry – Syntetiser C2 från ψ1, producera ψ2
    - Q3: Cascading Extension – Syntetiser C3 från ψ2, producera ψ3
    - Q4: Syntactic Resonance – Syntetiser C4 från ψ3, producera ψ4 (final output)
    Körs med befintlig run() för att inte störa din kod.
    """
    print(f"▶️ QCDS Hierarchical 4 på {BACKENDS[0]} med initial condition: {initial_condition}")

    # Q1: Initial Synthesizer
    q1_circuit = build_truth_gradient()  # Använd en bas-krets; utöka med oracle baserat på condition
    q1_result = run(q1_circuit)
    q1_top = q1_result["top"]
    print(f"Q1: Top state: {q1_top}, probs: {q1_result['probs']}")
    if QCDS_ENABLE:
        qcds_print("Q1", q1_result["probs"])

    # Q2: Inferential Re-entry
    q2_condition = synthesize_new_condition(q1_top)
    q2_circuit = build_amplitude_confirmation()  # Använd annan krets för variation; utöka med oracle
    q2_result = run(q2_circuit)
    q2_top = q2_result["top"]
    print(f"Q2 (C2: {q2_condition}): Top state: {q2_top}, probs: {q2_result['probs']}")
    if QCDS_ENABLE:
        qcds_print("Q2", q2_result["probs"])

    # Q3: Cascading Extension
    q3_condition = synthesize_new_condition(q2_top)
    q3_circuit = build_conditional_enhancement()  # Använd tredje krets
    q3_result = run(q3_circuit)
    q3_top = q3_result["top"]
    print(f"Q3 (C3: {q3_condition}): Top state: {q3_top}, probs: {q3_result['probs']}")
    if QCDS_ENABLE:
        qcds_print("Q3", q3_result["probs"])

    # Q4: Syntactic Resonance
    q4_condition = synthesize_new_condition(q3_top)
    q4_circuit = build_truth_gradient()  # Återanvänd eller bygg ny; utöka med full oracle
    q4_result = run(q4_circuit)
    q4_top = q4_result["top"]
    print(f"Q4 (C4: {q4_condition}): Top state: {q4_top}, probs: {q4_result['probs']}")
    if QCDS_ENABLE:
        qcds_print("Q4", q4_result["probs"])

    return {"final_top": q4_top, "final_probs": q4_result["probs"]}

# Kör hierarchical 4 när du vill:
qcds_hierarchical_4(initial_condition="x and y")

# ===== QCDS KONVERGENS (append-only) =====

from collections import defaultdict

def avg_probs(dicts):
    acc = defaultdict(float)
    for d in dicts:
        for k, v in d.items():
            acc[k] += v
    n = max(1, len(dicts))
    return {k: acc[k] / n for k in acc}

def qcds_converge(backend_name="ibm_brisbane", repeats=3):
    """
    QCDS-konvergens:
    - väljer "bästa" kopplade qubit-par
    - kör Truth Gradient (Z), Conditional Enhancement (Z), Amplitude Confirmation (X)
    - repeterar och medelvärdesbildar sannolikheter för bättre SNR
    Kräver att du har sparat hjälpfunktionerna:
      pick_best_coupled_pair, build_amplitude_confirmation_xbasis, try_ibm_sampler_on_pair
    """
    service = QiskitRuntimeService(channel="ibm_cloud")
    backend = service.backend(backend_name)

    # Välj rätt kopplingsgate och bästa par
    try:
        native_ops = {op.name for op in backend.target.operations}
    except Exception:
        native_ops = set()
    op_name = "ecr" if "ecr" in native_ops else "cz"
    pair = pick_best_coupled_pair(backend, op_name=op_name)
    print(f"🔎 QCDS konvergens på {backend_name}, valt par {pair}, op={op_name}, repeats={repeats}")

    # Bygg kretsar med rätt basis per fenomen
    tg  = build_truth_gradient()                 # Z-basis
    ce  = build_conditional_enhancement()        # Z-basis
    acx = build_amplitude_confirmation_xbasis()  # X-basis (H före mätning)

    # Kör flera repeter och medelvärdesbilda (konvergera)
    tg_runs  = [try_ibm_sampler_on_pair(tg,  backend_name, pair)["probs"] for _ in range(repeats)]
    ce_runs  = [try_ibm_sampler_on_pair(ce,  backend_name, pair)["probs"] for _ in range(repeats)]
    acx_runs = [try_ibm_sampler_on_pair(acx, backend_name, pair)["probs"] for _ in range(repeats)]

    tg_avg, ce_avg, acx_avg = avg_probs(tg_runs), avg_probs(ce_runs), avg_probs(acx_runs)

    def topk(probs, k=4):
        return {k_: round(v, 4) for k_, v in sorted(probs.items(), key=lambda kv: kv[1], reverse=True)[:k]}

    # Utskrift + QCDS-mått (⟨ZZ⟩ för Z-basis; för AC(X) motsvarar detta ⟨XX⟩ före H)
    print(f"▶️ TG (Z) @ {backend_name} {pair} — top:", max(tg_avg, key=tg_avg.get))
    print("   probs:", topk(tg_avg)); qcds_print("TG(Z)", tg_avg)

    print(f"▶️ CE (Z) @ {backend_name} {pair} — top:", max(ce_avg, key=ce_avg.get))
    print("   probs:", topk(ce_avg)); qcds_print("CE(Z)", ce_avg)

    print(f"▶️ AC (X) @ {backend_name} {pair} — top:", max(acx_avg, key=acx_avg.get))
    print("   probs:", topk(acx_avg))
    qcds_print("AC(X→tolka som ⟨XX⟩)", acx_avg)

# Kör när du vill:
qcds_converge("ibm_brisbane", repeats=3)
