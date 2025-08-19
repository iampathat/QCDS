# qcds_full_gpt_22.py
# QCDS – dynamiska lager + IBM Brisbane-läge + "weighted evidence" oracle
# Nytt i v22:
#  - Pilot p0-estimering (--pilot-p0, --pilot-p0-shots) för att uppskatta m* teoretiskt
#  - Dynamisk Grover-toppsökning med fönster runt m* (--pilot-window)
#  - Standard: maxprob + multi-hit + ASCII-graf över score(m)
#  - SANNINGS-GUARDS:
#       * visar p_true vid bekräftelse
#       * agreement (full & per aktiva index) mellan top-sträng och target_bits
#       * multi-hit: antal target-konsistenta (aktiva index) toppar
#  - NYTT: --print-hits för att skriva ut hela multi-hit-listan (state:prob), med ✓ för aktiva-index-konsistens

import os
import json
import math
import random
import argparse
from typing import List, Dict, Tuple, Optional

from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.circuit import AncillaRegister
from qiskit.circuit.library import MCXGate

# AER-simulator (fallback om qiskit-aer inte finns)
try:
    from qiskit_aer import AerSimulator
except Exception:
    from qiskit import Aer
    AerSimulator = None  # type: ignore


# ------------------------------- Utils ---------------------------------

def shannon_H_from_probs(pdict: Dict[str, float]) -> float:
    H = 0.0
    for p in pdict.values():
        if p > 0:
            H -= p * math.log2(p)
    return H


def probs_from_counts(counts: Dict[str, int], shots: int) -> Dict[str, float]:
    shots = max(1, shots)
    return {k: v / shots for k, v in counts.items()}


def ascii_spark(values, width=48) -> str:
    """Liten enradsgraf för score(m)."""
    if not values:
        return ""
    lo = min(values); hi = max(values)
    if hi - lo < 1e-12:
        return "█" * min(len(values), width)
    blocks = "▁▂▃▄▅▆▇█"
    out = []
    idxs = list(range(len(values)))
    if len(values) > width:
        step = len(values) / width
        idxs = [int(i * step) for i in range(width)]
    for i in idxs:
        v = values[i]
        t = (v - lo) / (hi - lo)
        out.append(blocks[min(len(blocks)-1, int(t * (len(blocks)-1) + 1e-9))])
    return "".join(out)


def agree_fraction(s: str, t: str) -> float:
    """Andel bitar som matchar mellan två bitsträngar av samma längd."""
    if not s or not t or len(s) != len(t):
        return 0.0
    return sum(1 for a, b in zip(s, t) if a == b) / max(1, len(t))


# ------------------------ Backend / Runner -----------------------------

class BackendBundle:
    def __init__(self,
                 mode: str,
                 backend=None,
                 backend_name: str = "qasm_simulator",
                 seed_sim: Optional[int] = None,
                 runner_kind: str = "aer",  # "aer" | "ibm_provider"
                 provider=None,
                 service=None):
        self.mode = mode
        self.backend = backend
        self.backend_name = backend_name
        self.seed_sim = seed_sim
        self.runner_kind = runner_kind
        self.provider = provider
        self.service = service


def _make_aer_backend(seed_sim: Optional[int]) -> Tuple[object, str]:
    if AerSimulator is not None:
        backend = AerSimulator()
        if seed_sim is not None:
            try:
                backend.set_options(seed_simulator=seed_sim)
            except Exception:
                pass
        return backend, "qasm_simulator"
    else:
        backend = Aer.get_backend("aer_simulator")  # type: ignore
        return backend, "aer_simulator"


def init_backend(mode: str = "aer",
                 ibm_backend_name: str = "ibm_brisbane",
                 seed_sim: Optional[int] = None) -> BackendBundle:
    mode = (mode or "aer").lower()
    if mode == "ibm":
        try:
            from qiskit_ibm_provider import IBMProvider
            provider = IBMProvider()
            backend = provider.get_backend(ibm_backend_name)
            return BackendBundle(mode="ibm",
                                 backend=backend,
                                 backend_name=ibm_backend_name,
                                 seed_sim=None,
                                 runner_kind="ibm_provider",
                                 provider=provider)
        except Exception as e1:
            print(f"⚠️ IBM-provider kunde inte initieras ({type(e1).__name__}). Faller tillbaka till AER.")
            backend, name = _make_aer_backend(seed_sim)
            return BackendBundle(mode="aer",
                                 backend=backend,
                                 backend_name=name,
                                 seed_sim=seed_sim,
                                 runner_kind="aer",
                                 provider=None)

    backend, name = _make_aer_backend(seed_sim)
    return BackendBundle(mode="aer",
                         backend=backend,
                         backend_name=name,
                         seed_sim=seed_sim,
                         runner_kind="aer")


def run_sampler_once(bundle: BackendBundle, qc: QuantumCircuit, shots: int) -> Dict[str, float]:
    if bundle.runner_kind == "aer":
        tqc = transpile(qc, bundle.backend)
        result = bundle.backend.run(tqc, shots=shots).result()
        counts = result.get_counts(tqc)
        return probs_from_counts(counts, shots)

    elif bundle.runner_kind == "ibm_provider":
        tqc = transpile(qc, bundle.backend)
        job = bundle.backend.run(tqc, shots=shots)
        result = job.result()
        counts = result.get_counts(tqc)
        return probs_from_counts(counts, shots)

    else:
        try:
            tqc = transpile(qc, bundle.backend)
            job = bundle.backend.run(tqc, shots=shots)
            result = job.result()
            counts = result.get_counts(tqc)
            return probs_from_counts(counts, shots)
        except Exception as e:
            print(f"⚠️ Okänt runner-läge ({bundle.runner_kind}): {e}. Faller tillbaka till AER lokalt.")
            aer_backend, _ = _make_aer_backend(bundle.seed_sim)
            tqc = transpile(qc, aer_backend)
            result = aer_backend.run(tqc, shots=shots).result()
            counts = result.get_counts(tqc)
            return probs_from_counts(counts, shots)


# --------------------- MCZ via MCX-tricket -----------------------------

def apply_mcz_on(qubits_list, circ: QuantumCircuit, anc_reg: Optional[AncillaRegister]):
    n = len(qubits_list)
    if n == 0:
        return
    if n == 1:
        circ.z(qubits_list[0]); return
    controls = qubits_list[:-1]
    target = qubits_list[-1]
    need_anc = max(0, len(controls) - 2)
    use_anc = []
    if anc_reg is not None and len(anc_reg) >= need_anc:
        use_anc = [anc_reg[i] for i in range(need_anc)]
    circ.h(target)
    try:
        circ.mcx(controls, target, ancillas=use_anc if use_anc else None)
    except Exception:
        gate = MCXGate(len(controls))
        circ.append(gate, [*controls, target])
    circ.h(target)


# -------------------- Oracle & Diffusion -------------------------------

def apply_weighted_bit_phases(circ: QuantumCircuit,
                              qreg: QuantumRegister,
                              active_indices: List[int],
                              target_bits: str,
                              angle: float):
    for i in active_indices:
        if target_bits[i] == '1':
            circ.p(angle, qreg[i])
        else:
            circ.x(qreg[i]); circ.p(angle, qreg[i]); circ.x(qreg[i])


def apply_oracle_phase_flip_on_state(circ: QuantumCircuit,
                                     qreg: QuantumRegister,
                                     active_indices: List[int],
                                     target_bits: str,
                                     anc_reg: Optional[AncillaRegister],
                                     oracle_type: str = "standard",
                                     weight_angle: float = math.pi,
                                     oracle_power: int = 1):
    active_qubits = [qreg[i] for i in active_indices]
    for _ in range(max(1, oracle_power)):
        if oracle_type == "weighted":
            apply_weighted_bit_phases(circ, qreg, active_indices, target_bits, weight_angle)
        flips = []
        for i in active_indices:
            if target_bits[i] == '0':
                circ.x(qreg[i]); flips.append(i)
        apply_mcz_on(active_qubits, circ, anc_reg)
        for i in flips:
            circ.x(qreg[i])


def apply_diffusion_on_active(circ: QuantumCircuit,
                              qreg: QuantumRegister,
                              active_indices: List[int],
                              anc_reg: Optional[AncillaRegister]):
    active_qubits = [qreg[i] for i in active_indices]
    if len(active_qubits) == 0:
        return
    for q in active_qubits:
        circ.h(q)
    for q in active_qubits:
        circ.x(q)
    apply_mcz_on(active_qubits, circ, anc_reg)
    for q in active_qubits:
        circ.x(q)
    for q in active_qubits:
        circ.h(q)


# ------------------ Grover-lager ---------------------------------------

def build_grover_layered_circuit(nbits: int,
                                 target_bits: str,
                                 active_indices: List[int],
                                 grover_iters: int = 1,
                                 oracle_type: str = "standard",
                                 weight_angle: float = math.pi,
                                 oracle_power: int = 1) -> QuantumCircuit:
    q = QuantumRegister(nbits, "q")
    anc = AncillaRegister(max(0, nbits - 2), "anc")
    c = ClassicalRegister(nbits, "c")
    qc = QuantumCircuit(q, anc, c, name="QCDS")

    for i in range(nbits):
        qc.h(q[i])

    for _ in range(max(1, grover_iters)):
        apply_oracle_phase_flip_on_state(qc, q, active_indices, target_bits, anc,
                                         oracle_type=oracle_type,
                                         weight_angle=weight_angle,
                                         oracle_power=oracle_power)
        apply_diffusion_on_active(qc, q, active_indices, anc)

    qc.barrier()
    qc.measure(q, c)
    return qc


# ------------------ Sched/Heuristik -----------------------------------

def rotating_exclusion(nbits: int, layer_index_1based: int) -> List[int]:
    if nbits <= 1:
        return [0] if nbits == 1 else []
    drop = (layer_index_1based - 1) % nbits
    return [i for i in range(nbits) if i != drop]


def parse_iters_sched(spec: str, L: int, max_m: int) -> List[int]:
    """
    - "auto": 1,1,3,5,7,...
    - "adaptive": hanteras i qcds_layers (dynamisk toppdetektor)
    - annars: kommaseparerad lista
    """
    spec = (spec or "").strip().lower()
    if spec == "adaptive":
        return [1] * L  # placeholder; används ej direkt
    if spec in ("auto", ""):
        seq = []
        m = 1
        step = 2
        while len(seq) < L:
            seq.append(min(max_m, m))
            if len(seq) % 2 == 0 and m < max_m:
                m = min(max_m, m + step)
        return seq
    raw = [s for s in spec.replace(" ", "").split(",") if s]
    vals = []
    for s in raw:
        try:
            vals.append(int(s))
        except:
            vals.append(1)
    if not vals:
        vals = [1]
    while len(vals) < L:
        vals.append(vals[-1])
    return [max(1, min(max_m, v)) for v in vals[:L]]


def optimal_grover_iters(p0: float, max_m: int) -> int:
    """Teoretiskt m* ≈ round(π/(4θ)-1/2), med sin^2 θ = p0."""
    p0 = max(1e-12, min(1.0 - 1e-12, p0))
    theta = math.asin(math.sqrt(p0))
    m_star = int(round((math.pi / (4.0 * theta)) - 0.5))
    return max(1, min(max_m, m_star))


# --------------- Pilot p0-estimering (frivillig) -----------------------

def pilot_estimate_m_star(bundle: BackendBundle,
                          nbits: int,
                          target_bits: str,
                          shots: int,
                          grover_max_iters: int,
                          oracle_type: str,
                          weight_angle: float,
                          oracle_power: int,
                          layman: bool = False) -> Tuple[int, float]:
    """
    Kör m=1 (roterande exkludering) med många shots för att skatta p_true,
    beräknar teoretiskt m*. Används som startfönster för dynamisk sökning.
    """
    active = rotating_exclusion(nbits, 1)
    qc = build_grover_layered_circuit(nbits, target_bits, active, 1,
                                      oracle_type=oracle_type,
                                      weight_angle=weight_angle,
                                      oracle_power=oracle_power)
    probs = run_sampler_once(bundle, qc, shots)
    p0 = probs.get(target_bits, 0.0)
    m_star = optimal_grover_iters(p0, grover_max_iters)
    if layman:
        print(f"   🧪 Pilot p₀≈{p0:.6f}  ⇒  m*≈{m_star} (shots={shots})")
    return m_star, p0


# --------------- Dynamisk Grover: maxprob + multi-hit ------------------

def dynamic_grover_peak_search_maxprob(bundle: BackendBundle,
                                       nbits: int,
                                       target_bits: str,
                                       m_max: int,
                                       oracle_type: str,
                                       weight_angle: float,
                                       oracle_power: int,
                                       shots_scout: int,
                                       shots_confirm: int,
                                       patience: int = 2,
                                       drop_tol: float = 0.05,
                                       multihit_ratio: float = 0.90,
                                       m_start: int = 1,
                                       m_end: Optional[int] = None,
                                       layman: bool = False,
                                       print_hits: bool = False) -> Tuple[int, Dict[str, float], Dict[str, float], list]:
    """
    Ökar m i intervallet [m_start, m_end] (standard 1..m_max) och använder max(prob) som score.
    Stoppar när:
      • score faller > drop_tol relativt föregående steg, eller
      • ingen förbättring på 'patience' steg.
    Bekräftar bästa m med fler shots; returnerar även multi-hit (states ≥ multihit_ratio * max).
    """
    if m_end is None:
        m_end = m_max
    m_start = max(1, min(m_start, m_max))
    m_end = max(m_start, min(m_end, m_max))

    if layman and (m_start != 1 or m_end != m_max):
        print(f"   🔎 Sökfönster: m ∈ [{m_start}, {m_end}] av max {m_max}")

    history = []
    best_score = -1.0
    best_m = m_start
    no_improve = 0
    prev_s = None

    for m in range(m_start, m_end + 1):
        active = rotating_exclusion(nbits, m)
        qc = build_grover_layered_circuit(nbits, target_bits, active, m,
                                          oracle_type=oracle_type,
                                          weight_angle=weight_angle,
                                          oracle_power=oracle_power)
        probs = run_sampler_once(bundle, qc, shots_scout)
        if not probs:
            history.append(0.0)
            continue

        top_str, top_p = max(probs.items(), key=lambda kv: kv[1])
        s = top_p
        history.append(s)

        if layman:
            print(f"   ↻ Scout m={m:>2}  maxprob={s:.4f}  top={top_str}")

        # förbättring?
        if s > best_score + 1e-12:
            best_score = s
            best_m = m
            no_improve = 0
        else:
            no_improve += 1

        # trendbrott? (fall > drop_tol)
        if prev_s is not None and s < (prev_s * (1.0 - drop_tol)):
            if layman:
                print(f"   ⛳ Stoppar: maxprob föll > {int(drop_tol*100)}% (m={m-1}→{m})")
            break

        # tålamod slut?
        if no_improve >= patience:
            if layman:
                print(f"   ⛳ Early stop (patience={patience}) vid m={m}")
            break

        prev_s = s

    # Bekräfta bästa m med fler shots
    active = rotating_exclusion(nbits, best_m)
    qc = build_grover_layered_circuit(nbits, target_bits, active, best_m,
                                      oracle_type=oracle_type,
                                      weight_angle=weight_angle,
                                      oracle_power=oracle_power)
    probs_confirm = run_sampler_once(bundle, qc, shots_confirm)

    # Multi-hit: alla strängar nära max
    hits = {}
    if probs_confirm:
        max_p = max(probs_confirm.values())
        thr = multihit_ratio * max_p
        hits = {k: v for k, v in probs_confirm.items() if v >= thr}

    if layman:
        maxprob_conf = max(probs_confirm.values()) if probs_confirm else 0.0
        p_true_conf = probs_confirm.get(target_bits, 0.0) if probs_confirm else 0.0
        print(f"   ✅ Bekräftat m*={best_m}  maxprob≈{maxprob_conf:.4f}  p_true≈{p_true_conf:.4f}  (confirm_shots={shots_confirm})")
        print(f"   ░ score(m): {ascii_spark(history)}")

        # Agreement för top-strängen
        top_str_conf = max(probs_confirm.items(), key=lambda kv: kv[1])[0] if probs_confirm else None
        agree_full = agree_fraction(top_str_conf, target_bits) if top_str_conf else 0.0
        active_idx = active
        agree_active = (
            sum(1 for i in active_idx if top_str_conf and top_str_conf[i] == target_bits[i]) / max(1, len(active_idx))
        ) if top_str_conf else 0.0

        print(f"   🔎 agreement: full={agree_full:.2%}, active={agree_active:.2%}")

        # Multi-hit: rapportera antal + target-konsistens (aktiva index)
        if hits:
            hits_consistent_active = {
                s: p for s, p in hits.items()
                if all(s[i] == target_bits[i] for i in active_idx)
            }
            top_str_display = top_str_conf if top_str_conf is not None else "-"
            print(f"   ⭐ multi-hit (≥{int(multihit_ratio*100)}% av max) — "
                  f"count={len(hits)}, target-konsistenta(aktiva)={len(hits_consistent_active)}, top={top_str_display}")

            # Full lista om begärt
            if print_hits:
                print("   📜 multi-hit lista (state: prob  [✓=aktivt konsistent]):")
                for s, p in sorted(hits.items(), key=lambda kv: kv[1], reverse=True):
                    ok_active = all(s[i] == target_bits[i] for i in active_idx)
                    mark = " ✓" if ok_active else ""
                    print(f"      {s}: {p:.6f}{mark}")

    return best_m, probs_confirm, hits, history


# ------------------ Dynamisk lager-körning -----------------------------

def qcds_layers(bundle: BackendBundle,
                nbits: int,
                target_bits: str,
                shots: int,
                L: int,
                iters_sched_spec: str,
                grover_max_iters: int,
                layman: bool = False,
                do_pilot_grover: bool = True,
                oracle_type: str = "standard",
                weight_angle: float = math.pi,
                oracle_power: int = 1,
                # adaptiva parametrar:
                adaptive_enabled: bool = False,
                adaptive_scout_shots: int = 256,
                adaptive_confirm_shots: int = 2048,
                adaptive_patience: int = 2,
                adaptive_drop_tol: float = 0.05,
                adaptive_multihit_ratio: float = 0.90,
                # pilot p0 / fönster:
                pilot_p0_enabled: bool = False,
                pilot_p0_shots: int = 4096,
                pilot_window: int = 3,
                # utskrifter:
                print_hits: bool = False) -> Tuple[Optional[str], Optional[float]]:
    """
    Kör L lager. Om adaptive_enabled=True eller iters_sched_spec='adaptive',
    körs dynamisk Grover-toppsökning (maxprob + multi-hit).
    Om pilot_p0_enabled=True: skattar p0 vid m=1 (många shots), beräknar m*,
    och begränsar den dynamiska sökningen till [m*-pilot_window, m*+pilot_window].
    """
    iters_sched = parse_iters_sched(iters_sched_spec, L, grover_max_iters)

    if layman:
        print(f"▶️ QCDS Hierarkisk {L}-lager (backend={bundle.backend_name})")

    last_probs = None
    last_top = None
    p_true_last = None

    for layer_idx in range(1, L + 1):
        active = rotating_exclusion(nbits, layer_idx)
        iters = iters_sched[layer_idx - 1]

        qc = build_grover_layered_circuit(nbits, target_bits, active, iters,
                                          oracle_type=oracle_type,
                                          weight_angle=weight_angle,
                                          oracle_power=oracle_power)
        probs = run_sampler_once(bundle, qc, shots)

        top = max(probs.items(), key=lambda kv: kv[1])[0] if probs else None
        p_true = probs.get(target_bits, 0.0)

        H = shannon_H_from_probs(probs)
        uniform = H / max(1, nbits)

        if layman:
            print(f"▶️ QCDS_L{layer_idx} (iters={iters}) — backend: {bundle.backend_name} — top: {top}  p_true≈{p_true:.4f}")
            print(f"   QCDS H({nbits} qubits): {H:.4f} | uniformitet≈{uniform:.4f}")

        last_probs = probs
        last_top = top
        p_true_last = p_true

    # Dynamisk toppjakt eller klassisk pilot
    if do_pilot_grover and last_probs is not None:
        if adaptive_enabled or (iters_sched_spec.strip().lower() == "adaptive"):
            # Ev. pilot p0 => fönster runt m*
            m_start = 1
            m_end = grover_max_iters
            if pilot_p0_enabled:
                m_star, p0 = pilot_estimate_m_star(
                    bundle, nbits, target_bits,
                    shots=pilot_p0_shots,
                    grover_max_iters=grover_max_iters,
                    oracle_type=oracle_type,
                    weight_angle=weight_angle,
                    oracle_power=oracle_power,
                    layman=layman
                )
                m_start = max(1, m_star - pilot_window)
                m_end = min(grover_max_iters, m_star + pilot_window)
                if layman:
                    print(f"   🎯 Pilotfönster: m ∈ [{m_start}, {m_end}] (window={pilot_window})")

            if layman:
                print("   🔄 Dynamisk Grover-toppsök (maxprob + multi-hit) startar…")
            m_best, probs_best, hits, history = dynamic_grover_peak_search_maxprob(
                bundle, nbits, target_bits,
                m_max=grover_max_iters,
                oracle_type=oracle_type,
                weight_angle=weight_angle,
                oracle_power=oracle_power,
                shots_scout=adaptive_scout_shots,
                shots_confirm=adaptive_confirm_shots,
                patience=adaptive_patience,
                drop_tol=adaptive_drop_tol,
                multihit_ratio=adaptive_multihit_ratio,
                m_start=m_start,
                m_end=m_end,
                layman=layman,
                print_hits=print_hits
            )
            if layman:
                maxprob_best = max(probs_best.values()) if probs_best else 0.0
                print(f"   ↳ Dynamisk Grover: m*={m_best}, maxprob≈{maxprob_best:.4f}, hits={list(hits.keys())}")
        else:
            # fallback: lokal pilot runt teoretisk m* baserat på p_true_last
            p_true_last = p_true_last if p_true_last is not None else 0.0
            m_star = optimal_grover_iters(p_true_last, grover_max_iters)
            m_list = sorted({m for m in [m_star-1, m_star, m_star+1] if 1 <= m <= grover_max_iters})
            best_p = -1.0
            best_m = None
            for m in m_list:
                active_m = rotating_exclusion(nbits, m)
                qc_m = build_grover_layered_circuit(nbits, target_bits, active_m, m,
                                                    oracle_type=oracle_type,
                                                    weight_angle=weight_angle,
                                                    oracle_power=oracle_power)
                probs_m = run_sampler_once(bundle, qc_m, shots)
                p_m = probs_m.get(target_bits, 0.0)
                if p_m > best_p:
                    best_p, best_m = p_m, m
            if layman:
                print(f"   ↳ Pilot Grover: bästa p_true≈{best_p:.4f}  (m={best_m})")

    return last_top, p_true_last


# ------------------ 2-qubit demo --------------------------------------

def two_qubit_demo(bundle: BackendBundle, shots: int, layman: bool, seed: Optional[int]):
    if not layman:
        return
    random.seed(seed or 1337)

    def do_once(label: str, circ_builder):
        q = QuantumRegister(2, "q")
        c = ClassicalRegister(2, "c")
        qc = QuantumCircuit(q, c, name=label)
        circ_builder(qc, q)
        qc.measure(q, c)
        probs = run_sampler_once(bundle, qc, shots)
        items = sorted(probs.items(), key=lambda kv: kv[0])
        top = max(probs.items(), key=lambda kv: kv[1])[0]
        p00 = probs.get("00", 0.0); p01 = probs.get("01", 0.0)
        p10 = probs.get("10", 0.0); p11 = probs.get("11", 0.0)
        exp_ZZ = (p00 + p11) - (p01 + p10)
        bias_q0 = (p00 + p01) - (p10 + p11)
        bias_q1 = (p00 + p10) - (p01 + p11)
        H = shannon_H_from_probs(probs)
        print(f"▶️ {label} @ backend — top: {top}")
        print("   probs: {")
        for k, v in items:
            print(f'  "{k}": {v:.12f},')
        print("}")
        print(f"   QCDS ⟨Z⊗Z⟩: {exp_ZZ:+0.4f} | bias(q0): {bias_q0:+0.4f} | bias(q1): {bias_q1:+0.4f} | H: {H:0.4f}")

    def b1(qc: QuantumCircuit, q):
        qc.h(q[0]); qc.rx(0.02, q[0])
        qc.h(q[1]); qc.rz(0.02, q[1])

    def b2(qc: QuantumCircuit, q):
        qc.h(q[0]); qc.h(q[1]); qc.cx(q[0], q[1])

    def b3(qc: QuantumCircuit, q):
        qc.ry(0.1, q[0]); qc.ry(0.1, q[1]); qc.cx(q[0], q[1]); qc.rx(0.05, q[1])

    print("🔎 2-qubit QCDS-kontroller init…")
    do_once("Truth Gradient", b1)
    do_once("Amplitude Confirmation", b2)
    do_once("Conditional Enhancement", b3)
    print("🔎 2-qubit QCDS-kontroller init…\n")


# ------------------- Syntetisk BRCA2-dataset ---------------------------

def make_brca2_synthetic_variants(nbits: int, max_variants: int = 60) -> List[Tuple[str, str]]:
    start = 32315474
    step = 1390
    out = []
    for i in range(max_variants):
        pos = start + i * step
        s = ["0"] * nbits
        s[i % nbits] = "1"
        target_bits = "".join(s)
        out.append((f"13:{pos:08d}", target_bits))
    return out


# ------------------------------ CLI/Main -------------------------------

def main():
    parser = argparse.ArgumentParser(description="QCDS – dynamiska lager (AER/IBM)")
    parser.add_argument("--mode", type=str, default="aer", help="aer | ibm")
    parser.add_argument("--ibm-backend", type=str, default="ibm_brisbane")
    parser.add_argument("--genome", type=str, default="GRCh38")
    parser.add_argument("--shots", type=int, default=1024)
    parser.add_argument("--grover-max-iters", type=int, default=7, dest="grover_max_iters")
    parser.add_argument("--signal-mode", type=str, default="rotate")
    parser.add_argument("--signal-seed", type=int, default=1337)
    parser.add_argument("--evidence-bits", type=int, default=12, dest="evidence_bits")
    parser.add_argument("--layers", type=int, default=4)
    parser.add_argument("--iters-sched", type=str, default="auto",
                        help='auto | adaptive | eller ex. "1,1,3,5,7"')
    parser.add_argument("--oracle-type", type=str, default="standard", choices=["standard", "weighted"])
    parser.add_argument("--oracle-weight-angle", type=float, default=1.0)
    parser.add_argument("--oracle-power", type=int, default=1)
    parser.add_argument("--layman", action="store_true")
    parser.add_argument("--out-dir", type=str, default="results",
                        help="Katalog att spara utdata i (skapas om den inte finns).")

    # Adaptiv/dynamisk Grover
    parser.add_argument("--adaptive", action="store_true",
                        help="Aktivera dynamisk Grover-toppsökning (maxprob + multi-hit).")
    parser.add_argument("--scout-shots", type=int, default=256,
                        help="Shots per scout-steg.")
    parser.add_argument("--confirm-shots", type=int, default=2048,
                        help="Shots för bekräftelse av bästa m.")
    parser.add_argument("--patience", type=int, default=2,
                        help="Avsluta efter så här många steg utan förbättring.")
    parser.add_argument("--drop-tol", type=float, default=0.05,
                        help="Stoppa om maxprob faller mer än denna andel mellan steg (t.ex. 0.05 = 5%).")
    parser.add_argument("--multihit-threshold", type=float, default=0.90,
                        help="Tröskel för multi-hit (andel av max, t.ex. 0.90).")

    # Pilot p0 (frivillig)
    parser.add_argument("--pilot-p0", action="store_true",
                        help="Kör pilot p0-estimering (m=1 med många shots) och begränsa sökfönstret runt beräknat m*.")
    parser.add_argument("--pilot-p0-shots", type=int, default=4096,
                        help="Shots för pilot p0-körning (m=1).")
    parser.add_argument("--pilot-window", type=int, default=3,
                        help="Sökfönstrets radie runt m* (m ∈ [m*-pilot-window, m*+pilot-window]).")

    # Utskrift av hela multi-hit-listan
    parser.add_argument("--print-hits", action="store_true",
                        help="Skriv ut hela multi-hit-listan (state: prob), markera ✓ om aktivt target-konsistent.")

    args = parser.parse_args()

    # Backend init
    bundle = init_backend(mode=args.mode, ibm_backend_name=args.ibm_backend, seed_sim=args.signal_seed)

    # Header
    iters_sched_preview = args.iters_sched
    print(f"⚙️  QCDS körläge: {bundle.mode}  |  LaymanExplain={bool(args.layman)}")
    print(f"     Profil för BRCA2: {args.genome}  |  default-region 13:32315474-32400266  | max_variants=AER:60")
    print(f"     iters_sched={iters_sched_preview}  |  shots={args.shots}  | adaptive={args.adaptive}  | pilot_p0={args.pilot_p0}  | print_hits={args.print_hits}")
    print(f"     Signal: --signal-mode {args.signal_mode} (seed={args.signal_seed})  |  Grover max iters: {args.grover_max_iters}  |  backend={bundle.backend_name}")

    # Liten 2-qubit-demo
    two_qubit_demo(bundle, args.shots, args.layman, args.signal_seed)

    # Dataset
    nbits = int(args.evidence_bits)
    L = int(args.layers)
    variants = make_brca2_synthetic_variants(nbits, max_variants=60)
    if args.layman:
        print(f"▶️ QCDS Hierarkisk {L}-lager (backend={bundle.backend_name})")
        print("📦 Dataset (syntetiskt): BRCA2 GRCh38 13:32315474-32400266")
        print(f"   Antal ‘varianter’ att köra: {len(variants)} (mode={bundle.mode}, shots={args.shots})")

    results = []
    for i, (pos, targ) in enumerate(variants):
        evs = []  # plats för framtida "evidence list"
        if args.layman:
            print(f"🎯 Variant {pos}  target_bits={targ}  evs={evs}")

        top, ptrue = qcds_layers(bundle,
                                 nbits=nbits,
                                 target_bits=targ,
                                 shots=args.shots,
                                 L=L,
                                 iters_sched_spec=args.iters_sched,
                                 grover_max_iters=args.grover_max_iters,
                                 layman=args.layman,
                                 do_pilot_grover=True,
                                 oracle_type=args.oracle_type,
                                 weight_angle=float(args.oracle_weight_angle),
                                 oracle_power=int(args.oracle_power),
                                 adaptive_enabled=bool(args.adaptive),
                                 adaptive_scout_shots=int(args.scout_shots),
                                 adaptive_confirm_shots=int(args.confirm_shots),
                                 adaptive_patience=int(args.patience),
                                 adaptive_drop_tol=float(args.drop_tol),
                                 adaptive_multihit_ratio=float(args.multihit_threshold),
                                 pilot_p0_enabled=bool(args.pilot_p0),
                                 pilot_p0_shots=int(args.pilot_p0_shots),
                                 pilot_window=int(args.pilot_window),
                                 print_hits=bool(args.print_hits))
        results.append({
            "pos": pos,
            "target_bits": targ,
            "top": top,
            "p_true": ptrue
        })

    # Sammanfattning + spara
    if args.layman:
        print("\n🧭 Lättförståelig sammanfattning (BRCA2 end-to-end):")
        print("• Region: 13:32315474-32400266 (GRCh38)  |  körda punkter: {}".format(len(variants)))
        print("• Dynamisk Grover: maxprob med tidigt stopp; multi-hit; pilot p₀ (valfritt) för att snäva in m.")
        if args.print_hits:
            print("• Utskrift: full multi-hit-lista är aktiverad (--print-hits).")
        print()

    os.makedirs(args.out_dir, exist_ok=True)
    out_path = os.path.join(
        args.out_dir,
        f"QCDS_end2end_{args.genome.lower()}_{bundle.mode}_{bundle.backend_name}.json"
    )
    payload = {
        "genome": args.genome,
        "mode": bundle.mode,
        "backend": bundle.backend_name,
        "shots": args.shots,
        "grover_max_iters": args.grover_max_iters,
        "signal_mode": args.signal_mode,
        "signal_seed": args.signal_seed,
        "evidence_bits": nbits,
        "layers": L,
        "iters_sched": args.iters_sched,
        "oracle_type": args.oracle_type,
        "oracle_weight_angle": float(args.oracle_weight_angle),
        "oracle_power": int(args.oracle_power),
        "adaptive": bool(args.adaptive),
        "scout_shots": int(args.scout_shots),
        "confirm_shots": int(args.confirm_shots),
        "patience": int(args.patience),
        "drop_tol": float(args.drop_tol),
        "multihit_threshold": float(args.multihit_threshold),
        "pilot_p0": bool(args.pilot_p0),
        "pilot_p0_shots": int(args.pilot_p0_shots),
        "pilot_window": int(args.pilot_window),
        "variants": results,
    }
    with open(out_path, "w") as f:
        json.dump(payload, f, indent=2)

    print(f"💾 QCDS end-to-end: sparade {out_path}")


if __name__ == "__main__":
    main()
