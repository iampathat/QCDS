#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
qcds_mutations_allinone.py

EN FIL SOM GÖR ALLT:
- Kör QCDS (AER eller IBM) för en eller flera varianter (batch).
- Har lekmanna + teknisk förklaring per variant.
- Samlar resultat (p_true=SANNING, m*, agreement) och skriver Markdown-tabell.
- Ger en icke-klinisk "trend/hypotes om sjukdomsfenotyp/diagnos" baserat på resultaten.
- Sparar per-variant JSON + sammanfattning till results/.

⚠️ VIKTIGT: Detta är forskning/experiment. Utskrifterna om "trend/diagnos"
är en icke-klinisk modellhypotes och ersätter inte medicinsk bedömning.

SNABBGUIDE
----------
Batch (AER, dina parametrar):
    python qcds_mutations_allinone.py --mode aer --batch \
      --layers 6 --shots 2048 \
      --iters-sched adaptive --adaptive \
      --pilot-p0 --pilot-window 10 --grover-max-iters 60 \
      --oracle-type standard --oracle-power 1 \
      --score-scout p-true --score-confirm p-true \
      --active-mode full \
      --confirm-shots 8192 --drop-tol 0.10 --patience 3 \
      --multihit-threshold 0.90 \
      --progress on --layman

Batch (IBM Brisbane):
    python qcds_mutations_allinone.py --mode ibm --ibm-backend ibm_brisbane --batch \
      --layers 6 --shots 2048 \
      --iters-sched adaptive --adaptive \
      --pilot-p0 --pilot-window 10 --grover-max-iters 60 \
      --oracle-type standard --oracle-power 1 \
      --score-scout p-true --score-confirm p-true \
      --active-mode full \
      --confirm-shots 8192 --drop-tol 0.10 --patience 3 \
      --multihit-threshold 0.90 \
      --progress on --layman

Enskild variant (t.ex. KRAS_G12D index 3):
    python qcds_mutations_allinone.py --mode aer \
      --single KRAS_G12D --target-index 3 \
      --layers 6 --shots 2048 --iters-sched adaptive --adaptive \
      --pilot-p0 --pilot-window 10 --grover-max-iters 60 \
      --oracle-type standard --oracle-power 1 \
      --score-scout p-true --score-confirm p-true \
      --active-mode full \
      --confirm-shots 8192 --drop-tol 0.10 --patience 3 \
      --multihit-threshold 0.90 \
      --progress on --layman
"""

import os
import sys
import json
import math
import random
import argparse
from typing import List, Dict, Tuple, Optional

# ------------------------- Qiskit Imports -------------------------
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.circuit import AncillaRegister
from qiskit.circuit.library import MCXGate

# AER-simulator (fallback om qiskit-aer inte finns)
try:
    from qiskit_aer import AerSimulator
except Exception:
    from qiskit import Aer
    AerSimulator = None  # type: ignore

# =============================== Utils ===============================

def shannon_H_from_probs(pdict: Dict[str, float]) -> float:
    """Shannon-entropi i bitar (log2)."""
    H = 0.0
    for p in pdict.values():
        if p > 0:
            H -= p * math.log2(p)
    return H

def probs_from_counts(counts: Dict[str, int], shots: int) -> Dict[str, float]:
    shots = max(1, shots)
    return {k: v / shots for k, v in counts.items()}

def ascii_spark(values, width=48) -> str:
    """Liten enradsgraf (sparklines) för t.ex. score(m) historik."""
    if not values:
        return ""
    lo = min(values); hi = max(values)
    if hi - lo < 1e-12:
        return "█" * min(len(values), width)
    blocks = "▁▂▃▄▅▆▇█"
    out, idxs = [], list(range(len(values)))
    if len(values) > width:
        step = len(values) / width
        idxs = [int(i * step) for i in range(width)]
    for i in idxs:
        v = values[i]
        t = (v - lo) / (hi - lo)
        out.append(blocks[min(len(blocks)-1, int(t*(len(blocks)-1)+1e-9))])
    return "".join(out)

def agree_fraction(s: str, t: str) -> float:
    """Andel bitar som matchar mellan två bitsträngar av samma längd."""
    if not s or not t or len(s) != len(t):
        return 0.0
    return sum(1 for a, b in zip(s, t) if a == b) / max(1, len(t))

def is_consistent_on_active(bitstr: str, target_bits: str, active_idx: List[int]) -> bool:
    return all(bitstr[i] == target_bits[i] for i in active_idx)

def first_one_index(target_bits: str) -> int:
    try:
        return target_bits.index('1')
    except ValueError:
        return 0

def make_onehot_bits(nbits: int, index: int) -> str:
    """Skapa en one-hot-sträng av längd nbits, '1' på given index (0-bas, vänster-justerad)."""
    if nbits <= 0:
        return ""
    i = max(0, min(nbits - 1, int(index)))
    return "".join('1' if k == i else '0' for k in range(nbits))

# ========================= Active-set choisir =========================

def choose_active_indices(target_bits: str,
                          nbits: int,
                          mode: str = "subset-k",
                          k: int = 8,
                          seed: Optional[int] = None,
                          drop_override: Optional[int] = None) -> Tuple[List[int], int]:
    """
    Välj aktivmängd enligt:
      - drop-one: uteslut 'drop'-index, aktiva = övriga
      - subset-k: uteslut 'drop'-index, välj K av resterande (deterministiskt på seed)
      - full:     alla index aktiva
    drop väljs som:
      - drop_override om satt
      - annars första '1' i target_bits (one-hot)
    Returnerar (active_indices, dropped_index)
    """
    drop = drop_override if (drop_override is not None) else first_one_index(target_bits)
    drop = max(0, min(nbits-1, drop)) if nbits > 0 else 0
    if mode == "drop-one":
        active = [i for i in range(nbits) if i != drop]
        return active, drop
    elif mode == "full":
        return list(range(nbits)), -1
    else:
        pool = [i for i in range(nbits) if i != drop]
        k_eff = max(1, min(len(pool), int(k)))
        tb_int = int(target_bits or "0", 2)
        seed_base = seed if (seed is not None) else 0
        mix = (seed_base * 1000003) ^ tb_int
        rng = random.Random(mix & 0xFFFFFFFF)
        active = sorted(rng.sample(pool, k_eff))
        return active, drop

# ========================= Status/Progress =========================

def print_progress_banner(pct: int):
    """Stor, tydlig banner + progress-bar."""
    pct = max(0, min(100, int(pct)))
    bar_len = 40
    filled = int(bar_len * pct / 100)
    bar = '█' * filled + '·' * (bar_len - filled)
    line = "=" * 64
    print("\n" + line)
    print(f"  🚀 PROGRESS: {pct:3d}%  [{bar}]")
    print(line + "\n")

# ============================= Backend =============================

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
    """Kompilerar och kör kretsen på valt backend och returnerar sannolikheter."""
    tqc = transpile(qc, bundle.backend)
    if bundle.runner_kind == "aer":
        result = bundle.backend.run(tqc, shots=shots).result()
    else:
        job = bundle.backend.run(tqc, shots=shots)
        result = job.result()
    counts = result.get_counts(tqc)
    return probs_from_counts(counts, shots)

# ==================== Oracle, Diffusion, Grover ====================

def apply_mcz_on(qubits_list, circ: QuantumCircuit, anc_reg: Optional[AncillaRegister]):
    """MCZ realiserad via MCX-tricket."""
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
    """Markerar måltillstånd (enligt aktiva index) genom fasflip."""
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
    for q in active_qubits: circ.h(q)
    for q in active_qubits: circ.x(q)
    apply_mcz_on(active_qubits, circ, anc_reg)
    for q in active_qubits: circ.x(q)
    for q in active_qubits: circ.h(q)

def build_grover_layered_circuit(nbits: int,
                                 target_bits: str,
                                 active_indices: List[int],
                                 grover_iters: int = 1,
                                 oracle_type: str = "standard",
                                 weight_angle: float = math.pi,
                                 oracle_power: int = 1) -> QuantumCircuit:
    """Bygger krets med initial H, oracle, diffusion, m gånger, och mäter (endian-fix q[::-1])."""
    q = QuantumRegister(nbits, "q")
    anc = AncillaRegister(max(0, nbits - 2), "anc")
    c = ClassicalRegister(nbits, "c")
    qc = QuantumCircuit(q, anc, c, name="QCDS")
    # Init superposition
    for i in range(nbits):
        qc.h(q[i])
    # Grover-itereringar
    for _ in range(max(1, grover_iters)):
        apply_oracle_phase_flip_on_state(qc, q, active_indices, target_bits, anc,
                                         oracle_type=oracle_type,
                                         weight_angle=weight_angle,
                                         oracle_power=oracle_power)
        apply_diffusion_on_active(qc, q, active_indices, anc)
    qc.barrier()
    # ENDIAN-FIX: mät i omvänd qubitordning så att bit 0 i target_bits motsvarar q[0]
    qc.measure(q[::-1], c)
    return qc

# ===================== Sched / Heuristik =====================

def rotating_exclusion(nbits: int, layer_index_1based: int) -> List[int]:
    if nbits <= 1:
        return [0] if nbits == 1 else []
    drop = (layer_index_1based - 1) % nbits
    return [i for i in range(nbits) if i != drop]

def parse_iters_sched(spec: str, L: int, max_m: int) -> List[int]:
    """
    - "adaptive": placeholder (används ej direkt i 'coarse' utan styr adaptiva delen)
    - "auto": 1,1,3,5,7,...
    - annars: kommaseparerad lista
    """
    spec = (spec or "").strip().lower()
    if spec == "adaptive":
        return [1] * L
    if spec in ("auto", ""):
        seq, m, step = [], 1, 2
        while len(seq) < L:
            seq.append(min(max_m, m))
            if len(seq) % 2 == 0 and m < max_m:
                m = min(max_m, m + step)
        return seq
    raw = [s for s in spec.replace(" ", "").split(",") if s]
    vals = []
    for s in raw:
        try: vals.append(int(s))
        except: vals.append(1)
    if not vals: vals = [1]
    while len(vals) < L: vals.append(vals[-1])
    return [max(1, min(max_m, v)) for v in vals[:L]]

def optimal_grover_iters(p0: float, max_m: int) -> int:
    """Teoretiskt m* ≈ round(π/(4θ)-1/2), med sin^2 θ = p0."""
    p0 = max(1e-12, min(1.0 - 1e-12, p0))
    theta = math.asin(math.sqrt(p0))
    m_star = int(round((math.pi / (4.0 * theta)) - 0.5))
    return max(1, min(max_m, m_star))

# ================= Pilot p0-estimering =================

def pilot_estimate_m_star(bundle: BackendBundle,
                          nbits: int,
                          target_bits: str,
                          active_indices: List[int],
                          shots: int,
                          grover_max_iters: int,
                          oracle_type: str,
                          weight_angle: float,
                          oracle_power: int,
                          score_mode: str = "p-true",
                          layman: bool = False) -> Tuple[int, float]:
    """Kör m=1 för att skatta p0 enligt valt scoreläge och deriv m*."""
    qc = build_grover_layered_circuit(nbits, target_bits, active_indices, 1,
                                      oracle_type=oracle_type,
                                      weight_angle=weight_angle,
                                      oracle_power=oracle_power)
    probs = run_sampler_once(bundle, qc, shots)
    active_sum = sum(p for s, p in probs.items() if is_consistent_on_active(s, target_bits, active_indices))
    maxprob = max(probs.values()) if probs else 0.0
    p_true = probs.get(target_bits, 0.0)
    if score_mode == "active-sum": p0 = active_sum
    elif score_mode == "maxprob":  p0 = maxprob
    elif score_mode == "hybrid":   p0 = 0.5 * maxprob + 0.5 * active_sum
    else:                          p0 = p_true
    m_star = optimal_grover_iters(p0, grover_max_iters)
    if layman:
        print(f"   🧪 Pilot p₀≈{p0:.6f}  ⇒  m*≈{m_star} (shots={shots}, score={score_mode})")
    return m_star, p0

# ===================== Score-motor =====================

def compute_components(probs: Dict[str, float],
                       target_bits: str,
                       active_indices: List[int]) -> Tuple[float, float, float, str]:
    """Returnerar (maxprob, active_sum, p_true, top_str)."""
    if not probs:
        return 0.0, 0.0, 0.0, ""
    top_str, top_p = max(probs.items(), key=lambda kv: kv[1])
    p_true = probs.get(target_bits, 0.0)
    active_sum = sum(p for s, p in probs.items() if is_consistent_on_active(s, target_bits, active_indices))
    return top_p, active_sum, p_true, top_str

def compute_score(probs: Dict[str, float],
                  score_mode: str,
                  target_bits: str,
                  active_indices: List[int],
                  multihit_ratio: float,
                  w_maxprob: float = 1.0,
                  w_ptrue: float = 1.0,
                  w_active: float = 0.5) -> Tuple[float, Dict[str, float], str, float]:
    """
    Returnerar (score, hits, top_str, p_true_exact).
    modes:
      - maxprob: score = max(prob). hits = ≥ r*max
      - active-sum: score = sum(prob över aktiva-konsistenta). hits = alla aktiva-konsistenta
      - p-true: score = p_true
      - hybrid: score = w_m*maxprob + w_t*p_true + w_a*active_sum
    """
    if not probs:
        return 0.0, {}, "", 0.0
    maxprob, active_sum, p_true_exact, top_str = compute_components(probs, target_bits, active_indices)
    if score_mode == "active-sum":
        hits = {s: p for s, p in probs.items() if is_consistent_on_active(s, target_bits, active_indices)}
        return active_sum, hits, top_str, p_true_exact
    thr = multihit_ratio * maxprob if maxprob > 0 else 1.0
    hits = {s: p for s, p in probs.items() if p >= thr}
    if score_mode == "p-true":
        return p_true_exact, hits, top_str, p_true_exact
    if score_mode == "hybrid":
        score = (w_maxprob * maxprob) + (w_ptrue * p_true_exact) + (w_active * active_sum)
        return score, hits, top_str, p_true_exact
    return maxprob, hits, top_str, p_true_exact  # default: maxprob

# =========== Dynamisk Grover (scout/confirm separata) ===========

def dynamic_grover_peak_search(bundle: BackendBundle,
                               nbits: int,
                               target_bits: str,
                               active_indices: List[int],
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
                               score_mode_scout: str = "maxprob",
                               score_mode_confirm: str = "hybrid",
                               w_maxprob: float = 1.0,
                               w_ptrue: float = 1.0,
                               w_active: float = 0.5,
                               layman: bool = False,
                               print_hits: bool = False,
                               confirm_span: int = 1) -> Tuple[int, Dict[str, float], Dict[str, float], list]:
    """
    Sveper m i [m_start, m_end] och hittar topp (scout). Bekräftar lokalt runt bästa m* (confirm).
    Early-stop via 'patience' och 'drop_tol'.
    """
    if m_end is None:
        m_end = m_max
    m_start = max(1, min(m_start, m_max))
    m_end = max(m_start, min(m_end, m_max))
    confirm_span = max(1, int(confirm_span))
    if layman and (m_start != 1 or m_end != m_max):
        print(f"   🔎 Sökfönster: m ∈ [{m_start}, {m_end}] av max {m_max}  |  scout={score_mode_scout}, confirm={score_mode_confirm}")

    history, best_score, best_m, no_improve, prev_s = [], -1.0, m_start, 0, None
    for m in range(m_start, m_end + 1):
        qc = build_grover_layered_circuit(nbits, target_bits, active_indices, m,
                                          oracle_type=oracle_type, weight_angle=weight_angle, oracle_power=oracle_power)
        probs = run_sampler_once(bundle, qc, shots_scout)
        score, _, top_str, _ = compute_score(
            probs, score_mode_scout, target_bits, active_indices, multihit_ratio,
            w_maxprob=w_maxprob, w_ptrue=w_ptrue, w_active=w_active
        )
        history.append(score)
        if layman: print(f"   ↻ Scout m={m:>2}  score={score:.4f}  top={top_str}")
        if score > best_score + 1e-12:
            best_score, best_m, no_improve = score, m, 0
        else:
            no_improve += 1
        if prev_s is not None and score < (prev_s * (1.0 - drop_tol)):
            if layman: print(f"   ⛳ Stoppar: score föll > {int(drop_tol*100)}% (m={m-1}→{m})")
            break
        if no_improve >= patience:
            if layman: print(f"   ⛳ Early stop (patience={patience}) vid m={m}")
            break
        prev_s = score

    # Bekräfta runt bästa m
    ms = sorted({x for x in range(best_m - confirm_span, best_m + confirm_span + 1)
                 if m_start <= x <= m_end})
    best_m_conf, best_score_conf, best_probs_conf = best_m, -1.0, {}
    for mm in ms:
        qc = build_grover_layered_circuit(nbits, target_bits, active_indices, mm,
                                          oracle_type=oracle_type, weight_angle=weight_angle, oracle_power=oracle_power)
        probs = run_sampler_once(bundle, qc, shots_confirm)
        score, _, _, _ = compute_score(
            probs, score_mode_confirm, target_bits, active_indices, multihit_ratio,
            w_maxprob=w_maxprob, w_ptrue=w_ptrue, w_active=w_active
        )
        if score > best_score_conf:
            best_score_conf, best_m_conf, best_probs_conf = score, mm, probs

    probs_confirm = best_probs_conf
    hits_max, hits_consistent_active = {}, {}
    if probs_confirm:
        max_p = max(probs_confirm.values())
        thr = multihit_ratio * max_p
        # ⚠️ FIX av tidigare NameError: använd v, inte 'p'
        hits_max = {k: v for k, v in probs_confirm.items() if v >= thr}
        hits_consistent_active = {s: p for s, p in probs_confirm.items()
                                  if is_consistent_on_active(s, target_bits, active_indices)}
    if layman:
        s_final, _, top_str_conf, p_true_conf = compute_score(
            probs_confirm, score_mode_confirm, target_bits, active_indices, multihit_ratio,
            w_maxprob=w_maxprob, w_ptrue=w_ptrue, w_active=w_active
        )
        print(f"   ✅ Bekräftat m*={best_m_conf}  score≈{s_final:.4f}  p_true≈{p_true_conf:.4f}  (confirm_shots={shots_confirm})")
        print(f"   ░ score(m): {ascii_spark(history)}")
        agree_full = agree_fraction(top_str_conf, target_bits) if top_str_conf else 0.0
        agree_active = (sum(1 for i in active_indices if top_str_conf and top_str_conf[i] == target_bits[i]) /
                        max(1, len(active_indices))) if top_str_conf else 0.0
        print(f"   🔎 agreement: full={agree_full:.2%}, active={agree_active:.2%}")
        top_disp = top_str_conf if top_str_conf else "-"
        print(f"   ⭐ multi-hit (≥{int(multihit_ratio*100)}% av max): count={len(hits_max)}, "
              f"target-konsistenta(aktiva)={len(hits_consistent_active)}, top={top_disp}")
    return best_m_conf, probs_confirm, hits_max, history

# ==================== Lagerkörning (coarse + adaptiv) ====================

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
                adaptive_enabled: bool = False,
                adaptive_scout_shots: int = 256,
                adaptive_confirm_shots: int = 4096,
                adaptive_patience: int = 2,
                adaptive_drop_tol: float = 0.05,
                adaptive_multihit_ratio: float = 0.90,
                pilot_p0_enabled: bool = False,
                pilot_p0_shots: int = 4096,
                pilot_window: int = 8,
                pilot_score_mode: str = "p-true",
                print_hits: bool = False,
                score_mode_scout: str = "maxprob",
                score_mode_confirm: str = "hybrid",
                w_maxprob: float = 1.0,
                w_ptrue: float = 1.0,
                w_active: float = 0.5,
                active_mode: str = "subset-k",
                active_k: int = 8,
                active_seed: Optional[int] = None,
                active_drop_index: Optional[int] = None,
                manual_search_start: Optional[int] = None,
                manual_search_end: Optional[int] = None,
                confirm_span: int = 1) -> Tuple[Optional[str], Optional[float], Optional[int]]:
    """
    Kör L "coarse" lager (bias/entropi-indikator), därefter adaptiv Grover med fast active-set.
    Returnerar (top, p_true, m*).
    """
    iters_sched = parse_iters_sched(iters_sched_spec, L, grover_max_iters)
    if layman:
        print(f"▶️ QCDS Hierarkisk {L}-lager (backend={bundle.backend_name})")
    last_probs = None
    last_top = None
    p_true_last = None
    for layer_idx in range(1, L + 1):
        active_layer = rotating_exclusion(nbits, layer_idx)
        iters = iters_sched[layer_idx - 1]
        qc = build_grover_layered_circuit(nbits, target_bits, active_layer, iters,
                                          oracle_type=oracle_type, weight_angle=weight_angle, oracle_power=oracle_power)
        probs = run_sampler_once(bundle, qc, shots)
        top = max(probs.items(), key=lambda kv: kv[1])[0] if probs else None
        p_true = probs.get(target_bits, 0.0)
        H = shannon_H_from_probs(probs)
        uniform = H / max(1, nbits)
        if layman:
            print(f"▶️ QCDS_L{layer_idx} (iters={iters}) — backend: {bundle.backend_name} — top: {top}  p_true≈{p_true:.4f}")
            print(f"   QCDS H({nbits} qubits): {H:.4f} | uniformitet≈{uniform:.4f}")
        last_probs, last_top, p_true_last = probs, top, p_true

    final_top = None
    final_p_true = None
    final_m = None

    if do_pilot_grover and last_probs is not None:
        active_fixed, dropped = choose_active_indices(
            target_bits, nbits, mode=active_mode, k=active_k, seed=active_seed, drop_override=active_drop_index
        )
        if layman:
            k_txt = f"K={active_k}" if active_mode == "subset-k" else "K=N/A"
            print(f"   🧭 Active-set mode={active_mode} ({k_txt}), drop={dropped}; scout={score_mode_scout}, confirm={score_mode_confirm}")

        # Bestäm m-intervall
        if manual_search_start is not None or manual_search_end is not None:
            m_start = max(1, manual_search_start if manual_search_start is not None else 1)
            m_end = min(grover_max_iters, manual_search_end if manual_search_end is not None else grover_max_iters)
            if layman: print(f"   🎛️ Manuellt sökfönster aktiverat: m ∈ [{m_start}, {m_end}]")
        elif adaptive_enabled or (iters_sched_spec.strip().lower() == "adaptive"):
            m_start = 1; m_end = grover_max_iters
            if pilot_p0_enabled:
                m_star, _ = pilot_estimate_m_star(
                    bundle, nbits, target_bits, active_indices=active_fixed, shots=pilot_p0_shots,
                    grover_max_iters=grover_max_iters, oracle_type=oracle_type,
                    weight_angle=weight_angle, oracle_power=oracle_power,
                    score_mode=pilot_score_mode, layman=layman
                )
                m_start = max(1, m_star - pilot_window)
                m_end = min(grover_max_iters, m_star + pilot_window)
                if layman: print(f"   🎯 Pilotfönster: m ∈ [{m_start}, {m_end}] (window={pilot_window})")
        else:
            # fallback: enkel pilot runt teoretisk m* baserat på p_true_last
            ptrue = p_true_last if p_true_last is not None else 0.0
            m_star = optimal_grover_iters(ptrue, grover_max_iters)
            m_start = max(1, m_star - pilot_window)
            m_end = min(grover_max_iters, m_star + pilot_window)
            if layman:
                print(f"   🎯 Fallback-pilotfönster: m ∈ [{m_start}, {m_end}] (window={pilot_window})")

        # Kör dynamisk peak-sök
        if layman: print("   🔄 Dynamisk Grover-toppsök (fast oracle) startar…")
        m_best, probs_best, hits, history = dynamic_grover_peak_search(
            bundle, nbits, target_bits, active_indices=active_fixed, m_max=grover_max_iters,
            oracle_type=oracle_type, weight_angle=weight_angle, oracle_power=oracle_power,
            shots_scout=adaptive_scout_shots, shots_confirm=adaptive_confirm_shots,
            patience=adaptive_patience, drop_tol=adaptive_drop_tol, multihit_ratio=adaptive_multihit_ratio,
            m_start=m_start, m_end=m_end,
            score_mode_scout=score_mode_scout, score_mode_confirm=score_mode_confirm,
            w_maxprob=w_maxprob, w_ptrue=w_ptrue, w_active=w_active,
            layman=layman, print_hits=print_hits, confirm_span=confirm_span
        )

        final_score, _, top_final, p_true_final = compute_score(
            probs_best, score_mode_confirm, target_bits, active_fixed, adaptive_multihit_ratio,
            w_maxprob=w_maxprob, w_ptrue=w_ptrue, w_active=w_active
        )
        final_top, final_p_true, final_m = top_final, p_true_final, m_best
        if layman:
            print(f"   ↳ Dynamisk Grover: m*={m_best}, score≈{final_score:.4f}, p_true≈{p_true_final:.4f}, top={top_final}")

    return (final_top if final_top is not None else last_top), \
           (final_p_true if final_p_true is not None else p_true_last), \
           final_m

# ===================== Förklaringar & Mutationslista =====================

MUTATION_INDEX = {
    "BRCA1_E23Vfs": 0,
    "BRCA2_5946delT": 1,
    "TP53_R175H": 2,
    "KRAS_G12D": 3,
    "EGFR_L858R": 4,
    "BRAF_V600E": 5,
    "PIK3CA_H1047R": 6,
    "PTEN_R233X": 7
}

EXPLAINS = {
    "BRCA1_E23Vfs": {
        "tech": "Tidig frameshift i BRCA1 → förlorad DNA-reparation (HR).",
        "lay":  "Som att säkerhetsbältet är avklippt redan innan bilen startar."
    },
    "BRCA2_5946delT": {
        "tech": "Deletion i BRCA2 → slår ut homolog rekombination.",
        "lay":  "Reparatören missar sista skruven – maskinen rasar till slut."
    },
    "TP53_R175H": {
        "tech": "Hotspot i p53 → tumörsuppressorn blir inaktiv.",
        "lay":  "Brandlarmet är urkopplat – bränder upptäcks inte."
    },
    "KRAS_G12D": {
        "tech": "Aktiverande KRAS → stark RAS-signalering (proliferation).",
        "lay":  "Gaspedalen fastnar i botten – cellen rusar."
    },
    "EGFR_L858R": {
        "tech": "Aktiverande EGFR → ihållande tillväxtsignal.",
        "lay":  "Telefonen ringer konstant utan att någon ringer."
    },
    "BRAF_V600E": {
        "tech": "Aktiverande BRAF-kinas → MAPK-vägen låst i 'på'.",
        "lay":  "Strömbrytaren fastnar i 'på' – lampan slocknar aldrig."
    },
    "PIK3CA_H1047R": {
        "tech": "Aktiverande PI3K → tillväxt/överlevnad upp.",
        "lay":  "Oändlig konstgödsel – växten fortsätter växa."
    },
    "PTEN_R233X": {
        "tech": "Stopkodon i PTEN → bromsen i PI3K/AKT försvinner.",
        "lay":  "Bromskabeln på cykeln är avklippt – inget saktar ner."
    },
}

# ===================== Sammanställning & Hypotes =====================

def agreement_str(top: Optional[str], target: str) -> str:
    return "100%" if (top == target and top not in (None, "")) else "—"

def phenotype_hypothesis(result_map: Dict[str, float]) -> Dict[str, str]:
    """
    Tar maxade p_true per nyckelväg och ger en icke-klinisk hypotes + “kommande trend”.
    Trösklar valda konservativt (heuristik).
    """
    thr_strong = 0.05
    thr_present = 0.02

    hrd = max(result_map.get("BRCA1_E23Vfs", 0.0), result_map.get("BRCA2_5946delT", 0.0))
    p53 = result_map.get("TP53_R175H", 0.0)
    mapk = max(result_map.get("KRAS_G12D", 0.0), result_map.get("EGFR_L858R", 0.0), result_map.get("BRAF_V600E", 0.0))
    pi3k = max(result_map.get("PIK3CA_H1047R", 0.0), result_map.get("PTEN_R233X", 0.0))

    lines = []
    if hrd >= thr_strong:
        lines.append("• **HRD-liknande profil** (BRCA1/2) – bristande DNA-reparation.")
    elif hrd >= thr_present:
        lines.append("• Tecken på HRD-likhet (BRCA1/2).")

    if p53 >= thr_strong:
        lines.append("• **p53-funktionsförlust** – ökad genomisk instabilitet.")
    elif p53 >= thr_present:
        lines.append("• Möjligen påverkan på p53-signalering.")

    if mapk >= thr_strong:
        lines.append("• **MAPK-driven tillväxt** (KRAS/EGFR/BRAF).")
    elif mapk >= thr_present:
        lines.append("• Tecken på MAPK-aktivering.")

    if pi3k >= thr_strong:
        lines.append("• **PI3K/AKT-aktivering** (PIK3CA/PTEN).")
    elif pi3k >= thr_present:
        lines.append("• Tecken på PI3K/AKT-aktivering.")

    if not lines:
        lines.append("• Ingen tydlig dominerande signal i modellen (under trösklarna).")

    trend = []
    ranked = sorted(result_map.items(), key=lambda kv: kv[1], reverse=True)
    top3 = ranked[:3]
    if top3 and top3[0][1] >= thr_present:
        trend.append("• **Kommande trend** (modell): " +
                     ", ".join([f"{k} (p_true {v:.3f})" for k, v in top3]) +
                     " – dessa tillstånd ser ut att dominera evidensen.")
    else:
        trend.append("• **Kommande trend** (modell): inga starka signaler (låga p_true).")

    # ”Diagnos”-språk (icke-kliniskt) – formulerat försiktigt
    diagnosis_lines = []
    if hrd >= thr_strong:
        diagnosis_lines.append("• Modellens diagnoslik tolkning: **HRD-präglad tumörbiologi**.")
    if mapk >= thr_strong:
        diagnosis_lines.append("• Modellens diagnoslik tolkning: **MAPK-beroende proliferation**.")
    if pi3k >= thr_strong:
        diagnosis_lines.append("• Modellens diagnoslik tolkning: **PI3K/AKT-aktiverad fenotyp**.")
    if p53 >= thr_strong:
        diagnosis_lines.append("• Modellens diagnoslik tolkning: **p53-dysfunktion** (instabilitet).")
    if not diagnosis_lines:
        diagnosis_lines.append("• Ingen enskild dominerande diagnoslik profil i modellen.")

    disclaimer = ("⚠️ Detta är en *icke-klinisk modellhypotes* baserad på QCDS-resultaten. "
                  "Tolkning/diagnos kräver klinisk kontext och medicinsk expertis.")
    return {
        "phenotype": "\n".join(lines),
        "trend": "\n".join(trend),
        "diagnosis_like": "\n".join(diagnosis_lines),
        "disclaimer": disclaimer
    }

def emit_table_markdown(rows: List[Tuple[str, str, str, str, str, str]]) -> str:
    header = "| Variant (index) | Teknisk beskrivning (kort) | Lekmannaförklaring | p_true (SANNING) | m* | Agreement |\n" \
             "|---|---|---|---:|--:|--:|"
    lines = [header]
    for r in rows:
        lines.append(f"| {r[0]} | {r[1]} | {r[2]} | {r[3]} | {r[4]} | {r[5]} |")
    return "\n".join(lines)

# ================================ Main ================================

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", type=str, default="aer", help="aer | ibm")
    ap.add_argument("--ibm-backend", type=str, default="ibm_brisbane")
    ap.add_argument("--genome", type=str, default=None, help="Tagg/namn för körning (används i filnamn).")
    ap.add_argument("--shots", type=int, default=2048)
    ap.add_argument("--grover-max-iters", type=int, default=60)
    ap.add_argument("--evidence-bits", type=int, default=12)
    ap.add_argument("--layers", type=int, default=6)
    ap.add_argument("--iters-sched", type=str, default="adaptive")
    ap.add_argument("--oracle-type", type=str, choices=["standard", "weighted"], default="standard")
    ap.add_argument("--oracle-weight-angle", type=float, default=math.pi)
    ap.add_argument("--oracle-power", type=int, default=1)
    ap.add_argument("--layman", action="store_true")
    ap.add_argument("--out-dir", type=str, default="results")
    ap.add_argument("--adaptive", action="store_true")
    ap.add_argument("--scout-shots", type=int, default=256)
    ap.add_argument("--confirm-shots", type=int, default=4096)
    ap.add_argument("--patience", type=int, default=2)
    ap.add_argument("--drop-tol", type=float, default=0.05)
    ap.add_argument("--multihit-threshold", type=float, default=0.90)
    ap.add_argument("--pilot-p0", action="store_true")
    ap.add_argument("--pilot-p0-shots", type=int, default=4096)
    ap.add_argument("--pilot-window", type=int, default=10)
    ap.add_argument("--pilot-score", type=str, choices=["maxprob","active-sum","p-true","hybrid"], default="p-true")
    ap.add_argument("--search-start", type=int, default=None)
    ap.add_argument("--search-end", type=int, default=None)
    ap.add_argument("--score-scout", type=str, choices=["maxprob","active-sum","p-true","hybrid"], default="p-true")
    ap.add_argument("--score-confirm", type=str, choices=["maxprob","active-sum","p-true","hybrid"], default="p-true")
    ap.add_argument("--hybrid-w-maxprob", type=float, default=1.0)
    ap.add_argument("--hybrid-w-ptrue", type=float, default=1.0)
    ap.add_argument("--hybrid-w-active", type=float, default=0.5)
    ap.add_argument("--active-mode", type=str, choices=["drop-one","subset-k","full"], default="full")
    ap.add_argument("--active-k", type=int, default=8)
    ap.add_argument("--active-seed", type=int, default=None)
    ap.add_argument("--active-drop-index", type=int, default=None)
    ap.add_argument("--print-hits", action="store_true")
    ap.add_argument("--progress", type=str, choices=["on","off"], default="on")
    ap.add_argument("--progress-steps", type=int, default=10)
    ap.add_argument("--target-index", type=int, default=0, help="Vilken qubit ska vara '1' i one-hot-target.")
    ap.add_argument("--target-bits", type=str, default=None, help="Fri bitmask, override index.")
    # Batch/Single
    ap.add_argument("--batch", action="store_true", help="Kör alla fördefinierade mutationer (8 st).")
    ap.add_argument("--single", type=str, default=None, help="Kör en enda variant med namn i MUTATION_INDEX.")
    ap.add_argument("--confirm-span", type=int, default=1, help="Bekräftelse-span för m: {m*−s .. m*+s}")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    bundle = init_backend(mode=args.mode, ibm_backend_name=args.ibm_backend, seed_sim=None)

    # ---- Hjälpfunktion: kör EN variant och returnera resultat ----
    def run_variant(name: str, index: int) -> Dict[str, any]:
        # Sätt target_bits/nbits
        if args.target_bits:
            tb = "".join('1' if ch == '1' else '0' for ch in args.target_bits.strip())
            target_bits = tb
            nbits = len(tb)
        else:
            nbits = args.evidence_bits
            target_bits = make_onehot_bits(nbits, index)

        if args.layman:
            src = "target-bits override" if args.target_bits else f"target-index={index}"
            print("\n" + "—"*40)
            print(f"🔬 {name} (index {index})")
            print(f"▶️ Start QCDS-run (nbits={nbits}, target={target_bits}, {src})")

            expl = EXPLAINS.get(name, {})
            print(f"  • Teknisk: {expl.get('tech','(okänd)')}")
            print(f"  • Lekman:  {expl.get('lay','(okänd)')}")

        top, p_true, m_star = qcds_layers(
            bundle=bundle,
            nbits=nbits, target_bits=target_bits,
            shots=args.shots, L=args.layers,
            iters_sched_spec=args.iters_sched, grover_max_iters=args.grover_max_iters,
            layman=args.layman, do_pilot_grover=True,
            oracle_type=args.oracle_type, weight_angle=args.oracle_weight_angle, oracle_power=args.oracle_power,
            adaptive_enabled=args.adaptive, adaptive_scout_shots=args.scout_shots, adaptive_confirm_shots=args.confirm_shots,
            adaptive_patience=args.patience, adaptive_drop_tol=args.drop_tol, adaptive_multihit_ratio=args.multihit_threshold,
            pilot_p0_enabled=args.pilot_p0, pilot_p0_shots=args.pilot_p0_shots, pilot_window=args.pilot_window,
            pilot_score_mode=args.pilot_score,
            print_hits=args.print_hits, score_mode_scout=args.score_scout, score_mode_confirm=args.score_confirm,
            w_maxprob=args.hybrid_w_maxprob, w_ptrue=args.hybrid_w_ptrue, w_active=args.hybrid_w_active,
            active_mode=args.active_mode, active_k=args.active_k, active_seed=args.active_seed,
            active_drop_index=args.active_drop_index,
            manual_search_start=args.search_start, manual_search_end=args.search_end,
            confirm_span=args.confirm_span
        )

        # Spara JSON
        genome_tag = args.genome if args.genome else name
        fname = os.path.join(args.out_dir, f"QCDS_end2end_{genome_tag}_{bundle.backend_name}.json")
        data = {
            "variant": name,
            "index": index,
            "top": top,
            "p_true": p_true,
            "target": target_bits,
            "nbits": len(target_bits),
            "backend": bundle.backend_name,
            "m_star": m_star,
            "score_scout": args.score_scout,
            "score_confirm": args.score_confirm,
            "pilot_score": args.pilot_score,
            "confirm_span": args.confirm_span,
            "multihit_threshold": args.multihit_threshold
        }
        with open(fname, "w") as f:
            json.dump(data, f, indent=2)
        if args.layman:
            print(f"💾 Sparade: {fname}")
        return data

    if args.progress == "on":
        print_progress_banner(0)

    results: List[Dict[str, any]] = []

    if args.batch:
        # Kör alla åtta
        for name, idx in MUTATION_INDEX.items():
            res = run_variant(name, idx)
            results.append(res)
    else:
        # Enskild
        if args.single:
            if args.single not in MUTATION_INDEX:
                print(f"❌ --single {args.single} okänt, välj en av: {', '.join(MUTATION_INDEX.keys())}")
                sys.exit(2)
            idx = MUTATION_INDEX[args.single]
            results.append(run_variant(args.single, idx))
        else:
            # Om inget valt: kör BRCA1 som demo
            name = "BRCA1_E23Vfs"; idx = MUTATION_INDEX[name]
            results.append(run_variant(name, idx))

    # =================== Sammanställning ===================
    rows = []
    ptrue_map = {}
    for r in results:
        name = r["variant"]
        idx = r["index"]
        tech = EXPLAINS.get(name, {}).get("tech", "(okänd)")
        lay  = EXPLAINS.get(name, {}).get("lay", "(okänd)")
        p_true = r.get("p_true", None)
        m_star = r.get("m_star", None)
        top, target = r.get("top",""), r.get("target","")
        agree = agreement_str(top, target)
        p_disp = f"{p_true:.4f}" if isinstance(p_true, (int,float)) else "—"
        m_disp = f"{m_star}" if isinstance(m_star, int) else "—"
        rows.append((f"{name} ({idx})", tech, lay, p_disp, m_disp, agree))
        if isinstance(p_true, (int,float)):
            ptrue_map[name] = p_true

    md = emit_table_markdown(rows)
    summary_path = os.path.join(args.out_dir, "summary_mutations.md")
    with open(summary_path, "w") as f:
        f.write(md + "\n")
    print("\n" + md)
    print(f"\n📄 Skrev tabell till {summary_path}")

    # =================== Trend/Hypotes (icke-klinisk) ===================
    hyp = phenotype_hypothesis(ptrue_map)
    hyp_md = (
        "\n## Modellhypotes (icke-klinisk)\n"
        + hyp["phenotype"] + "\n\n"
        + hyp["trend"] + "\n\n"
        + hyp["diagnosis_like"] + "\n\n"
        + hyp["disclaimer"] + "\n"
    )
    with open(summary_path, "a") as f:
        f.write(hyp_md)
    print(hyp_md)

    if args.progress == "on":
        print_progress_banner(100)

if __name__ == "__main__":
    main()
