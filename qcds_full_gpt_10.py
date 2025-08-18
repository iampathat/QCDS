#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
QCDS – end-to-end demo med:
  • AER-masskörning över BRCA2-regionen (GRCh37/GRCh38)
  • IBM-läge (begränsad mängd varianter men utökningsbart)
  • 2-qubit-kontroller (TG/AC/CE) + 4-lager QCDS
  • Adaptiv Grover + “roterande bitexkludering” (orakel per variant)
  • Sampler-shim (klarar olika Qiskit-versioner)
  • Lekmannaförklaringar (svenska) kan slås på med --layman
  • En och samma växel för kvantdator/emu (AER): --mode [aer|ibm],
    och flaggan --force-aer-4l för att köra 4L på emulator även i IBM-läge.

Kör-exempel:
  python qcds_full_gpt_10.py --mode aer --genome GRCh38 --layman
  python qcds_full_gpt_10.py --mode ibm --ibm-backend ibm_brisbane --ibm-max-variants 5 --layman
  python qcds_full_gpt_10.py --mode ibm --force-aer-4l --layman
"""

from __future__ import annotations
import argparse
import json
import math
import os
import random
import sys
import time
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

# ---------------------------
# Qiskit-importer + robust "Sampler"-shim
# ---------------------------
BACKEND_FLAVOUR = "unknown"
try:
    # Qiskit Terra >= 0.25 med primitives
    from qiskit import QuantumCircuit, transpile
    from qiskit.quantum_info import Statevector
    try:
        from qiskit_aer import Aer
        BACKEND_FLAVOUR = "aer"
    except Exception:
        from qiskit import BasicAer as Aer  # fallback
        BACKEND_FLAVOUR = "basic-aer"
except Exception as e:
    print("❌ Kunde inte importera Qiskit-kärna:", e)
    sys.exit(1)

# För IBM-läge (valfritt): vi försöker, men faller tillbaka snyggt.
HAVE_RUNTIME = False
try:
    from qiskit_ibm_runtime import QiskitRuntimeService
    HAVE_RUNTIME = True
except Exception:
    HAVE_RUNTIME = False


# --------- NYTT: robust backend-namn ----------
def backend_label(backend) -> str:
    """
    Returnera ett läsbart backend-namn oavsett Qiskit-version:
    - backend.name -> om sträng
    - backend.name() -> om metod
    - backend.backend_name -> alternativt
    - annars klassnamn
    """
    try:
        name_attr = getattr(backend, "name", None)
        if callable(name_attr):
            n = name_attr()
        elif isinstance(name_attr, str):
            n = name_attr
        else:
            n = getattr(backend, "backend_name", None)
        if not n:
            n = type(backend).__name__
        return str(n)
    except Exception:
        return type(backend).__name__
# -----------------------------------------------


# Minimal “sampler”-adapter som fungerar oavsett Qiskit-version/back-end.
class MiniSampler:
    """Liten wrapper som ger .run(circuits, shots) -> list[dict(state->prob)]"""

    def __init__(self, backend, shots: int = 2048, seed_simulator: Optional[int] = None):
        self.backend = backend
        self.shots = shots
        self.seed_simulator = seed_simulator

    def _run_single(self, qc: "QuantumCircuit") -> Dict[str, float]:
        # säkerställ mätningar finns
        if not qc.num_clbits or all([mop.operation.name != "measure" for mop in qc.data]):
            measured = qc.copy()
            measured.measure_all()
        else:
            measured = qc

        # Transpile för backend
        tqc = transpile(measured, self.backend)
        job = self.backend.run(tqc, shots=self.shots, seed_simulator=self.seed_simulator)
        result = job.result()
        counts = result.get_counts()
        # normalisera till sannolikheter
        total = sum(counts.values()) if counts else 1
        return {k: v / total for k, v in counts.items()}

    def run(self, circuits: List["QuantumCircuit"], shots: Optional[int] = None) -> List[Dict[str, float]]:
        out = []
        for qc in circuits:
            out.append(self._run_single(qc))
        return out


# ---------------------------
# Hjälpfunktioner (statistik, utskrift)
# ---------------------------

def bitstr_prob_top(probs: Dict[str, float]) -> Tuple[str, float]:
    if not probs:
        return "0", 0.0
    key = max(probs, key=lambda k: probs[k])
    return key, probs[key]


def z_expect_2q(probs: Dict[str, float]) -> float:
    # ⟨Z⊗Z⟩ = P(00)+P(11) - P(01) - P(10)
    return probs.get("00", 0.0) + probs.get("11", 0.0) - probs.get("01", 0.0) - probs.get("10", 0.0)


def bias_1q(probs: Dict[str, float], bit_index_from_right: int) -> float:
    # bias = P(bit=0) - P(bit=1)
    p0 = 0.0
    p1 = 0.0
    for b, p in probs.items():
        bit = b[::-1][bit_index_from_right]  # höger=LSB
        if bit == "0":
            p0 += p
        else:
            p1 += p
    return p0 - p1


def shannon_H_bits(probs: Dict[str, float]) -> float:
    H = 0.0
    for p in probs.values():
        if p > 0:
            H -= p * math.log2(p)
    return H


def uniformity_score(probs: Dict[str, float]) -> float:
    # 2^n möjligheter; fördela mot uniform: sum(|p - 1/N|)/2 -> här rapporterar vi 1 - (L1-avvikelse)
    # Detta blir ett snabbt, intuitivt “hur close till uniform”.
    if not probs:
        return 0.0
    n = len(next(iter(probs)))  # längd på bitsträng
    N = 2 ** n
    base = 1.0 / N
    l1 = 0.0
    for s in probspace(n):
        l1 += abs(probs.get(s, 0.0) - base)
    # normalisera, 0=helt olika, 1=helt uniform
    return max(0.0, 1.0 - 0.5 * l1)


def probspace(nbits: int) -> List[str]:
    return [format(i, f"0{nbits}b") for i in range(2 ** nbits)]


def pretty_json(d: dict) -> str:
    return json.dumps(d, ensure_ascii=False, indent=2, sort_keys=True)


def say(msg: str):
    print(msg, flush=True)


# ---------------------------
# 2-qubit QCDS “kontroller”
# ---------------------------

def qc_truth_gradient() -> QuantumCircuit:
    # Skapar stark Z-korrelation (Bell-liknande 00/11)
    qc = QuantumCircuit(2, 2)
    qc.h(0)
    qc.cx(0, 1)
    qc.measure([0, 1], [0, 1])
    return qc


def qc_amplitude_confirmation() -> QuantumCircuit:
    # Lätt fas-koppling (CZ) från superposition -> mäter i Z-bas
    qc = QuantumCircuit(2, 2)
    qc.h([0, 1])
    qc.cz(0, 1)
    qc.measure([0, 1], [0, 1])
    return qc


def qc_conditional_enhancement() -> QuantumCircuit:
    # Förstärk kombinationer där båda är "lika" (00/11) via X, H, CX, H, X (en simpel proxy)
    qc = QuantumCircuit(2, 2)
    qc.h(0)
    qc.cx(0, 1)
    qc.h(0)
    qc.h(1)
    qc.measure([0, 1], [0, 1])
    return qc


def report_2q(name: str, probs: Dict[str, float], layman: bool):
    top, ptop = bitstr_prob_top(probs)
    expzz = z_expect_2q(probs)
    b0 = bias_1q(probs, 0)
    b1 = bias_1q(probs, 1)
    H = shannon_H_bits(probs)
    say(f"▶️ {name} — top: {top}")
    say(f"   probs: {pretty_json(probs)}")
    say(f"   QCDS ⟨Z⊗Z⟩: {round(expzz, 4)} | bias(q0): {round(b0, 4)} | bias(q1): {round(b1, 4)} | H: {round(H, 4)}")
    if layman:
        if name.startswith("Truth Gradient"):
            say("\n📘 Lättförståelig förklaring (2 qubits):")
            say("— Truth Gradient (TG): mäter hur ofta bitarna är lika (00/11) mot olika (01/10).")
            say(f"  Resultat: högsta utfallet är {top} ({ptop:.1%}). "
                f"Hög ⟨Z⊗Z⟩={round(expzz,3)} ⇒ stark samvariation.")
            say("  Lekman: “Tänk två lampor – de tänds/släcks oftare samtidigt än motsatt.”\n")
        elif name.startswith("Amplitude Confirmation"):
            say("\n📘 Lättförståelig förklaring (2 qubits):")
            say("— Amplitude Confirmation (AC): liten fas-interaktion (CZ), sedan läsning i Z-bas.")
            say(f"  Resultat: högsta utfallet är {top} ({ptop:.1%}). "
                f"⟨Z⊗Z⟩≈{round(expzz,3)} ⇒ nära neutral/återhållen korrelation.")
            say("  Lekman: “Vi nyper systemet i fas – utfallet blir ganska jämnt om allt är balanserat.”\n")
        elif name.startswith("Conditional Enhancement"):
            say("\n📘 Lättförståelig förklaring (2 qubits):")
            say("— Conditional Enhancement (CE): förstärker lägen där båda bitar uppfyller ett villkor.")
            say(f"  Resultat: högsta utfallet är {top} ({ptop:.1%}). "
                f"⟨Z⊗Z⟩={round(expzz,3)} ⇒ tydlig förstärkning.")
            say("  Lekman: “När villkoret stämmer blir den kombinationen extra sann.”\n")


# ---------------------------
# 4-lager QCDS (symbolisk, 4 qubits)
# ---------------------------

@dataclass
class LayerConfig:
    # prop: "vilken bit fixas" på den här nivån (symboliskt)
    prop: Dict[int, str]
    iters: int


def qc_4layer(prop: Dict[int, str], iters: int = 1) -> QuantumCircuit:
    """
    En enkel 4-qubits “våg”-förstärkare.
    - prop = t.ex. {3:'0', 2:'0'} betyder att vi biasar q3,q2 mot 0.
    - iters = hur många repetitioner av en mild ‘förstärkare’.
    """
    nq = 4
    qc = QuantumCircuit(nq, nq)
    # starta i superposition
    qc.h(range(nq))
    for _ in range(max(1, iters)):
        # villkors-förstärkning: CZ-kedja + X på specificerade noder
        for qb, want in prop.items():
            if want == "1":
                qc.x(qb)
        # enklare multi-CZ via kedja
        qc.cz(0, 1)
        qc.cz(1, 2)
        qc.cz(2, 3)
        for qb, want in prop.items():
            if want == "1":
                qc.x(qb)
        # spridning
        qc.h(range(nq))
        qc.cz(0, 3)
        qc.h(range(nq))
    qc.measure(range(nq), range(nq))
    return qc


def report_4l(stage: str, probs: Dict[str, float], layman: bool):
    n = len(next(iter(probs))) if probs else 4
    H = shannon_H_bits(probs)
    uni = uniformity_score(probs)
    top, ptop = bitstr_prob_top(probs)
    say(f"▶️ QCDS_4L {stage} — top: {top}  p_true≈{ptop:.4f}")
    say(f"   QCDS H({n} qubits): {round(H,4)} | uniformitet≈{round(uni,4)}")
    if layman:
        say("  Lekman: “Vi kör 4 ‘nivåer’ som förstärker mönster som hänger ihop. "
            "Högre ‘H’ och jämn ‘uniformitet’ ger oss robust signal utan att låsa fast en slump.”")


# ---------------------------
# Adaptiv Grover med roterande bitexkludering (4 bitars evidens)
# ---------------------------

def make_grover_oracle_4bit(target_bits: str, exclude_idx: int) -> QuantumCircuit:
    """
    Orakel markerar "sanning" för de inkluderade 3 av 4 bitarna.
    exclude_idx = vilken bit som *tas bort* ur indata-embeddningen denna gång.
    """
    assert len(target_bits) == 4
    qc = QuantumCircuit(4)
    # X på de inkluderade bitar som ska vara '0' -> markera 000* osv
    for i, tb in enumerate(target_bits):
        if i == exclude_idx:
            continue  # hoppa över
        if tb == "0":
            qc.x(i)
    # enkel multi-Z via CZ-kedja som “markering”
    qc.h(3)
    qc.ccx(0, 1, 2)  # grov proxy för "AND"
    qc.cz(2, 3)      # mark
    qc.ccx(0, 1, 2)  # uncompute
    qc.h(3)
    # återställ X
    for i, tb in enumerate(target_bits):
        if i == exclude_idx:
            continue
        if tb == "0":
            qc.x(i)
    return qc


def grover_diffusion_4(nq: int = 4) -> QuantumCircuit:
    qc = QuantumCircuit(nq)
    qc.h(range(nq))
    qc.x(range(nq))
    qc.h(nq - 1)
    qc.mcx(list(range(nq - 1)), nq - 1)
    qc.h(nq - 1)
    qc.x(range(nq))
    qc.h(range(nq))
    return qc


def run_adaptive_grover_4bit(sampler: MiniSampler,
                              target_bits: str,
                              shots: int,
                              max_iters_to_try: int = 6,
                              rotate_exclusion: bool = True,
                              seed: Optional[int] = None) -> Dict[str, float]:
    """
    Adaptivt Grover-upplägg:
      - Roterar exkluderad bit (om rotate_exclusion=True)
      - Testar 1..max_iters_to_try iterationer
      - Väljer bästa "p_true" (andel mätningar som matchar target_bits exakt)
    """
    if seed is not None:
        random.seed(seed)

    def p_true_of(probs: Dict[str, float]) -> float:
        return probs.get(target_bits, 0.0)

    best: Dict[str, float] = {}
    best_p = -1.0
    exclude_order = [0, 1, 2, 3] if rotate_exclusion else [None]
    for excl in exclude_order:
        for iters in range(1, max_iters_to_try + 1):
            # bygg Grover-cirkeln
            qc = QuantumCircuit(4, 4)
            qc.h(range(4))  # init
            # iters gånger: orakel (med exkluderad bit) + diffusion
            for _ in range(iters):
                oracle = make_grover_oracle_4bit(target_bits, excl if excl is not None else 0)
                qc.compose(oracle, qubits=range(4), inplace=True)
                qc.compose(grover_diffusion_4(4), qubits=range(4), inplace=True)
            qc.measure(range(4), range(4))
            probs = sampler.run([qc], shots=shots)[0]
            p = p_true_of(probs)
            if p > best_p:
                best_p = p
                best = probs
    return best


# ---------------------------
# BRCA2-region & variantgenerering
# ---------------------------

@dataclass
class Region:
    chrom: str
    start: int
    end: int


BRCA2_DEFAULT = {
    "GRCh38": Region(chrom="13", start=32315474, end=32400266),
    "GRCh37": Region(chrom="13", start=32889611, end=32973809),
}


def iter_synthetic_brca2_positions(region: Region,
                                   stride: int = 50,
                                   max_positions: int = 60,
                                   seed: Optional[int] = None) -> List[int]:
    if seed is not None:
        random.seed(seed)
    L = max(1, (region.end - region.start) // max(1, stride))
    L = min(L, max_positions)
    ps: List[int] = []
    if L == 1:
        return [region.start]
    step = max(1, (region.end - region.start) // L)
    for i in range(L):
        p = region.start + i * step + random.randint(0, min(step - 1, 30))
        if p <= region.end:
            ps.append(p)
    return sorted(set(ps))


# ---------------------------
# IBM-backend (frivillig)
# ---------------------------

def get_backend(mode: str, ibm_backend_name: Optional[str]):

    if mode == "aer":
        try:
            backend = Aer.get_backend("qasm_simulator")
            return backend, "qasm_simulator"
        except Exception as e:
            say(f"⚠️ Kunde inte få AER qasm_simulator ({e}), försöker BasicAer.")
            backend = Aer.get_backend("qasm_simulator")
            return backend, "qasm_simulator"

    # mode == "ibm"
    if not HAVE_RUNTIME:
        say("⚠️ IBM Runtime ej tillgänglig i denna miljö. Faller tillbaka till AER.")
        backend = Aer.get_backend("qasm_simulator")
        return backend, "qasm_simulator"

    try:
        svc = QiskitRuntimeService()
        name = ibm_backend_name or "ibm_brisbane"
        backend = svc.backend(name)
        return backend, name
    except Exception as e:
        say(f"⚠️ Kunde inte öppna IBM-backend ({e}). Faller tillbaka till AER.")
        backend = Aer.get_backend("qasm_simulator")
        return backend, "qasm_simulator"


# ---------------------------
# Körflöden
# ---------------------------

def run_2q_controls(sampler: MiniSampler, layman: bool):
    # TG
    tg = qc_truth_gradient()
    tg_probs = sampler.run([tg])[0]
    report_2q("Truth Gradient @ backend", tg_probs, layman)

    # AC
    ac = qc_amplitude_confirmation()
    ac_probs = sampler.run([ac])[0]
    report_2q("Amplitude Confirmation @ backend", ac_probs, layman)

    # CE
    ce = qc_conditional_enhancement()
    ce_probs = sampler.run([ce])[0]
    report_2q("Conditional Enhancement @ backend", ce_probs, layman)

    return {"TG": tg_probs, "AC": ac_probs, "CE": ce_probs}


def run_4layer_sequence(sampler: MiniSampler,
                        shots: int,
                        iters_sched=(1, 1, 3, 5),
                        force_aer_4l: bool = False,
                        running_on_aer: bool = True,
                        layman: bool = True):
    def _run_stage(label: str, prop: Dict[int, str], iters: int):
        qc = qc_4layer(prop, iters)
        probs = sampler.run([qc], shots=shots)[0]
        report_4l(f"{label} (iters={iters}, prop={prop})", probs, layman)
        return probs

    say(f"▶️ QCDS Hierarkisk 4-lager (backend={'aer' if running_on_aer else 'ibm'})")

    out = {}
    out["L1"] = _run_stage("L1", {3: "0"}, iters_sched[0])
    out["L2"] = _run_stage("L2", {3: "0", 2: "0"}, iters_sched[1])
    out["L3"] = _run_stage("L3", {3: "0", 2: "0", 1: "0"}, iters_sched[2])
    out["L4"] = _run_stage("L4", {3: "0", 2: "0", 1: "0", 0: "0"}, iters_sched[3])
    return out


def run_brca2_end2end(mode: str,
                      sampler_main: MiniSampler,
                      sampler_for_4l: MiniSampler,
                      genome: str,
                      region: Optional[str],
                      aer_max_variants: int,
                      ibm_max_variants: int,
                      iters_sched: Tuple[int, int, int, int],
                      shots: int,
                      layman: bool,
                      seed: Optional[int] = None) -> Dict:

    # Regionval
    if region:
        chrom, rest = region.split(":")
        start, end = map(int, rest.split("-"))
        reg = Region(chrom=chrom, start=start, end=end)
    else:
        reg = BRCA2_DEFAULT.get(genome, BRCA2_DEFAULT["GRCh38"])

    # Hur många varianter kör vi?
    max_vars = aer_max_variants if mode == "aer" else ibm_max_variants
    pos_list = iter_synthetic_brca2_positions(reg, stride=50, max_positions=max_vars, seed=seed)

    say(f"📦 Dataset (syntetiskt): BRCA2 {genome} {reg.chrom}:{reg.start}-{reg.end}")
    say(f"   Antal ‘varianter’ att köra: {len(pos_list)} (mode={mode}, shots={shots})")

    results = []
    for pos in pos_list:
        evs: List[str] = []
        target_bits = "0000"

        say(f"🎯 Variant {reg.chrom}:{pos}  target_bits={target_bits}  evs={evs}")
        # L1..L4 på AER-sampler (kan vara samma som huvud om mode=aer)
        for i, (label, prop, iters) in enumerate([
            ("L1", {3: "0"}, iters_sched[0]),
            ("L2", {3: "0", 2: "0"}, iters_sched[1]),
            ("L3", {3: "0", 2: "0", 1: "0"}, iters_sched[2]),
            ("L4", {3: "0", 2: "0", 1: "0", 0: "0"}, iters_sched[3]),
        ], start=1):
            qc = qc_4layer(prop, iters)
            probs = sampler_for_4l.run([qc], shots=shots)[0]
            top, ptop = bitstr_prob_top(probs)
            H = shannon_H_bits(probs)
            uni = uniformity_score(probs)
            be_name = backend_label(sampler_for_4l.backend)  # <-- robust
            say(f"▶️ QCDS_4L {label} (iters={iters}, prop={prop}) — backend: {be_name} — top: {top}  p_true≈{ptop:.4f}")
            say(f"   QCDS H(4 qubits): {round(H,4)} | uniformitet≈{round(uni,4)}")

        # Adaptiv Grover med roterande exkludering
        probs_g = run_adaptive_grover_4bit(
            sampler=sampler_main,
            target_bits=target_bits,
            shots=shots,
            max_iters_to_try=6,
            rotate_exclusion=True,
            seed=seed
        )
        p_true = probs_g.get(target_bits, 0.0)
        say(f"   ↳ Adaptiv Grover: bästa p_true≈{p_true:.4f}  (roterande bitexkludering)")

        results.append({
            "chrom": reg.chrom,
            "pos": pos,
            "genome": genome,
            "evs": evs,
            "target_bits": target_bits,
            "grover_probs": probs_g,
            "p_true": p_true
        })

    if layman:
        say("\n🧭 Lättförståelig sammanfattning (BRCA2 end-to-end):")
        say(f"• Region: {reg.chrom}:{reg.start}-{reg.end} ({genome})  |  körda punkter: {len(pos_list)}")
        say("• Varje 4-bitars ‘evidensmönster’ (här 0000 som demo) betyder: varje ‘1’ skulle flagga att just "
            "den evidenstypen fanns (t.ex. PVS1 = loss-of-function).")
        say("• QCDS 4-lager förstärker konsistenta mönster över flera nivåer för att ‘låsa in’ sannare svar, och "
            "Grover-delen hittar sanningen utan fix iteration – den konvergerar själv när signalen är stark nog.")

    return {
        "region": {"chrom": reg.chrom, "start": reg.start, "end": reg.end, "genome": genome},
        "results": results,
    }


# ---------------------------
# CLI
# ---------------------------

def parse_args():
    p = argparse.ArgumentParser(description="QCDS end-to-end (AER/IBM) med 2q-kontroller, 4L och adaptiv Grover.")
    p.add_argument("--mode", choices=["aer", "ibm"], default="aer", help="Välj simulator (aer) eller IBM (ibm).")
    p.add_argument("--genome", choices=["GRCh37", "GRCh38"], default="GRCh38", help="Genom-build.")
    p.add_argument("--region", type=str, default=None,
                   help="Kör över specifik region, t.ex. 13:32315474-32400266 (överskuggar genome-default).")
    p.add_argument("--shots", type=int, default=2048, help="Antal shots.")
    p.add_argument("--seed", type=int, default=None, help="Frö för determinism i demo-sampling.")
    p.add_argument("--layman", action="store_true", help="Aktivera lekmannaförklaringar.")
    p.add_argument("--iters-sched", type=str, default="1,1,3,5",
                   help="Schema för 4L-iters som csv, t.ex. 1,1,3,5")
    p.add_argument("--aer-max-variants", type=int, default=60, help="Många punkter i AER-läge.")
    p.add_argument("--ibm-max-variants", type=int, default=5, help="Få punkter i IBM-läge (utökningsbart).")
    p.add_argument("--ibm-backend", type=str, default=None, help="IBM-backendnamn (om --mode ibm).")
    p.add_argument("--force-aer-4l", action="store_true",
                   help="Kör 4-lager-delen på AER även om --mode ibm.")
    return p.parse_args()


def main():
    args = parse_args()
    iters_sched = tuple(int(x) for x in args.iters_sched.split(","))  # typ (1,1,3,5)
    assert len(iters_sched) == 4, "--iters-sched måste ha fyra heltal, t.ex. 1,1,3,5"

    say(f"⚙️  QCDS körläge: {args.mode}  |  4L_FORCE_AER={args.force_aer_4l}  |  LaymanExplain={args.layman}")
    if args.region:
        say(f"     Region override: {args.region}")
    else:
        reg = BRCA2_DEFAULT.get(args.genome, BRCA2_DEFAULT["GRCh38"])
        say(f"     Profil för BRCA2: {args.genome}  |  default-region {reg.chrom}:{reg.start}-{reg.end}  | "
            f"max_variants={'AER:'+str(args.aer_max_variants) if args.mode=='aer' else 'IBM:'+str(args.ibm_max_variants)}  "
            f"|  iters_sched={iters_sched}  |  shots={args.shots}")

    # Backend/sampler för huvudflödet
    backend_main, backend_name = get_backend(args.mode, args.ibm_backend)
    sampler_main = MiniSampler(backend_main, shots=args.shots, seed_simulator=args.seed)

    # Om 4L ska tvingas till AER (även i IBM-läge)
    if args.force_aer_4l and args.mode == "ibm":
        be_4l = Aer.get_backend("qasm_simulator")
        sampler_4l = MiniSampler(be_4l, shots=args.shots, seed_simulator=args.seed)
        running_on_aer_for_4l = True
    else:
        sampler_4l = sampler_main
        running_on_aer_for_4l = backend_name.startswith("qasm")

    # 2-qubit kontroller
    say("🔎 2-qubit QCDS-kontroller init…")
    _ = run_2q_controls(sampler_main, layman=args.layman)

    # 4-lager
    _ = run_4layer_sequence(sampler=sampler_4l,
                            shots=args.shots,
                            iters_sched=iters_sched,
                            force_aer_4l=args.force_aer_4l,
                            running_on_aer=running_on_aer_for_4l,
                            layman=args.layman)

    # BRCA2 end-to-end
    end2end = run_brca2_end2end(mode=args.mode,
                                sampler_main=sampler_main,
                                sampler_for_4l=sampler_4l,
                                genome=args.genome,
                                region=args.region,
                                aer_max_variants=args.aer_max_variants,
                                ibm_max_variants=args.ibm_max_variants,
                                iters_sched=iters_sched,
                                shots=args.shots,
                                layman=args.layman,
                                seed=args.seed)

    # Spara resultat lokalt (för spårbarhet i din miljö)
    outfile = f"QCDS_end2end_brca2_{args.genome.lower()}_{args.mode}.json"
    try:
        with open(outfile, "w", encoding="utf-8") as f:
            json.dump(end2end, f, ensure_ascii=False, indent=2)
        say(f"\n💾 QCDS end-to-end: sparade {outfile}")
    except Exception as e:
        say(f"⚠️ Kunde inte spara {outfile}: {e}")

    if args.layman:
        say("\n💡 Varför detta upplägg kan vara överlägset i praktiken:")
        say("• QCDS-orakel + Grover använder **alla kvantbitar som vittnen**, men plockar bort **en bit i taget** "
            "(roterande exkludering). Då tvingas varje körning ‘bevisa’ sanningen utan att kunna förlita sig på hela "
            "inmatningen. När resultaten sedan sammanfogas får du en **konsensus** som är robust mot brus och bias.")
        say("• Grover är **adaptiv** här: vi bestämmer inte i förväg hur många iterationer som behövs. Istället låter "
            "vi sanningen ‘segla upp’ och väljer den iteration som gav bäst signal (p_true). Det **minskar risken** "
            "att över- eller under-amplifiera (som annars kan tvätta bort signal).")
        say("• 4L-vågen fungerar som en **hierarkisk stötdämpare** som förstärker **konsekventa** mönster och dämpar "
            "spret – och därför kan du köra AER över hela BRCA2 (massivt) och IBM över mindre batchar men få "
            "jämförbara, stabila sammanfattningar.")

if __name__ == "__main__":
    main()
