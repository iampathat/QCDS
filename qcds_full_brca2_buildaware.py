#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
qcds_full_brca2_buildaware.py
--------------------------------
En fristående QCDS-körning för BRCA2 med:
- Stöd för både GRCh37 och GRCh38 (eller "auto" via VCF/koordinater)
- AER-masskörning över hela BRCA2-regionen
- IBM-bakände som kan slås på — samma variantpolicy men med rimligt tak
- En master-flagga som tvingar *emulering* (AER) för alla kvantdelar
- Spara JSON-rapporter (4L wave + end-to-end)
- Massor av lekmannaförklaringar i text
- Hierarkiskt 4 för cancer-caset med brus-simulering
"""

from __future__ import annotations
import os, sys, json, math, gzip, io, random, time
from dataclasses import dataclass, asdict
from typing import Dict, List, Tuple, Optional
from collections import Counter
from qiskit_aer.noise import NoiseModel, depolarizing_error

# ---------- Konfiguration via miljövariabler ----------

def env_bool(name: str, default: bool=False) -> bool:
    v = os.getenv(name, "")
    if v == "": return default
    return v.strip().lower() in ("1","true","yes","y","on")

def env_int(name: str, default: int) -> int:
    v = os.getenv(name, "")
    if v == "": return default
    try:
        return int(v)
    except:
        return default

def env_float(name: str, default: float) -> float:
    v = os.getenv(name, "")
    if v == "": return default
    try:
        return float(v)
    except:
        return default

def env_list_int(name: str, default: List[int]) -> List[int]:
    v = os.getenv(name, "")
    if v == "": return default
    try:
        parts = [int(x.strip()) for x in v.split(",") if x.strip()!=""]
        return parts if parts else default
    except:
        return default

QCDS_MODE           = os.getenv("QCDS_MODE", "aer").strip().lower()   # 'aer' eller 'ibm'
QCDS_FORCE_AER      = env_bool("QCDS_FORCE_AER", False)               # Master-switch: kör ALLT i AER
QCDS_ENABLE         = env_bool("QCDS_ENABLE", True)                   # Aktivera QCDS-utskrifter
QCDS_BUILD          = os.getenv("QCDS_BUILD", "auto").strip().upper() # 'GRCH37' | 'GRCH38' | 'AUTO'
QCDS_VCF            = os.getenv("QCDS_VCF", "").strip()               # valfri VCF/VCF.GZ
QCDS_LAYMAN         = env_bool("QCDS_LAYMAN", True)                   # lägg med lekmannatexter
MAX_VARIANTS_AER    = env_int("QCDS_MAX_VARIANTS_AER", 200)           # hur många punkter AER kan testa
MAX_VARIANTS_IBM    = env_int("QCDS_MAX_VARIANTS_IBM", 20)            # tak för IBM
REPEATS             = env_int("QCDS_REPEATS", 3)                      # hur många upprepningar i 2Q/konvergens
ITERS_SCHED         = env_list_int("QCDS_ITERS_SCHED", [1,1,3,5])     # lager-iter: L1,L2,L3,L4
SHOTS               = env_int("QCDS_SHOTS", 2048)                     # shots per körning
IBM_BACKEND         = os.getenv("QCDS_IBM_BACKEND", "ibm_brisbane").strip()
PRINT_LIMIT         = env_int("QCDS_PRINT_LIMIT", 60)                 # hur många variantloggar som skrivs

# Lägg till brus-nivå (0.0 = ingen brus, 0.01 = 1% depolarisering)
NOISE_LEVEL = env_float("QCDS_NOISE_LEVEL", 0.01)  # Standard 1% brus för Aer-simulering

# ---------- Build & BRCA2 region ----------

BRCA2_REGIONS = {
    "GRCH38": {"contig": "13",    "start": 32315086, "end": 32400268, "alt_contig": "chr13"},
    "GRCH37": {"contig": "13",    "start": 32889611, "end": 32973809, "alt_contig": "chr13"},
}

def pick_brca2_region(build: str) -> Tuple[str,int,int,str]:
    b = build.upper()
    if b not in BRCA2_REGIONS:
        raise ValueError(f"Okänd build: {build}. Tillåtna: GRCh37/GRCh38/auto")
    info = BRCA2_REGIONS[b]
    return info["contig"], info["start"], info["end"], info["alt_contig"]

def guess_build_from_positions(pos_values: List[int]) -> str:
    if not pos_values:
        return "GRCh38" # default
    avg = sum(pos_values)/len(pos_values)
    return "GRCh38" if abs(avg - 32350000) < abs(avg - 32930000) else "GRCh37"

def detect_build_from_vcf_header(vcf_path: str) -> Optional[str]:
    try:
        opener = gzip.open if vcf_path.endswith(".gz") else open
        with opener(vcf_path, "rt", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if not line.startswith("##"):
                    break
                low = line.lower()
                if "grch38" in low: return "GRCh38"
                if "grch37" in low or "hg19" in low: return "GRCh37"
        return None
    except Exception:
        return None

def open_vcf_lines(vcf_path: str):
    opener = gzip.open if vcf_path.endswith(".gz") else open
    return opener(vcf_path, "rt", encoding="utf-8", errors="ignore")

def load_variants_from_vcf(vcf_path: str, contig: str, start: int, end: int, alt_contig: str, max_keep: int) -> List[Tuple[str,int,str,str]]:
    hits: List[Tuple[str,int,str,str]] = []
    wanted = {contig, alt_contig}
    try:
        with open_vcf_lines(vcf_path) as f:
            for line in f:
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 5:
                    continue
                c, p, _, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
                if c not in wanted:
                    continue
                try:
                    pos = int(p)
                except:
                    continue
                if start <= pos <= end:
                    alts = alt.split(",")
                    alt0 = alts[0]
                    hits.append((c, pos, ref, alt0))
                    if len(hits) >= max_keep:
                        break
    except FileNotFoundError:
        print(f"[!] VCF saknas: {vcf_path}")
    except Exception as e:
        print(f"[!] VCF-läsfel: {e}")
    return hits

def synthesize_positions(contig: str, start: int, end: int, step: int, max_keep: int) -> List[Tuple[str,int,str,str]]:
    res: List[Tuple[str,int,str,str]] = []
    length = end - start + 1
    n = max(1, min(max_keep, length // step))
    if n <= 0: n = 1
    stride = max(1, length // n)
    pos = start
    while pos <= end and len(res) < max_keep:
        res.append((contig, pos, "A", "T"))
        pos += stride
    return res

# ---------- QCDS mätningar (2 qubits + 4-lager på 4 qubits) ----------

def shannon_H(prob: Dict[str,float]) -> float:
    eps = 1e-12
    return -sum(p*math.log2(p+eps) for p in prob.values() if p > 0)

def expect_zz_2q(prob: Dict[str,float]) -> float:
    return (prob.get("00",0)+prob.get("11",0)) - (prob.get("01",0)+prob.get("10",0))

def bias_bit(prob: Dict[str,float], bit_index: int) -> float:
    s1 = 0.0
    for s,p in prob.items():
        if len(s) <= bit_index:
            continue
        if s[-1-bit_index] == "1":
            s1 += p
    return s1 - (1.0 - s1)

def top_state(prob: Dict[str,float]) -> Tuple[str,float]:
    if not prob: return ("", 0.0)
    k = max(prob.items(), key=lambda kv: kv[1])
    return k[0], k[1]

# --- Qiskit set-up ---
def get_aer_backend():
    try:
        from qiskit_aer import Aer
        return Aer.get_backend("qasm_simulator")
    except Exception:
        try:
            from qiskit import BasicAer
            return BasicAer.get_backend("qasm_simulator")
        except Exception as e:
            raise RuntimeError("Ingen AER/BasicAer tillgänglig för simulering.") from e

def run_circuits_aer(circuits, shots: int, noise_model=None) -> List[Dict[str,int]]:
    from qiskit import transpile
    backend = get_aer_backend()
    if noise_model:
        backend.set_options(noise_model=noise_model)
    tc = transpile(circuits, backend=backend, optimization_level=1)
    job = backend.run(tc, shots=shots)
    result = job.result()
    counts_list = []
    for i in range(len(circuits)):
        counts = result.get_counts(i)
        counts_norm = {}
        for k,v in counts.items():
            s = k.replace(" ", "")
            counts_norm[s] = v
        counts_list.append(counts_norm)
    return counts_list

def run_circuits_ibm(circuits, shots: int, backend_name: str) -> List[Dict[str,int]]:
    if QCDS_FORCE_AER:
        noise_model = NoiseModel()
        error = depolarizing_error(NOISE_LEVEL, 1)
        noise_model.add_all_qubit_quantum_error(error, ['h', 'cx', 'measure'])
        return run_circuits_aer(circuits, shots, noise_model)
    try:
        from qiskit_ibm_runtime import QiskitRuntimeService, Sampler, Options
    except Exception as e:
        print("[!] qiskit-ibm-runtime saknas eller fel. Faller tillbaka till AER med brus.")
        noise_model = NoiseModel()
        error = depolarizing_error(NOISE_LEVEL, 1)
        noise_model.add_all_qubit_quantum_error(error, ['h', 'cx', 'measure'])
        return run_circuits_aer(circuits, shots, noise_model)
    try:
        service = QiskitRuntimeService()
        backend = service.backend(backend_name)
        options = Options(resilience_level=0)
        sampler = Sampler(session=None, options=options)
        job = sampler.run(circuits=circuits, shots=shots, backend=backend)
        result = job.result()
        counts_list = []
        for quasi in result.quasi_dists:
            counts = {k: int(round(v*shots)) for k,v in quasi.items()}
            counts_list.append(counts)
        return counts_list
    except Exception as e:
        print(f"[!] IBM-körning misslyckades ({backend_name}): {e}. Faller tillbaka till AER med brus.")
        noise_model = NoiseModel()
        error = depolarizing_error(NOISE_LEVEL, 1)
        noise_model.add_all_qubit_quantum_error(error, ['h', 'cx', 'measure'])
        return run_circuits_aer(circuits, shots, noise_model)

# --- Konstruera experiment ---
from qiskit import QuantumCircuit

def two_qubit_tg_circuit() -> QuantumCircuit:
    qc = QuantumCircuit(2,2)
    qc.h(0); qc.cx(0,1)
    qc.measure([0,1],[0,1])
    return qc

def two_qubit_ac_circuit() -> QuantumCircuit:
    qc = QuantumCircuit(2,2)
    qc.h([0,1])
    qc.cx(0,1); qc.rz(0.2, 1); qc.cx(0,1)
    qc.h([0,1])
    qc.measure([0,1],[0,1])
    return qc

def two_qubit_ce_circuit() -> QuantumCircuit:
    qc = QuantumCircuit(2,2)
    qc.h(0); qc.cx(0,1)
    qc.cz(0,1)
    qc.measure([0,1],[0,1])
    return qc

def four_qubit_layered(prop_zero: Dict[int,str]) -> QuantumCircuit:
    qc = QuantumCircuit(4,4)
    qc.h([0,1,2,3])
    qc.cx(0,1); qc.cx(1,2); qc.cx(2,3)
    for q,s in prop_zero.items():
        if s == '0':
            qc.rz(0.35, q)
            qc.rz(0.1, q)
    qc.rx(0.15, [0,1,2,3])
    qc.barrier()
    qc.measure([0,1,2,3],[0,1,2,3])
    return qc

def counts_to_probs(counts: Dict[str,int], shots: int) -> Dict[str,float]:
    total = sum(counts.values()) if counts else 0
    if total <= 0: return {}
    return {k: v/total for k,v in counts.items()}

# ---------- Högre nivå: körningar och utskrifter ----------

def pretty_probs(order: List[str], probs: Dict[str,float]) -> Dict[str,float]:
    out = {}
    for s in order:
        if s in probs: out[s] = round(probs[s], 4)
    for s,p in sorted(probs.items(), key=lambda kv: kv[1], reverse=True):
        if s not in out and len(out) < 4:
            out[s] = round(p,4)
    return out

def layman_explain_2q(name: str, top: str, ptop: float, expzz: float) -> str:
    if not QCDS_LAYMAN: return ""
    lines = []
    if name == "TG":
        lines.append("— Truth Gradient (TG): mäter hur ofta bitarna är lika (00/11) jämfört med olika (01/10).")
        lines.append(f"  Resultat: högsta utfallet är {top} ({ptop:.1%}). Hög ⟨Z⊗Z⟩≈{expzz:.2f} ⇒ stark samvariation.")
        lines.append("  Lekman: “Tänk två lampor – de tänds/släcks oftare samtidigt än motsatt.”")
    elif name == "AC":
        lines.append("— Amplitude Confirmation (AC): vi ger en liten fas-knuff och läser av i vanlig mätning.")
        lines.append(f"  Resultat: högsta utfallet är {top} ({ptop:.1%}). ⟨Z⊗Z⟩≈{expzz:.2f} ⇒ nära neutral/balanserad korrelation.")
        lines.append("  Lekman: “Vi nyper systemet lite – blir utfallet fortsatt jämnt är fasen balanserad.”")
    elif name == "CE":
        lines.append("— Conditional Enhancement (CE): förstärker kombinationer som uppfyller ett villkor (t.ex. att bitarna är lika).")
        lines.append(f"  Resultat: högsta utfallet är {top} ({ptop:.1%}). ⟨Z⊗Z⟩≈{expzz:.2f} ⇒ tydlig förstärkning av lika-utfall.")
        lines.append("  Lekman: “När villkoret stämmer blir just den kombinationen extra sann.”")
    return "\n".join(lines)

def layman_explain_4q(variant_tag: str, final_top: str) -> str:
    if not QCDS_LAYMAN: return ""
    return (
        f"• {variant_tag}: 4-bitars evidens-mönster = {final_top} "
        "[none] om alla bitar är 0.\n"
        "  Lekman: Varje ‘1’ betyder att just den evidenstypen flaggats (t.ex. “förlust av funktion”).\n"
        "  Metoden jobbar i fyra steg för att förstärka svar som är *konsekventa* mellan nivåerna—"
        "lite som att fråga samma sak på olika sätt och låta svaret “låsa sig” om det håller."
    )

def run_two_qubit_suite(mode: str, shots: int, repeats: int):
    circuits = [two_qubit_tg_circuit(), two_qubit_ac_circuit(), two_qubit_ce_circuit()]
    if mode == "ibm" and not QCDS_FORCE_AER:
        counts_all = run_circuits_ibm(circuits, shots, IBM_BACKEND)
        backend_label = IBM_BACKEND
    else:
        counts_all = run_circuits_aer(circuits, shots)
        backend_label = "qasm_simulator"
    names = ["TG", "AC", "CE"]
    for name,counts in zip(names, counts_all):
        probs = counts_to_probs(counts, shots)
        top, ptop = top_state(probs)
        expzz = expect_zz_2q(probs)
        H = shannon_H(probs)
        order = ["00","11","01","10"]
        shown = pretty_probs(order, probs)
        print(f"▶️ { {'TG':'Truth Gradient','AC':'Amplitude Confirmation','CE':'Conditional Enhancement'}[name] } @ {backend_label} — top: {top}")
        print(f"   probs: {shown}")
        print(f"   QCDS ⟨Z⊗Z⟩: {round(expzz,4)} | bias(q0): {round(bias_bit(probs,0),4)} | bias(q1): {round(bias_bit(probs,1),4)} | H: {round(H,4)}\n")
        if QCDS_LAYMAN:
            print("📘 Lättförståelig förklaring (2 qubits):")
            print(layman_explain_2q(name, top, ptop, expzz))
            print()

def run_qcds_convergence_pair(mode: str, shots: int, repeats: int):
    if mode == "ibm" and not QCDS_FORCE_AER:
        backend_label = IBM_BACKEND
    else:
        backend_label = "qasm_simulator"
    print(f"🔎 QCDS konvergens (mode={'ibm' if backend_label!= 'qasm_simulator' else 'aer'}) init… (repeats={repeats})")
    tg, ce, ac = two_qubit_tg_circuit(), two_qubit_ce_circuit(), two_qubit_ac_circuit()
    if backend_label == "qasm_simulator":
        counts = run_circuits_aer([tg,ce,ac], shots)
    else:
        counts = run_circuits_ibm([tg,ce,ac], shots, backend_label)
    names = ["TG (Z)", "CE (Z)", "AC (X)"]
    for name,c in zip(names, counts):
        probs = counts_to_probs(c, shots)
        top, ptop = top_state(probs)
        expzz = round(expect_zz_2q(probs),4)
        H = round(shannon_H(probs),4)
        order = ["00","11","01","10"]
        shown = pretty_probs(order, probs)
        print(f"▶️ {name} @ {backend_label} (1, 0) — top: {top}")
        print(f"   probs: {shown}")
        print(f"   QCDS ⟨Z⊗Z⟩: {expzz} | bias(q0): {round(bias_bit(probs,0),4)} | bias(q1): {round(bias_bit(probs,1),4)} | H: {H}")

def run_qcds_hier4(mode: str, init_cond: str="x and y", shots: int=2048) -> List[Dict]:
    if mode == "ibm" and not QCDS_FORCE_AER:
        backend_label = IBM_BACKEND
    else:
        backend_label = "qasm_simulator"
    print(f"▶️ QCDS Hierarchical 4 (mode={'ibm' if backend_label!= 'qasm_simulator' else 'aer'}) med initial condition: {init_cond}")
    out = []
    for qi in range(1,5):
        circ = two_qubit_tg_circuit() if qi%2==1 else two_qubit_ce_circuit()
        counts = (run_circuits_ibm([circ], shots, backend_label) if backend_label!= "qasm_simulator"
                  else run_circuits_aer([circ], shots))[0]
        probs = counts_to_probs(counts, shots)
        top, _ = top_state(probs)
        expzz = round(expect_zz_2q(probs),4)
        H = round(shannon_H(probs),4)
        order = ["00","11","01","10"]
        shown = pretty_probs(order, probs)
        print(f"Q{qi}: Top state: {top}, probs: {shown}")
        print(f"   QCDS ⟨Z⊗Z⟩: {expzz} | bias(q0): {round(bias_bit(probs,0),4)} | bias(q1): {round(bias_bit(probs,1),4)} | H: {H}")
        out.append({"Qi":qi,"top":top,"probs":probs,"H":H,"expZZ":expzz})
    return out

# ---------- 4-qubit "4L" körning per BRCA2-variant ----------

@dataclass
class FourLStep:
    layer: int
    iters: int
    prop: Dict[int,str]
    backend: str
    top_state: str
    p_true: float
    H: float

@dataclass
class VariantResult:
    contig: str
    pos: int
    ref: str
    alt: str
    target_bits: str   # slutligt "evidensmönster" (stub)
    steps: List[FourLStep]
    H4: float
    uniformity: float
    final_top: str

def run_4L_for_variant(variant: Tuple[str,int,str,str], mode: str, shots: int, iters_sched: List[int]) -> VariantResult:
    contig, pos, ref, alt = variant
    target_bits = "0000" if ref == alt else "0001"
    props = [
        {3:'0'},
        {3:'0', 2:'0'},
        {3:'0', 2:'0', 1:'0'},
        {3:'0', 2:'0', 1:'0', 0: target_bits[0]}
    ]
    backend_label = (IBM_BACKEND if (mode=="ibm" and not QCDS_FORCE_AER) else "qasm_simulator")
    steps: List[FourLStep] = []

    H4 = 0.0; final_top = "0000"
    for L, iters, prop in zip(range(1,5), iters_sched, props):
        circ = four_qubit_layered(prop)
        if backend_label == "qasm_simulator":
            counts = run_circuits_aer([circ], shots)[0]
        else:
            counts = run_circuits_ibm([circ], shots, backend_label)[0]
        probs = counts_to_probs(counts, shots)
        for s in [f"{i:04b}" for i in range(16)]:
            probs.setdefault(s, 0.0)
        H_now = shannon_H(probs)
        top, ptop = top_state(probs)
        p_true = probs.get(target_bits, 0.0) if target_bits in probs else 0.0
        if L == 1:
            H4 = shannon_H({s:1/16 for s in [f"{i:04b}" for i in range(16)]})
        steps.append(FourLStep(L, iters, prop, backend_label, top, p_true, H_now))
        final_top = top

    uniformity = min(1.0, shannon_H({s:1/16 for s in [f"{i:04b}" for i in range(16)]}) / 4.0)
    return VariantResult(contig, pos, ref, alt, target_bits, steps, H4, uniformity, final_top)

# ---------- Hämta/Skapa BRCA2-variantlista ----------

def build_and_variants(vcf_path: str, build_opt: str, max_keep: int) -> Tuple[str, List[Tuple[str,int,str,str]], Tuple[str,int,int,str]]:
    if build_opt != "AUTO":
        build = "GRCh38" if build_opt == "GRCH38" else "GRCh37"
    else:
        build = detect_build_from_vcf_header(vcf_path) if vcf_path else None
        if build is None and vcf_path:
            try:
                with open_vcf_lines(vcf_path) as f:
                    pos_samp = []
                    for line in f:
                        if line.startswith("#"): continue
                        parts = line.split("\t")
                        if len(parts) < 2: continue
                        try: pos_samp.append(int(parts[1]))
                        except: continue
                        if len(pos_samp) >= 50: break
                build = guess_build_from_positions(pos_samp)
            except: build = "GRCh38"
        if build is None: build = "GRCh38"
    contig, start, end, alt_contig = pick_brca2_region(build)
    if vcf_path and os.path.exists(vcf_path):
        vars_list = load_variants_from_vcf(vcf_path, contig, start, end, alt_contig, max_keep)
        if not vars_list:
            step = max(10, (end-start)//max_keep if max_keep>0 else 100)
            vars_list = synthesize_positions(contig, start, end, step, max_keep)
    else:
        step = max(10, (end-start)//max_keep if max_keep>0 else 100)
        vars_list = synthesize_positions(contig, start, end, step, max_keep)
    return build, vars_list, (contig, start, end, alt_contig)

# ---------- Huvudkörning ----------

def main():
    mode = QCDS_MODE
    if QCDS_FORCE_AER:
        mode_label = "aer"
    else:
        mode_label = mode

    fourL_force_aer = QCDS_FORCE_AER

    print(f"⚙️  QCDS körläge: {mode_label}  |  4L_FORCE_AER={fourL_force_aer}  |  LaymanExplain={QCDS_LAYMAN}")
    run_two_qubit_suite(mode, SHOTS, REPEATS)
    run_qcds_convergence_pair(mode, SHOTS, REPEATS)
    run_qcds_hier4(mode, init_cond="x and y", shots=SHOTS)

    max_keep = MAX_VARIANTS_AER if (mode_label=="aer") else MAX_VARIANTS_IBM
    build, variants, (contig, start, end, alt_contig) = build_and_variants(QCDS_VCF if QCDS_VCF else "", QCDS_BUILD, max_keep)
    print(f"✅ BRCA2-region ({build}): {contig}:{start}-{end}  |  profil: {mode_label.upper()}  |  max_variants={max_keep}  |  repeats={REPEATS}  |  iters_sched={tuple(ITERS_SCHED)}  |  shots={SHOTS}")

    if QCDS_VCF:
        base = os.path.basename(QCDS_VCF)
        print(f"📦 Dataset: {base}")
    else:
        print("📦 Dataset: [syntetisk positionslista — ingen VCF inskickad]")

    wave_steps_all = []
    end2end = {
        "dataset": os.path.basename(QCDS_VCF) if QCDS_VCF else "synthetic",
        "build": build,
        "region": f"{contig}:{start}-{end}",
        "ran_variants": 0,
        "variants": []
    }

    shown = 0
    for (c, p, ref, alt) in variants:
        tag = f"{c}:{p}"
        print(f"🎯 Variant {tag} → evs=[]")  # Ta bort target_bits från utskrift
        vr = run_4L_for_variant((c,p,ref,alt), mode, SHOTS, ITERS_SCHED)
        for st in vr.steps:
            if shown < PRINT_LIMIT:
                print(f"▶️ QCDS_4L L{st.layer} (iters={st.iters}, prop={st.prop}) — backend: {st.backend} — top: {st.top_state}  p_true≈{st.p_true:.4f}")
            shown += 1
        if shown < PRINT_LIMIT:
            print(f"   QCDS H(4 qubits): {vr.steps[-1].H:.4f} | uniformitet≈{vr.uniformity:.4f}")
        wave_steps_all.append({"variant": tag, "steps": [asdict(s) for s in vr.steps]})
        end2end["ran_variants"] += 1
        end2end["variants"].append({
            "variant": tag,
            "target_bits": vr.target_bits,
            "final_top": vr.final_top,
            "H4": vr.steps[-1].H
        })

    wave_path = "QCDS_4L_wave_report.json"
    end_path = "QCDS_end2end_brca2_quantum.json"
    with open(wave_path, "w", encoding="utf-8") as f:
        json.dump(wave_steps_all, f, ensure_ascii=False, indent=2)
    with open(end_path, "w", encoding="utf-8") as f:
        json.dump(end2end, f, ensure_ascii=False, indent=2)
    print(f"💾 QCDS-4L: sparade {wave_path}")
    print(f"💾 QCDS end-to-end: sparade {end_path}")

    if QCDS_LAYMAN:
        print("\n🧭 Lättförståelig sammanfattning (BRCA2 end-to-end):")
        print(f"• Region: {contig}:{start}-{end} ({build})  |  varianter körda: {end2end['ran_variants']}")
        if end2end["ran_variants"] == 0:
            print("• Inga BRCA2-varianter kördes i denna omgång.")
        else:
            for v in end2end["variants"][:10]:
                print(f"• {v['variant']}: 4-bitars evidens-mönster = {v['final_top']} (just nu en prototyp).")
            if end2end["ran_variants"] > 10:
                print(f"  …och {end2end['ran_variants']-10} till (se JSON-rapporten).")
            print("  Lekman: Varje bit motsvarar en grov evidenstyp. 0 betyder “inget påträffat” i denna prototyp.")
            print("  För klinisk tolkning krävs riktiga annots (ex. konsekvenser i BRCA2-proteinet).")

    print(f"\n▶️ QCDS Hierarchical 4 Cancer Case (mode={mode_label}) med initial condition: tp53 and brca1")
    cancer_result = run_qcds_hierarchical_4_cancer(mode_label, initial_condition="tp53 and brca1", shots=SHOTS, enable=QCDS_ENABLE)
    final_risk_state = cancer_result["final_risk_state"]
    final_probs = cancer_result["final_probs"]
    print(f"🏁 Final risk-state: {final_risk_state}, probs: {pretty_probs(['00', '01', '10', '11'], final_probs)}")
    if QCDS_ENABLE:
        qcds_print("Final Cancer Risk", final_probs)

# ---------- Ny funktion: Hierarkiskt 4 för cancer-caset med brus ----------

def synthesize_cancer_condition(previous_top: str) -> str:
    if previous_top == '11':
        return "tp53 or brca1"
    elif previous_top == '00':
        return "tp53 and brca1"
    else:
        return "tp53 and brca1"

def build_cancer_oracle(condition: str, num_qubits: int = 2) -> QuantumCircuit:
    qc = QuantumCircuit(num_qubits, num_qubits)
    qc.h(range(num_qubits))
    if condition.lower() == "tp53 and brca1":
        qc.cz(0, 1)  # Markera |11> som högrisk
    elif condition.lower() == "tp53 or brca1":
        qc.x([0, 1])
        qc.cz(0, 1)
        qc.x([0, 1])
    qc.measure(range(num_qubits), range(num_qubits))
    return qc

def run_qcds_hierarchical_4_cancer(mode: str, initial_condition: str, shots: int, enable: bool) -> Dict:
    backend_label = IBM_BACKEND if (mode == "ibm" and not QCDS_FORCE_AER) else "qasm_simulator"
    print(f"▶️ QCDS Hierarchical 4 Cancer Case (mode={mode}) med initial condition: {initial_condition}")

    # Q1: Initial Synthesizer
    q1_circuit = build_cancer_oracle(initial_condition)
    q1_counts = (run_circuits_ibm([q1_circuit], shots, backend_label) if backend_label != "qasm_simulator"
                 else run_circuits_aer([q1_circuit], shots))[0]
    q1_probs = counts_to_probs(q1_counts, shots)
    q1_top, _ = top_state(q1_probs)
    print(f"Q1: Top risk-state: {q1_top}, probs: {pretty_probs(['00', '01', '10', '11'], q1_probs)}")
    if enable:
        qcds_print("Q1 Cancer", q1_probs)

    # Q2: Inferential Re-entry
    q2_condition = synthesize_cancer_condition(q1_top)
    q2_circuit = build_cancer_oracle(q2_condition)
    q2_counts = (run_circuits_ibm([q2_circuit], shots, backend_label) if backend_label != "qasm_simulator"
                 else run_circuits_aer([q2_circuit], shots))[0]
    q2_probs = counts_to_probs(q2_counts, shots)
    q2_top, _ = top_state(q2_probs)
    print(f"Q2 (C2: {q2_condition}): Top risk-state: {q2_top}, probs: {pretty_probs(['00', '01', '10', '11'], q2_probs)}")
    if enable:
        qcds_print("Q2 Cancer", q2_probs)

    # Q3: Cascading Extension
    q3_condition = synthesize_cancer_condition(q2_top)
    q3_circuit = build_cancer_oracle(q3_condition)
    q3_counts = (run_circuits_ibm([q3_circuit], shots, backend_label) if backend_label != "qasm_simulator"
                 else run_circuits_aer([q3_circuit], shots))[0]
    q3_probs = counts_to_probs(q3_counts, shots)
    q3_top, _ = top_state(q3_probs)
    print(f"Q3 (C3: {q3_condition}): Top risk-state: {q3_top}, probs: {pretty_probs(['00', '01', '10', '11'], q3_probs)}")
    if enable:
        qcds_print("Q3 Cancer", q3_probs)

    # Q4: Syntactic Resonance
    q4_condition = synthesize_cancer_condition(q3_top)
    q4_circuit = build_cancer_oracle(q4_condition)
    q4_counts = (run_circuits_ibm([q4_circuit], shots, backend_label) if backend_label != "qasm_simulator"
                 else run_circuits_aer([q4_circuit], shots))[0]
    q4_probs = counts_to_probs(q4_counts, shots)
    q4_top, _ = top_state(q4_probs)
    print(f"Q4 (C4: {q4_condition}): Final risk-state: {q4_top}, probs: {pretty_probs(['00', '01', '10', '11'], q4_probs)}")
    if enable:
        qcds_print("Q4 Cancer", q4_probs)

    return {"final_risk_state": q4_top, "final_probs": q4_probs}

# ---------- Huvudkörning ----------

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n[Avbrutet]")
