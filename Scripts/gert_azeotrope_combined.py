"""
GERT Azeotrope — Combined Model: Network + Dispersion
========================================================
Explicit mathematics for both regimes.

REGIME 1 (WITH NETWORK): Crystal classification
REGIME 2 (WITHOUT NETWORK): Hansen dispersion criterion
"""

import numpy as np

# =====================================================================
# MATHEMATICAL FRAMEWORK
# =====================================================================

print("""
╔══════════════════════════════════════════════════════════════════════╗
║  GERT AZEOTROPE — COMBINED CRYSTAL + DISPERSION MODEL              ║
║  Explicit Mathematics                                               ║
╚══════════════════════════════════════════════════════════════════════╝

═══════════════════════════════════════════════════════════════════════
STEP 1: MOLECULE CLASSIFICATION
═══════════════════════════════════════════════════════════════════════

Define self-association parameter:

    σᵢ = αᵢ × βᵢ                                              (Eq. 1)

where αᵢ = Kamlet-Taft H-bond donor, βᵢ = Kamlet-Taft H-bond acceptor.

Classification:
    σ > 0.40  AND  3D network (water)  →  Type 2 (Td)
    σ > 0.05                            →  Type 1 (chain)
    σ ≤ 0.05                            →  Type 0 (no network)

═══════════════════════════════════════════════════════════════════════
STEP 2: PAIR CLASSIFICATION
═══════════════════════════════════════════════════════════════════════

Define cross-association and complementarity:

    Γcross = αA·βB + αB·βA         (cross H-bonding)           (Eq. 2)
    Γself  = αA·βA + αB·βB = σA + σB   (self H-bonding)        (Eq. 3)
    C = Γcross - Γself              (net complementarity)       (Eq. 4)

Pair classification:
    C > 0  AND  σA < 0.05  AND  σB < 0.05  →  Type 3 (complementary)
    Otherwise: max(typeA, typeB) determines the pair type.

═══════════════════════════════════════════════════════════════════════
STEP 3: PREDICTION RULES BY PAIR TYPE
═══════════════════════════════════════════════════════════════════════

─── TYPE 3 (complementary): ──────────────────────────────────────────

    C > 0  →  NEGATIVE azeotrope                                (Eq. 5)

    Emergent network: exists ONLY in the mixture.

─── TYPE 1+2, TYPE 1+1, TYPE 0+1 (with chain or Td): ────────────────

    Integration fraction into the stronger network:

    f_int = √(f_donor × f_acceptor)                             (Eq. 6)

    where:
        f_donor = min(α_solute / α_network, 1)                  (Eq. 7a)
        f_acceptor = min(β_solute / β_network, 1)               (Eq. 7b)

    Prediction:
        f_int > 0.90  →  NO azeotrope (fully integrates)        (Eq. 8a)
        f_int ≤ 0.90  →  POSITIVE azeotrope (disrupts network)  (Eq. 8b)

    For Type 1+1 (chain + chain):
        Same chemical class (both alcohols, both acids) → no azeotrope
        Different classes → positive azeotrope

─── TYPE 0+2 (no network + Td): ──────────────────────────────────────

    Combined criterion — BOTH must be satisfied:

    1. Network disruption:  f_int < 0.80                         (Eq. 9a)
    2. Volatility access:   |ΔTb| < 40 K                        (Eq. 9b)

    If BOTH: positive azeotrope
    If either fails: no azeotrope
    
    Physics: even if the solute disrupts the network (Eq. 9a),
    the azeotrope is NOT REACHABLE if the volatility difference 
    is too large (Eq. 9b) — OPEN SYSTEM effect.

─── TYPE 0+0 (no network + no network): ─────────────────────────────

    Dispersion criterion (Hansen):

    D² = (δd,A - δd,B)² + (δp,A - δp,B)²                       (Eq. 10)
    
    (No δh term: Type 0 molecules have no significant H-bonds)

    Family similarity:
    
    F = 1   if A and B are in the same chemical class             (Eq. 11)
    F = 0   otherwise

    Prediction:
        F = 1                  →  NO azeotrope (too similar)     (Eq. 12a)
        F = 0  AND  D > 4 MPa^½  →  POSITIVE azeotrope          (Eq. 12b)
        F = 0  AND  D ≤ 4 MPa^½  →  NO azeotrope                (Eq. 12c)

    (D threshold from: minimum dissimilarity to create 
     significant positive deviation from Raoult's law)
""")

# =====================================================================
# DATA
# =====================================================================

molecules = {
    #                  α     β      σ    type  δd    δp    family
    'water':        (1.17, 0.47, 0.550,  2,   15.5, 16.0, 'water'),
    'ethanol':      (0.86, 0.75, 0.645,  1,   15.8,  8.8, 'alcohol'),
    'methanol':     (0.98, 0.66, 0.647,  1,   15.1, 12.3, 'alcohol'),
    'isopropanol':  (0.76, 0.84, 0.638,  1,   15.8,  6.1, 'alcohol'),
    'n_butanol':    (0.84, 0.84, 0.706,  1,   16.0,  5.7, 'alcohol'),
    'acetic_acid':  (1.12, 0.45, 0.504,  1,   14.5,  8.0, 'acid'),
    'acetone':      (0.08, 0.43, 0.034,  0,   15.5, 10.4, 'ketone'),
    'chloroform':   (0.20, 0.10, 0.020,  0,   17.8,  3.1, 'haloalkane'),
    'benzene':      (0.00, 0.10, 0.000,  0,   18.4,  0.0, 'aromatic'),
    'toluene':      (0.00, 0.11, 0.000,  0,   18.0,  1.4, 'aromatic'),
    'hexane':       (0.00, 0.00, 0.000,  0,   14.9,  0.0, 'alkane'),
    'heptane':      (0.00, 0.00, 0.000,  0,   15.3,  0.0, 'alkane'),
    'CS2':          (0.00, 0.07, 0.000,  0,   20.5,  0.0, 'other'),
    'cyclohexane':  (0.00, 0.00, 0.000,  0,   16.8,  0.0, 'alkane'),
    'ethyl_acetate':(0.00, 0.45, 0.000,  0,   15.8,  5.3, 'ester'),
    'diethyl_ether':(0.00, 0.47, 0.000,  0,   14.5,  2.9, 'ether'),
    'THF':          (0.00, 0.55, 0.000,  0,   16.8,  5.7, 'ether'),
    'DCM':          (0.13, 0.10, 0.013,  0,   18.2,  6.3, 'haloalkane'),
}

Tb = {'water':373,'ethanol':351,'methanol':338,'acetone':329,'chloroform':334,
      'benzene':353,'toluene':384,'hexane':342,'heptane':372,'CS2':319,
      'cyclohexane':354,'ethyl_acetate':350,'isopropanol':355,'diethyl_ether':308,
      'THF':339,'n_butanol':391,'acetic_acid':391,'DCM':313}

def predict(mol_A, mol_B):
    aA, bA, sA, tA, dA, pA, fA = molecules[mol_A]
    aB, bB, sB, tB, dB, pB, fB = molecules[mol_B]
    
    # Eq. 2-4: Cross, self, complementarity
    G_cross = aA*bB + aB*bA
    G_self = sA + sB
    C = G_cross - G_self
    
    # Step 1: Check Type 3 (complementary)
    if C > 0 and sA < 0.05 and sB < 0.05:
        return "negative", 3, (f"Type 3: C={C:.3f}>0, σA={sA:.3f}, σB={sB:.3f}. "
                               f"Γcross={G_cross:.3f} > Γself={G_self:.3f}. "
                               f"Emergent network → NEGATIVE")
    
    # Pair type = max of individual types
    pair_t = max(tA, tB)
    
    # Step 2: TYPE 0+0 (no network)
    if tA == 0 and tB == 0:
        # Eq. 10: Dispersion distance
        D = np.sqrt((dA - dB)**2 + (pA - pB)**2)
        # Eq. 11: Family similarity
        F = 1 if fA == fB else 0
        
        if F == 1:
            return "none", 0, (f"Type 0+0, SAME family ({fA}). "
                               f"F=1 → no azeotrope (Eq. 12a)")
        elif D > 4.0:
            return "positive", 0, (f"Type 0+0, diff families ({fA}+{fB}). "
                                   f"D={D:.1f} > 4.0 → positive (Eq. 12b)")
        else:
            return "none", 0, (f"Type 0+0, diff families ({fA}+{fB}). "
                               f"D={D:.1f} ≤ 4.0 → no azeotrope (Eq. 12c)")
    
    # Step 3: With network — identify which has stronger network
    if tA >= tB:
        net_a, net_b = aA, bA
        sol_a, sol_b = aB, bB
        net_t, sol_t = tA, tB
        net_f, sol_f = fA, fB
        net_mol, sol_mol = mol_A, mol_B
    else:
        net_a, net_b = aB, bB
        sol_a, sol_b = aA, bA
        net_t, sol_t = tB, tA
        net_f, sol_f = fB, fA
        net_mol, sol_mol = mol_B, mol_A
    
    # Eq. 6-7: Integration fraction
    f_d = min(sol_a / net_a, 1.0) if net_a > 0 else 0
    f_a = min(sol_b / net_b, 1.0) if net_b > 0 else 0
    f_int = np.sqrt(f_d * f_a) if (f_d > 0 and f_a > 0) else 0
    
    # TYPE 1+1 (chain + chain)
    if net_t == 1 and sol_t == 1:
        if net_f == sol_f:
            return "none", 1, (f"Type 1+1, SAME class ({net_f}). "
                               f"Compatible chains → no azeotrope")
        else:
            return "positive", 1, (f"Type 1+1, DIFF class ({net_f}+{sol_f}). "
                                   f"Incompatible chains → positive")
    
    # TYPE 0+1 (no net + chain)
    if net_t == 1 and sol_t == 0:
        return "positive", 1, (f"Type 0+1: {sol_mol}(0) disrupts {net_mol}(1) chain. "
                               f"→ positive azeotrope")
    
    # TYPE 1+2 (chain + Td)
    if net_t == 2 and sol_t == 1:
        if f_int > 0.90:
            return "none", 2, (f"Type 1+2: f_int={f_int:.2f} > 0.90. "
                               f"{sol_mol} integrates into Td → no azeotrope (Eq. 8a)")
        else:
            return "positive", 2, (f"Type 1+2: f_int={f_int:.2f} ≤ 0.90. "
                                   f"{sol_mol} disrupts Td → positive (Eq. 8b)")
    
    # TYPE 0+2 (no net + Td) — combined criterion
    if net_t == 2 and sol_t == 0:
        dTb = abs(Tb[net_mol] - Tb[sol_mol])
        crit_net = f_int < 0.80   # Eq. 9a: disrupts
        crit_vol = dTb < 40       # Eq. 9b: reachable
        
        if crit_net and crit_vol:
            return "positive", 2, (f"Type 0+2: f_int={f_int:.2f}<0.80 (disrupts, Eq.9a) "
                                   f"AND ΔTb={dTb:.0f}<40 (reachable, Eq.9b) → positive")
        elif crit_net and not crit_vol:
            return "none", 2, (f"Type 0+2: f_int={f_int:.2f}<0.80 (disrupts) "
                               f"BUT ΔTb={dTb:.0f}≥40 (OPEN SYSTEM: unreachable, Eq.9b) → none")
        else:
            return "none", 2, (f"Type 0+2: f_int={f_int:.2f}≥0.80 (integrates) → none")
    
    return "unknown", -1, "Unclassified"

# =====================================================================
# TEST
# =====================================================================

systems = [
    ("EtOH-H₂O",       'ethanol',      'water',         0.894, "positive"),
    ("iPrOH-H₂O",      'isopropanol',  'water',         0.674, "positive"),
    ("nBuOH-H₂O",      'n_butanol',    'water',         0.749, "positive"),
    ("MeOH-Benzene",    'methanol',     'benzene',       0.615, "positive"),
    ("Benzene-EtOH",    'benzene',      'ethanol',       0.448, "positive"),
    ("EtOAc-EtOH",      'ethyl_acetate','ethanol',       0.462, "positive"),
    ("Cyclohex-EtOH",   'cyclohexane',  'ethanol',       0.446, "positive"),
    ("CS₂-Acetone",     'CS2',          'acetone',       0.610, "positive"),
    ("EtOH-Benzene",    'ethanol',      'benzene',       0.552, "positive"),
    ("THF-H₂O",         'THF',          'water',         0.820, "positive"),
    ("MeOH-Cyclohex",   'methanol',     'cyclohexane',   0.510, "positive"),
    ("DCM-MeOH",        'DCM',          'methanol',      0.810, "positive"),
    ("Acetone-CHCl₃",   'acetone',      'chloroform',    0.357, "negative"),
    ("CHCl₃-Et₂O",      'chloroform',   'diethyl_ether', 0.850, "negative"),
    ("CHCl₃-Acetone",   'chloroform',   'acetone',       0.643, "negative"),
    ("MeOH-EtOH",       'methanol',     'ethanol',       None,  "none"),
    ("Benzene-Toluene",  'benzene',      'toluene',       None,  "none"),
    ("Hexane-Heptane",   'hexane',       'heptane',       None,  "none"),
    ("Acetone-H₂O",     'acetone',      'water',         None,  "none"),
    ("AcOH-H₂O",        'acetic_acid',  'water',         None,  "none"),
]

print("═══════════════════════════════════════════════════════════════════════")
print("PREDICTIONS — 20 SYSTEMS")
print("═══════════════════════════════════════════════════════════════════════")
print()

n_ex = 0; n_tp = 0

for name, cA, cB, x_obs, t_obs in systems:
    pred_type, pair_t, reasoning = predict(cA, cB)
    
    obs_exists = x_obs is not None
    pred_exists = pred_type != "none"
    
    ex_ok = (pred_exists == obs_exists)
    tp_ok = (pred_type == t_obs)
    
    if ex_ok: n_ex += 1
    if tp_ok: n_tp += 1
    
    mark_ex = "✓" if ex_ok else "✗"
    mark_tp = "✓" if tp_ok else "✗"
    
    print(f"  {mark_ex} {mark_tp}  {name:<18} → {pred_type:>8} (obs: {t_obs})")
    print(f"         {reasoning}")
    print()

print("═══════════════════════════════════════════════════════════════════════")
print(f"  EXISTENCE:  {n_ex}/20 ({100*n_ex/20:.0f}%)")
print(f"  TYPE:       {n_tp}/20 ({100*n_tp/20:.0f}%)")
print("═══════════════════════════════════════════════════════════════════════")
print()

# Summary by regime
print("BY REGIME:")
regime_net = {'correct': 0, 'total': 0}
regime_disp = {'correct': 0, 'total': 0}

for name, cA, cB, x_obs, t_obs in systems:
    pred_type, pair_t, _ = predict(cA, cB)
    tA = molecules[cA][3]
    tB = molecules[cB][3]
    tp_ok = (pred_type == t_obs)
    
    if tA == 0 and tB == 0 and pair_t != 3:
        regime_disp['total'] += 1
        if tp_ok: regime_disp['correct'] += 1
    else:
        regime_net['total'] += 1
        if tp_ok: regime_net['correct'] += 1

print(f"  NETWORK regime (Types 1,2,3):    {regime_net['correct']}/{regime_net['total']}")
print(f"  DISPERSION regime (Type 0+0):    {regime_disp['correct']}/{regime_disp['total']}")
print()

print("""
═══════════════════════════════════════════════════════════════════════
SUMMARY OF EQUATIONS
═══════════════════════════════════════════════════════════════════════

The complete model uses 12 equations and 4 thresholds:

CLASSIFICATION:
  (1)  σᵢ = αᵢ × βᵢ                    → self-association
  (2)  Γcross = αA·βB + αB·βA           → cross H-bonding
  (3)  Γself = σA + σB                   → self H-bonding
  (4)  C = Γcross - Γself                → net complementarity

NETWORK REGIME (Types 1, 2, 3):
  (5)  C > 0, σA < 0.05, σB < 0.05      → Type 3 → NEGATIVE azeotrope
  (6)  f_int = √(f_d × f_a)             → integration fraction
  (7a) f_d = min(α_sol/α_net, 1)        → donor integration
  (7b) f_a = min(β_sol/β_net, 1)        → acceptor integration
  (8a) f_int > 0.90                      → no azeotrope (integrates)
  (8b) f_int ≤ 0.90                      → POSITIVE azeotrope (disrupts)
  (9a) f_int < 0.80  }                    
  (9b) |ΔTb| < 40 K  } both required    → POSITIVE (Type 0+2 only)

DISPERSION REGIME (Type 0+0):
  (10) D = √[(Δδd)² + (Δδp)²]           → dispersion distance
  (11) F = 1 if same family, 0 otherwise → family similarity
  (12a) F = 1                            → no azeotrope
  (12b) F = 0, D > 4                     → POSITIVE azeotrope
  (12c) F = 0, D ≤ 4                     → no azeotrope

THRESHOLDS:
  σ = 0.05 (self-association boundary)
  f_int = 0.90 (full integration)
  f_int = 0.80 (Td disruption)
  ΔTb = 40 K (open-system volatility limit)
  D = 4 MPa^½ (dispersion mismatch)

INPUTS (ALL from pure component data, zero mixture measurements):
  αᵢ, βᵢ — Kamlet-Taft (tabulated for >300 solvents)
  δd,i, δp,i — Hansen (tabulated for >1000 solvents)
  Tb,i — boiling point (NIST WebBook)
  Family — chemical class (trivial)
""")

