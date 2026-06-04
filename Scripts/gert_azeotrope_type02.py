"""
GERT Azeotrope — Type 0+2 fix: polarity-corrected integration
================================================================
The "k_gas" of azeotropes: for Type 0 acceptor in water Td,
integration depends on POLARITY MATCH, not just β.

Physics: to sit in the polar Td network without donating (α≈0),
the solute needs polarity compatible with water. β alone says
"can accept H-bonds" but δ_p says "fits in the polar environment."

f_int(Type 0+2) = (β_sol/β_water) × (δ_p,sol/δ_p,water)
"""

import numpy as np

data = {
    #                  α     β     σ    type  δd    δp    family  Vm
    'water':         (1.17, 0.47, 0.550, 2,  15.5, 16.0, 'water', 18.1),
    'ethanol':       (0.86, 0.75, 0.645, 1,  15.8,  8.8, 'alcohol', 58.4),
    'methanol':      (0.98, 0.66, 0.647, 1,  15.1, 12.3, 'alcohol', 40.7),
    'isopropanol':   (0.76, 0.84, 0.638, 1,  15.8,  6.1, 'alcohol', 76.9),
    'n_butanol':     (0.84, 0.84, 0.706, 1,  16.0,  5.7, 'alcohol', 91.5),
    'acetic_acid':   (1.12, 0.45, 0.504, 1,  14.5,  8.0, 'acid', 57.5),
    'acetone':       (0.08, 0.43, 0.034, 0,  15.5, 10.4, 'ketone', 73.5),
    'chloroform':    (0.20, 0.10, 0.020, 0,  17.8,  3.1, 'haloalkane', 80.7),
    'benzene':       (0.00, 0.10, 0.000, 0,  18.4,  0.0, 'aromatic', 89.4),
    'toluene':       (0.00, 0.11, 0.000, 0,  18.0,  1.4, 'aromatic', 106.9),
    'hexane':        (0.00, 0.00, 0.000, 0,  14.9,  0.0, 'alkane', 131.6),
    'heptane':       (0.00, 0.00, 0.000, 0,  15.3,  0.0, 'alkane', 147.5),
    'CS2':           (0.00, 0.07, 0.000, 0,  20.5,  0.0, 'other', 60.0),
    'cyclohexane':   (0.00, 0.00, 0.000, 0,  16.8,  0.0, 'alkane', 108.8),
    'ethyl_acetate': (0.00, 0.45, 0.000, 0,  15.8,  5.3, 'ester', 98.5),
    'diethyl_ether': (0.00, 0.47, 0.000, 0,  14.5,  2.9, 'ether', 104.7),
    'THF':           (0.00, 0.55, 0.000, 0,  16.8,  5.7, 'ether', 81.7),
    'DCM':           (0.13, 0.10, 0.013, 0,  18.2,  6.3, 'haloalkane', 63.9),
}

# =====================================================================
# TYPE 0+2 ANALYSIS
# =====================================================================

print("=" * 80)
print("TYPE 0+2 ANALYSIS — Polarity-corrected integration")
print("=" * 80)
print()
print("For Type 0 solute in Type 2 water network:")
print("  f_int = (β_sol/β_water) × (δ_p,sol/δ_p,water)")
print("  Threshold: f_int > 0.50 → integrates → no azeotrope")
print("             f_int ≤ 0.50 → disrupts → positive azeotrope")
print()

type02_systems = [
    # (solute, observed_azeotrope, x_az_observed)
    ('acetone',       False, None,  "β=0.43, δp=10.4 → polar acceptor"),
    ('THF',           True,  0.820, "β=0.55, δp=5.7 → moderate polar"),
    ('diethyl_ether', True,  0.953, "β=0.47, δp=2.9 → weak polar"),
    ('ethyl_acetate', True,  0.700, "β=0.45, δp=5.3 → moderate polar"),
    ('hexane',        True,  None,  "β=0.00, δp=0.0 → no integration"),
    ('cyclohexane',   True,  None,  "β=0.00, δp=0.0 → no integration"),
    ('benzene',       True,  None,  "β=0.10, δp=0.0 → negligible"),
    ('DCM',           True,  None,  "β=0.10, δp=6.3 → weak"),
    ('CS2',           True,  None,  "β=0.07, δp=0.0 → no integration"),
    ('chloroform',    True,  None,  "β=0.10, δp=3.1 → negligible"),
]

print(f"{'Solute':<16} {'β':>5} {'δp':>5} {'β/βw':>6} {'δp/δpw':>7} {'f_int':>6} {'Pred':>8} {'Obs':>8} {'':>3}")
print("-" * 78)

b_water = data['water'][1]  # 0.47
dp_water = data['water'][5]  # 16.0
n_correct = 0

for sol, has_az, x_obs, note in type02_systems:
    a_sol, b_sol = data[sol][0], data[sol][1]
    dp_sol = data[sol][5]
    
    # Polarity-corrected f_int
    b_ratio = b_sol / b_water if b_water > 0 else 0
    dp_ratio = dp_sol / dp_water if dp_water > 0 else 0
    f_int = b_ratio * dp_ratio
    
    pred = "none" if f_int > 0.50 else "positive"
    obs = "positive" if has_az else "none"
    ok = pred == obs
    if ok: n_correct += 1
    
    print(f"{sol:<16} {b_sol:>5.2f} {dp_sol:>5.1f} {b_ratio:>6.2f} {dp_ratio:>7.3f} "
          f"{f_int:>6.3f} {pred:>8} {obs:>8} {'✓' if ok else '✗':>3}")

print()
print(f"Score: {n_correct}/{len(type02_systems)} correct")
print()

# =====================================================================
# PHYSICAL INTERPRETATION
# =====================================================================

print("=" * 80)
print("PHYSICAL INTERPRETATION")
print("=" * 80)
print("""
The polarity-corrected integration fraction:

  f_int(Type 0+2) = (β_sol / β_water) × (δ_p,sol / δ_p,water)    (Eq. NEW)

Physics: for a Type 0 molecule to sit in the Td water network
WITHOUT donating (α ≈ 0), it needs TWO things:

  1. β > 0: can ACCEPT H-bonds from water's O-H donors
     (measured by β_sol/β_water)
     
  2. δ_p match: is POLAR enough to be compatible with the 
     polar environment of the network
     (measured by δ_p,sol/δ_p,water)

Why acetone integrates and Et₂O doesn't:
  Acetone: β/βw = 0.91, δp/δpw = 0.65 → f_int = 0.594 > 0.50 → INTEGRATES
  Et₂O:   β/βw = 1.00, δp/δpw = 0.18 → f_int = 0.181 < 0.50 → DISRUPTS

  Both are good acceptors (β ≈ 0.45).
  But acetone is POLAR (δp = 10.4) and Et₂O is NOT (δp = 2.9).
  The Td water network is highly polar.
  A non-polar acceptor sits in it poorly → disrupts → azeotrope.

This is the "k_gas" of the azeotrope framework:
  One correction (polarity) fixes the entire Type 0+2 class.
  Like k_gas, it is NOT a free parameter — it uses δ_p already tabulated.
  Like k_gas, it captures physics that the simpler model missed.
""")

# =====================================================================
# FULL TEST: ALL 24 SYSTEMS (original 20 + 4 Type 0+2)
# =====================================================================

def predict_full(mol_A, mol_B):
    aA,bA,sA,tA,dA,pA,fA,vA = data[mol_A]
    aB,bB,sB,tB,dB,pB,fB,vB = data[mol_B]
    
    # Type 3
    G_cross = aA*bB + aB*bA
    G_self = sA + sB
    C = G_cross - G_self
    if C > 0 and sA < 0.05 and sB < 0.05:
        return "negative", 3, f"C={C:.3f}>0 → emergent"
    
    pair_t = max(tA, tB)
    
    # Type 0+0
    if tA == 0 and tB == 0:
        D = np.sqrt((dA-dB)**2 + (pA-pB)**2)
        F = 1 if fA == fB else 0
        if F == 1: return "none", 0, f"Same family ({fA})"
        elif D > 4.0: return "positive", 0, f"D={D:.1f}>4"
        else: return "none", 0, f"D={D:.1f}≤4"
    
    # Identify network holder
    if tA >= tB:
        net_a,net_b,net_t,net_f,net_p = aA,bA,tA,fA,pA
        sol_a,sol_b,sol_t,sol_f,sol_p = aB,bB,tB,fB,pB
        net_mol, sol_mol = mol_A, mol_B
    else:
        net_a,net_b,net_t,net_f,net_p = aB,bB,tB,fB,pB
        sol_a,sol_b,sol_t,sol_f,sol_p = aA,bA,tA,fA,pA
        net_mol, sol_mol = mol_B, mol_A
    
    # Type 0+2 — POLARITY-CORRECTED f_int
    if net_t == 2 and sol_t == 0:
        b_ratio = sol_b / net_b if net_b > 0 else 0
        dp_ratio = sol_p / net_p if net_p > 0 else 0
        f_int = b_ratio * dp_ratio
        
        if f_int > 0.50:
            return "none", 2, f"Type 0+2: f_int={f_int:.3f}>0.50 (polarity OK)"
        else:
            return "positive", 2, f"Type 0+2: f_int={f_int:.3f}≤0.50 (disrupts)"
    
    # Type 1+2 — standard f_int
    if net_t == 2 and sol_t == 1:
        f_d = min(sol_a/net_a, 1.0) if net_a > 0 else 0
        f_a = min(sol_b/net_b, 1.0) if net_b > 0 else 0
        f_int = np.sqrt(f_d * f_a)
        if f_int > 0.90:
            return "none", 2, f"Type 1+2: f_int={f_int:.2f}>0.90 → integrates"
        else:
            return "positive", 2, f"Type 1+2: f_int={f_int:.2f}≤0.90 → disrupts"
    
    # Type 1+1
    if net_t == 1 and sol_t == 1:
        if net_f == sol_f: return "none", 1, f"Same class ({net_f})"
        else: return "positive", 1, f"Diff class ({net_f}+{sol_f})"
    
    # Type 0+1
    if net_t == 1 and sol_t == 0:
        return "positive", 1, f"{sol_mol}(0) disrupts {net_mol}(1)"
    
    return "unknown", -1, "?"

# All systems including Type 0+2
all_systems = [
    ("MeOH–Benzene",     'methanol',     'benzene',       0.615, "positive"),
    ("Benzene–EtOH",     'benzene',      'ethanol',       0.448, "positive"),
    ("EtOAc–EtOH",       'ethyl_acetate','ethanol',       0.462, "positive"),
    ("Cyclohex–EtOH",    'cyclohexane',  'ethanol',       0.446, "positive"),
    ("MeOH–Cyclohex",    'methanol',     'cyclohexane',   0.510, "positive"),
    ("DCM–MeOH",         'DCM',          'methanol',      0.810, "positive"),
    ("THF–EtOH",         'THF',          'ethanol',       0.650, "positive"),
    ("Hexane–EtOH",      'hexane',       'ethanol',       0.345, "positive"),
    ("Acetone–MeOH",     'acetone',      'methanol',      0.800, "positive"),
    ("CS₂–Acetone",      'CS2',          'acetone',       0.610, "positive"),
    ("EtOH–H₂O",        'ethanol',      'water',         0.894, "positive"),
    ("iPrOH–H₂O",       'isopropanol',  'water',         0.674, "positive"),
    ("nBuOH–H₂O",       'n_butanol',    'water',         0.749, "positive"),
    ("Acetone–CHCl₃",    'acetone',      'chloroform',    0.357, "negative"),
    ("CHCl₃–Et₂O",       'chloroform',   'diethyl_ether', 0.850, "negative"),
    ("MeOH–EtOH",        'methanol',     'ethanol',       None,  "none"),
    ("EtOH–iPrOH",       'ethanol',      'isopropanol',   None,  "none"),
    ("Benzene–Toluene",   'benzene',      'toluene',       None,  "none"),
    ("Hexane–Heptane",    'hexane',       'heptane',       None,  "none"),
    ("AcOH–H₂O",         'acetic_acid',  'water',         None,  "none"),
    # TYPE 0+2 — previously excluded
    ("Acetone–H₂O",      'acetone',      'water',         None,  "none"),
    ("THF–H₂O",          'THF',          'water',         0.820, "positive"),
    ("Et₂O–H₂O",         'diethyl_ether','water',         0.953, "positive"),
    ("EtOAc–H₂O",        'ethyl_acetate','water',         0.700, "positive"),
]

print()
print("=" * 80)
print(f"FULL TEST — {len(all_systems)} SYSTEMS (INCLUDING TYPE 0+2)")
print("=" * 80)
print()

n_ex=0; n_tp=0
for i,(name,cA,cB,x_obs,t_obs) in enumerate(all_systems,1):
    pred,pt,reason = predict_full(cA,cB)
    tA,tB = data[cA][3], data[cB][3]
    obs_ex = x_obs is not None; pred_ex = pred!="none"
    ex_ok = pred_ex==obs_ex; tp_ok = pred==t_obs
    if ex_ok: n_ex+=1
    if tp_ok: n_tp+=1
    m = "✓✓" if (ex_ok and tp_ok) else "✗✗"
    marker = " ← NEW" if i > 20 else ""
    print(f"{i:>2} {m} {name:<18} {tA}+{tB}→{pt} {pred:>8} {t_obs:>8} {reason}{marker}")

print()
print("=" * 80)
print(f"EXISTENCE: {n_ex}/{len(all_systems)} ({100*n_ex/len(all_systems):.0f}%)")
print(f"TYPE:      {n_tp}/{len(all_systems)} ({100*n_tp/len(all_systems):.0f}%)")
print("=" * 80)

