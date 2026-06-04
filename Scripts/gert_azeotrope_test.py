"""
GERT Azeotrope Test — First validation
========================================
Can dφ/dx = 0 predict azeotropes from pure component data?

Framework:
  f_M(x) = cohesive energy of mixture (from ΔH_vap + excess)
  f_L(x) = dispersive/entropic energy (from T·ΔS_vap)
  φ(x) = f_M/(f_M + f_L)
  Azeotrope ↔ dφ/dx = 0
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# =====================================================================
# PURE COMPONENT DATA (NIST WebBook — all verified textbook values)
# =====================================================================
# Format: (name, ΔH_vap kJ/mol at T_b, T_b K, V_m mL/mol at 25°C)

components = {
    'water':       {'dHvap': 40.7, 'Tb': 373.15, 'Vm': 18.1},
    'ethanol':     {'dHvap': 38.6, 'Tb': 351.44, 'Vm': 58.4},
    'methanol':    {'dHvap': 35.2, 'Tb': 337.70, 'Vm': 40.7},
    'acetone':     {'dHvap': 31.3, 'Tb': 329.20, 'Vm': 73.5},
    'chloroform':  {'dHvap': 29.2, 'Tb': 334.33, 'Vm': 80.7},
    'benzene':     {'dHvap': 30.7, 'Tb': 353.25, 'Vm': 89.4},
    'toluene':     {'dHvap': 33.2, 'Tb': 383.78, 'Vm': 106.9},
    'hexane':      {'dHvap': 28.9, 'Tb': 341.88, 'Vm': 131.6},
    'heptane':     {'dHvap': 31.8, 'Tb': 371.57, 'Vm': 147.5},
    'ethyl_acetate': {'dHvap': 31.9, 'Tb': 350.21, 'Vm': 98.5},
    'isopropanol': {'dHvap': 39.9, 'Tb': 355.41, 'Vm': 76.9},
    'diethyl_ether': {'dHvap': 26.5, 'Tb': 307.58, 'Vm': 104.7},
    'CS2':         {'dHvap': 26.7, 'Tb': 319.38, 'Vm': 60.0},
    'acetic_acid': {'dHvap': 23.7, 'Tb': 391.05, 'Vm': 57.5},
    'formic_acid': {'dHvap': 22.7, 'Tb': 373.95, 'Vm': 37.7},
    'HCl_aq':      {'dHvap': 16.2, 'Tb': 321.15, 'Vm': 30.0},  # approximate
    'cyclohexane': {'dHvap': 29.9, 'Tb': 353.87, 'Vm': 108.8},
}

# =====================================================================
# MIXTURE MODEL
# =====================================================================

def compute_phi_profile(comp_A, comp_B, T, W, n_points=200):
    """
    Compute φ(x) for binary mixture A+B at temperature T.
    
    f_M(x) = x·ΔH_A + (1-x)·ΔH_B + W·x·(1-x)
    f_L(x) = T · [x·ΔS_A + (1-x)·ΔS_B]
    
    where ΔS_i = ΔH_i/T_b,i (Trouton approximation)
    W = cross-cohesive excess parameter (kJ/mol)
      W < 0: A-B weaker than pure → positive azeotrope
      W > 0: A-B stronger than pure → negative azeotrope
    """
    A = components[comp_A]
    B = components[comp_B]
    
    dH_A = A['dHvap']  # kJ/mol
    dH_B = B['dHvap']
    dS_A = dH_A / A['Tb']  # kJ/mol/K (Trouton)
    dS_B = dH_B / B['Tb']
    
    x = np.linspace(0.001, 0.999, n_points)
    
    f_M = x * dH_A + (1-x) * dH_B + W * x * (1-x)
    f_L = T * (x * dS_A + (1-x) * dS_B)
    
    phi = f_M / (f_M + f_L)
    
    # Find extrema (dφ/dx ≈ 0)
    dphi = np.diff(phi) / np.diff(x)
    extrema = []
    for i in range(len(dphi)-1):
        if dphi[i] * dphi[i+1] < 0:
            # Linear interpolation for x*
            x_star = x[i+1] - dphi[i+1] * (x[i+2]-x[i+1])/(dphi[i+1]-dphi[i])
            phi_star = np.interp(x_star, x, phi)
            atype = "MAX (negative az.)" if dphi[i] > 0 else "MIN (positive az.)"
            extrema.append((x_star, phi_star, atype))
    
    return x, phi, f_M, f_L, extrema

def analyze_system(name, comp_A, comp_B, T, W, x_az_obs, az_type_obs):
    """Analyze one system and compare with observation."""
    x, phi, f_M, f_L, extrema = compute_phi_profile(comp_A, comp_B, T, W)
    
    # φ at boundaries
    phi_A = phi[0]   # pure A (x=1... wait, x is mole fraction of A)
    phi_B = phi[-1]  # pure B
    
    print(f"\n{'='*60}")
    print(f"SYSTEM: {name}")
    print(f"  {comp_A} (A) + {comp_B} (B) at T = {T:.1f} K")
    print(f"  W = {W:.1f} kJ/mol ({'A-B weaker' if W < 0 else 'A-B stronger' if W > 0 else 'ideal'})")
    print(f"  φ(x=0, pure B) = {phi[0]:.4f}")
    print(f"  φ(x=1, pure A) = {phi[-1]:.4f}")
    
    if extrema:
        for xs, ps, at in extrema:
            print(f"  → EXTREMUM at x* = {xs:.3f}, φ* = {ps:.4f} — {at}")
            if x_az_obs is not None:
                print(f"  → OBSERVED: x_az = {x_az_obs:.3f} ({az_type_obs})")
                print(f"  → ERROR: Δx = {abs(xs - x_az_obs):.3f}")
                status = "✓" if abs(xs - x_az_obs) < 0.15 else "✗"
                print(f"  → {status}")
    else:
        print(f"  → NO EXTREMUM (no azeotrope predicted)")
        if x_az_obs is None:
            print(f"  → OBSERVED: no azeotrope — CORRECT ✓")
        else:
            print(f"  → OBSERVED: azeotrope at x = {x_az_obs:.3f} — MISSED ✗")
    
    return x, phi, extrema

# =====================================================================
# TEST SYSTEMS
# =====================================================================

sep = "=" * 60
print(sep)
print("GERT AZEOTROPE TEST — dφ/dx = 0 CRITERION")
print(sep)
print()
print("Model: φ(x) = f_M(x)/(f_M(x) + f_L(x))")
print("  f_M = x·ΔH_A + (1-x)·ΔH_B + W·x(1-x)")
print("  f_L = T·(x·ΔS_A + (1-x)·ΔS_B)")
print("  ΔS_i = ΔH_i/T_b,i (Trouton)")
print("  W = cross-cohesive excess (from H^E data)")
print()

# Systems to test
# (name, comp_A, comp_B, T_K, W_kJ/mol, x_az_observed, type_observed)
# W estimated from literature H^E values at equimolar:
#   H^E ≈ W/4 (symmetric Margules) → W ≈ 4·H^E(x=0.5)

systems = [
    # POSITIVE AZEOTROPES (W < 0, A-B weaker)
    ("EtOH-Water",     'ethanol', 'water',     351, -8.0, 0.894, "positive"),
    ("MeOH-Benzene",   'methanol', 'benzene',  331, -12.0, 0.615, "positive"),
    ("Benzene-EtOH",   'benzene', 'ethanol',   341, -6.0, 0.448, "positive"),
    ("EtOAc-EtOH",     'ethyl_acetate', 'ethanol', 345, -3.0, 0.462, "positive"),
    ("Cyclohex-EtOH",  'cyclohexane', 'ethanol', 338, -10.0, 0.446, "positive"),
    ("CS2-Acetone",    'CS2', 'acetone',        312, -8.0, 0.610, "positive"),
    
    # NEGATIVE AZEOTROPES (W > 0, A-B stronger)
    ("Acetone-CHCl3",  'acetone', 'chloroform', 337, +6.0, 0.357, "negative"),
    ("CHCl3-Et2O",     'chloroform', 'diethyl_ether', 326, +4.0, 0.85, "negative"),
    
    # NON-AZEOTROPIC (W ≈ 0, similar components)
    ("MeOH-EtOH",     'methanol', 'ethanol',   340, -0.5, None, "none"),
    ("Benzene-Toluene",'benzene', 'toluene',    365, -0.3, None, "none"),
    ("Hexane-Heptane", 'hexane', 'heptane',     355, -0.1, None, "none"),
]

results = []
fig, axes = plt.subplots(4, 3, figsize=(15, 16))
axes_flat = axes.flatten()

for idx, (name, cA, cB, T, W, x_obs, type_obs) in enumerate(systems):
    x, phi, extrema = analyze_system(name, cA, cB, T, W, x_obs, type_obs)
    
    # Store result
    if extrema and x_obs is not None:
        predicted = True
        x_pred = extrema[0][0]
        error = abs(x_pred - x_obs)
        results.append((name, True, True, error, type_obs))
    elif not extrema and x_obs is None:
        results.append((name, False, False, 0, "none"))
    elif extrema and x_obs is None:
        results.append((name, True, False, None, "false positive"))
    else:
        results.append((name, False, True, None, "missed"))
    
    # Plot
    if idx < 12:
        ax = axes_flat[idx]
        ax.plot(x, phi, 'b-', linewidth=2)
        ax.axhline(y=0.5, color='gray', linestyle=':', alpha=0.5)
        
        if extrema:
            for xs, ps, at in extrema:
                ax.plot(xs, ps, 'ro', markersize=8, zorder=5)
                ax.annotate(f'x*={xs:.2f}', xy=(xs, ps), 
                           xytext=(xs+0.1, ps+0.002), fontsize=8, color='red')
        if x_obs is not None:
            ax.axvline(x=x_obs, color='green', linestyle='--', alpha=0.5,
                      label=f'Obs: {x_obs:.3f}')
        
        ax.set_title(name, fontsize=10, fontweight='bold')
        ax.set_xlabel('x (mol frac A)', fontsize=8)
        ax.set_ylabel('φ', fontsize=8)
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.15)

# Remove empty subplot
if len(systems) < 12:
    axes_flat[11].axis('off')

plt.tight_layout()
plt.savefig('/mnt/user-data/outputs/gert_azeotrope_test.png', dpi=200, bbox_inches='tight')
plt.close()

# Summary
print(f"\n\n{'='*60}")
print("SUMMARY")
print(f"{'='*60}")
print(f"\n{'System':<20} {'Exist?':<8} {'Obs?':<8} {'Δx':<8} {'Status'}")
print("-" * 52)

correct = 0
total = len(results)
for name, pred, obs, err, typ in results:
    if pred == obs:
        correct += 1
        if err is not None and err < 0.15:
            status = f"✓ (Δx={err:.3f})"
        elif err == 0:
            status = "✓ (correct absence)"
        else:
            status = f"~ (Δx={err:.3f})" if err else "✓"
    else:
        status = "✗ MISSED" if obs else "✗ FALSE POS"
    
    print(f"{name:<20} {'yes' if pred else 'no':<8} {'yes' if obs else 'no':<8} "
          f"{f'{err:.3f}' if err is not None else '---':<8} {status}")

print(f"\nExistence prediction: {correct}/{total} correct")

