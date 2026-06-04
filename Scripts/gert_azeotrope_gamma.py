#!/usr/bin/env python3
"""
GERT Azeotrope — Complete Formula with γ∞

φ_int = f_M / (f_M + f_L)    with threshold 1/2

f_M = (β_sol/β_w)(δp_sol/δp_w)        — integration capacity
f_L = max(ln(γ∞)/C, logP/C')          — disruption (dominant mechanism)

C  = 3.73  (calibrated on acetone/MEK gap)
C' = 0.48  (calibrated on cyclopentanone/cyclohexanone)

Author: V. P. Dutra
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import matplotlib.patches as mpatches

# ============================================================
# CONSTANTS
# ============================================================
BETA_W = 0.47       # β water (Marcus 1998, Table 4.3)
DP_W   = 16.0       # δp water (Hansen 2007, Table A.1)
C_GAMMA = 3.729     # calibrated on acetone/MEK gap, max min margin
C_LOGP  = 0.48      # calibrated on cyclopentanone/cyclohexanone

# ============================================================
# DATABASE — ALL VERIFIED SYSTEMS
# ============================================================
# (name, family, β, δp, γ∞, logP, observed, az_type,
#  β_source, γ∞_source, logP_source)

TYPE_02_SYSTEMS = [
    # KETONES — homologous series
    ("Acetone",         "Ketone",  0.48, 10.4,   7.65, -0.24, "none",     "homo", "K83","Brouwer","PubChem"),
    ("MEK",             "Ketone",  0.48,  9.0,  27.6,   0.29, "positive", "homo", "K83","Brouwer","PubChem"),
    ("2-Pentanone",     "Ketone",  0.50,  7.5, 102.0,   0.91, "positive", "homo", "K83","Brouwer","PubChem"),
    ("3-Pentanone",     "Ketone",  0.45,  7.6, 113.0,   0.82, "positive", "homo", "K83","Brouwer","PubChem"),
    ("Cyclopentanone",  "Ketone",  0.52, 11.9,  16.1,   0.40, "positive", "hetero","K83","Brouwer","XLogP3"),
    ("Cyclohexanone",   "Ketone",  0.53,  8.4,   5.14,  0.80, "positive", "hetero","K83","Brouwer","XLogP3"),

    # ESTERS — homologous series
    ("Methyl formate",  "Ester",   0.37,  8.4,  15.75,  0.03, "positive", "homo", "K83","Brouwer","PubChem"),
    ("Ethyl formate",   "Ester",   0.36,  7.2,  46.65,  0.23, "positive", "homo", "K83","Brouwer","PubChem"),
    ("Methyl acetate",  "Ester",   0.42,  7.6,  22.5,   0.18, "positive", "homo", "K83","Brouwer","PubChem"),
    ("Ethyl acetate",   "Ester",   0.45,  5.3,  75.6,   0.73, "positive", "homo", "K83","Brouwer","PubChem"),
    ("Propyl acetate",  "Ester",   0.45,  3.3, 274.5,   1.24, "positive", "homo", "K83","Brouwer","PubChem"),
    ("Butyl acetate",   "Ester",   0.45,  3.7,1058.0,   1.78, "positive", "hetero","K83","Brouwer","PubChem"),

    # ETHERS — homologous series
    ("THF",             "Ether",   0.55,  5.7,  17.0,   0.46, "positive", "homo", "K83","Brouwer","PubChem"),
    ("Tetrahydropyran", "Ether",   0.54,  4.5,  70.5,   0.82, "positive", "homo", "K83","Brouwer","PubChem"),
    ("1,3-Dioxolane",   "Ether",   0.45,  6.6,   9.71, -0.37, "positive", "homo", "K83","Brouwer","PubChem"),
    ("1,4-Dioxane",     "Ether",   0.37,  1.8,   5.44, -0.27, "positive", "homo", "K83","Brouwer","PubChem"),
    ("Diethyl ether",   "Ether",   0.47,  2.9,  77.1,   0.89, "positive", "homo", "K83","Brouwer","PubChem"),
    ("MTBE",            "Ether",   0.55,  4.3, 113.0,   0.94, "positive", "homo", "est","Brouwer","PubChem"),
    ("DiisoPr ether",   "Ether",   0.49,  3.4, 628.0,   1.52, "positive", "homo", "est","Brouwer","PubChem"),
    ("Di-n-propyl ether","Ether",  0.46,  2.4,2313.0,   1.21, "positive", "homo", "est","Brouwer","PubChem"),
    ("Di-n-butyl ether","Ether",   0.46,  1.6,47180.,   3.21, "positive", "homo", "est","Brouwer","PubChem"),
    ("Anisole",         "Ether",   0.22,  4.1,4000.0,   2.11, "positive", "homo", "K83","Brouwer","PubChem"),

    # AMIDES — strong integrators
    ("NMP",             "Amide",   0.77, 12.3,   0.37, -0.38, "none",     "homo", "K83","Brouwer","PubChem"),
    ("DMF",             "Amide",   0.69, 13.7,   0.70, -1.01, "none",     "homo", "K83","Brouwer","PubChem"),
    ("DMAc",            "Amide",   0.76, 11.5,   1.04, -0.77, "none",     "homo", "est","Brouwer","PubChem"),

    # NITRILES
    ("Propionitrile",   "Nitrile", 0.37, 16.1,  36.4,   0.16, "positive", "homo", "K83","Brouwer","PubChem"),
    ("Benzonitrile",    "Nitrile", 0.37,  8.4,1741.0,   1.56, "positive", "homo", "est","Brouwer","PubChem"),

    # N-AROMATICS
    ("Pyridine",        "N-arom",  0.64,  8.8,  31.8,   0.65, "positive", "homo", "K83","Brouwer","PubChem"),
    ("3-MePyridine",    "N-arom",  0.68,  7.8,  49.1,   1.20, "positive", "homo", "K83","Brouwer","PubChem"),
    ("4-MePyridine",    "N-arom",  0.67,  7.8,  42.3,   1.22, "positive", "homo", "K83","Brouwer","PubChem"),

    # ALDEHYDES
    ("Acetaldehyde",    "Aldehyde",0.45, 11.3,   4.15, -0.34, "none",     "homo", "est","Brouwer","Ullmann"),
    ("Furfural",        "Aldehyde",0.34, 14.9,  97.4,   0.41, "positive", "homo", "est","Brouwer","PubChem"),

    # HALOGENATED (donors, β≈0)
    ("DCM",             "Haloalk", 0.00,  6.3, 250.0,   1.25, "positive", "homo", "K83","Brouwer","PubChem"),
    ("Chloroform",      "Haloalk", 0.00,  3.1, 835.0,   1.97, "positive", "homo", "K83","Brouwer","PubChem"),
    ("1,2-DCE",         "Haloalk", 0.00, 10.4, 641.0,   1.48, "positive", "homo", "K83","Brouwer","PubChem"),

    # OTHER
    ("DMSO",            "Sulfox",  0.76, 16.4,   1.3,  -1.35, "none",     "homo", "M98","NIST",  "PubChem"),
    ("Triethylamine",   "Amine",   0.71,  0.4,  67.5,   1.45, "positive", "homo", "K83","Brouwer","PubChem"),
    ("Nitromethane",    "Nitro",   0.06, 18.8,  31.6,  -0.33, "positive", "homo", "K83","Brouwer","PubChem"),
    ("Nitroethane",     "Nitro",   0.25, 15.5,  75.0,   0.18, "positive", "homo", "est","Brouwer","PubChem"),
]


def compute_phi(beta, dp, gamma, logP):
    """Compute φ_int = f_M / (f_M + f_L) with f_L = max(ln(γ∞)/C, logP/C')"""
    f_M = (beta / BETA_W) * (dp / DP_W)
    f_L_gamma = np.log(gamma) / C_GAMMA if gamma > 0 else 0
    f_L_logp  = logP / C_LOGP
    f_L = max(f_L_gamma, f_L_logp)
    if (f_M + f_L) == 0:
        return f_M, f_L_gamma, f_L_logp, f_L, 0.0
    phi = f_M / (f_M + f_L)
    return f_M, f_L_gamma, f_L_logp, f_L, phi


def predict(phi):
    return "none" if phi > 0.5 else "positive"


# ============================================================
# MAIN CALCULATION
# ============================================================
print("=" * 110)
print(f"φ = f_M / (f_M + max(ln(γ∞)/{C_GAMMA:.3f}, logP/{C_LOGP:.2f})),  threshold = 1/2")
print("=" * 110)
print(f"\n{'Solute':<20} {'Fam':<8} {'f_M':>6} {'γ∞':>8} {'logP':>6} {'f_L(γ)':>7} {'f_L(P)':>7} "
      f"{'f_L':>7} {'φ':>6} {'Pred':<9} {'Obs':<9} {'Dom':<5}")
print("-" * 110)

ok = 0; total = 0
results = []
for entry in TYPE_02_SYSTEMS:
    name, fam, beta, dp, gamma, logP, obs, az_type = entry[:8]
    f_M, f_Lg, f_Lp, f_L, phi = compute_phi(beta, dp, gamma, logP)
    pred = predict(phi)
    dom = "logP" if f_Lp > f_Lg else "γ∞"
    match = "✓" if pred == obs else "✗"
    if pred == obs: ok += 1
    total += 1
    results.append((name, fam, beta, dp, gamma, logP, obs, az_type,
                     f_M, f_Lg, f_Lp, f_L, phi, pred, dom))
    print(f"{name:<20} {fam:<8} {f_M:6.3f} {gamma:8.1f} {logP:6.2f} {f_Lg:7.3f} {f_Lp:7.3f} "
          f"{f_L:7.3f} {phi:6.3f} {pred:<9} {obs:<9} {dom:<5} {match}")

print(f"\nSCORE: {ok}/{total}")

# ============================================================
# R ANALYSIS
# ============================================================
print(f"\n{'=' * 80}")
print("R = ln(γ∞)/f_M — score de separação")
print(f"{'=' * 80}")

R_data = []
for r in results:
    name, fam, beta, dp, gamma, logP, obs = r[:7]
    f_M = r[8]
    R = np.log(gamma) / f_M if f_M > 0.001 else 999
    R_data.append((name, fam, R, f_M, gamma, obs))

nones = sorted([x for x in R_data if x[5] == "none"], key=lambda x: -x[2])
positives = sorted([x for x in R_data if x[5] == "positive"], key=lambda x: x[2])

R_max_none = max(x[2] for x in nones if x[2] < 100)
R_min_pos = min(x[2] for x in positives if x[2] > -100)

print(f"\n  max(R_none)     = {R_max_none:.3f}  (frontier: {[x[0] for x in nones if x[2]==R_max_none][0]})")
print(f"  min(R_positive) = {R_min_pos:.3f}  (frontier: {[x[0] for x in positives if x[2]==R_min_pos][0]})")
print(f"  GAP = [{R_max_none:.3f}, {R_min_pos:.3f}], width = {R_min_pos - R_max_none:.3f}")
print(f"  C = {C_GAMMA:.3f} (optimal midpoint)")

# ============================================================
# HOMOLOGOUS SERIES ANALYSIS
# ============================================================
print(f"\n{'=' * 80}")
print("HOMOLOGOUS SERIES")
print(f"{'=' * 80}")

series_data = {
    "Ketones (R-CO-R)": [
        ("Acetone (C3)",    3,   7.65),
        ("MEK (C4)",        4,  27.6),
        ("2-Pentanone (C5)",5, 102.0),
        ("3-Pentanone (C5)",5, 113.0),
    ],
    "Acetate esters (R-OAc)": [
        ("Methyl acetate (C3)",  3,  22.5),
        ("Ethyl acetate (C4)",   4,  75.6),
        ("Propyl acetate (C5)",  5, 274.5),
        ("Butyl acetate (C6)",   6,1058.0),
    ],
    "Ethers (R-O-R)": [
        ("THF (cyc-C4)",         4,  17.0),
        ("THP (cyc-C5)",         5,  70.5),
        ("Et₂O (C4)",           4,  77.1),
        ("DiisoPr (C6)",        6, 628.0),
    ],
    "Amides": [
        ("NMP",    5,  0.37),
        ("DMF",    3,  0.70),
        ("DMAc",   4,  1.04),
    ],
}

for title, members in series_data.items():
    print(f"\n  {title}:")
    for name, nc, g in members:
        print(f"    {name:<25} C{nc}  γ∞ = {g:>8.1f}  ln(γ∞) = {np.log(g):>6.2f}")


# ============================================================
# FIGURES
# ============================================================
plt.rcParams.update({
    'font.size': 11, 'axes.labelsize': 13, 'axes.titlesize': 13,
    'figure.dpi': 150, 'savefig.dpi': 300
})

# --- FIGURE 1: φ for all systems ---
fig, ax = plt.subplots(figsize=(14, 7))

# Sort by φ
sorted_results = sorted(results, key=lambda x: x[12])  # phi
names_sorted = [r[0] for r in sorted_results]
phis_sorted = [r[12] for r in sorted_results]
obs_sorted = [r[6] for r in sorted_results]
dom_sorted = [r[14] for r in sorted_results]

colors = []
for obs, dom in zip(obs_sorted, dom_sorted):
    if obs == "none":
        colors.append('#2ecc71')  # green
    elif dom == "logP":
        colors.append('#e74c3c')  # red (logP-driven)
    else:
        colors.append('#3498db')  # blue (γ∞-driven)

bars = ax.barh(range(len(names_sorted)), phis_sorted, color=colors, edgecolor='white', linewidth=0.5)
ax.axvline(x=0.5, color='black', linestyle='--', linewidth=2, label='φ = 1/2 (GERT threshold)')
ax.set_yticks(range(len(names_sorted)))
ax.set_yticklabels(names_sorted, fontsize=8)
ax.set_xlabel('φ_int = f_M / (f_M + f_L)')
ax.set_title('Integration fraction φ for Type 0+2 systems in water\n'
             'Green: integrates (none) | Blue: disrupts via γ∞ | Red: disrupts via logP')
ax.set_xlim(0, 1.4)
ax.legend(loc='lower right', fontsize=10)
plt.tight_layout()
plt.savefig('fig_gamma_phi_all_systems.png')
plt.close()
print("\n  → fig_gamma_phi_all_systems.png")


# --- FIGURE 2: Homologous series γ∞ vs carbon number ---
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Ketones
ax = axes[0]
carbons_k = [3, 4, 5]
gammas_k = [7.65, 27.6, 107.5]  # avg of 2-pent + 3-pent
ax.semilogy(carbons_k, gammas_k, 'o-', color='#e74c3c', markersize=10, linewidth=2)
ax.axhline(y=np.exp(C_GAMMA * 0.5), color='grey', linestyle=':', alpha=0.5)
for c, g, n in zip(carbons_k, gammas_k, ['Acetone', 'MEK', '2/3-Pentanone']):
    ax.annotate(n, (c, g), textcoords="offset points", xytext=(8, 5), fontsize=9)
ax.set_xlabel('Carbon number')
ax.set_ylabel('γ∞ in water')
ax.set_title('Ketones: R-CO-R')
ax.set_xticks([3, 4, 5])
ax.grid(True, alpha=0.3)

# Esters
ax = axes[1]
carbons_e = [3, 4, 5, 6]
gammas_e = [22.5, 75.6, 274.5, 1058]
ax.semilogy(carbons_e, gammas_e, 's-', color='#3498db', markersize=10, linewidth=2)
for c, g, n in zip(carbons_e, gammas_e, ['MeOAc', 'EtOAc', 'PrOAc', 'BuOAc']):
    ax.annotate(n, (c, g), textcoords="offset points", xytext=(8, 5), fontsize=9)
ax.set_xlabel('Carbon number')
ax.set_title('Acetate esters: R-OAc')
ax.set_xticks([3, 4, 5, 6])
ax.grid(True, alpha=0.3)

# Ethers
ax = axes[2]
carbons_et = [4, 4, 5, 6]
gammas_et = [17.0, 77.1, 70.5, 628]
labels_et = ['THF', 'Et₂O', 'THP', 'DiisoPr']
ax.semilogy(carbons_et, gammas_et, 'D-', color='#9b59b6', markersize=10, linewidth=0)
for c, g, n in zip(carbons_et, gammas_et, labels_et):
    ax.annotate(n, (c, g), textcoords="offset points", xytext=(8, 5), fontsize=9)
    ax.plot(c, g, 'D', color='#9b59b6', markersize=10)
ax.set_xlabel('Carbon number')
ax.set_title('Ethers: R-O-R')
ax.set_xticks([3, 4, 5, 6, 7])
ax.grid(True, alpha=0.3)

plt.suptitle('Each CH₂ group multiplies γ∞ by ~3–4×\nThe hydrophobic tail dominates the dispersive penalty',
             fontsize=12, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig('fig_gamma_homologous_series.png')
plt.close()
print("  → fig_gamma_homologous_series.png")


# --- FIGURE 3: f_M vs f_L scatter with threshold ---
fig, ax = plt.subplots(figsize=(10, 8))

for r in results:
    name, fam, beta, dp, gamma, logP, obs, az_type = r[:8]
    f_M, f_Lg, f_Lp, f_L, phi = r[8:13]

    if obs == "none":
        color = '#2ecc71'; marker = 's'; zorder = 5
    elif az_type == "hetero":
        color = '#e74c3c'; marker = '^'; zorder = 6
    else:
        color = '#3498db'; marker = 'o'; zorder = 4

    ax.plot(f_L, f_M, marker, color=color, markersize=10, markeredgecolor='white',
            markeredgewidth=0.5, zorder=zorder)
    # Label key systems
    if name in ['Acetone', 'MEK', 'Pyridine', 'DMF', 'DMSO', 'NMP',
                'Cyclopentanone', 'Cyclohexanone', 'Acetaldehyde',
                'Diethyl ether', 'Anisole', 'Triethylamine']:
        ax.annotate(name, (f_L, f_M), textcoords="offset points",
                    xytext=(6, 4), fontsize=7, alpha=0.8)

# Threshold line: f_M = f_L (φ = 1/2)
x_line = np.linspace(-0.5, 3.5, 100)
ax.plot(x_line, x_line, 'k--', linewidth=2, label='φ = 1/2  (f_M = f_L)')
ax.fill_between(x_line, x_line, 2.0, alpha=0.05, color='green')
ax.fill_between(x_line, 0, x_line, alpha=0.05, color='red')

ax.set_xlabel('f_L = max(ln(γ∞)/C, logP/C\')')
ax.set_ylabel('f_M = (β/β_w)(δp/δp_w)')
ax.set_title('GERT integration diagram: f_M vs f_L\nAbove diagonal: integrates (none) | Below: disrupts (positive)')
ax.set_xlim(-0.5, 3.5)
ax.set_ylim(-0.05, 1.8)

# Legend
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='s', color='w', markerfacecolor='#2ecc71', markersize=10, label='None (integrates)'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#3498db', markersize=10, label='Positive (homo, γ∞-driven)'),
    Line2D([0], [0], marker='^', color='w', markerfacecolor='#e74c3c', markersize=10, label='Positive (hetero, logP-driven)'),
    Line2D([0], [0], color='black', linestyle='--', linewidth=2, label='φ = 1/2 (GERT threshold)'),
]
ax.legend(handles=legend_elements, loc='upper right', fontsize=10)
ax.grid(True, alpha=0.2)
plt.tight_layout()
plt.savefig('fig_gamma_fM_vs_fL.png')
plt.close()
print("  → fig_gamma_fM_vs_fL.png")


# --- FIGURE 4: R = ln(γ∞)/f_M — the gap ---
fig, ax = plt.subplots(figsize=(12, 6))

R_none = [(x[0], x[2]) for x in R_data if x[5] == "none" and x[2] < 100]
R_pos  = [(x[0], x[2]) for x in R_data if x[5] == "positive" and x[2] < 100]

# Sort
R_none_sorted = sorted(R_none, key=lambda x: x[1])
R_pos_sorted = sorted(R_pos, key=lambda x: x[1])

# Plot none
for i, (name, R) in enumerate(R_none_sorted):
    ax.barh(i, R, color='#2ecc71', height=0.7, edgecolor='white')
    ax.text(R + 0.1, i, f'{name} (R={R:.1f})', va='center', fontsize=8)

offset = len(R_none_sorted) + 1
# Gap indicator
ax.axhline(y=offset - 0.5, color='black', linestyle=':', alpha=0.3)
ax.text(C_GAMMA, offset - 0.5, '← GAP →', ha='center', va='center',
        fontsize=12, fontweight='bold', color='red',
        bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.8))

# Plot positive (first 15)
for i, (name, R) in enumerate(R_pos_sorted[:15]):
    ax.barh(i + offset, R, color='#3498db', height=0.7, edgecolor='white')
    if R < 25:
        ax.text(R + 0.1, i + offset, f'{name} (R={R:.1f})', va='center', fontsize=8)

# Threshold line
ax.axvline(x=C_GAMMA, color='red', linestyle='--', linewidth=2,
           label=f'C = {C_GAMMA:.2f} (threshold)')

ax.set_xlabel('R = ln(γ∞) / f_M')
ax.set_title('Separation score R: clear gap between "none" and "positive"\n'
             f'max(R_none) = {R_max_none:.2f} | min(R_positive) = {R_min_pos:.2f} | '
             f'Gap width = {R_min_pos - R_max_none:.2f}')
ax.set_yticks([])
ax.set_xlim(-1.5, 20)
ax.legend(loc='lower right', fontsize=11)
plt.tight_layout()
plt.savefig('fig_gamma_R_gap.png')
plt.close()
print("  → fig_gamma_R_gap.png")


# --- FIGURE 5: Two mechanisms diagram ---
fig, ax = plt.subplots(figsize=(10, 6))

for r in results:
    name, fam = r[0], r[1]
    f_Lg, f_Lp = r[9], r[10]
    obs = r[6]
    dom = r[14]

    if obs == "none":
        color = '#2ecc71'; marker = 's'
    else:
        color = '#e74c3c' if dom == "logP" else '#3498db'; marker = 'o'

    ax.plot(f_Lg, f_Lp, marker, color=color, markersize=9,
            markeredgecolor='white', markeredgewidth=0.5)
    if name in ['Acetone', 'MEK', 'Cyclopentanone', 'Cyclohexanone',
                'DMF', 'NMP', 'Pyridine', 'Dioxane', 'Nitromethane']:
        ax.annotate(name, (f_Lg, f_Lp), textcoords="offset points",
                    xytext=(6, 4), fontsize=7)

# Diagonal: where both agree
diag = np.linspace(-1, 3, 100)
ax.plot(diag, diag, 'k:', alpha=0.3, label='f_L(γ∞) = f_L(logP)')

ax.axhline(y=0, color='grey', alpha=0.3)
ax.axvline(x=0, color='grey', alpha=0.3)
ax.set_xlabel('f_L(γ∞) = ln(γ∞) / C')
ax.set_ylabel('f_L(logP) = logP / C\'')
ax.set_title('Two disruption mechanisms\n'
             'Above diagonal: logP dominates | Below: γ∞ dominates')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.2)
plt.tight_layout()
plt.savefig('fig_gamma_two_mechanisms.png')
plt.close()
print("  → fig_gamma_two_mechanisms.png")


# ============================================================
# SUMMARY
# ============================================================
n_homo = sum(1 for r in results if r[7] == "homo")
n_hetero = sum(1 for r in results if r[7] == "hetero")
n_none = sum(1 for r in results if r[6] == "none")
n_pos = sum(1 for r in results if r[6] == "positive")
n_logp_needed = sum(1 for r in results if r[14] == "logP" and r[6] == "positive"
                    and r[9] < r[10])  # cases where logP > γ∞ AND matters

print(f"\n{'=' * 80}")
print(f"SUMMARY")
print(f"{'=' * 80}")
print(f"  Total systems:        {total}")
print(f"  Homogeneous:          {n_homo}")
print(f"  Heterogeneous:        {n_hetero}")
print(f"  Predicted 'none':     {n_none}")
print(f"  Predicted 'positive': {n_pos}")
print(f"  Score:                {ok}/{total} ({100*ok/total:.1f}%)")
print(f"  Cases needing logP:   {n_hetero} (cyclopentanone, cyclohexanone)")
print(f"  C (γ∞):               {C_GAMMA:.3f}")
print(f"  C' (logP):            {C_LOGP:.2f}")
print(f"  R gap:                [{R_max_none:.3f}, {R_min_pos:.3f}]")
print(f"  R gap width:          {R_min_pos - R_max_none:.3f}")
