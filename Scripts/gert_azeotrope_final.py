"""
GERT Azeotrope — Final consolidated script
=============================================
24 systems, 6 rules, 3 figures, corrected f_L (Trouton)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# =====================================================================
# DATA
# =====================================================================
# (α, β, σ, type, δd, δp, family, ΔHvap kJ/mol, Tb K, Vm mL/mol)

D = {
 'water':        (1.17,0.47,0.550,2, 15.5,16.0,'water',   40.7,373.15,18.1),
 'ethanol':      (0.86,0.75,0.645,1, 15.8, 8.8,'alcohol', 38.6,351.44,58.4),
 'methanol':     (0.93,0.62,0.577,1, 15.1,12.3,'alcohol', 35.2,337.70,40.7),
 'isopropanol':  (0.76,0.84,0.638,1, 15.8, 6.1,'alcohol', 39.9,355.41,76.9),
 'n_butanol':    (0.84,0.84,0.706,1, 16.0, 5.7,'alcohol', 43.3,390.88,91.5),
 'acetic_acid':  (1.12,0.45,0.504,1, 14.5, 8.0,'acid',    23.7,391.05,57.5),
 'acetone':      (0.08,0.48,0.038,0, 15.5,10.4,'ketone',  29.1,329.20,73.5),
 'chloroform':   (0.44,0.00,0.000,0, 17.8, 3.1,'haloalk', 29.2,334.33,80.7),
 'benzene':      (0.00,0.10,0.000,0, 18.4, 0.0,'aromatic',30.7,353.25,89.4),
 'toluene':      (0.00,0.11,0.000,0, 18.0, 1.4,'aromatic',33.2,383.78,106.9),
 'hexane':       (0.00,0.00,0.000,0, 14.9, 0.0,'alkane',  28.9,341.88,131.6),
 'heptane':      (0.00,0.00,0.000,0, 15.3, 0.0,'alkane',  31.8,371.57,147.5),
 'CS2':          (0.00,0.07,0.000,0, 20.5, 0.0,'other',   26.7,319.38,60.0),
 'cyclohexane':  (0.00,0.00,0.000,0, 16.8, 0.0,'alkane',  29.9,353.87,108.8),
 'ethyl_acetate':(0.00,0.45,0.000,0, 15.8, 5.3,'ester',   31.9,350.21,98.5),
 'diethyl_ether':(0.00,0.47,0.000,0, 14.5, 2.9,'ether',   26.5,307.58,104.7),
 'THF':          (0.00,0.55,0.000,0, 16.8, 5.7,'ether',   29.81,339.12,81.7),
 'DCM':          (0.30,0.00,0.000,0, 18.2, 6.3,'haloalk', 28.1,312.90,63.9),
}

def predict(mA, mB):
    """Full 6-rule classification."""
    aA,bA,sA,tA,dA,pA,fA,hA,TbA,vA = D[mA]
    aB,bB,sB,tB,dB,pB,fB,hB,TbB,vB = D[mB]
    # Type 3
    Gc = aA*bB + aB*bA; Gs = sA + sB; C = Gc - Gs
    if C > 0 and sA < 0.05 and sB < 0.05:
        return "negative", 3, f"C={C:.3f}>0 → emergent"
    # Type 0+0
    if tA == 0 and tB == 0:
        Dd = np.sqrt((dA-dB)**2+(pA-pB)**2)
        F = 1 if fA == fB else 0
        if F == 1: return "none", 0, f"Same family"
        elif Dd > 4: return "positive", 0, f"D={Dd:.1f}>4"
        else: return "none", 0, f"D={Dd:.1f}≤4"
    # Network holder
    if tA >= tB:
        na,nb,nt,nf,np_,nm = aA,bA,tA,fA,pA,mA
        sa,sb,st,sf,sp,sm = aB,bB,tB,fB,pB,mB
    else:
        na,nb,nt,nf,np_,nm = aB,bB,tB,fB,pB,mB
        sa,sb,st,sf,sp,sm = aA,bA,tA,fA,pA,mA
    # Type 0+2 (Rule 5)
    if nt == 2 and st == 0:
        br = sb/nb if nb > 0 else 0
        pr = sp/np_ if np_ > 0 else 0
        fi = br * pr
        if fi > 0.50: return "none", 2, f"f_int={fi:.3f}>0.50"
        else: return "positive", 2, f"f_int={fi:.3f}≤0.50"
    # Type 1+2 (Rule 4)
    if nt == 2 and st == 1:
        fd = min(sa/na,1) if na > 0 else 0
        fa = min(sb/nb,1) if nb > 0 else 0
        fi = np.sqrt(fd*fa) if fd > 0 and fa > 0 else 0
        if fi > 0.90: return "none", 2, f"f_int={fi:.2f}>0.90"
        else: return "positive", 2, f"f_int={fi:.2f}≤0.90"
    # Type 1+1
    if nt == 1 and st == 1:
        if nf == sf: return "none", 1, f"Same class ({nf})"
        else: return "positive", 1, f"Diff ({nf}+{sf})"
    # Type 0+1
    if nt == 1 and st == 0:
        return "positive", 1, f"{sm}(0) disrupts {nm}(1)"
    return "unknown", -1, "?"

def phi_profile(mA, mB, T, W=None, n=500):
    """Compute φ(x) with corrected f_L = T·ΔS_vap (Trouton)."""
    A, B = D[mA], D[mB]
    hA, TbA = A[7], A[8]
    hB, TbB = B[7], B[8]
    dSA = hA / TbA  # Trouton
    dSB = hB / TbB
    x = np.linspace(0.001, 0.999, n)
    if W is None: W = 0
    fM = x*hA + (1-x)*hB + W*x*(1-x)
    fL = T * (x*dSA + (1-x)*dSB)
    phi = fM / (fM + fL)
    # Find extrema
    dp = np.diff(phi)/np.diff(x)
    ext = []
    for i in range(len(dp)-1):
        if dp[i]*dp[i+1] < 0:
            xs = x[i+1] - dp[i+1]*(x[i+2]-x[i+1])/(dp[i+1]-dp[i])
            ps = np.interp(xs, x, phi)
            tp = "negative" if dp[i] > 0 else "positive"
            ext.append((xs, ps, tp))
    return x, phi, ext

# =====================================================================
# 24 SYSTEMS (corrected compositions)
# =====================================================================
systems = [
    ("MeOH–Benzene",    'methanol',     'benzene',       0.615,"positive"),
    ("Benzene–EtOH",    'benzene',      'ethanol',       0.448,"positive"),
    ("EtOAc–EtOH",      'ethyl_acetate','ethanol',       0.460,"positive"),
    ("Cyclohex–EtOH",   'cyclohexane',  'ethanol',       0.430,"positive"),
    ("MeOH–Cyclohex",   'methanol',     'cyclohexane',   0.510,"positive"),
    ("DCM–MeOH",        'DCM',          'methanol',      0.850,"positive"),
    ("THF–EtOH",        'THF',          'ethanol',       0.908,"positive"),
    ("Hexane–EtOH",     'hexane',       'ethanol',       0.332,"positive"),
    ("Acetone–MeOH",    'acetone',      'methanol',      0.790,"positive"),
    ("CS₂–Acetone",     'CS2',          'acetone',       0.610,"positive"),
    ("EtOH–H₂O",        'ethanol',      'water',         0.894,"positive"),
    ("iPrOH–H₂O",       'isopropanol',  'water',         0.684,"positive"),
    ("nBuOH–H₂O",       'n_butanol',    'water',         0.237,"positive"),
    ("Acetone–CHCl₃",   'acetone',      'chloroform',    0.357,"negative"),
    ("CHCl₃–Et₂O",      'chloroform',   'diethyl_ether', 0.850,"negative"),
    ("MeOH–EtOH",       'methanol',     'ethanol',       None, "none"),
    ("EtOH–iPrOH",      'ethanol',      'isopropanol',   None, "none"),
    ("Benzene–Toluene",  'benzene',      'toluene',       None, "none"),
    ("Hexane–Heptane",   'hexane',       'heptane',       None, "none"),
    ("AcOH–H₂O",        'acetic_acid',  'water',         None, "none"),
    ("Acetone–H₂O",     'acetone',      'water',         None, "none"),
    ("THF–H₂O",         'THF',          'water',         0.830,"positive"),
    ("Et₂O–H₂O",        'diethyl_ether','water',         0.950,"positive"),
    ("EtOAc–H₂O",       'ethyl_acetate','water',         0.677,"positive"),
]

# =====================================================================
# RUN VALIDATION
# =====================================================================
print("="*80)
print("GERT AZEOTROPE — FINAL VALIDATION (24 systems)")
print("="*80)
print()
n_ex=0; n_tp=0
for i,(nm,cA,cB,xo,to) in enumerate(systems,1):
    pr,pt,rs = predict(cA,cB)
    ex_ok = (pr!="none") == (xo is not None)
    tp_ok = pr == to
    if ex_ok: n_ex+=1
    if tp_ok: n_tp+=1
    m = "✓✓" if (ex_ok and tp_ok) else "✗✗"
    print(f"{i:>2} {m} {nm:<18} {pr:>8} {to:>8}  {rs}")
print(f"\nEXISTENCE: {n_ex}/24  TYPE: {n_tp}/24")

# =====================================================================
# FIGURE 1: φ(x) profiles — 4 representative systems
# =====================================================================
fig1, axes = plt.subplots(2, 2, figsize=(12, 10))

# Panel (a): Benzene-EtOH (Type 0+1, positive)
ax = axes[0,0]
x, phi, ext = phi_profile('benzene','ethanol', 341, W=-6.0)
ax.plot(x, phi, 'b-', linewidth=2.5)
if ext: ax.plot(ext[0][0], ext[0][1], 'ro', markersize=10, zorder=5)
ax.axvline(x=0.448, color='green', linestyle='--', linewidth=1.5, alpha=0.6, label='Observed')
ax.set_title('(a) Benzene–EtOH (Type 0+1, positive)', fontweight='bold', fontsize=11)
ax.set_xlabel('$x$ (benzene mol fraction)', fontsize=10)
ax.set_ylabel('$\\varphi$', fontsize=12)
ax.legend(fontsize=9)
ax.annotate('minimum → positive\nazeotrope', xy=(ext[0][0],ext[0][1]),
            xytext=(0.15,ext[0][1]+0.004), fontsize=9, color='red',
            arrowprops=dict(arrowstyle='->',color='red'))
ax.grid(True, alpha=0.15)

# Panel (b): Acetone-CHCl₃ (Type 3, negative)
ax = axes[0,1]
x, phi, ext = phi_profile('acetone','chloroform', 337, W=+6.0)
ax.plot(x, phi, 'b-', linewidth=2.5)
if ext: ax.plot(ext[0][0], ext[0][1], 'ro', markersize=10, zorder=5)
ax.axvline(x=0.357, color='green', linestyle='--', linewidth=1.5, alpha=0.6, label='Observed')
ax.set_title('(b) Acetone–CHCl₃ (Type 3, negative)', fontweight='bold', fontsize=11)
ax.set_xlabel('$x$ (acetone mol fraction)', fontsize=10)
ax.set_ylabel('$\\varphi$', fontsize=12)
ax.legend(fontsize=9)
ax.annotate('maximum → negative\nazeotrope', xy=(ext[0][0],ext[0][1]),
            xytext=(0.55,ext[0][1]-0.004), fontsize=9, color='red',
            arrowprops=dict(arrowstyle='->',color='red'))
ax.grid(True, alpha=0.15)

# Panel (c): EtOH-H₂O (Type 1+2, positive, asymmetric)
ax = axes[1,0]
x, phi, ext = phi_profile('ethanol','water', 351, W=-3.5)
# Add asymmetric W₂ manually
A, B = D['ethanol'], D['water']
W1, W2 = -3.5, 1.0
hA,TbA,hB,TbB = A[7],A[8],B[7],B[8]
dSA,dSB = hA/TbA, hB/TbB
x2 = np.linspace(0.001,0.999,500)
fM2 = x2*hA+(1-x2)*hB+x2*(1-x2)*(W1+W2*(2*x2-1))
fL2 = 351*(x2*dSA+(1-x2)*dSB)
phi2 = fM2/(fM2+fL2)
ax.plot(x2, phi2, 'b-', linewidth=2.5)
# Find min
imin = np.argmin(phi2[10:-10]) + 10
ax.plot(x2[imin], phi2[imin], 'ro', markersize=10, zorder=5)
ax.axvline(x=0.894, color='green', linestyle='--', linewidth=1.5, alpha=0.6, label='Observed (0.894)')
ax.set_title('(c) EtOH–H₂O (Type 1+2, asymmetric)', fontweight='bold', fontsize=11)
ax.set_xlabel('$x$ (ethanol mol fraction)', fontsize=10)
ax.set_ylabel('$\\varphi$', fontsize=12)
ax.legend(fontsize=9)
ax.annotate(f'$x^*$={x2[imin]:.3f}', xy=(x2[imin],phi2[imin]),
            xytext=(0.5,phi2[imin]+0.004), fontsize=9, color='red',
            arrowprops=dict(arrowstyle='->',color='red'))
ax.grid(True, alpha=0.15)

# Panel (d): Hexane-Heptane (Type 0+0, no azeotrope)
ax = axes[1,1]
x, phi, ext = phi_profile('hexane','heptane', 355, W=-0.1)
ax.plot(x, phi, 'b-', linewidth=2.5)
ax.set_title('(d) Hexane–Heptane (Type 0+0, no azeotrope)', fontweight='bold', fontsize=11)
ax.set_xlabel('$x$ (hexane mol fraction)', fontsize=10)
ax.set_ylabel('$\\varphi$', fontsize=12)
ax.annotate('monotonic → no\nextremum → no azeotrope', xy=(0.5,phi[250]),
            xytext=(0.15,phi[250]+0.003), fontsize=9, color='blue')
ax.grid(True, alpha=0.15)

plt.tight_layout()
plt.savefig('/mnt/user-data/outputs/fig1_phi_profiles.png', dpi=300, bbox_inches='tight')
plt.close()
print("\nFig. 1 saved: fig1_phi_profiles.png")

# =====================================================================
# FIGURE 2: Network classification diagram
# =====================================================================
fig2, ax = plt.subplots(1, 1, figsize=(14, 7))
ax.set_xlim(-0.5, 4.5)
ax.set_ylim(-0.5, 3.5)
ax.axis('off')

colors = {'0':'#E8E8E8', '1':'#FFD700', '2':'#4169E1', '3':'#FF6347'}
titles = ['Type 0\nNo network', 'Type 1\nChain', 'Type 2\nTetrahedral (Tₐ)', 'Type 3\nComplementary']
symms = ['—', 'C∞ᵥ', 'Tₐ ⊂ Oₕ', 'Emergent']
examples = [
    'hexane, benzene,\nacetone, CHCl₃,\nTHF, Et₂O, EtOAc',
    'ethanol, methanol,\nisopropanol, butanol,\nacetic acid',
    'water',
    'acetone + CHCl₃\n(donor meets acceptor)\n→ exists ONLY in mixture'
]
criteria = ['σ = α×β < 0.05', 'σ = α×β > 0.05', 'σ > 0.40 + 3D\ncooperative', 'Γcross > Γself\nboth σ < 0.05']
azeotype = ['Depends on D\n(Hansen distance)', 'Disrupted by\nType 0 → positive', 'f_int determines\n(Rules 4, 5)', 'NEGATIVE\n(emergent network)']

for i in range(4):
    x0 = i * 1.1
    # Box
    rect = plt.Rectangle((x0, 0.3), 0.95, 2.9, facecolor=colors[str(i)],
                          edgecolor='black', linewidth=2, alpha=0.7)
    ax.add_patch(rect)
    # Title
    ax.text(x0+0.475, 3.0, titles[i], ha='center', va='top',
            fontsize=11, fontweight='bold')
    # Symmetry
    ax.text(x0+0.475, 2.45, symms[i], ha='center', va='center',
            fontsize=10, style='italic', color='#333')
    # Criterion
    ax.text(x0+0.475, 2.0, criteria[i], ha='center', va='center',
            fontsize=8, color='#555',
            bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8))
    # Examples
    ax.text(x0+0.475, 1.35, examples[i], ha='center', va='center', fontsize=8)
    # Azeotrope
    ax.text(x0+0.475, 0.6, azeotype[i], ha='center', va='center',
            fontsize=8, fontweight='bold', color='darkred')

ax.set_title('Network Classification of Liquid Molecules', fontsize=14, fontweight='bold', pad=20)

plt.tight_layout()
plt.savefig('/mnt/user-data/outputs/fig2_network_classification.png', dpi=300, bbox_inches='tight')
plt.close()
print("Fig. 2 saved: fig2_network_classification.png")

# =====================================================================
# FIGURE 3: Type 0+2 — f_int vs polarity
# =====================================================================
fig3, ax = plt.subplots(1, 1, figsize=(9, 6))

# All Type 0 molecules paired with water
type0_in_water = [
    ('acetone',      0.08, 0.48, 10.4, False, 'Acetone'),
    ('THF',          0.00, 0.55,  5.7, True,  'THF'),
    ('diethyl_ether',0.00, 0.47,  2.9, True,  'Et₂O'),
    ('ethyl_acetate',0.00, 0.45,  5.3, True,  'EtOAc'),
    ('hexane',       0.00, 0.00,  0.0, True,  'Hexane'),
    ('cyclohexane',  0.00, 0.00,  0.0, True,  'Cyclohex'),
    ('benzene',      0.00, 0.10,  0.0, True,  'Benzene'),
    ('DCM',          0.30, 0.00,  6.3, True,  'DCM'),
    ('CS2',          0.00, 0.07,  0.0, True,  'CS₂'),
    ('chloroform',   0.44, 0.00,  3.1, True,  'CHCl₃'),
]

b_w = D['water'][1]  # 0.47
p_w = D['water'][5]  # 16.0

for mol, a, b, dp, has_az, label in type0_in_water:
    fi = (b/b_w) * (dp/p_w) if b_w > 0 and p_w > 0 else 0
    color = 'red' if has_az else 'green'
    marker = 'o' if has_az else 's'
    ax.plot(dp, fi, marker, color=color, markersize=12, zorder=5)
    offset = (0.3, 0.02) if label != 'Acetone' else (0.3, -0.04)
    if label in ('Hexane','Cyclohex','Benzene','CS₂','CHCl₃'):
        offset = (0.3, 0.015)
    ax.annotate(label, xy=(dp, fi), xytext=(dp+offset[0], fi+offset[1]),
                fontsize=9, color=color, fontweight='bold')

# Threshold line
ax.axhline(y=0.50, color='black', linestyle='--', linewidth=2, alpha=0.7, label='Threshold $f_{int}^{(0+2)}$ = 0.50')
ax.fill_between([0,18], 0.50, 1.0, color='green', alpha=0.08, label='Integrates (no azeotrope)')
ax.fill_between([0,18], 0, 0.50, color='red', alpha=0.08, label='Disrupts (positive azeotrope)')

ax.set_xlabel('$\\delta_p$ of solute (MPa$^{1/2}$)', fontsize=12)
ax.set_ylabel('$f_{int}^{(0+2)} = (\\beta_{sol}/\\beta_w) \\times (\\delta_{p,sol}/\\delta_{p,w})$', fontsize=11)
ax.set_title('Type 0+2: Polarity-corrected integration into water Tₐ network', fontsize=12, fontweight='bold')
ax.set_xlim(-0.5, 14)
ax.set_ylim(-0.05, 0.75)
ax.legend(fontsize=9, loc='upper left')
ax.grid(True, alpha=0.15)

plt.tight_layout()
plt.savefig('/mnt/user-data/outputs/fig3_type02_polarity.png', dpi=300, bbox_inches='tight')
plt.close()
print("Fig. 3 saved: fig3_type02_polarity.png")

print("\n" + "="*60)
print("ALL OUTPUTS:")
print("  fig1_phi_profiles.png         — φ(x) for 4 representative systems")
print("  fig2_network_classification.png — Types 0/1/2/3 diagram")
print("  fig3_type02_polarity.png      — f_int vs δ_p for Type 0+2")
print("  gert_azeotrope_final.py       — this script (consolidated)")
print("="*60)

