# GERT-Azeotrope
Scripts relacionados ao paper  Liquids Predicts Azeotrope 
**Title:** Crystalline Order in Liquids Predicts Azeotrope Formation and Emulsion Stability: A GERT Network Classification

**Target:** PCCP (Physical Chemistry Chemical Physics)

**Author:** V. P. Dutra

---

## Summary

The paper applies the GERT cohesive fraction φ = f_M/(f_M + f_L) = 1/2 to predict azeotrope formation in binary liquid mixtures. Two levels of prediction:

- **Level 1 (Screening):** Six rules using pure-component parameters only (α, β, δd, δp, Tb, chemical family). Score: 24/24.
- **Level 2 (Complete formula):** φ_int = f_M/(f_M + f_L) with f_L = max(ln(γ∞)/C, logP/C'). Score: 44/44.

Constants: C = 3.73 (calibrated on acetone/MEK gap), C' = 0.48 (heterogeneous correction). Threshold: φ = 1/2 (universal GERT).

---

## File Inventory

### LaTeX

| File                         | Description                                                |
| ---------------------------- | ---------------------------------------------------------- |
| `GERT_Azeotrope_UNIFIED.tex` | Complete unified paper (12 sections, ~840 lines)           |
| `gert_azeotrope.bib`         | Bibliography (TO DO: add Brouwer, Sander, Ullmann entries) |

### Python Scripts

| Script                       | Function                                      |
| ---------------------------- | --------------------------------------------- |
| `gert_azeotrope_final.py`    | Level 1: 24-system screening rules + Figs 1-3 |
| `gert_azeotrope_gamma.py`    | Level 2: 44-system γ∞ formula + Figs 4-8      |
| `gert_azeotrope_test.py`     | dφ/dx = 0 with W from H^E (11 systems)        |
| `gert_azeotrope_type02.py`   | Type 0+2 polarity-corrected f_int             |
| `gert_azeotrope_combined.py` | Network classification + dispersion           |

### Figures (8 total)

| Figure | File                              | Section | Description                                     |
| ------ | --------------------------------- | ------- | ----------------------------------------------- |
| Fig 1  | `fig1_phi_profiles.png`           | §2      | φ(x) profiles for 4 representative systems      |
| Fig 2  | `fig2_network_classification.png` | §3      | Network types 0-3 diagram                       |
| Fig 3  | `fig3_type02_polarity.png`        | §4      | Type 0+2 f_int vs δp                            |
| Fig 4  | `fig_gamma_two_mechanisms.png`    | §6      | f_L(γ∞) vs f_L(logP): two disruption mechanisms |
| Fig 5  | `fig_gamma_homologous_series.png` | §7      | CH₂ increment: γ∞ vs carbon number (3 series)   |
| Fig 6  | `fig_gamma_fM_vs_fL.png`          | §8      | GERT integration diagram: f_M vs f_L            |
| Fig 7  | `fig_gamma_phi_all_systems.png`   | §8      | Bar chart of φ for all 44 systems               |
| Fig 8  | `fig_gamma_R_gap.png`             | §8      | Separation score R gap                          |

### Data

| File                            | Description                                     |
| ------------------------------- | ----------------------------------------------- |
| `gert_phi_complete_type02.xlsx` | All 39 Type 0+2 systems with f_M, γ∞, φ         |
| `gert_water_gamma_series.xlsx`  | Brouwer-derived γ∞ data (Veronica's extraction) |

---

## Paper Structure

```
§1  Introduction (two readers, methodological inversion)
§2  The Framework (φ = f_M/(f_M+f_L), Trouton, azeotrope criterion)
§3  Network Classification (Types 0-3, σ = αβ)
§4  Screening Rules 1-6 (zero mixture data)
§5  Validation I: 24/24 (Level 1, screening)
§6  Complete Water-Integration Criterion (Level 2)
    §6.1  Screening rules as one-sided estimates
    §6.2  γ∞ as f_L (Eq. f_L = ln(γ∞)/C)
    §6.3  Two mechanisms (γ∞ local + logP bulk)
    §6.4  Complete formula (Eq. φ = f_M/(f_M+f_L))
    §6.5  Calibration (C = 3.73, C' = 0.48)
    §6.6  Unification of Rules 4 and 5
§7  Homologous Series
    §7.1  Ketones (C3→C4→C5, transition at C3/C4)
    §7.2  Esters (MeOAc→BuOAc, 4 members)
    §7.3  Ethers (THF→di-n-Bu, geometry + volume)
    §7.4  Amides (NMP/DMF/DMAc, γ∞ < 1, integration)
    §7.5  CH₂ increment as quantitative law
§8  Validation II: 44/44 (Level 2, complete)
§9  Emulsion Extension (5/5, R_crit)
§10 Discussion (two levels, scope, limitations)
§11 Bridge: Cosmos to Flask
§12 Conclusions (10 items)
```

---

## Key Results

| Level     | Systems | Score        | Data Required       |
| --------- | ------- | ------------ | ------------------- |
| Screening | 24      | 24/24 (100%) | Pure-component only |
| Complete  | 44      | 44/44 (100%) | + γ∞ and logP       |
| Emulsion  | 5       | 5/5 (100%)   | Interfacial tension |

## Key Equations

```
φ_int = f_M / (f_M + f_L) = 1/2

f_M (Type 0+2) = (β_sol/β_w)(δp_sol/δp_w)
f_M (Type 1+2) = √(min(α_sol/α_w,1) × min(β_sol/β_w,1))

f_L = max(ln(γ∞)/C, logP/C')
C = 3.73    C' = 0.48

R_crit = 3γV_m/(R_gas T)  [emulsion]
```

## Parameter Sources

| Parameter                | Source                                             |
| ------------------------ | -------------------------------------------------- |
| α, β (non-amphiprotic)   | Kamlet et al. 1983, J. Org. Chem. 48, 2877         |
| α, β (amphiprotic + CS₂) | Marcus 1998, The Properties of Solvents, Table 4.3 |
| δd, δp                   | Hansen 2007, HSP Handbook, Table A.1               |
| γ∞                       | Brouwer et al. (DDB open compilation), 298.15 K    |
| γ∞ (DMSO)                | NIST/Sander Henry constant, 55.51/(k_H × P_sat)    |
| γ∞ (pyridine)            | NIST/Sander Henry + Brouwer (31.8, unstable fit)   |
| logP                     | PubChem/XLogP3                                     |
| Azeotrope data           | DECHEMA, CRC Handbook, Ullmann's Encyclopedia      |

---

## 
