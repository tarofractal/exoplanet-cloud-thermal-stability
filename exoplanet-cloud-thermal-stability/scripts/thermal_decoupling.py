# -*- coding: utf-8 -*-
"""
Thermal-relaxation decoupling figure (recovered 2026-04-06 version)
See docs/RECOVERY_STATUS.md for differences from the published paper.
------------------------------------
修正点:
- P_MIN_BAR を 1e-12 に拡張 (プロットの 1e-10 をカバー)
- TP_MAX を 20000K に拡張 (Fe等の高温対策)
- 分割数を増やして線のガタつきを抑制
"""

import os
import math
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.constants import h, c, k as k_B, pi
from scipy.optimize import brentq
from scipy.interpolate import PchipInterpolator
from scipy.integrate import trapezoid

# ==========================================================
# File names (base names only)
# ==========================================================
QABS_BASENAMES = {
    "Al2O3":   "Al2O3_logR700_Qabs.csv",
    "Fe":      "Fe_logR700_Qabs.csv",
    "KCl":     "KCl_logR700_Qabs.csv",
    "Mg2SiO4": "Mg2SiO4_logR700_Qabs.csv",
    "MnS":     "MnS_logR700_Qabs.csv",
    "Na2S":    "Na2S_logR700_Qabs.csv",
    "SiO2":    "SiO2_logR700_Qabs.csv",
    "ZnS":     "ZnS_logR700_Qabs.csv",
}

ALL_MATERIALS_8 = ["Al2O3", "Fe", "KCl", "Mg2SiO4", "MnS", "Na2S", "SiO2", "ZnS"]
SAVE_DIR = "./thermal_relaxation_decoupling_outputs"

# ==========================================================
# User settings
# ==========================================================
# Physics
ALPHA_ACC = 0.5
MU = 2.3
AMU = 1.66053906660e-27
M_GAS = MU * AMU

# Decoupling criterion
DELTA_FRAC = 0.05              # 5%
DECOUPLING_REFERENCE = "Tgas"  # "Tgas" or "Trad"
TGAS_FACTOR = 1.0              # Tgas = TGAS_FACTOR * Teq

# Pressure grid [bar] (重要: プロット範囲より広く設定)
P_MIN_BAR = 1e-12              # 1e-8から1e-12へ拡張
P_MAX_BAR = 1e2                # 1e1から1e2へ拡張
N_PRESSURE = 300               # 滑らかさのために増加

# Temperature solving
TP_MIN = 10.0
TP_MAX = 20000.0               # 10000.0から20000.0へ拡張 (Fe対策)
N_T_SOLVE_SCAN = 300           # 探索精度向上のために増加
EMISSION_T_GRID = np.logspace(np.log10(TP_MIN), np.log10(TP_MAX), 500)

# Panel 1: P_dec vs Teq
TEFF_PANEL1 = 5000.0
A_REF_UM = 0.1
TEQ_GRID = np.linspace(300.0, 2200.0, 60)

# Panel 2: P_dec vs a
TEFF_PANEL2 = 5000.0
TEQ_REF_K = 1000.0

# Axis controls
PANEL1_XLIM = (300.0, 2200.0)
PANEL1_YLIM = (1e-10, 1e1)       # プロット表示範囲
PANEL2_XLIM = (1e-2, 1e2)
PANEL2_YLIM = (1e-8, 1e1)

MATERIAL_COLORS = {
    "Al2O3": "tab:blue", "KCl": "tab:orange", "Mg2SiO4": "tab:green",
    "Fe": "tab:red", "ZnS": "tab:purple", "MnS": "#8c564b",
    "SiO2": "#e377c2", "Na2S": "tab:gray",
}

# ==========================================================
# Functions
# ==========================================================

def find_file_anywhere(basename, data_dir="."):
    path = os.path.join(data_dir, basename)
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Could not find {path}; provide --data-dir with the Qabs CSV files.")
    return os.path.abspath(path)

def planck_lambda(lam_m, T):
    x = h * c / (lam_m * k_B * T)
    x = np.clip(x, 1e-12, 700.0)
    return (2.0 * h * c**2) / (lam_m**5) / np.expm1(x)

def fill_nonfinite_qabs(lam_um, q):
    finite = np.isfinite(lam_um) & (lam_um > 0) & np.isfinite(q)
    if finite.sum() < 2: return np.clip(q, 0.0, None)
    x_good, y_good = np.log10(lam_um[finite]), q[finite]
    q_filled = np.interp(np.log10(lam_um), x_good, y_good, left=y_good[0], right=y_good[-1])
    return np.clip(q_filled, 0.0, None)

def load_qabs_table_exact(filepath):
    df = pd.read_csv(filepath)
    col_lam = "Wavelength (μm)"
    lam_um = df[col_lam].to_numpy(dtype=float)
    mask = (lam_um > 0) & np.isfinite(lam_um)
    lam_um = lam_um[mask]
    qabs_dict = {}
    for col in df.columns:
        if "Qabs_" in col:
            r_str = col.replace("Qabs_", "").replace("μm", "").replace("um", "")
            r = float(r_str)
            q = fill_nonfinite_qabs(lam_um, df[col].to_numpy()[mask])
            qabs_dict[r] = q
    order = np.argsort(lam_um)
    return lam_um[order] * 1e-6, {r: q[order] for r, q in qabs_dict.items()}

class EmissionInterpolator:
    def __init__(self, lam_m, qabs, t_grid):
        e_grid = [4.0 * trapezoid(qabs * planck_lambda(lam_m, T), lam_m) for T in t_grid]
        self.t_grid, self.e_grid = t_grid, np.asarray(e_grid)
        mask = (self.e_grid > 0) & np.isfinite(self.e_grid)
        self._interp = PchipInterpolator(np.log(self.t_grid[mask]), np.log(self.e_grid[mask]), extrapolate=True)
    def __call__(self, T):
        return np.exp(self._interp(np.log(np.clip(T, self.t_grid.min(), self.t_grid.max()))))

def solve_tp(Hstar, Tgas, Pbar, emission_func, collision=True):
    coeff = 0.0
    if collision:
        Ppa = Pbar * 1.0e5
        vth = math.sqrt(8.0 * k_B * Tgas / (pi * M_GAS))
        coeff = 2.0 * ALPHA_ACC * vth * (Ppa / (k_B * Tgas)) * k_B
    def f(Tp): return Hstar - emission_func(Tp) + coeff * (Tgas - Tp)
    Ts = np.logspace(np.log10(TP_MIN), np.log10(TP_MAX), N_T_SOLVE_SCAN)
    fvals = np.array([f(T) for T in Ts])
    sgn = np.sign(fvals)
    for i in range(len(Ts)-1):
        if sgn[i]*sgn[i+1] <= 0:
            return brentq(f, Ts[i], Ts[i+1])
    return np.nan

def find_threshold_pressure_for_case(lam_m, qabs, Teff, Teq, Pscan):
    Tgas = TGAS_FACTOR * Teq
    Hstar = 4.0 * (Teq/Teff)**4 * trapezoid(qabs * planck_lambda(lam_m, Teff), lam_m)
    em_func = EmissionInterpolator(lam_m, qabs, EMISSION_T_GRID)

    # Check if decoupling is even possible at low-P limit
    Trad = solve_tp(Hstar, Tgas, 0.0, em_func, collision=False)
    if not np.isfinite(Trad) or abs(Trad - Tgas)/Tgas < DELTA_FRAC:
        return np.nan # Too transparent to decouple 5%

    metrics = []
    for P in Pscan:
        Tp = solve_tp(Hstar, Tgas, P, em_func)
        metrics.append(abs(Tp - Tgas)/Tgas)

    metrics = np.asarray(metrics)
    cond = np.isfinite(metrics) & (metrics >= DELTA_FRAC)
    if not np.any(cond): return np.nan

    idx = np.where(cond)[0][0]
    if idx == 0: return Pscan[0]

    # Log-linear interpolation for smooth lines
    logP1, logP2 = np.log10(Pscan[idx-1]), np.log10(Pscan[idx])
    m1, m2 = metrics[idx-1], metrics[idx]
    frac = (DELTA_FRAC - m1) / (m2 - m1)
    return 10.0**(logP1 + frac * (logP2 - logP1))

# ==========================================================
# Main Calculations
# ==========================================================

def main(argv=None):
    parser = argparse.ArgumentParser(description='Recovered thermal-decoupling calculation')
    parser.add_argument('--data-dir', default='.', help='Directory containing the eight Qabs CSV files')
    parser.add_argument('--output-dir', default=SAVE_DIR)
    args = parser.parse_args(argv)
    os.makedirs(args.output_dir, exist_ok=True)
    species_data = {}
    for mat in ALL_MATERIALS_8:
        path = find_file_anywhere(QABS_BASENAMES[mat], args.data_dir)
        lam_m, q_dict = load_qabs_table_exact(path)
        species_data[mat] = {"lam_m": lam_m, "qabs": q_dict, "radii": sorted(q_dict.keys())}

    Pscan = np.logspace(np.log10(P_MAX_BAR), np.log10(P_MIN_BAR), N_PRESSURE)

    # Panel 1
    p1_rows = []
    for mat in ALL_MATERIALS_8:
        print(f"Calculating Panel 1: {mat}...")
        r_use = species_data[mat]["radii"][np.argmin(np.abs(np.array(species_data[mat]["radii"]) - A_REF_UM))]
        qabs = species_data[mat]["qabs"][r_use]
        for Teq in TEQ_GRID:
            Pdec = find_threshold_pressure_for_case(species_data[mat]["lam_m"], qabs, TEFF_PANEL1, Teq, Pscan)
            p1_rows.append({"Material": mat, "Teq": Teq, "Pdec": Pdec})

    # Panel 2
    p2_rows = []
    for mat in ALL_MATERIALS_8:
        print(f"Calculating Panel 2: {mat}...")
        for r in species_data[mat]["radii"]:
            Pdec = find_threshold_pressure_for_case(species_data[mat]["lam_m"], species_data[mat]["qabs"][r], TEFF_PANEL2, TEQ_REF_K, Pscan)
            p2_rows.append({"Material": mat, "Radius": r, "Pdec": Pdec})

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.6), constrained_layout=True)

    for mat in ALL_MATERIALS_8:
        d1 = pd.DataFrame(p1_rows); d1 = d1[d1["Material"]==mat]
        axes[0].plot(d1["Teq"], d1["Pdec"], color=MATERIAL_COLORS[mat], lw=2, label=mat)

        d2 = pd.DataFrame(p2_rows); d2 = d2[d2["Material"]==mat]
        axes[1].plot(d2["Radius"], d2["Pdec"], color=MATERIAL_COLORS[mat], lw=2)

    for ax, xlim, ylim, xscale in zip(axes, [PANEL1_XLIM, PANEL2_XLIM], [PANEL1_YLIM, PANEL2_YLIM], ["linear", "log"]):
        ax.set_xscale(xscale); ax.set_yscale("log")
        ax.set_xlim(xlim); ax.set_ylim(ylim)
        ax.invert_yaxis(); ax.set_ylabel("Decoupling Pressure [bar]")

    axes[0].set_xlabel("Equilibrium Temperature [K]")
    axes[1].set_xlabel("Particle Radius [μm]")
    fig.legend(loc="center right", bbox_to_anchor=(1.12, 0.5), title="Material")

    plt.savefig(os.path.join(args.output_dir, "decoupling_full.pdf"), bbox_inches="tight")
    print("Done. Saved to 'decoupling_full.pdf'")

if __name__ == "__main__":
    main()