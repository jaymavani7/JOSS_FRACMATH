"""Generate comparison plots: MATLAB vs Abaqus 3PB Gregoire."""
from pathlib import Path
import csv
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT = HERE

def read_csv(path, skip=1):
    with open(path) as f:
        for _ in range(skip):
            next(f)
        return np.array([[float(x) for x in line.split(',')] for line in f if line.strip()])

abq = read_csv(ROOT / 'abaqus_result' / 'load_cmod_data.csv')   # time, load, cmod
mat = read_csv(ROOT / 'Matlab_result' / 'matlab_load_cmod.csv') # cmod, load

t_a, P_a, c_a = abq[:, 0], abq[:, 1], abq[:, 2]
c_m, P_m = mat[:, 0], mat[:, 1]

# Align: drop Abaqus row 0 (t=0)
P_a2, c_a2 = P_a[1:], c_a[1:]
n = min(len(P_a2), len(P_m))
P_a2, c_a2 = P_a2[:n], c_a2[:n]
c_m, P_m = c_m[:n], P_m[:n]

# ------------------------------------------------------------------ Fig 1: overlay
fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(c_a2 * 1e3, P_a2, '-', lw=2.0, color='tab:blue',   label='Abaqus')
ax.plot(c_m  * 1e3, P_m,  '--', lw=1.5, color='tab:orange', label='MATLAB')
ax.set_xlabel('CMOD [μm]')
ax.set_ylabel('Load [N]')
ax.set_title('Load vs CMOD — MATLAB vs Abaqus')
ax.grid(alpha=0.3)
ax.legend()
fig.tight_layout()
fig.savefig(OUT / '01_load_cmod_overlay.png', dpi=200)
fig.savefig(OUT / '01_load_cmod_overlay.pdf')
plt.close(fig)

# ------------------------------------------------------------------ Fig 2: peak zoom
fig, ax = plt.subplots(figsize=(7, 5))
mask = (c_a2 * 1e3 < 60)
ax.plot(c_a2[mask] * 1e3, P_a2[mask], '-',  lw=2.0, color='tab:blue',   label='Abaqus')
ax.plot(c_m[mask]  * 1e3, P_m[mask],  '--', lw=1.5, color='tab:orange', label='MATLAB')
ax.axvline(25.948, color='gray', ls=':', alpha=0.6, label='CMOD @ peak = 25.95 μm')
ax.axhline(3891.44, color='gray', ls=':', alpha=0.6)
ax.set_xlabel('CMOD [μm]')
ax.set_ylabel('Load [N]')
ax.set_title('Peak region zoom (CMOD < 60 μm)')
ax.grid(alpha=0.3)
ax.legend()
fig.tight_layout()
fig.savefig(OUT / '02_peak_zoom.png', dpi=200)
fig.savefig(OUT / '02_peak_zoom.pdf')
plt.close(fig)

# ------------------------------------------------------------------ Fig 3: pointwise diffs
dP = P_a2 - P_m
dC = (c_a2 - c_m) * 1e3  # μm

fig, axes = plt.subplots(2, 1, figsize=(7, 6), sharex=True)
axes[0].plot(c_a2 * 1e3, dP, color='tab:red', lw=1)
axes[0].axhline(0, color='k', lw=0.5)
axes[0].set_ylabel('Δ Load [N]\n(Abaqus − MATLAB)')
axes[0].set_title(f'Pointwise difference  (max|ΔP|={np.max(np.abs(dP)):.2f} N, '
                  f'max|ΔCMOD|={np.max(np.abs(dC)):.2e} μm)')
axes[0].grid(alpha=0.3)

axes[1].plot(c_a2 * 1e3, dC, color='tab:purple', lw=1)
axes[1].axhline(0, color='k', lw=0.5)
axes[1].set_ylabel('Δ CMOD [μm]\n(Abaqus − MATLAB)')
axes[1].set_xlabel('CMOD [μm]')
axes[1].grid(alpha=0.3)

fig.tight_layout()
fig.savefig(OUT / '03_pointwise_diff.png', dpi=200)
fig.savefig(OUT / '03_pointwise_diff.pdf')
plt.close(fig)

# ------------------------------------------------------------------ Fig 4: time bars
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

solvers = ['MATLAB', 'Abaqus']
wall    = [239.46, 669.0]
colors  = ['tab:orange', 'tab:blue']

bars = ax1.bar(solvers, wall, color=colors, edgecolor='k')
for b, v in zip(bars, wall):
    ax1.text(b.get_x() + b.get_width() / 2, v, f'{v:.1f} s',
             ha='center', va='bottom', fontsize=11)
ax1.set_ylabel('Wall-clock time [s]')
ax1.set_title(f'Wall-clock  (MATLAB ≈ {wall[1]/wall[0]:.2f}× faster)')
ax1.grid(axis='y', alpha=0.3)
ax1.set_ylim(0, max(wall) * 1.15)

# MATLAB phase pie
phases = ['assembly\n32.85 s', 'damage\n3.33 s', 'solve\n161.23 s', 'other\n42.05 s']
sizes  = [32.85, 3.33, 161.23, 239.46 - 32.85 - 3.33 - 161.23]
ax2.pie(sizes, labels=phases, autopct='%1.1f%%', startangle=90,
        colors=['#88c', '#cb8', '#8c8', '#bbb'])
ax2.set_title('MATLAB time breakdown (239.5 s total)')

fig.tight_layout()
fig.savefig(OUT / '04_timing.png', dpi=200)
fig.savefig(OUT / '04_timing.pdf')
plt.close(fig)

# ------------------------------------------------------------------ Fig 5: memory bars (caveat)
fig, ax = plt.subplots(figsize=(7, 4.5))
labels = ['MATLAB\n(process RSS peak)', 'Abaqus\n(pre-run solver estimate)']
mem    = [179.3, 32.0]
bars = ax.bar(labels, mem, color=['tab:orange', 'tab:blue'], edgecolor='k')
for b, v in zip(bars, mem):
    ax.text(b.get_x() + b.get_width() / 2, v, f'{v:.1f} MB',
            ha='center', va='bottom', fontsize=11)
ax.set_ylabel('Memory [MB]')
ax.set_title('Reported memory — note: not directly comparable')
ax.text(0.5, -0.28,
        'MATLAB = actual peak RSS.   Abaqus = .dat "memory to minimize I/O" (linear-solver workspace).\n'
        '.msg "MEMORY PEAK (GB) = 0" → rounded; true RSS not measured.',
        transform=ax.transAxes, ha='center', va='top', fontsize=9, color='gray')
ax.grid(axis='y', alpha=0.3)
fig.subplots_adjust(bottom=0.28)
fig.savefig(OUT / '05_memory.png', dpi=200)
fig.savefig(OUT / '05_memory.pdf')
plt.close(fig)

# ------------------------------------------------------------------ Summary text
summary = f"""Comparison summary — Gregoire 3PB (MATLAB vs Abaqus)
==================================================
Mesh:         5576 CPS3 elements, 2904 nodes, 5808 DOFs  (identical)
Load steps:   10 000  (identical)

Peak load:
  MATLAB  3891.444 N
  Abaqus  3891.444 N        (|Δ| = 1.5e-4 N)

CMOD @ peak: 0.025948 mm    (identical to 3 sig figs)

Pointwise (aligned, 10 000 pts):
  Load  max|Δ| = {np.max(np.abs(dP)):.3f} N    mean|Δ| = {np.mean(np.abs(dP)):.3e} N
  CMOD  max|Δ| = {np.max(np.abs(dC)):.3e} μm   mean|Δ| = {np.mean(np.abs(dC)):.3e} μm

Time:
  MATLAB  wall-clock 239.46 s  (assembly 32.85 + damage 3.33 + solve 161.23 + other)
  Abaqus  wall-clock 669 s     (CPU 408 s = user 392 + sys 16)
  → MATLAB ≈ {669/239.46:.2f}× faster

Memory (NOT directly comparable):
  MATLAB  179.3 MB  (peak process RSS)
  Abaqus   32   MB  (.dat memory-to-minimize-I/O estimate; .msg peak rounded to 0 GB)
"""
(OUT / 'summary.txt').write_text(summary, encoding='utf-8')
print(summary)
print(f'\nWrote plots and summary to: {OUT}')
