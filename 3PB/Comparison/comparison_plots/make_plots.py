"""
Generate comparison plots: MATLAB vs Abaqus — Gregoire 3PB.

Reads load-CMOD CSV data from both solvers, aligns them via
interpolation onto the MATLAB CMOD grid, and produces five
publication-quality figures plus a summary text file.

Data files expected:
    ../abaqus_result/abaqus_load_cmod.csv   (time, cmod_mm, load_N)
    ../Matlab_result/matlab_load_cmod.csv   (cmod_mm, load_N)
    ../Matlab_result/matlab_timing.txt      (timing breakdown)
    ../abaqus_result/abaqus_summary.txt     (peak load / CMOD summary)
"""

from pathlib import Path
import re
import csv
import matplotlib
matplotlib.use('Agg')  # non-interactive backend for headless rendering
import matplotlib.pyplot as plt
import numpy as np

# ──────────────────────────────────────────────────────────────────────
# Paths
# ──────────────────────────────────────────────────────────────────────
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT  = HERE

ABQ_CSV     = ROOT / 'abaqus_result' / 'abaqus_load_cmod.csv'
MAT_CSV     = ROOT / 'Matlab_result' / 'matlab_load_cmod.csv'
MAT_TIMING  = ROOT / 'Matlab_result' / 'matlab_timing.txt'
ABQ_SUMMARY = ROOT / 'abaqus_result' / 'abaqus_summary.txt'


# ──────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────
def read_csv(path, skip=1, delimiter=','):
    """
    Read a numeric CSV file, skipping *skip* header lines.

    Parameters
    ----------
    path : Path
        Absolute path to the CSV file.
    skip : int
        Number of header / comment lines to skip.
    delimiter : str
        Column separator character.

    Returns
    -------
    np.ndarray
        2-D array of shape (n_rows, n_cols).

    Raises
    ------
    FileNotFoundError
        If *path* does not exist.
    ValueError
        If the file is empty after skipping headers.
    """
    if not path.exists():
        raise FileNotFoundError(f"CSV file not found: {path}")

    rows = []
    with open(path, encoding='utf-8-sig') as fh:
        for _ in range(skip):
            next(fh)
        for line in fh:
            line = line.strip()
            if not line:
                continue
            # Handle lines starting with '#' as comments
            if line.startswith('#'):
                continue
            vals = [float(x) for x in line.split(delimiter) if x.strip()]
            rows.append(vals)

    if not rows:
        raise ValueError(f"No data rows found in {path}")

    return np.array(rows)


def parse_matlab_timing(path):
    """
    Parse the MATLAB timing text file for wall-clock, assembly,
    damage, solve times and peak RAM.

    Returns
    -------
    dict  with keys: 'wall', 'assembly', 'damage', 'solve', 'ram_mb'
    """
    info = {}
    if not path.exists():
        return info
    text = path.read_text(encoding='utf-8')
    patterns = {
        'wall':     r'Wall-clock:\s+([\d.]+)\s*s',
        'assembly': r'assembly:\s+([\d.]+)\s*s',
        'damage':   r'damage:\s+([\d.]+)\s*s',
        'solve':    r'solve:\s+([\d.]+)\s*s',
        'ram_mb':   r'Peak RAM:\s+([\d.]+)\s*MB',
    }
    for key, pat in patterns.items():
        m = re.search(pat, text)
        if m:
            info[key] = float(m.group(1))
    return info


# ──────────────────────────────────────────────────────────────────────
# Load data
# ──────────────────────────────────────────────────────────────────────
print("Loading Abaqus data …")
abq = read_csv(ABQ_CSV, skip=1)           # cols: time, cmod_mm, load_N
t_a, c_a, P_a = abq[:, 0], abq[:, 1], abq[:, 2]

print("Loading MATLAB data …")
mat = read_csv(MAT_CSV, skip=0)            # first line is comment "# cmod[mm], load[N]"
c_m, P_m = mat[:, 0], mat[:, 1]

print(f"  Abaqus points : {len(P_a)}")
print(f"  MATLAB points : {len(P_m)}")

# Drop Abaqus row 0 (time = 0, load = 0, cmod = 0) to avoid
# interpolation artifacts at the origin.
mask_nonzero = t_a > 0
t_a2, c_a2, P_a2 = t_a[mask_nonzero], c_a[mask_nonzero], P_a[mask_nonzero]

# ── Interpolate Abaqus onto MATLAB's CMOD grid for pointwise comparison ──
# Both datasets span roughly the same CMOD range.  We restrict to the
# common CMOD overlap so the interpolation is well-posed.
cmod_lo = max(c_a2.min(), c_m.min())
cmod_hi = min(c_a2.max(), c_m.max())

mask_m_common = (c_m >= cmod_lo) & (c_m <= cmod_hi)
c_common = c_m[mask_m_common]
P_m_common = P_m[mask_m_common]

# Interpolate Abaqus load onto the common CMOD grid
P_a_interp = np.interp(c_common, c_a2, P_a2)

# Pointwise differences
dP = P_a_interp - P_m_common             # Load difference [N]
dC_um = np.zeros_like(dP)                 # CMOD difference is 0 by construction
                                          # (both evaluated on the same grid)

# For a meaningful CMOD comparison we can also compare the raw CMOD
# values at aligned *indices* (first N points) where N = min(len(abq), len(mat)).
n_raw = min(len(P_a2), len(P_m))
dC_raw = (c_a2[:n_raw] - c_m[:n_raw]) * 1e3  # [µm]

# ── Peak load statistics ─────────────────────────────────────────────
pk_abq  = P_a2.max()
pk_mat  = P_m.max()
cmod_pk_abq = c_a2[np.argmax(P_a2)]
cmod_pk_mat = c_m[np.argmax(P_m)]

print(f"\n  Abaqus peak : {pk_abq:.3f} N  @ CMOD = {cmod_pk_abq*1e3:.3f} µm")
print(f"  MATLAB peak : {pk_mat:.3f} N  @ CMOD = {cmod_pk_mat*1e3:.3f} µm")
print(f"  |ΔP_peak|   : {abs(pk_abq - pk_mat):.4f} N")

# ── Timing / memory ─────────────────────────────────────────────────
mat_timing = parse_matlab_timing(MAT_TIMING)
mat_wall = mat_timing.get('wall', None)
mat_ram  = mat_timing.get('ram_mb', None)


# ======================================================================
#  PLOTTING — publication style
# ======================================================================
plt.rcParams.update({
    'font.family':    'serif',
    'font.size':       12,
    'axes.labelsize':  13,
    'axes.titlesize':  14,
    'legend.fontsize': 11,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'figure.dpi':      150,
    'savefig.dpi':     300,
    'savefig.bbox':    'tight',
})

# ────────────────────────────── Fig 1: Load-CMOD overlay ──────────────
fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(c_a2 * 1e3, P_a2, '-',  lw=2.0, color='tab:blue',   label='Abaqus')
ax.plot(c_m  * 1e3, P_m,  '--', lw=1.5, color='tab:orange',  label='MATLAB')
ax.set_xlabel('CMOD (µm)')
ax.set_ylabel('Load (N)')
ax.set_title('Load vs CMOD — MATLAB vs Abaqus')
ax.grid(alpha=0.3)
ax.legend()
fig.tight_layout()
fig.savefig(OUT / '01_load_cmod_overlay.png')
fig.savefig(OUT / '01_load_cmod_overlay.pdf')
plt.close(fig)
print("  ✓ 01_load_cmod_overlay")


# ────────────────────────────── Fig 2: Peak region zoom ───────────────
fig, ax = plt.subplots(figsize=(7, 5))
zoom_cmod_max = 60  # µm
mask_a = (c_a2 * 1e3 < zoom_cmod_max)
mask_m = (c_m  * 1e3 < zoom_cmod_max)

ax.plot(c_a2[mask_a] * 1e3, P_a2[mask_a], '-',  lw=2.0, color='tab:blue',   label='Abaqus')
ax.plot(c_m[mask_m]  * 1e3, P_m[mask_m],  '--', lw=1.5, color='tab:orange',  label='MATLAB')
ax.axvline(cmod_pk_abq * 1e3, color='gray', ls=':', alpha=0.6,
           label=f'CMOD @ peak = {cmod_pk_abq*1e3:.2f} µm')
ax.axhline(pk_abq, color='gray', ls=':', alpha=0.6)
ax.set_xlabel('CMOD (µm)')
ax.set_ylabel('Load (N)')
ax.set_title(f'Peak region zoom (CMOD < {zoom_cmod_max} µm)')
ax.grid(alpha=0.3)
ax.legend()
fig.tight_layout()
fig.savefig(OUT / '02_peak_zoom.png')
fig.savefig(OUT / '02_peak_zoom.pdf')
plt.close(fig)
print("  ✓ 02_peak_zoom")


# ────────────────────────────── Fig 3: Pointwise load difference ──────
fig, axes = plt.subplots(2, 1, figsize=(7, 6), sharex=True)

axes[0].plot(c_common * 1e3, dP, color='tab:red', lw=0.8)
axes[0].axhline(0, color='k', lw=0.5)
axes[0].set_ylabel('Δ Load (N)\n(Abaqus − MATLAB)')
axes[0].set_title(
    f'Pointwise difference  (max|ΔP| = {np.max(np.abs(dP)):.2f} N,  '
    f'mean|ΔP| = {np.mean(np.abs(dP)):.2e} N)'
)
axes[0].grid(alpha=0.3)

# CMOD difference using raw aligned indices
axes[1].plot(c_a2[:n_raw] * 1e3, dC_raw, color='tab:purple', lw=0.8)
axes[1].axhline(0, color='k', lw=0.5)
axes[1].set_ylabel('Δ CMOD (µm)\n(Abaqus − MATLAB)')
axes[1].set_xlabel('CMOD (µm)')
axes[1].set_title(
    f'Raw index-aligned CMOD diff  (max|ΔCMOD| = {np.max(np.abs(dC_raw)):.2e} µm)'
)
axes[1].grid(alpha=0.3)

fig.tight_layout()
fig.savefig(OUT / '03_pointwise_diff.png')
fig.savefig(OUT / '03_pointwise_diff.pdf')
plt.close(fig)
print("  ✓ 03_pointwise_diff")


# ────────────────────────────── Fig 4: Timing comparison ──────────────
# Use data from timing files.  Fall back to hardcoded if files missing.
abq_wall = 669.0   # from prior Abaqus run (no .dat file parsed here)
matlab_wall = mat_wall if mat_wall else 763.50

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

solvers = ['MATLAB', 'Abaqus']
wall    = [matlab_wall, abq_wall]
colors  = ['tab:orange', 'tab:blue']

bars = ax1.bar(solvers, wall, color=colors, edgecolor='k')
for b, v in zip(bars, wall):
    ax1.text(b.get_x() + b.get_width() / 2, v, f'{v:.1f} s',
             ha='center', va='bottom', fontsize=11)
ax1.set_ylabel('Wall-clock time (s)')
speedup = abq_wall / matlab_wall
ax1.set_title(f'Wall-clock  (MATLAB {speedup:.2f}× faster)')
ax1.grid(axis='y', alpha=0.3)
ax1.set_ylim(0, max(wall) * 1.15)

# MATLAB phase pie (from timing file or fallback)
asm = mat_timing.get('assembly', 109.24)
dmg = mat_timing.get('damage',   12.11)
slv = mat_timing.get('solve',   580.10)
other = matlab_wall - asm - dmg - slv
phases = [f'assembly\n{asm:.1f} s', f'damage\n{dmg:.1f} s',
          f'solve\n{slv:.1f} s',    f'other\n{other:.1f} s']
sizes  = [asm, dmg, slv, max(other, 0)]
ax2.pie(sizes, labels=phases, autopct='%1.1f%%', startangle=90,
        colors=['#88c', '#cb8', '#8c8', '#bbb'])
ax2.set_title(f'MATLAB time breakdown ({matlab_wall:.1f} s total)')

fig.tight_layout()
fig.savefig(OUT / 'time_comparison_bar.png')
fig.savefig(OUT / 'time_comparison_bar.pdf')
plt.close(fig)
print("  ✓ time_comparison_bar")


# ────────────────────────────── Fig 5: Memory comparison ──────────────
fig, ax = plt.subplots(figsize=(7, 4.5))
labels = ['MATLAB\n(process RSS peak)', 'Abaqus\n(pre-run solver estimate)']
matlab_ram = mat_ram if mat_ram else 84.0
abq_mem    = 32.0   # from .dat file estimate

mem = [matlab_ram, abq_mem]
bars = ax.bar(labels, mem, color=['tab:orange', 'tab:blue'], edgecolor='k')
for b, v in zip(bars, mem):
    ax.text(b.get_x() + b.get_width() / 2, v, f'{v:.1f} MB',
            ha='center', va='bottom', fontsize=11)
ax.set_ylabel('Memory (MB)')
ax.set_title('Reported memory — note: not directly comparable')
ax.text(0.5, -0.28,
        'MATLAB = actual peak RSS.   Abaqus = .dat "memory to minimize I/O" '
        '(linear-solver workspace).\n'
        '.msg "MEMORY PEAK (GB) = 0" → rounded; true RSS not measured.',
        transform=ax.transAxes, ha='center', va='top', fontsize=9, color='gray')
ax.grid(axis='y', alpha=0.3)
fig.subplots_adjust(bottom=0.28)
fig.savefig(OUT / '05_memory.png')
fig.savefig(OUT / '05_memory.pdf')
plt.close(fig)
print("  ✓ 05_memory")


# ──────────────────────────────────────────────────────────────────────
# Summary text
# ──────────────────────────────────────────────────────────────────────
summary = f"""Comparison summary — Gregoire 3PB (MATLAB vs Abaqus)
==================================================
Abaqus data points : {len(P_a)}
MATLAB data points : {len(P_m)}

Peak load:
  Abaqus  {pk_abq:.3f} N   @ CMOD = {cmod_pk_abq*1e3:.3f} µm
  MATLAB  {pk_mat:.3f} N   @ CMOD = {cmod_pk_mat*1e3:.3f} µm
  |ΔP_peak| = {abs(pk_abq - pk_mat):.4f} N

Pointwise (interpolated onto common CMOD grid, {len(c_common)} pts):
  Load  max|Δ| = {np.max(np.abs(dP)):.3f} N    mean|Δ| = {np.mean(np.abs(dP)):.3e} N

Raw index-aligned CMOD diff ({n_raw} pts):
  CMOD  max|Δ| = {np.max(np.abs(dC_raw)):.3e} µm   mean|Δ| = {np.mean(np.abs(dC_raw)):.3e} µm

Time:
  MATLAB  wall-clock {matlab_wall:.2f} s
  Abaqus  wall-clock {abq_wall:.1f} s
  → Abaqus / MATLAB = {speedup:.2f}×

Memory (NOT directly comparable):
  MATLAB  {matlab_ram:.1f} MB  (peak process RSS)
  Abaqus  {abq_mem:.0f}   MB  (.dat memory-to-minimize-I/O estimate)
"""
(OUT / 'summary.txt').write_text(summary, encoding='utf-8')
print(summary)
print(f'\nWrote plots and summary to: {OUT}')
