# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
plot_results.py
===============
Standalone Python script (run with a normal CPython / Anaconda, NOT Abaqus).
Reads the CSVs already in results/ and produces:

    results/load_cmod_plot.png / .pdf          -- Load vs CMOD (publication quality)
    results/damage_postpeak.png / .pdf         -- Damage map at post-peak frame
    results/damage_last_step.png / .pdf        -- Damage map at last step frame
    results/damage_peak.png / .pdf             -- Damage map at peak load (if data exists)

RUN:
    python plot_results.py
    python plot_results.py "path/to/results"

The script auto-locates the results folder relative to itself if not supplied.

NOTE on omega_peak.csv being empty:
    The Abaqus field output was written every FIELD_FREQ increments (default 100).
    The peak load frame often falls between two field frames.  The postpeak and
    last frames are used instead.  If you need peak-frame damage, re-run with
    ABQ_FIELD_FREQ=1 or extract again with extract_damage.py (abaqus python).
"""

from __future__ import print_function
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches

# ── tuneable ──────────────────────────────────────────────────────────────────
FULL_DAMAGE_THRESHOLD = 0.99   # omega >= this → element is "fully damaged"
CMOD_XLIM             = None   # set e.g. (0, 0.35) to clip x-axis; None = auto
LOAD_YLIM             = None   # set e.g. (0, 5.0)  to clip y-axis (kN); None = auto
DPI                   = 300    # output resolution
# ─────────────────────────────────────────────────────────────────────────────


def _find_results(arg):
    if arg and os.path.isdir(arg):
        return os.path.abspath(arg)
    # search relative to this script
    base = os.path.dirname(os.path.abspath(__file__))
    cands = [
        os.path.join(base, 'results'),
        os.path.join(base, 'Gregoire_3PB', 'results'),
    ]
    for c in cands:
        if os.path.isdir(c):
            return c
    raise FileNotFoundError(
        'Cannot find results/ folder.  Pass it as an argument:\n'
        '  python plot_results.py "path/to/results"')


def read_csv2(path):
    """Read a 2-column CSV (skip # comments).  Returns (col0, col1) as float arrays."""
    rows = []
    with open(path, 'r') as f:
        for ln in f:
            s = ln.strip()
            if not s or s.startswith('#'):
                continue
            parts = s.replace(',', ' ').split()
            try:
                rows.append((float(parts[0]), float(parts[1])))
            except (ValueError, IndexError):
                continue
    if not rows:
        return np.zeros(0), np.zeros(0)
    arr = np.array(rows)
    return arr[:, 0], arr[:, 1]


def read_omega_csv(path):
    """Read omega CSV -> dict {element_label: max_omega_over_IPs}."""
    d = {}
    if not os.path.exists(path):
        return d
    with open(path, 'r') as f:
        for ln in f:
            s = ln.strip()
            if not s or s.startswith('#'):
                continue
            parts = s.replace(',', ' ').split()
            try:
                eid, omega = int(parts[0]), float(parts[2])
            except (ValueError, IndexError):
                continue
            d[eid] = max(d.get(eid, 0.0), omega)
    return d


def load_mesh(pd_dir):
    """
    Returns:
        nodes  – (N,2) float array of [x, y]
        node_map – {label: index}
        elems  – (M,3) int array of node indices (into nodes)
        elabels – (M,) int array of element labels
        elem_map – {label: index}
    """
    nraw = []
    with open(os.path.join(pd_dir, 'mesh_nodes.csv'), 'r') as f:
        for ln in f:
            s = ln.strip()
            if not s or s.startswith('#'):
                continue
            p = s.replace(',', ' ').split()
            nraw.append((int(p[0]), float(p[1]), float(p[2])))
    nlabels = [r[0] for r in nraw]
    nodes = np.array([[r[1], r[2]] for r in nraw])
    node_map = {lbl: i for i, lbl in enumerate(nlabels)}

    eraw = []
    with open(os.path.join(pd_dir, 'mesh_elements.csv'), 'r') as f:
        for ln in f:
            s = ln.strip()
            if not s or s.startswith('#'):
                continue
            p = s.replace(',', ' ').split()
            eraw.append((int(p[0]), int(p[1]), int(p[2]), int(p[3])))
    elabels = np.array([r[0] for r in eraw], dtype=int)
    elems = np.array([[node_map[r[1]], node_map[r[2]], node_map[r[3]]]
                      for r in eraw], dtype=int)
    elem_map = {lbl: i for i, lbl in enumerate(elabels)}
    return nodes, node_map, elems, elabels, elem_map


def omega_array(omega_dict, elabels, elem_map):
    """Map omega dict to per-element array aligned with elabels."""
    w = np.zeros(len(elabels))
    for lbl, val in omega_dict.items():
        if lbl in elem_map:
            w[elem_map[lbl]] = val
    return w


# ─────────────────────────────────── FIGURE 1 ─────────────────────────────────
def fig_load_cmod(CMOD, F_N, res_dir):
    """
    Publication-quality Load vs CMOD.
    F_N  – load in Newtons.
    """
    F_kN = F_N / 1000.0

    ip = int(np.argmax(F_kN))
    pk_cmod, pk_load = CMOD[ip], F_kN[ip]

    fig, ax = plt.subplots(figsize=(8.8 / 2.54, 7.0 / 2.54), facecolor='w')
    ax.set_facecolor('white')

    col_fill = (0.08, 0.30, 0.72)
    col_line = (0.04, 0.18, 0.56)

    # shaded area under curve
    ax.fill_between(CMOD, F_kN, 0,
                    color=col_fill, alpha=0.12, linewidth=0)

    # main curve
    ax.plot(CMOD, F_kN, '-', color=col_line, linewidth=1.8, zorder=3)

    # peak marker
    ax.plot(pk_cmod, pk_load, marker='*', markersize=14,
            markerfacecolor=(0.98, 0.80, 0.0),
            markeredgecolor=(0.38, 0.24, 0.0),
            markeredgewidth=0.8,
            linewidth=0, zorder=5)

    # annotation
    x_max = CMOD.max()
    ha = 'left' if pk_cmod < 0.6 * x_max else 'right'
    offset = 0.02 * x_max if ha == 'left' else -0.02 * x_max
    ax.annotate(
        '$P_{\\rm peak} = %.2f$ kN\nCMOD $= %.4f$ mm' % (pk_load, pk_cmod),
        xy=(pk_cmod, pk_load),
        xytext=(pk_cmod + offset, pk_load * 1.05),
        fontsize=8.0, fontweight='bold',
        color=(0.28, 0.18, 0.0), ha=ha, va='bottom',
        bbox=dict(facecolor='w', edgecolor=(0.60, 0.48, 0.18),
                  linewidth=0.5, pad=3),
        arrowprops=dict(arrowstyle='->', color=(0.50, 0.38, 0.10),
                        lw=0.8))

    # grid
    ax.grid(True, linestyle=':', linewidth=0.5, color=(0.80, 0.80, 0.80), zorder=0)

    # axes formatting
    for sp in ax.spines.values():
        sp.set_linewidth(0.7)
    ax.tick_params(labelsize=9, direction='out', top=False, right=False)
    ax.set_xlabel('CMOD  [mm]', fontsize=10)
    ax.set_ylabel('Load  [kN]', fontsize=10)
    ax.set_title('Load – CMOD response\n(Gregoire 3PB, CDM UMAT)',
                 fontsize=10, fontweight='bold', pad=6)

    if CMOD_XLIM:
        ax.set_xlim(*CMOD_XLIM)
    else:
        ax.set_xlim(0, x_max * 1.04)

    if LOAD_YLIM:
        ax.set_ylim(*LOAD_YLIM)
    else:
        ax.set_ylim(0, pk_load * 1.28)

    fig.tight_layout()
    base = os.path.join(res_dir, 'load_cmod_plot')
    fig.savefig(base + '.png', dpi=DPI, facecolor='w', bbox_inches='tight')
    fig.savefig(base + '.pdf', dpi=DPI, facecolor='w', bbox_inches='tight')
    plt.close(fig)
    print('  Wrote %s.png / .pdf' % base)
    return pk_load, pk_cmod


# ─────────────────────────────────── FIGURE 2 ─────────────────────────────────
def fig_damage(nodes, elems, omega, label, base):
    """
    Damage map (MATLAB style):
      • Grey full mesh
      • Only elements with omega >= threshold shown in solid black
      • No partial-damage colormap — binary: mesh grey OR crack black
    """
    th = FULL_DAMAGE_THRESHOLD
    w_max = omega.max() if omega.size else 0.0

    verts = nodes[elems]                       # (M, 3, 2)
    x_min, x_max = nodes[:, 0].min(), nodes[:, 0].max()
    y_min, y_max = nodes[:, 1].min(), nodes[:, 1].max()
    aspect = (x_max - x_min) / max(y_max - y_min, 1e-9)

    # figure size scales with specimen aspect ratio
    fw = min(max(20.0 / 2.54, aspect * 8.5 / 2.54), 34.0 / 2.54)
    fh_fig = 9.0 / 2.54

    # layout: main axes | text panel
    left_pad = 0.07
    panel_w  = 0.22
    gap      = 0.02
    ax_w     = 1.0 - left_pad - panel_w - gap
    bottom   = 0.15
    top_     = 0.88
    ax_h     = top_ - bottom

    fig = plt.figure(figsize=(fw, fh_fig), facecolor='w')
    ax  = fig.add_axes([left_pad, bottom, ax_w, ax_h])

    mesh_fc  = (0.94, 0.94, 0.95)
    mesh_ec  = (0.80, 0.82, 0.84)
    crack_fc = (0.00, 0.00, 0.00)

    # Layer 0: full grey mesh
    ax.add_collection(PolyCollection(verts,
                                     facecolors=mesh_fc,
                                     edgecolors=mesh_ec,
                                     linewidths=0.06))

    # Layer 1: ONLY fully damaged elements (omega >= threshold) in black
    full_mask = omega >= th
    n_cracked = int(full_mask.sum())
    if full_mask.any():
        ax.add_collection(PolyCollection(verts[full_mask],
                                         facecolors=crack_fc,
                                         edgecolors=crack_fc,
                                         linewidths=0.18))

    ax.set_aspect('equal')
    ax.set_xlim(x_min - 3, x_max + 3)
    ax.set_ylim(y_min - 3, y_max + 3)
    ax.set_xlabel(r'$x$  [mm]', fontsize=11)
    ax.set_ylabel(r'$y$  [mm]', fontsize=11)
    ax.set_title('Fully damaged elements  ($\\omega \\geq %.2f$)' % th,
                 fontsize=12, fontweight='bold', pad=4)

    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(labelsize=9, direction='out', top=False, right=False)

    # ── right-side legend panel ────────────────────────────────────────────
    axp = fig.add_axes([left_pad + ax_w + gap, bottom, panel_w - 0.02, ax_h])
    axp.set_xlim(0, 1); axp.set_ylim(0, 1); axp.axis('off')

    axp.text(0.0, 1.00, 'Plot key', fontsize=10.5, fontweight='bold',
             ha='left', va='top')

    legend_handles = [
        Line2D([0], [0], color=mesh_ec, lw=7, solid_capstyle='butt',
               label='FE mesh'),
        Line2D([0], [0], color=crack_fc, lw=7, solid_capstyle='butt',
               label='$\\omega \\geq %.2f$' % th),
    ]
    axp.legend(handles=legend_handles,
               loc='upper left',
               bbox_to_anchor=(0.0, 0.85),
               frameon=False,
               fontsize=9.5,
               handlelength=1.8,
               handletextpad=0.7,
               borderaxespad=0.0,
               labelspacing=0.9)

    info = [
        label,
        'Threshold: %.2f' % th,
        'Cracked: %d elems' % n_cracked,
        'Max $\\omega$: %.4f' % w_max,
    ]
    y_txt = 0.38
    for txt in info:
        axp.text(0.0, y_txt, txt, fontsize=9.0, fontweight='bold',
                 ha='left', va='top')
        y_txt -= 0.14

    fig.savefig(base + '.png', dpi=DPI, facecolor='w', bbox_inches='tight')
    fig.savefig(base + '.pdf', dpi=DPI, facecolor='w', bbox_inches='tight')
    plt.close(fig)
    print('  Wrote %s.png / .pdf  [cracked: %d, max omega: %.4f]'
          % (os.path.basename(base), n_cracked, w_max))


# ──────────────────────────────────────────────────────────────────────────────
def main():
    arg = sys.argv[1] if len(sys.argv) > 1 else None
    res_dir = _find_results(arg)
    pd_dir = os.path.join(res_dir, 'plotdata')

    print('Results dir : %s' % res_dir)
    print('Plotdata dir: %s' % pd_dir)
    print()

    # ── Load-CMOD ────────────────────────────────────────────────────────────
    csv_path = os.path.join(res_dir, 'abaqus_load_cmod.csv')
    if not os.path.exists(csv_path):
        print('ERROR: %s not found.  Cannot plot Load-CMOD.' % csv_path)
    else:
        CMOD, F = read_csv2(csv_path)
        if CMOD.size == 0:
            print('ERROR: Load-CMOD CSV is empty.')
        else:
            print('Load-CMOD: %d data points' % len(CMOD))
            print('  Peak load  = %.2f kN' % (F.max() / 1000.0))
            print('  CMOD range = [%.4f, %.4f] mm' % (CMOD.min(), CMOD.max()))
            print('Plotting Load vs CMOD ...')
            pk_load, pk_cmod = fig_load_cmod(CMOD, F, res_dir)
            print()

    # ── mesh ─────────────────────────────────────────────────────────────────
    nodes_csv = os.path.join(pd_dir, 'mesh_nodes.csv')
    elems_csv = os.path.join(pd_dir, 'mesh_elements.csv')
    if not (os.path.exists(nodes_csv) and os.path.exists(elems_csv)):
        print('ERROR: mesh CSVs not found in %s.  Cannot plot damage.' % pd_dir)
        return

    print('Loading mesh ...')
    nodes, node_map, elems, elabels, elem_map = load_mesh(pd_dir)
    print('  Nodes   : %d' % len(nodes))
    print('  Elements: %d' % len(elems))
    print()

    # ── damage figures ────────────────────────────────────────────────────────
    cases = [
        ('omega_peak.csv',     'Peak load',  'damage_peak'),
        ('omega_postpeak.csv', 'Post-peak',  'damage_postpeak'),
        ('omega_last.csv',     'Last step',  'damage_last_step'),
    ]

    for csv_name, label, out_name in cases:
        csv_path = os.path.join(pd_dir, csv_name)
        if not os.path.exists(csv_path):
            print('SKIP: %s not found.' % csv_name)
            continue
        omap = read_omega_csv(csv_path)
        if not omap:
            print('SKIP: %s is empty (no damage data at this frame).' % csv_name)
            continue
        w = omega_array(omap, elabels, elem_map)
        w_max = w.max()
        print('Plotting damage (%s): max omega = %.4f ...' % (label, w_max))
        base = os.path.join(res_dir, out_name)
        fig_damage(nodes, elems, w, label, base)

    print()
    print('Done.  All figures saved to:')
    print('  %s' % res_dir)


if __name__ == '__main__':
    main()
