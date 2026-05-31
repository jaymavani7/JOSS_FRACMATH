# -*- coding: utf-8 -*-
# =====================================================================
#  plot_matlab_style.py
#
#  Standalone plotter — run with NORMAL python (python3 + matplotlib),
#  NOT inside Abaqus. Reads the data files written by
#  RUN_3PB_ABAQUS_FULL_MATCH_MATLAB.py and produces plots IDENTICAL in
#  style to the MATLAB solver figures.
#
#  Use this if matplotlib was not available inside Abaqus python.
#
#  RUN:
#    python plot_matlab_style.py
#
#  Reads from abaqus_results/:
#    abaqus_load_cmod.csv
#    damage_peak.dat
#    damage_postpeak.dat
#  Writes to abaqus_results/:
#    abaqus_load_cmod.png
#    abaqus_damage_peak.png
#    abaqus_damage_postpeak.png
# =====================================================================
import os
import csv

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize
import numpy as np

RESULT_DIR = 'abaqus_results'


# ---- EXACT MATLAB crack_cmap ----------------------------------------
def crack_cmap_rgb(v):
    stops = [
        (0.00, (0.84, 0.88, 0.95)),
        (0.10, (0.18, 0.42, 0.86)),
        (0.30, (0.05, 0.72, 0.88)),
        (0.50, (0.18, 0.80, 0.32)),
        (0.70, (0.96, 0.90, 0.08)),
        (0.85, (0.98, 0.44, 0.04)),
        (1.00, (0.82, 0.04, 0.04)),
    ]
    v = max(0.0, min(1.0, float(v)))
    for i in range(len(stops) - 1):
        if stops[i][0] <= v <= stops[i + 1][0]:
            t = (v - stops[i][0]) / max(stops[i + 1][0] - stops[i][0], 1e-12)
            c0 = stops[i][1]; c1 = stops[i + 1][1]
            return tuple(c0[j] * (1 - t) + c1[j] * t for j in range(3))
    return stops[-1][1]


def crack_colormap():
    grid = np.linspace(0, 1, 256)
    return LinearSegmentedColormap.from_list(
        'crack', [crack_cmap_rgb(t) for t in grid], N=256)


# ---- Load-CMOD plot (MATLAB-exact) ----------------------------------
def plot_load_cmod(csv_path, png_path):
    xs = []; ys = []
    with open(csv_path) as f:
        r = csv.reader(f); next(r, None)
        for row in r:
            if len(row) >= 3:
                try:
                    xs.append(float(row[1])); ys.append(float(row[2]) / 1000.0)
                except Exception:
                    pass
    if len(xs) < 2:
        print('too few points'); return

    pk = max(range(len(ys)), key=lambda i: ys[i])
    pk_load = ys[pk]; pk_cmod = xs[pk]

    fig = plt.figure(figsize=(5.0, 3.9), dpi=300)
    ax = fig.add_axes([0.14, 0.14, 0.82, 0.78])
    ax.grid(True, linestyle=':', color=(0.80, 0.80, 0.80), linewidth=0.7)
    ax.set_axisbelow(True)

    ax.fill_between(xs, ys, 0, color=(0.08, 0.30, 0.72), alpha=0.10,
                    linewidth=0)
    ax.plot(xs, ys, '-', color=(0.08, 0.30, 0.72), linewidth=2.0)
    ax.plot([pk_cmod], [pk_load], marker='*', markersize=15,
            markerfacecolor=(0.98, 0.82, 0.0),
            markeredgecolor=(0.40, 0.28, 0.0), markeredgewidth=1.0,
            linestyle='None')
    ax.annotate('$P_{\\rm peak}=%.2f$ kN\nCMOD$=%.4f$ mm'
                % (pk_load, pk_cmod),
                xy=(pk_cmod, pk_load),
                xytext=(pk_cmod + 0.012, pk_load * 1.02),
                fontsize=8.5, color=(0.30, 0.20, 0.0),
                bbox=dict(boxstyle='round,pad=0.3', fc='white',
                          ec=(0.60, 0.50, 0.20), lw=0.5))

    ax.set_xlabel('CMOD [mm]', fontsize=10)
    ax.set_ylabel('Load [kN]', fontsize=10)
    ax.set_title('Load--CMOD response', fontsize=10, fontweight='bold')
    ax.set_xlim(0.0, 0.35)
    ax.set_ylim(0.0, pk_load * 1.28)
    ax.tick_params(labelsize=9)
    for s in ax.spines.values():
        s.set_linewidth(0.7)

    fig.savefig(png_path)
    fig.savefig(png_path.replace('.png', '.pdf'))
    plt.close(fig)
    print('saved ' + png_path)


# ---- Damage contour (MATLAB-exact) ----------------------------------
def plot_damage(dat_path, png_path, title_label):
    nodes = {}; disp = {}; elems = []
    with open(dat_path) as f:
        for line in f:
            p = line.split()
            if not p:
                continue
            if p[0] == 'N':
                lb = int(p[1])
                nodes[lb] = (float(p[2]), float(p[3]))
                disp[lb] = (float(p[4]), float(p[5]))
            elif p[0] == 'E':
                elems.append((int(p[1]), int(p[2]), int(p[3]),
                              int(p[4]), float(p[5])))
    if not nodes or not elems:
        print('empty geometry ' + dat_path); return

    umax = 0.0
    for lb in nodes:
        ux, uy = disp[lb]
        umax = max(umax, abs(ux), abs(uy))
    scale = min(1.0, max(30.0, 2.0 / max(umax, 1e-12)))

    defxy = {}
    for lb in nodes:
        x, y = nodes[lb]; ux, uy = disp[lb]
        defxy[lb] = (x + scale * ux, y + scale * uy)

    xs = [p[0] for p in defxy.values()]; ys = [p[1] for p in defxy.values()]
    xmin, xmax = min(xs), max(xs); ymin, ymax = min(ys), max(ys)

    fig = plt.figure(figsize=(9.0, 4.0), dpi=300)
    ax = fig.add_axes([0.10, 0.15, 0.75, 0.75])

    polys_all = []
    for eid, n1, n2, n3, om in elems:
        if n1 in defxy and n2 in defxy and n3 in defxy:
            polys_all.append([defxy[n1], defxy[n2], defxy[n3]])
    ax.add_collection(PolyCollection(
        polys_all, facecolors=(0.92, 0.92, 0.93),
        edgecolors=(0.75, 0.77, 0.80), linewidths=0.15))

    polys_cr = []; cols_cr = []
    for eid, n1, n2, n3, om in elems:
        if om > 0.01 and n1 in defxy and n2 in defxy and n3 in defxy:
            polys_cr.append([defxy[n1], defxy[n2], defxy[n3]])
            cols_cr.append(crack_cmap_rgb(om))
    if polys_cr:
        ax.add_collection(PolyCollection(
            polys_cr, facecolors=cols_cr, edgecolors='none'))

    sm = plt.cm.ScalarMappable(cmap=crack_colormap(),
                               norm=Normalize(vmin=0, vmax=1))
    sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label(r'$\omega$  (damage)', fontsize=10)

    ax.set_xlim(xmin - 5, xmax + 5)
    ax.set_ylim(ymin - 5, ymax + 5)
    ax.set_aspect('equal')
    ax.set_xlabel('$x$ [mm]', fontsize=10)
    ax.set_ylabel('$y$ [mm]', fontsize=10)
    for s in ax.spines.values():
        s.set_linewidth(0.8)

    ax.text(0.03, 0.94, title_label, transform=ax.transAxes,
            fontsize=11, fontweight='bold', va='top', ha='left',
            bbox=dict(boxstyle='round,pad=0.4', fc='white',
                      ec=(0.2, 0.2, 0.2), lw=0.8))

    fig.savefig(png_path)
    fig.savefig(png_path.replace('.png', '.pdf'))
    plt.close(fig)
    print('saved ' + png_path)


def main():
    csv_path = os.path.join(RESULT_DIR, 'abaqus_load_cmod.csv')
    if os.path.exists(csv_path):
        plot_load_cmod(csv_path,
                       os.path.join(RESULT_DIR, 'abaqus_load_cmod.png'))
    else:
        print('missing ' + csv_path)

    pk = os.path.join(RESULT_DIR, 'damage_peak.dat')
    if os.path.exists(pk):
        # read peak load from summary for the label
        label = 'Peak load'
        sp = os.path.join(RESULT_DIR, 'abaqus_summary.txt')
        if os.path.exists(sp):
            for line in open(sp):
                if 'kN' in line and 'Peak' in line:
                    label = 'Peak load: %s' % line.split(':')[1].strip()
                    break
        plot_damage(pk, os.path.join(RESULT_DIR, 'abaqus_damage_peak.png'),
                    label)

    pp = os.path.join(RESULT_DIR, 'damage_postpeak.dat')
    if os.path.exists(pp):
        plot_damage(pp,
                    os.path.join(RESULT_DIR, 'abaqus_damage_postpeak.png'),
                    'Post-peak')


if __name__ == '__main__':
    main()
