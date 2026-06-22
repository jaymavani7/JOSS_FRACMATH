
from __future__ import print_function
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

BASE   = r"C:\Users\jmavani\OneDrive - University of New Mexico\Desktop\JOSS\3pb_2"
MAT_R  = os.path.join(BASE, "matlab",  "Gregoire_3PB", "results")
ABQ_R  = os.path.join(BASE, "abaqus",  "Gregoire_3PB", "results")
OUT    = os.path.join(BASE, "comparison")
os.makedirs(OUT, exist_ok=True)

DPI = 200

def side_by_side(left_path, right_path, left_label, right_label,
                 title, out_name, subtitle=""):

    imgs = []
    for p in (left_path, right_path):
        if p and os.path.exists(p):
            imgs.append(mpimg.imread(p))
        else:
            imgs.append(None)

    fig = plt.figure(figsize=(18, 7), facecolor='white')

    fig.text(0.5, 0.97, title,
             ha='center', va='top', fontsize=15, fontweight='bold', color='#222222')
    if subtitle:
        fig.text(0.5, 0.92, subtitle,
                 ha='center', va='top', fontsize=10, color='#555555', style='italic')

    top_start = 0.88 if subtitle else 0.91

    ax_l = fig.add_axes([0.01, 0.06, 0.47, top_start - 0.06])
    ax_l.axis('off')
    if imgs[0] is not None:
        ax_l.imshow(imgs[0], aspect='equal')
    else:
        ax_l.text(0.5, 0.5, 'Image not found\n%s' % left_path,
                  ha='center', va='center', fontsize=10, color='red',
                  transform=ax_l.transAxes)

    bbox_l = ax_l.get_position()
    bar_l  = fig.add_axes([bbox_l.x0, top_start, bbox_l.width, 0.035])
    bar_l.set_facecolor('#0F4C81')
    bar_l.axis('off')
    bar_l.text(0.5, 0.45, left_label,
               ha='center', va='center', fontsize=12, fontweight='bold',
               color='white', transform=bar_l.transAxes)

    fig.add_axes([0.487, 0.04, 0.003, top_start - 0.04 + 0.035]).set_visible(False)
    fig.lines.extend([
        plt.Line2D([0.488, 0.488], [0.04, top_start + 0.035],
                   transform=fig.transFigure,
                   color='#cccccc', linewidth=1.5)
    ])

    ax_r = fig.add_axes([0.52, 0.06, 0.47, top_start - 0.06])
    ax_r.axis('off')
    if imgs[1] is not None:
        ax_r.imshow(imgs[1], aspect='equal')
    else:
        ax_r.text(0.5, 0.5, 'Image not found\n%s' % right_path,
                  ha='center', va='center', fontsize=10, color='red',
                  transform=ax_r.transAxes)

    bbox_r = ax_r.get_position()
    bar_r  = fig.add_axes([bbox_r.x0, top_start, bbox_r.width, 0.035])
    bar_r.set_facecolor('#C0392B')
    bar_r.axis('off')
    bar_r.text(0.5, 0.45, right_label,
               ha='center', va='center', fontsize=12, fontweight='bold',
               color='white', transform=bar_r.transAxes)

    out_png = os.path.join(OUT, out_name + '.png')
    out_pdf = os.path.join(OUT, out_name + '.pdf')
    fig.savefig(out_png, dpi=DPI, facecolor='white', bbox_inches='tight')
    fig.savefig(out_pdf,           facecolor='white', bbox_inches='tight')
    plt.close(fig)
    print('  Saved: %s' % out_name)

def read_csv2(path):
    rows = []
    with open(path, 'r') as f:
        for ln in f:
            s = ln.strip()
            if not s or s.startswith('#'):
                continue
            p = s.replace(',', ' ').split()
            try:
                rows.append((float(p[0]), float(p[1])))
            except Exception:
                continue
    if not rows:
        return np.zeros(0), np.zeros(0)
    arr = np.array(rows)
    return arr[:, 0], arr[:, 1]

print('Damage post-peak comparison ...')
side_by_side(
    left_path   = os.path.join(MAT_R, 'fig_damage_postpeak.png'),
    right_path  = os.path.join(ABQ_R, 'damage_postpeak.png'),
    left_label  = 'MATLAB — FRACMATH',
    right_label = 'Abaqus — UMAT',
    title       = 'Crack geometry at Post-peak load  (ω ≥ 0.99)',
    subtitle    = 'Three-Point Bending · Grégoire specimen · CDM model',
    out_name    = 'comparison_damage_postpeak',
)

print('Damage last-step comparison ...')
side_by_side(
    left_path   = os.path.join(MAT_R, 'fig_damage_last_step.png'),
    right_path  = os.path.join(ABQ_R, 'damage_last_step.png'),
    left_label  = 'MATLAB — FRACMATH',
    right_label = 'Abaqus — UMAT',
    title       = 'Crack geometry at Last step  (ω ≥ 0.99)',
    subtitle    = 'Three-Point Bending · Grégoire specimen · CDM model',
    out_name    = 'comparison_damage_last_step',
)

print('Damage peak comparison ...')
abq_peak = os.path.join(ABQ_R, 'damage_peak.png')
side_by_side(
    left_path   = os.path.join(MAT_R, 'fig_damage_peak.png'),
    right_path  = abq_peak if os.path.exists(abq_peak) else None,
    left_label  = 'MATLAB — FRACMATH',
    right_label = 'Abaqus — UMAT  (run abaqus python extract_peak_omega.py to generate)',
    title       = 'Crack geometry at Peak load  (ω ≥ 0.99)',
    subtitle    = 'Three-Point Bending · Grégoire specimen · CDM model',
    out_name    = 'comparison_damage_peak',
)

print('Load-CMOD panel comparison ...')
mat_csv = os.path.join(MAT_R, 'matlab_load_cmod.csv')
abq_csv = os.path.join(ABQ_R, 'abaqus_load_cmod.csv')

cmod_m, F_m = read_csv2(mat_csv)
cmod_a, F_a = read_csv2(abq_csv)
F_m_kN = F_m / 1000.0
F_a_kN = F_a / 1000.0

fig, axes = plt.subplots(1, 2, figsize=(16, 6), facecolor='white',
                         sharey=False)
fig.suptitle('Load vs. CMOD — Individual Panels\n'
             'Three-Point Bending · Grégoire specimen · CDM model',
             fontsize=13, fontweight='bold', y=0.98, color='#222222')

COL_M = '#0F4C81'
COL_A = '#C0392B'

for ax, cmod, F_kN, col, solver in [
        (axes[0], cmod_m, F_m_kN, COL_M, 'MATLAB — FRACMATH'),
        (axes[1], cmod_a, F_a_kN, COL_A, 'Abaqus — UMAT'),
]:
    ip = int(np.argmax(F_kN))
    pk_c, pk_f = cmod[ip], F_kN[ip]

    ax.fill_between(cmod, F_kN, 0, color=col, alpha=0.10, linewidth=0)
    ax.plot(cmod, F_kN, '-', color=col, linewidth=2.2)
    ax.plot(pk_c, pk_f, '*', markersize=14,
            markerfacecolor='#FFD700', markeredgecolor='#5a4000',
            markeredgewidth=0.8, zorder=5)

    x_max = cmod.max()
    ha = 'left' if pk_c < 0.55 * x_max else 'right'
    off = 0.015 * x_max * (1 if ha == 'left' else -1)
    ax.annotate('Peak: %.2f kN\nCMOD: %.4f mm' % (pk_f, pk_c),
                xy=(pk_c, pk_f),
                xytext=(pk_c + off, pk_f * 1.08),
                fontsize=9, fontweight='bold', color=col,
                ha=ha, va='bottom',
                bbox=dict(facecolor='white', edgecolor=col, linewidth=0.7, pad=3),
                arrowprops=dict(arrowstyle='->', color=col, lw=0.9))

    ax.set_xlim(0, x_max * 1.05)
    ax.set_ylim(0, F_kN.max() * 1.30)
    ax.set_xlabel('CMOD  [mm]', fontsize=11)
    ax.set_ylabel('Load  [kN]', fontsize=11)
    ax.set_title(solver, fontsize=12, fontweight='bold', color=col, pad=8)
    ax.grid(True, linestyle=':', linewidth=0.5, color='#bbbbbb')
    for sp in ax.spines.values():
        sp.set_linewidth(0.6)
    ax.tick_params(labelsize=9, direction='out', top=False, right=False)

fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(os.path.join(OUT, 'comparison_load_cmod_panels.png'),
            dpi=300, facecolor='white', bbox_inches='tight')
fig.savefig(os.path.join(OUT, 'comparison_load_cmod_panels.pdf'),
            facecolor='white', bbox_inches='tight')
plt.close(fig)
print('  Saved: comparison_load_cmod_panels')

print('Mesh comparison ...')
side_by_side(
    left_path   = os.path.join(MAT_R, 'fig_mesh.png'),
    right_path  = None,
    left_label  = 'MATLAB — FRACMATH',
    right_label = 'Abaqus — (same mesh, no standalone figure)',
    title       = 'Finite element mesh  (shared between both solvers)',
    subtitle    = 'Grégoire 3PB · 7319 nodes · 14268 T3 elements',
    out_name    = 'comparison_mesh',
)

print()
print('Done. All comparison images saved to:')
print('  %s' % OUT)
