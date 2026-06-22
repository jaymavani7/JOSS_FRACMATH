import csv
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
MAT_R = os.path.join(BASE, "matlab", "Gregoire_3PB", "results")
ABQ_R = os.path.join(BASE, "abaqus", "Gregoire_3PB", "results")
OUT = os.path.join(BASE, "comparison")
os.makedirs(OUT, exist_ok=True)

def read_xy(path):
    x_vals, y_vals = [], []
    with open(path, "r") as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue
            first = row[0].strip()
            if not first or first.startswith("#") or first.lower().startswith("cmod"):
                continue
            try:
                x_vals.append(float(row[0]))
                y_vals.append(float(row[1]))
            except (ValueError, IndexError):
                continue
    if not x_vals:
        raise ValueError("No numeric data found in %s" % path)
    return np.array(x_vals), np.array(y_vals)

def first_float(pattern, text, default=None):
    match = re.search(pattern, text, flags=re.IGNORECASE)
    return float(match.group(1)) if match else default

def read_matlab_timing(path):
    with open(path, "r") as f:
        text = f.read()
    solver = first_float(r"Solver wall-clock:\s*([0-9.]+)\s*s", text)
    if solver is None:
        raise ValueError("Could not parse MATLAB solver time from %s" % path)
    total = first_float(r"End-to-end:\s*([0-9.]+)\s*s", text, solver)
    return {
        "solver": solver,
        "total": total,
        "assembly": first_float(r"assembly:\s*([0-9.]+)\s*s", text, 0.0),
        "damage": first_float(r"damage:\s*([0-9.]+)\s*s", text, 0.0),
        "solve": first_float(r"solve:\s*([0-9.]+)\s*s", text, 0.0),
    }

def read_abaqus_timing(path):
    with open(path, "r") as f:
        text = f.read()
    submit = first_float(r"Solver wall-clock\s*\(submit->done\):\s*([0-9.]+)\s*s", text)
    internal = first_float(r"Abaqus WALLCLOCK\s*\([^)]*\):\s*([0-9.]+)\s*s", text, submit)
    rows = first_float(r"Rows:\s*([0-9.]+)", text, None)
    if submit is None and internal is None:
        raise ValueError("Could not parse Abaqus timing from %s" % path)
    return {
        "submit": submit if submit is not None else internal,
        "internal": internal if internal is not None else submit,
        "rows": int(rows) if rows is not None else None,
    }

def save(fig, stem, dpi=300):
    fig.savefig(os.path.join(OUT, stem + ".png"), dpi=dpi, facecolor="white", bbox_inches="tight")
    fig.savefig(os.path.join(OUT, stem + ".pdf"), facecolor="white", bbox_inches="tight")
    plt.close(fig)

def plot_load_cmod(cmod_m, load_m, cmod_a, load_a):
    peak_m_i = int(np.argmax(load_m))
    peak_a_i = int(np.argmax(load_a))
    peak_m = (cmod_m[peak_m_i], load_m[peak_m_i])
    peak_a = (cmod_a[peak_a_i], load_a[peak_a_i])

    fig, ax = plt.subplots(figsize=(9, 6.5), facecolor="white")
    ax.plot(cmod_m, load_m, label="MATLAB (FRACMATH)", color="#0F4C81", linewidth=2.5)
    ax.plot(cmod_a, load_a, label="Abaqus (UMAT)", color="#F25C54", linewidth=2.0, linestyle="--")
    ax.scatter([peak_m[0]], [peak_m[1]], color="#092f52", s=70, zorder=5, edgecolor="white", linewidth=1.5)
    ax.scatter([peak_a[0]], [peak_a[1]], color="#c2332b", s=70, zorder=5, edgecolor="white", linewidth=1.5)
    x_max = max(cmod_m.max(), cmod_a.max())
    y_max = max(load_m.max(), load_a.max())
    ax.annotate(
        "MATLAB Peak: %.1f N\nCMOD: %.4f mm" % (peak_m[1], peak_m[0]),
        xy=peak_m,
        xytext=(0.16 * x_max, 1.05 * y_max),
        arrowprops=dict(facecolor="#0F4C81", shrink=0.08, width=1.5, headwidth=6, headlength=6),
        fontsize=10,
        fontweight="bold",
        color="#0F4C81",
        bbox=dict(boxstyle="round,pad=0.3", fc="#f4f7f9", ec="#0F4C81", lw=0.5),
    )
    ax.annotate(
        "Abaqus Peak: %.1f N\nCMOD: %.4f mm" % (peak_a[1], peak_a[0]),
        xy=peak_a,
        xytext=(0.16 * x_max, 0.92 * y_max),
        arrowprops=dict(facecolor="#F25C54", shrink=0.08, width=1.5, headwidth=6, headlength=6),
        fontsize=10,
        fontweight="bold",
        color="#F25C54",
        bbox=dict(boxstyle="round,pad=0.3", fc="#fff5f5", ec="#F25C54", lw=0.5),
    )
    ax.set_title("Load vs. CMOD Comparison (Three-Point Bending)", fontsize=14, fontweight="bold", pad=15)
    ax.set_xlabel("Crack Mouth Opening Displacement, CMOD (mm)", fontsize=12, labelpad=10)
    ax.set_ylabel("Reaction Load (N)", fontsize=12, labelpad=10)
    ax.set_xlim(0, x_max * 1.05)
    ax.set_ylim(0, y_max * 1.15)
    ax.grid(True, linestyle=":", alpha=0.6, color="#888888")
    ax.legend(loc="upper right", frameon=True, facecolor="white", edgecolor="#e5e5e5", fontsize=11)
    fig.tight_layout()
    save(fig, "comparison_load_cmod")

    fig, ax = plt.subplots(figsize=(9, 6.5), facecolor="white")
    ax.plot(cmod_m, load_m, label="MATLAB (FRACMATH)", color="#0F4C81", linewidth=2.5)
    ax.plot(cmod_a, load_a, label="Abaqus (UMAT)", color="#F25C54", linewidth=2.0, linestyle="--")
    ax.set_title("Load vs. CMOD Comparison (Three-Point Bending)", fontsize=14, fontweight="bold", pad=15)
    ax.set_xlabel("CMOD (mm)")
    ax.set_ylabel("Reaction Load (N)")
    ax.grid(True, linestyle=":", alpha=0.6, color="#888888")
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "load_cmod_comparison.png"), dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)

def plot_performance(mt, at):
    matlab_speed_ratio = at["submit"] / mt["solver"]
    if matlab_speed_ratio >= 1.0:
        speed_label = "MATLAB approx %.1fx faster" % matlab_speed_ratio
    else:
        speed_label = "Abaqus approx %.1fx faster" % (1.0 / matlab_speed_ratio)
    other = max(mt["total"] - (mt["assembly"] + mt["damage"] + mt["solve"]), 0.0)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 5), facecolor="white")
    ax1.set_title("Wall-clock (%s)" % speed_label, fontsize=12, pad=12)
    times = [mt["solver"], at["submit"]]
    bars = ax1.bar(["MATLAB", "Abaqus"], times, color=["#f97316", "#1f77b4"], width=0.55, edgecolor="black", linewidth=0.8)
    ax1.set_ylabel("Wall-clock time [s]", fontsize=11)
    ax1.grid(axis="y", linestyle=":", alpha=0.5, color="#888888")
    ax1.set_axisbelow(True)
    ax1.set_ylim(0, max(times) * 1.14)
    for bar in bars:
        h = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width() / 2.0, h + max(times) * 0.025, "%.2f s" % h,
                 ha="center", va="bottom", fontsize=10, fontweight="bold")

    ax2.set_title("MATLAB time breakdown (%.1f s total)" % mt["total"], fontsize=12, pad=12)
    pie_labels = [
        "assembly\n%.2f s" % mt["assembly"],
        "damage\n%.2f s" % mt["damage"],
        "solve\n%.2f s" % mt["solve"],
        "other\n%.2f s" % other,
    ]
    wedges, texts, autotexts = ax2.pie(
        [mt["assembly"], mt["damage"], mt["solve"], other],
        labels=pie_labels,
        autopct="%1.1f%%",
        startangle=120,
        colors=["#8f94d4", "#d6c59b", "#81c784", "#b3b3b3"],
        textprops=dict(fontsize=10),
        pctdistance=0.6,
        labeldistance=1.1,
    )
    for atxt in autotexts:
        atxt.set_fontsize(9.5)
    fig.tight_layout()
    save(fig, "performance_breakdown")

    src = os.path.join(OUT, "performance_breakdown.png")
    dst = os.path.join(OUT, "time_comparison_bar.png")
    with open(src, "rb") as fsrc, open(dst, "wb") as fdst:
        fdst.write(fsrc.read())

    fig, ax = plt.subplots(figsize=(10, 6.5), facecolor="white")
    labels = [
        "Abaqus (UMAT)\nSubmit-to-Done",
        "MATLAB (FRACMATH)\nSolver",
        "Abaqus\n.msg/.dat WALLCLOCK",
        "MATLAB\nEnd-to-End",
    ]
    values = [at["submit"], mt["solver"], at["internal"], mt["total"]]
    bars = ax.bar(labels, values, color=["#E87070", "#5B8DB8", "#E87070", "#1B3A5C"], width=0.55)
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2.0, bar.get_height() + max(values) * 0.025,
                "%.2f s" % val, ha="center", va="bottom", fontsize=11, fontweight="bold")
    ax.set_ylabel("Runtime (seconds)", fontsize=12)
    ax.set_title("Computational Performance Comparison\n%s" % speed_label,
                 fontsize=13, fontweight="bold", pad=14)
    ax.set_ylim(0, max(values) * 1.15)
    ax.grid(axis="y", linestyle=":", alpha=0.55, color="#aaaaaa")
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(labelsize=10, bottom=False)
    fig.tight_layout()
    save(fig, "runtime_comparison")

def side_by_side(left_path, right_path, left_label, right_label, title, subtitle, stem):
    imgs = []
    for path in (left_path, right_path):
        imgs.append(mpimg.imread(path) if path and os.path.exists(path) else None)

    fig = plt.figure(figsize=(18, 7), facecolor="white")
    fig.text(0.5, 0.97, title, ha="center", va="top", fontsize=15, fontweight="bold", color="#222222")
    fig.text(0.5, 0.92, subtitle, ha="center", va="top", fontsize=10, color="#555555", style="italic")

    top = 0.88
    specs = [
        (0.01, 0.47, imgs[0], left_path, left_label, "#0F4C81"),
        (0.52, 0.47, imgs[1], right_path, right_label, "#C0392B"),
    ]
    for left, width, img, path, label, color in specs:
        ax = fig.add_axes([left, 0.06, width, top - 0.06])
        ax.axis("off")
        if img is not None:
            ax.imshow(img, aspect="equal")
        else:
            ax.text(0.5, 0.5, "Image not found\n%s" % path, ha="center", va="center",
                    fontsize=10, color="red", transform=ax.transAxes)
        bar = fig.add_axes([left, top, width, 0.035])
        bar.set_facecolor(color)
        bar.axis("off")
        bar.text(0.5, 0.45, label, ha="center", va="center", fontsize=12,
                 fontweight="bold", color="white", transform=bar.transAxes)

    fig.lines.extend([plt.Line2D([0.488, 0.488], [0.04, top + 0.035],
                                 transform=fig.transFigure, color="#cccccc", linewidth=1.5)])
    save(fig, stem, dpi=200)

def plot_image_comparisons(cmod_m, load_m, cmod_a, load_a):
    sub = "Three-Point Bending - Gregoire specimen - CDM model"
    side_by_side(
        os.path.join(MAT_R, "fig_damage_postpeak.png"),
        os.path.join(ABQ_R, "damage_postpeak.png"),
        "MATLAB - FRACMATH",
        "Abaqus - UMAT",
        "Crack geometry at Post-peak load (omega >= 0.99)",
        sub,
        "comparison_damage_postpeak",
    )
    side_by_side(
        os.path.join(MAT_R, "fig_damage_last_step.png"),
        os.path.join(ABQ_R, "damage_last_step.png"),
        "MATLAB - FRACMATH",
        "Abaqus - UMAT",
        "Crack geometry at Last step (omega >= 0.99)",
        sub,
        "comparison_damage_last_step",
    )
    abq_peak = os.path.join(ABQ_R, "damage_peak.png")
    side_by_side(
        os.path.join(MAT_R, "fig_damage_peak.png"),
        abq_peak if os.path.exists(abq_peak) else None,
        "MATLAB - FRACMATH",
        "Abaqus - UMAT",
        "Crack geometry at Peak load (omega >= 0.99)",
        sub,
        "comparison_damage_peak",
    )
    side_by_side(
        os.path.join(MAT_R, "fig_mesh.png"),
        None,
        "MATLAB - FRACMATH",
        "Abaqus - same mesh",
        "Finite element mesh (shared between both solvers)",
        "Gregoire 3PB - 7319 nodes - 14268 T3 elements",
        "comparison_mesh",
    )

    fig, axes = plt.subplots(1, 2, figsize=(16, 6), facecolor="white")
    fig.suptitle("Load vs. CMOD - Individual Panels\n" + sub, fontsize=13, fontweight="bold", y=0.98)
    for ax, cmod, load, col, solver in [
        (axes[0], cmod_m, load_m / 1000.0, "#0F4C81", "MATLAB - FRACMATH"),
        (axes[1], cmod_a, load_a / 1000.0, "#C0392B", "Abaqus - UMAT"),
    ]:
        ip = int(np.argmax(load))
        pk_c, pk_f = cmod[ip], load[ip]
        ax.fill_between(cmod, load, 0, color=col, alpha=0.10, linewidth=0)
        ax.plot(cmod, load, "-", color=col, linewidth=2.2)
        ax.plot(pk_c, pk_f, "*", markersize=14, markerfacecolor="#FFD700",
                markeredgecolor="#5a4000", markeredgewidth=0.8, zorder=5)
        x_max = float(cmod.max())
        ha = "left" if pk_c < 0.55 * x_max else "right"
        off = 0.015 * x_max * (1 if ha == "left" else -1)
        ax.annotate("Peak: %.2f kN\nCMOD: %.4f mm" % (pk_f, pk_c),
                    xy=(pk_c, pk_f), xytext=(pk_c + off, pk_f * 1.08),
                    fontsize=9, fontweight="bold", color=col, ha=ha, va="bottom",
                    bbox=dict(facecolor="white", edgecolor=col, linewidth=0.7, pad=3),
                    arrowprops=dict(arrowstyle="->", color=col, lw=0.9))
        ax.set_xlim(0, x_max * 1.05)
        ax.set_ylim(0, float(load.max()) * 1.30)
        ax.set_xlabel("CMOD [mm]", fontsize=11)
        ax.set_ylabel("Load [kN]", fontsize=11)
        ax.set_title(solver, fontsize=12, fontweight="bold", color=col, pad=8)
        ax.grid(True, linestyle=":", linewidth=0.5, color="#bbbbbb")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    save(fig, "comparison_load_cmod_panels")

def write_summary(cmod_m, load_m, cmod_a, load_a, mt, at):
    pm_i = int(np.argmax(load_m))
    pa_i = int(np.argmax(load_a))
    peak_m_kn = load_m[pm_i] / 1000.0
    peak_a_kn = load_a[pa_i] / 1000.0
    cmod_pm = cmod_m[pm_i]
    cmod_pa = cmod_a[pa_i]
    ratio_load = peak_m_kn / peak_a_kn
    ratio_cmod = cmod_pm / cmod_pa
    ratio_time = mt["solver"] / at["submit"]
    ratio_internal = mt["total"] / at["internal"]
    rows = at["rows"] if at["rows"] is not None else len(cmod_a)

    latex = r"""\begin{table}[htbp]
\centering
\caption{MATLAB vs. Abaqus on the 2D 3PB case (identical mesh, material, and tolerances). Abaqus ran on 4 threads, MATLAB on 1.}
\label{tab:comparison}
\begin{tabular}{lrrr}
\hline
Quantity & MATLAB & Abaqus + UMAT & Ratio \\
\hline
Peak load (kN) & %.2f & %.2f & %.3f \\
CMOD at peak (mm) & %.3f & %.3f & %.3f \\
Solver/submit wall-clock (s) & %.2f & %.2f & %.3f \\
End-to-end / internal time (s) & %.2f & %.2f & %.3f \\
\hline
\end{tabular}
\end{table}""" % (
        peak_m_kn, peak_a_kn, ratio_load,
        cmod_pm, cmod_pa, ratio_cmod,
        mt["solver"], at["submit"], ratio_time,
        mt["total"], at["internal"], ratio_internal,
    )

    summary = """# Abaqus vs MATLAB 3PB Comparison Summary

This report compares the computational results and performance of the three-point bending (3PB) simulation using Abaqus (UMAT) and the MATLAB solver (FRACMATH).

| Quantity | MATLAB | Abaqus + UMAT | Ratio |
| :--- | :---: | :---: | :---: |
| **Peak load (kN)** | %.2f | %.2f | %.3f |
| **CMOD at peak (mm)** | %.3f | %.3f | %.3f |
| **Solver/submit wall-clock (s)** | %.2f | %.2f | %.3f |
| **End-to-end / internal time (s)** | %.2f | %.2f | %.3f |
| **Load-CMOD rows** | %d | %d | - |

```latex
%s
```

- Combined bar and breakdown pie chart: [performance_breakdown.png](performance_breakdown.png)
- Load vs CMOD comparison plot: [comparison_load_cmod.png](comparison_load_cmod.png)
- Side-by-side damage/mesh panels: [comparison_damage_last_step.png](comparison_damage_last_step.png)
""" % (
        peak_m_kn, peak_a_kn, ratio_load,
        cmod_pm, cmod_pa, ratio_cmod,
        mt["solver"], at["submit"], ratio_time,
        mt["total"], at["internal"], ratio_internal,
        len(cmod_m), rows,
        latex,
    )
    with open(os.path.join(OUT, "comparison_summary.md"), "w") as f:
        f.write(summary)

    readme = """# 3PB comparison figures

This folder contains regenerated MATLAB and Abaqus/UMAT 3PB comparison figures.

Run from this folder:

```bash
python plot_comparison.py
```

- `../matlab/Gregoire_3PB/results/`
- `../abaqus/Gregoire_3PB/results/`

- MATLAB peak load: %.2f kN.
- Abaqus + UMAT peak load: %.2f kN.
- CMOD at peak: MATLAB %.6f mm; Abaqus %.6f mm.
- MATLAB solver wall-clock time: %.2f s.
- Abaqus submit-to-done wall-clock time: %.2f s.
- Abaqus internal `.msg/.dat` wall-clock time: %.2f s.
- Abaqus Load-CMOD rows: %d.
""" % (peak_m_kn, peak_a_kn, cmod_pm, cmod_pa, mt["solver"], at["submit"], at["internal"], rows)
    with open(os.path.join(OUT, "README.md"), "w") as f:
        f.write(readme)

def main():
    cmod_m, load_m = read_xy(os.path.join(MAT_R, "matlab_load_cmod.csv"))
    cmod_a, load_a = read_xy(os.path.join(ABQ_R, "abaqus_load_cmod.csv"))
    mt = read_matlab_timing(os.path.join(MAT_R, "matlab_timing.txt"))
    at = read_abaqus_timing(os.path.join(ABQ_R, "abaqus_timing.txt"))

    plot_load_cmod(cmod_m, load_m, cmod_a, load_a)
    plot_performance(mt, at)
    plot_image_comparisons(cmod_m, load_m, cmod_a, load_a)
    write_summary(cmod_m, load_m, cmod_a, load_a, mt, at)

    print("Saved comparison outputs to:")
    print("  %s" % OUT)
    print("Abaqus rows: %d" % (at["rows"] if at["rows"] is not None else len(cmod_a)))
    print("Abaqus submit-to-done: %.2f s" % at["submit"])

if __name__ == "__main__":
    main()
