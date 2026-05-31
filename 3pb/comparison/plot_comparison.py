import os
import csv
import matplotlib.pyplot as plt

# Define absolute paths
matlab_csv = r"C:\Users\jmavani\OneDrive - University of New Mexico\Desktop\JOSS\3pb_2\matlab\Gregoire_3PB\results\matlab_load_cmod.csv"
abaqus_csv = r"C:\Users\jmavani\OneDrive - University of New Mexico\Desktop\JOSS\3pb_2\abaqus\Gregoire_3PB\results\abaqus_load_cmod.csv"
output_dir = r"C:\Users\jmavani\OneDrive - University of New Mexico\Desktop\JOSS\3pb_2\comparison"

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Helper function to read CSV files
def read_csv(file_path):
    cmods = []
    loads = []
    if not os.path.exists(file_path):
        print(f"Error: File not found at {file_path}")
        return cmods, loads
    with open(file_path, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue
            first_cell = row[0].strip()
            if first_cell.startswith('#') or first_cell.startswith('cmod'):
                continue
            try:
                cmods.append(float(row[0]))
                loads.append(float(row[1]))
            except (ValueError, IndexError):
                continue
    return cmods, loads

# Load the data
print("Reading MATLAB data...")
cmods_matlab, loads_matlab = read_csv(matlab_csv)
print("Reading Abaqus data...")
cmods_abaqus, loads_abaqus = read_csv(abaqus_csv)

if not cmods_matlab or not cmods_abaqus:
    print("Error: Could not load data from one or both files.")
    exit(1)

# Extract peak values
peak_load_mat = max(loads_matlab)
peak_cmod_mat = cmods_matlab[loads_matlab.index(peak_load_mat)]

peak_load_abq = max(loads_abaqus)
peak_cmod_abq = cmods_abaqus[loads_abaqus.index(peak_load_abq)]

# Convert loads to kN for the table
peak_load_mat_kn = peak_load_mat / 1000.0
peak_load_abq_kn = peak_load_abq / 1000.0

# Performance runtimes  ── updated values from matlab_timing.txt
runtime_matlab_wallclock = 97.18    # solver wall-clock (compare to Abaqus)
runtime_matlab_total     = 103.90   # end-to-end incl. visualization
runtime_abaqus_wallclock = 1336.00  # Abaqus solver wall-clock
runtime_abaqus_total     = 1356.20  # Abaqus total process

# MATLAB breakdown components (from matlab_timing.txt)
assembly_time = 86.42
damage_time   = 2.69
solve_time    = 7.92
other_time    = runtime_matlab_total - (assembly_time + damage_time + solve_time)  # 6.87 s

# Ratio calculations (wall-clock based)
ratio_load = peak_load_mat_kn / peak_load_abq_kn
ratio_cmod = peak_cmod_mat / peak_cmod_abq
ratio_time = runtime_matlab_wallclock / runtime_abaqus_wallclock
speedup    = runtime_abaqus_wallclock / runtime_matlab_wallclock

# ----------------- PLOT 1: Load vs CMOD -----------------
print("Generating Load vs CMOD Comparison Plot...")
plt.figure(figsize=(9, 6.5))

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Liberation Sans']
plt.rcParams['axes.edgecolor'] = '#cccccc'
plt.rcParams['axes.linewidth'] = 0.8

plt.plot(cmods_matlab, loads_matlab, label='MATLAB (FRACMATH)', color='#0F4C81', linewidth=2.5, linestyle='-')
plt.plot(cmods_abaqus, loads_abaqus, label='Abaqus (UMAT)', color='#F25C54', linewidth=2.0, linestyle='--')

plt.scatter([peak_cmod_mat], [peak_load_mat], color='#092f52', s=70, zorder=5, edgecolor='white', linewidth=1.5)
plt.scatter([peak_cmod_abq], [peak_load_abq], color='#c2332b', s=70, zorder=5, edgecolor='white', linewidth=1.5)

plt.annotate(f"MATLAB Peak: {peak_load_mat:.1f} N\nCMOD: {peak_cmod_mat:.4f} mm",
             xy=(peak_cmod_mat, peak_load_mat),
             xytext=(peak_cmod_mat + 0.005, peak_load_mat - 400),
             arrowprops=dict(facecolor='#0F4C81', shrink=0.08, width=1.5, headwidth=6, headlength=6),
             fontsize=10, fontweight='bold', color='#0F4C81',
             bbox=dict(boxstyle="round,pad=0.3", fc="#f4f7f9", ec="#0F4C81", lw=0.5))

plt.annotate(f"Abaqus Peak: {peak_load_abq:.1f} N\nCMOD: {peak_cmod_abq:.4f} mm",
             xy=(peak_cmod_abq, peak_load_abq),
             xytext=(peak_cmod_abq + 0.015, peak_load_abq + 150),
             arrowprops=dict(facecolor='#F25C54', shrink=0.08, width=1.5, headwidth=6, headlength=6),
             fontsize=10, fontweight='bold', color='#F25C54',
             bbox=dict(boxstyle="round,pad=0.3", fc="#fff5f5", ec="#F25C54", lw=0.5))

plt.title('Load vs. CMOD Comparison (Three-Point Bending)', fontsize=14, fontweight='bold', pad=15, color='#333333')
plt.xlabel('Crack Mouth Opening Displacement, CMOD (mm)', fontsize=12, labelpad=10, color='#333333')
plt.ylabel('Reaction Load (N)', fontsize=12, labelpad=10, color='#333333')

plt.xlim(0, max(max(cmods_matlab), max(cmods_abaqus)) * 1.05)
plt.ylim(0, max(max(loads_matlab), max(loads_abaqus)) * 1.15)

plt.grid(True, linestyle=':', alpha=0.6, color='#888888')
plt.legend(loc='upper right', frameon=True, facecolor='white', edgecolor='#e5e5e5', fontsize=11)

plt.tight_layout()
plt.savefig(os.path.join(output_dir, "comparison_load_cmod.png"), dpi=300)
plt.savefig(os.path.join(output_dir, "comparison_load_cmod.pdf"), format='pdf')
plt.close()


# ----------------- PLOT 2: Bar + Pie (wall-clock vs MATLAB breakdown) -----------------
print("Generating Performance breakdown plot layout...")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 5))

# Left subplot: wall-clock bar chart
ax1.set_title(f'Wall-clock  (MATLAB $\\approx$ {speedup:.1f}$\\times$ faster)', fontsize=12, pad=12)

categories = ['MATLAB', 'Abaqus']
times = [runtime_matlab_wallclock, runtime_abaqus_wallclock]
bar_colors = ['#f97316', '#1f77b4']

bars = ax1.bar(categories, times, color=bar_colors, width=0.55, edgecolor='black', linewidth=0.8)
ax1.set_ylabel('Wall-clock time [s]', fontsize=11)
ax1.grid(axis='y', linestyle=':', alpha=0.5, color='#888888')
ax1.set_axisbelow(True)
ax1.set_ylim(0, max(times) * 1.12)
for bar in bars:
    h = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2.0, h + 20,
             f"{h:.2f} s", ha='center', va='bottom', fontsize=10, fontweight='bold')

# Right subplot: MATLAB breakdown pie
ax2.set_title(f'MATLAB time breakdown ({runtime_matlab_total:.1f} s total)', fontsize=12, pad=12)
pie_labels = [
    f'assembly\n{assembly_time:.2f} s',
    f'damage\n{damage_time:.2f} s',
    f'solve\n{solve_time:.2f} s',
    f'other\n{other_time:.2f} s'
]
pie_sizes  = [assembly_time, damage_time, solve_time, other_time]
pie_colors = ['#8f94d4', '#d6c59b', '#81c784', '#b3b3b3']
wedges, texts, autotexts = ax2.pie(
    pie_sizes, labels=pie_labels, autopct='%1.1f%%',
    startangle=120, colors=pie_colors,
    textprops=dict(fontsize=10),
    pctdistance=0.6, labeldistance=1.1)
for at in autotexts:
    at.set_fontsize(9.5)

plt.tight_layout()
for fname in ['performance_breakdown.png', 'time_comparison_bar.png']:
    plt.savefig(os.path.join(output_dir, fname), dpi=300)
for fname in ['performance_breakdown.pdf']:
    plt.savefig(os.path.join(output_dir, fname), format='pdf')
plt.close()
print(f"Saved: performance_breakdown.png / time_comparison_bar.png")

# ----------------- PLOT 3: 4-bar chart (wall-clock + total process) -----------------
print("Generating 4-bar runtime comparison plot...")
fig, ax = plt.subplots(figsize=(10, 6.5), facecolor='white')
ax.set_facecolor('white')

labels = [
    'Abaqus (UMAT)\nWall-Clock',
    'MATLAB (FRACMATH)\nWall-Clock',
    'Abaqus\nTotal Process',
    'MATLAB\nTotal Process'
]
values = [runtime_abaqus_wallclock, runtime_matlab_wallclock,
          runtime_abaqus_total,     runtime_matlab_total]
colors = ['#E87070', '#5B8DB8', '#E87070', '#1B3A5C']

bars4 = ax.bar(labels, values, color=colors, width=0.55,
               edgecolor='white', linewidth=0)
for bar, val in zip(bars4, values):
    ax.text(bar.get_x() + bar.get_width()/2.0,
            bar.get_height() + 18,
            f"{val:.2f} s",
            ha='center', va='bottom', fontsize=11, fontweight='bold', color='#222222')

ax.set_ylabel('Runtime (seconds)', fontsize=12, labelpad=8)
ax.set_title(
    f'Computational Performance Comparison\n(MATLAB is {speedup:.1f}x Faster)',
    fontsize=13, fontweight='bold', pad=14, color='#222222')
ax.set_ylim(0, max(values) * 1.13)
ax.grid(axis='y', linestyle=':', alpha=0.55, color='#aaaaaa')
ax.set_axisbelow(True)
for sp in ax.spines.values():
    sp.set_visible(False)
ax.tick_params(labelsize=10, bottom=False)

plt.tight_layout()
for fname in ['runtime_comparison.png']:
    plt.savefig(os.path.join(output_dir, fname), dpi=300, facecolor='white')
for fname in ['runtime_comparison.pdf']:
    plt.savefig(os.path.join(output_dir, fname), format='pdf', facecolor='white')
plt.close()
print(f"Saved: runtime_comparison.png")


# ----------------- Generate Markdown & LaTeX Comparison Summary -----------------
print("Generating Markdown and LaTeX Comparison Summary...")
summary_path = os.path.join(output_dir, "comparison_summary.md")

latex_table = f"""\\begin{{table}}[htbp]
\\centering
\\caption{{MATLAB vs. Abaqus on the 2D 3PB case (identical mesh, material, and tolerances). Abaqus ran on 4 threads, MATLAB on 1.}}
\\label{{tab:comparison}}
\\begin{{tabular}}{{lrrr}}
\\hline
Quantity & MATLAB & Abaqus + UMAT & Ratio \\\\
\\hline
Peak load (kN) & {peak_load_mat_kn:.2f} & {peak_load_abq_kn:.2f} & {ratio_load:.3f} \\\\
CMOD at peak (mm) & {peak_cmod_mat:.3f} & {peak_cmod_abq:.3f} & {ratio_cmod:.3f} \\\\
Wall-clock time (s) & {runtime_matlab_wallclock:.2f} & {runtime_abaqus_wallclock:.2f} & {ratio_time:.3f} \\\\
Total process (s) & {runtime_matlab_total:.2f} & {runtime_abaqus_total:.2f} & {runtime_matlab_total/runtime_abaqus_total:.3f} \\\\
\\hline
\\end{{tabular}}
\\end{{table}}"""

with open(summary_path, 'w') as sf:
    sf.write(f"""# Abaqus vs MATLAB 3PB Comparison Summary

This report compares the computational results and performance of the three-point bending (3PB) simulation using Abaqus (UMAT) and the MATLAB solver (FRACMATH).

## Summary Table (Markdown)

| Quantity | MATLAB | Abaqus + UMAT | Ratio |
| :--- | :---: | :---: | :---: |
| **Peak load (kN)** | {peak_load_mat_kn:.2f} | {peak_load_abq_kn:.2f} | {ratio_load:.3f} |
| **CMOD at peak (mm)** | {peak_cmod_mat:.3f} | {peak_cmod_abq:.3f} | {ratio_cmod:.3f} |
| **Wall-clock time (s)** | {runtime_matlab_wallclock:.2f} | {runtime_abaqus_wallclock:.2f} | {ratio_time:.3f} |
| **Total process time (s)** | {runtime_matlab_total:.2f} | {runtime_abaqus_total:.2f} | {runtime_matlab_total/runtime_abaqus_total:.3f} |

## LaTeX Table Format
```latex
{latex_table}
```

## Visualizations
- Combined bar and breakdown pie chart: [performance_breakdown.png](performance_breakdown.png)
- Load vs CMOD comparison plot: [comparison_load_cmod.png](comparison_load_cmod.png)
""")

print(f"Saved: {summary_path}")
print("All tasks completed successfully!")
