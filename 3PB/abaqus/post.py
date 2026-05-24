# -*- coding: utf-8 -*-
# =============================================================================
# post.py  -- Robust post-processing for Gregoire_3PB Abaqus ODB
#
# FIX FOR YOUR ERROR:
#   Your run_3pb_abaqus.py deletes all field output requests and keeps only
#   HISTORY output. Therefore step.frames can be 0 even when the job completes.
#   This script extracts Load-CMOD from step.historyRegions, not step.frames.
#
# RUN:
#   cd Gregoire_3PB
#   abaqus python ..\post.py
#
# OUTPUT:
#   ./results/load_cmod_data.csv
#   ./results/load_vs_cmod.png
#   ./results/load_vs_cmod.pdf
#
# NOTE:
#   Damage contour needs field output SDV/U/RF frames. Your current driver
#   intentionally does NOT write field outputs, so damage plotting is skipped
#   unless those outputs exist in the ODB.
# =============================================================================

from __future__ import print_function

import os
import sys
import math
import numpy as np

CONFIG = {
    # Job / paths
    "odb_name"        : "Gregoire_3PB",
    "odb_path"        : None,          # None = auto-search current folder and ./Gregoire_3PB
    "output_dir"      : "./results",

    # Step
    "step_name"       : "Loading",

    # Node set names from run_3pb_abaqus.py
    "load_nset"       : "LOAD_NODES",
    "cmod_node_left"  : "CMOD1",
    "cmod_node_right" : "CMOD2",

    # History variable names
    "load_hist_key"   : "RF2",
    "cmod_hist_key"   : "U1",

    # Sign and units
    "load_sign"       : -1.0,          # downward reaction becomes positive load
    "load_scale_plot" : 1.0e-3,        # N -> kN
    "cmod_abs"        : True,

    # Plot
    "dpi"             : 200,
    "fig_size_curve"  : (7, 5),

    # Optional damage plot. This only works if the ODB has field frames + SDV2.
    "try_damage_plot" : True,
    "damage_var"      : "SDV2",
    "instance_name"   : "BEAMINST",
    "frame_index"     : -1,
    "fig_size_damage" : (8, 5),
    "colormap"        : "hot_r",
}

try:
    from odbAccess import openOdb
    import abaqusConstants
    ABAQUS_ENV = True
except Exception:
    ABAQUS_ENV = False

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable


# =============================================================================
# Utilities
# =============================================================================

def _upper(s):
    return str(s).upper()

def find_odb_path(cfg):
    """Find ODB robustly whether post.py is run from parent folder or case folder."""
    if cfg.get("odb_path"):
        p = cfg["odb_path"]
        if os.path.exists(p):
            return p
        raise IOError("ODB not found at configured odb_path: " + p)

    name = cfg["odb_name"]
    candidates = [
        name + ".odb",
        os.path.join(name, name + ".odb"),
        os.path.join(os.getcwd(), name + ".odb"),
        os.path.join(os.getcwd(), name, name + ".odb"),
    ]

    # Also try script folder, useful when post.py is outside the case folder.
    try:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        candidates += [
            os.path.join(script_dir, name + ".odb"),
            os.path.join(script_dir, name, name + ".odb"),
        ]
    except Exception:
        pass

    for p in candidates:
        if os.path.exists(p):
            return os.path.abspath(p)

    raise IOError(
        "ODB not found. Tried:\n  " + "\n  ".join(candidates) +
        "\nRun from the folder containing Gregoire_3PB.odb or set CONFIG['odb_path']."
    )

def resolve_output_dir(cfg, odb_path):
    """Save results next to the ODB by default."""
    out_dir = cfg["output_dir"]
    if not os.path.isabs(out_dir):
        out_dir = os.path.join(os.path.dirname(os.path.abspath(odb_path)), out_dir)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    return out_dir

def get_node_labels_from_set(assembly, set_name):
    """
    Abaqus assembly nodeSets[NAME].nodes is usually a tuple of node arrays,
    one per instance. This returns all labels in that set.
    """
    set_name = _upper(set_name)
    if set_name not in assembly.nodeSets.keys():
        raise KeyError("Node set not found in ODB rootAssembly: " + set_name)

    labels = set()
    nset = assembly.nodeSets[set_name]

    for node_array in nset.nodes:
        for nd in node_array:
            labels.add(int(nd.label))

    if not labels:
        raise ValueError("Node set exists but contains no nodes: " + set_name)

    return labels

def parse_history_node_label(history_region_key):
    """
    Common Abaqus history region keys look like:
      'Node BEAMINST.123'
      'Node PART-1-1.45'
    Return the final integer label, or None.
    """
    if not history_region_key.startswith("Node "):
        return None
    try:
        return int(history_region_key.split(".")[-1])
    except Exception:
        return None

def add_history_series(series, data):
    """Add data values into an existing [(time,value)] series."""
    data = list(data)
    if series is None:
        return [(float(t), float(v)) for t, v in data]

    n = min(len(series), len(data))
    out = []
    for i in range(n):
        out.append((float(series[i][0]), float(series[i][1]) + float(data[i][1])))
    return out

def sum_history_for_nodes(step, labels, variable_name):
    """Sum a history output over all node labels."""
    labels = set([int(x) for x in labels])
    series = None
    matched = 0

    for hkey, hreg in step.historyRegions.items():
        lbl = parse_history_node_label(hkey)
        if lbl is None or lbl not in labels:
            continue
        if variable_name not in hreg.historyOutputs.keys():
            continue

        series = add_history_series(series, hreg.historyOutputs[variable_name].data)
        matched += 1

    return series, matched

def avg_history_for_nodes(step, labels, variable_name):
    """Average a history output over all node labels."""
    series, matched = sum_history_for_nodes(step, labels, variable_name)
    if series is None or matched == 0:
        return None, 0
    return [(t, v / float(matched)) for (t, v) in series], matched

def print_history_debug(step):
    print("[DEBUG] Available history regions / outputs:")
    count = 0
    for hkey, hreg in step.historyRegions.items():
        outs = list(hreg.historyOutputs.keys())
        print("  " + str(hkey) + " -> " + ", ".join(outs))
        count += 1
        if count >= 25:
            print("  ... more history regions not printed")
            break


# =============================================================================
# Load-CMOD extraction from HISTORY OUTPUT
# =============================================================================

def extract_load_cmod_from_history(cfg):
    odb_path = find_odb_path(cfg)
    print("[INFO] Opening ODB: " + odb_path)

    odb = openOdb(odb_path, readOnly=True)

    try:
        step_name = cfg["step_name"]
        if step_name is None:
            step_name = list(odb.steps.keys())[-1]
        if step_name not in odb.steps.keys():
            raise KeyError("Step not found: " + str(step_name) +
                           ". Available steps: " + ", ".join(odb.steps.keys()))

        step = odb.steps[step_name]
        print("[INFO] Step: '" + step_name + "'")
        print("[INFO] Field frames in step      : " + str(len(step.frames)))
        print("[INFO] History regions in step  : " + str(len(step.historyRegions)))

        asm = odb.rootAssembly

        load_labels  = get_node_labels_from_set(asm, cfg["load_nset"])
        cmod1_labels = get_node_labels_from_set(asm, cfg["cmod_node_left"])
        cmod2_labels = get_node_labels_from_set(asm, cfg["cmod_node_right"])

        print("[INFO] LOAD_NODES labels: " + str(sorted(load_labels)))
        print("[INFO] CMOD1 labels     : " + str(sorted(cmod1_labels)))
        print("[INFO] CMOD2 labels     : " + str(sorted(cmod2_labels)))

        rf, n_rf = sum_history_for_nodes(step, load_labels, cfg["load_hist_key"])
        u1, n_u1 = avg_history_for_nodes(step, cmod1_labels, cfg["cmod_hist_key"])
        u2, n_u2 = avg_history_for_nodes(step, cmod2_labels, cfg["cmod_hist_key"])

        print("[INFO] Matched RF2 node histories : " + str(n_rf))
        print("[INFO] Matched CMOD1 U1 histories : " + str(n_u1))
        print("[INFO] Matched CMOD2 U1 histories : " + str(n_u2))

        if rf is None or u1 is None or u2 is None:
            print_history_debug(step)
            raise RuntimeError(
                "Missing required history output. Your ODB must contain RF2 at LOAD_NODES "
                "and U1 at CMOD1/CMOD2. Re-run the analysis after confirming "
                "HistoryOutputRequest exists in run_3pb_abaqus.py."
            )

        n = min(len(rf), len(u1), len(u2))
        if n == 0:
            raise RuntimeError("History outputs were found, but contain zero data points.")

        times = np.zeros(n, dtype=float)
        loads = np.zeros(n, dtype=float)
        cmods = np.zeros(n, dtype=float)

        for i in range(n):
            times[i] = float(rf[i][0])
            loads[i] = float(cfg["load_sign"]) * float(rf[i][1])
            cmod = float(u2[i][1]) - float(u1[i][1])
            cmods[i] = abs(cmod) if cfg["cmod_abs"] else cmod

        print("[INFO] Extracted " + str(n) + " Load-CMOD data points.")
        return odb_path, times, loads, cmods

    finally:
        odb.close()


# =============================================================================
# Optional damage extraction from FIELD OUTPUT
# =============================================================================

def extract_damage_field_if_available(cfg, odb_path):
    """
    Damage plot needs field frames and SDV2. Your current driver deletes
    field output, so this will usually return None and skip safely.
    """
    odb = openOdb(odb_path, readOnly=True)
    try:
        step_name = cfg["step_name"]
        if step_name is None:
            step_name = list(odb.steps.keys())[-1]
        step = odb.steps[step_name]

        if len(step.frames) == 0:
            print("[INFO] Damage plot skipped: ODB has 0 field frames.")
            return None

        frame = step.frames[cfg["frame_index"]]
        if cfg["damage_var"] not in frame.fieldOutputs.keys():
            print("[INFO] Damage plot skipped: " + cfg["damage_var"] +
                  " is not in fieldOutputs for this frame.")
            print("[INFO] Available field outputs: " +
                  ", ".join(frame.fieldOutputs.keys()))
            return None

        inst_name = _upper(cfg["instance_name"])
        if inst_name not in odb.rootAssembly.instances.keys():
            print("[INFO] Damage plot skipped: instance not found: " + inst_name)
            print("[INFO] Available instances: " +
                  ", ".join(odb.rootAssembly.instances.keys()))
            return None

        inst = odb.rootAssembly.instances[inst_name]

        node_map = {}
        for node in inst.nodes:
            node_map[int(node.label)] = (float(node.coordinates[0]), float(node.coordinates[1]))

        connectivity = []
        elem_labels = []
        for elem in inst.elements:
            connectivity.append(tuple([int(x) for x in elem.connectivity]))
            elem_labels.append(int(elem.label))

        dmg_field = frame.fieldOutputs[cfg["damage_var"]]
        try:
            dmg_subset = dmg_field.getSubset(region=inst, position=abaqusConstants.CENTROID)
        except Exception:
            dmg_subset = dmg_field.getSubset(region=inst)

        elem_damage = {}
        for val in dmg_subset.values:
            if not hasattr(val, "elementLabel"):
                continue
            d = val.data
            try:
                d_val = float(np.mean(d))
            except Exception:
                d_val = float(d)
            elem_damage[int(val.elementLabel)] = d_val

        coords = np.array([node_map[k] for k in sorted(node_map.keys())], dtype=float)
        node_idx = {}
        for i, k in enumerate(sorted(node_map.keys())):
            node_idx[k] = i

        ex, ey, ed = [], [], []
        for lbl, conn in zip(elem_labels, connectivity):
            pts = [node_map[n] for n in conn if n in node_map]
            if not pts:
                continue
            xs = [p[0] for p in pts]
            ys = [p[1] for p in pts]
            ex.append(sum(xs) / float(len(xs)))
            ey.append(sum(ys) / float(len(ys)))
            ed.append(elem_damage.get(lbl, 0.0))

        return np.array(ex), np.array(ey), np.array(ed), coords, connectivity, node_idx

    finally:
        odb.close()


# =============================================================================
# Save + plot
# =============================================================================

def save_csv(times, loads, cmods, out_dir):
    path = os.path.join(out_dir, "load_cmod_data.csv")
    data = np.column_stack([times, loads, cmods])
    np.savetxt(path, data, delimiter=",", header="time,load_N,cmod_mm", comments="")
    print("[SAVED] " + path)

def plot_load_cmod(times, loads, cmods, cfg, out_dir):
    if len(loads) == 0:
        raise RuntimeError("No Load-CMOD rows to plot.")

    load_kN = loads * float(cfg["load_scale_plot"])

    peak_idx = int(np.argmax(load_kN))
    peak_load = float(load_kN[peak_idx])
    peak_cmod = float(cmods[peak_idx])

    fig, ax = plt.subplots(figsize=cfg["fig_size_curve"])
    ax.plot(cmods, load_kN, "-", linewidth=2.0, label="Abaqus + SCM UMAT")
    ax.plot(peak_cmod, peak_load, marker="*", markersize=16, linestyle="None",
            label="Peak: %.3f kN @ %.4f mm" % (peak_load, peak_cmod))

    ax.set_xlabel("CMOD [mm]", fontsize=12)
    ax.set_ylabel("Load [kN]", fontsize=12)
    ax.set_title("Abaqus 3PB - Load vs CMOD", fontsize=12)
    ax.grid(True, alpha=0.3, linestyle=":")
    ax.legend(loc="best", frameon=True, framealpha=0.95, fontsize=10)

    xmax = float(np.max(cmods))
    ymax = float(np.max(load_kN))
    ax.set_xlim(left=0.0, right=xmax * 1.05 if xmax > 0.0 else 1.0)
    ax.set_ylim(bottom=0.0, top=ymax * 1.20 if ymax > 0.0 else 1.0)

    fig.tight_layout()

    png_path = os.path.join(out_dir, "load_vs_cmod.png")
    pdf_path = os.path.join(out_dir, "load_vs_cmod.pdf")
    fig.savefig(png_path, dpi=cfg["dpi"], bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    print("[SAVED] " + png_path)
    print("[SAVED] " + pdf_path)
    print("[INFO] Peak load: %.6f kN" % peak_load)
    print("[INFO] CMOD@peak: %.6f mm" % peak_cmod)

def plot_damage(damage_data, cfg, out_dir):
    if damage_data is None:
        return

    ex, ey, ed, coords, connectivity, node_idx = damage_data

    fig, ax = plt.subplots(figsize=cfg["fig_size_damage"])

    for conn in connectivity:
        pts = [coords[node_idx[n]] for n in conn if n in node_idx]
        if len(pts) < 2:
            continue
        poly = plt.Polygon(pts, edgecolor="0.70", facecolor="none", linewidth=0.25, zorder=1)
        ax.add_patch(poly)

    norm = Normalize(vmin=0.0, vmax=1.0)
    ax.scatter(ex, ey, c=ed, cmap=cfg["colormap"], norm=norm, s=5, zorder=2, linewidths=0)

    crack_mask = np.array(ed) >= 0.95
    if crack_mask.any():
        ax.scatter(ex[crack_mask], ey[crack_mask], s=12, marker="x", zorder=3,
                   label="Crack zone D >= 0.95")
        ax.legend(fontsize=8, loc="upper right")

    cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cfg["colormap"]), ax=ax)
    cbar.set_label("Damage (" + str(cfg["damage_var"]) + ")", fontsize=10)

    ax.set_aspect("equal")
    ax.autoscale()
    ax.set_xlabel("X [mm]", fontsize=10)
    ax.set_ylabel("Y [mm]", fontsize=10)
    ax.set_title("Damage contour - final available field frame", fontsize=11)

    fig.tight_layout()
    png_path = os.path.join(out_dir, "damage_plot.png")
    pdf_path = os.path.join(out_dir, "damage_plot.pdf")
    fig.savefig(png_path, dpi=cfg["dpi"], bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    print("[SAVED] " + png_path)
    print("[SAVED] " + pdf_path)


# =============================================================================
# Plain Python replot mode
# =============================================================================

def replot_from_csv(cfg):
    out_dir = cfg["output_dir"]
    csv_path = os.path.join(out_dir, "load_cmod_data.csv")
    if not os.path.exists(csv_path):
        raise IOError("CSV not found: " + csv_path)

    data = np.loadtxt(csv_path, delimiter=",", skiprows=1)
    if data.ndim == 1:
        data = data.reshape((1, -1))

    times = data[:, 0]
    loads = data[:, 1]
    cmods = data[:, 2]
    plot_load_cmod(times, loads, cmods, cfg, out_dir)


# =============================================================================
# Main
# =============================================================================

def main():
    cfg = CONFIG

    if not ABAQUS_ENV:
        if "--no-odb" in sys.argv or "--replot" in sys.argv:
            print("[INFO] Running plain-Python replot mode from saved CSV.")
            replot_from_csv(cfg)
            print("[INFO] Done.")
            return

        print("[ERROR] odbAccess is not available.")
        print("        Run ODB extraction with:")
        print("        abaqus python post.py")
        print("        Or replot existing CSV with:")
        print("        python post.py --replot")
        sys.exit(1)

    print("=" * 60)
    print("  STEP 1 / 3 - Extracting Load & CMOD from HISTORY output")
    print("=" * 60)
    odb_path, times, loads, cmods = extract_load_cmod_from_history(cfg)
    out_dir = resolve_output_dir(cfg, odb_path)
    save_csv(times, loads, cmods, out_dir)

    print("")
    print("=" * 60)
    print("  STEP 2 / 3 - Plotting Load vs CMOD")
    print("=" * 60)
    plot_load_cmod(times, loads, cmods, cfg, out_dir)

    print("")
    print("=" * 60)
    print("  STEP 3 / 3 - Optional damage contour")
    print("=" * 60)
    if cfg.get("try_damage_plot", True):
        damage_data = extract_damage_field_if_available(cfg, odb_path)
        plot_damage(damage_data, cfg, out_dir)
    else:
        print("[INFO] Damage plot disabled by CONFIG['try_damage_plot'].")

    print("")
    print("OK Done. Output files written to: " + os.path.abspath(out_dir))


if __name__ == "__main__":
    main()
