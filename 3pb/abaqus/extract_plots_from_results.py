
from __future__ import print_function

import glob
import os
import re
import runpy
import sys

MODEL = "Gregoire_3PB"
POSTPEAK_CMOD = 0.30
DEFAULT_THRESHOLD = 0.99

def script_dir():
    return os.path.abspath(os.path.dirname(__file__) if "__file__" in globals() else os.getcwd())

def ensure_dir(path):
    if not os.path.isdir(path):
        os.makedirs(path)

def norm_name(name):
    return str(name).strip().upper().replace("-", "_").replace(" ", "_")

def get_by_name_ci(container, wanted):
    if wanted in container.keys():
        return container[wanted]
    target = norm_name(wanted)
    for key in container.keys():
        if norm_name(key) == target:
            return container[key]
    return None

def find_odb(arg=None):
    if arg and os.path.exists(arg):
        return os.path.abspath(arg)

    base = script_dir()
    candidates = []
    candidates.extend(glob.glob(os.path.join(base, MODEL, "*.odb")))
    candidates.extend(glob.glob(os.path.join(base, "*.odb")))
    candidates = [p for p in candidates if os.path.exists(p)]
    if not candidates:
        raise RuntimeError("No .odb found. Pass the ODB path explicitly.")
    candidates.sort(key=os.path.getmtime)
    return os.path.abspath(candidates[-1])

def _flatten_nodes(obj, out):
    if hasattr(obj, "label"):
        out.add(int(obj.label))
        return
    try:
        for item in obj:
            _flatten_nodes(item, out)
    except TypeError:
        pass

def labels_from_nodeset(nodeset):
    labels = set()
    _flatten_nodes(nodeset.nodes, labels)
    return labels

def _last_int(text):
    values = re.findall(r"\d+", str(text))
    if not values:
        return None
    try:
        return int(values[-1])
    except Exception:
        return None

def history_node_label(hkey, hreg):
    label = _last_int(hkey)
    if label is not None:
        return label
    for attr in ("point", "position"):
        try:
            label = _last_int(str(getattr(hreg, attr)))
            if label is not None:
                return label
        except Exception:
            pass
    return None

def collect_history(step, labels, component, mode):
    labels = set(int(x) for x in labels)
    accum, count, nregion = {}, {}, 0
    for hkey, hreg in step.historyRegions.items():
        label = history_node_label(hkey, hreg)
        if label is None or label not in labels:
            continue
        if component not in hreg.historyOutputs.keys():
            continue
        nregion += 1
        for time_value, value in hreg.historyOutputs[component].data:
            t = round(float(time_value), 12)
            accum[t] = accum.get(t, 0.0) + float(value)
            count[t] = count.get(t, 0) + 1

    if nregion == 0 or not accum:
        return None

    rows = []
    for t in sorted(accum.keys()):
        value = accum[t]
        if mode == "avg":
            value = value / float(max(count.get(t, 1), 1))
        rows.append((t, value))
    return rows

def force_positive_load(rows):
    if not rows:
        return rows
    loads = [load for _, load in rows]
    if max(loads) <= 0.0 and min(loads) < 0.0:
        return [(cmod, -load) for cmod, load in rows]
    return rows

def make_load_cmod_rows(rf_hist, u1_hist, u2_hist):
    rf = dict(rf_hist)
    u1 = dict(u1_hist)
    u2 = dict(u2_hist)
    common_times = sorted(set(rf).intersection(u1).intersection(u2))
    rows, times = [], []

    if common_times:
        for t in common_times:
            rows.append((abs(u2[t] - u1[t]), -rf[t]))
            times.append(t)
        return times, force_positive_load(rows)

    n = min(len(rf_hist), len(u1_hist), len(u2_hist))
    for i in range(n):
        rows.append((abs(u2_hist[i][1] - u1_hist[i][1]), -rf_hist[i][1]))
        times.append(rf_hist[i][0])
    return times, force_positive_load(rows)

def float_component(data, idx):
    try:
        return float(data[idx])
    except Exception:
        return float(data)

def get_damage_field(frame):
    keys = frame.fieldOutputs.keys()
    for key in ("SDV2", "SDV_2", "STATEV2", "STATEV_2"):
        if key in keys:
            return frame.fieldOutputs[key], 0, key
    for key in keys:
        if norm_name(key) in ("SDV2", "SDV_2", "STATEV2", "STATEV_2"):
            return frame.fieldOutputs[key], 0, key
    for key in ("SDV", "STATEV"):
        if key in keys:
            return frame.fieldOutputs[key], 1, key
    return None, None, None

def frame_has_damage(frame):
    field, _, _ = get_damage_field(frame)
    if field is None:
        return False
    try:
        return len(field.values) > 0
    except Exception:
        return False

def closest_damage_frame_index(step, target_time, prefer_at_or_after=False):
    candidates = []
    for idx, frame in enumerate(step.frames):
        if not frame_has_damage(frame):
            continue
        t = float(frame.frameValue)
        if prefer_at_or_after and t < float(target_time):
            continue
        candidates.append((abs(t - float(target_time)), idx))

    if not candidates and prefer_at_or_after:
        return closest_damage_frame_index(step, target_time, prefer_at_or_after=False)
    if not candidates:
        return None
    candidates.sort()
    return candidates[0][1]

def last_damage_frame_index(step):
    for idx in range(len(step.frames) - 1, -1, -1):
        if frame_has_damage(step.frames[idx]):
            return idx
    return len(step.frames) - 1 if len(step.frames) else None

def write_omega_csv(frame, path):
    with open(path, "w") as handle:
        handle.write("# element_id, ip, omega\n")
        if frame is None:
            return 0, -1.0
        field, comp, _ = get_damage_field(frame)
        if field is None:
            return 0, -1.0

        count, max_omega = 0, -1.0e99
        for value in field.values:
            try:
                omega = float_component(value.data, comp)
            except Exception:
                continue
            ip = getattr(value, "integrationPoint", 1)
            handle.write("%s, %s, %.6e\n" % (str(value.elementLabel), str(ip), omega))
            count += 1
            if omega > max_omega:
                max_omega = omega
    return count, max_omega if max_omega > -1.0e98 else -1.0

def dump_mesh(odb, plotdata_dir):
    assembly = odb.rootAssembly
    instance = list(assembly.instances.values())[0]
    if len(assembly.instances.keys()) > 1:
        for name in assembly.instances.keys():
            if norm_name(name) != "ASSEMBLY":
                instance = assembly.instances[name]
                break

    with open(os.path.join(plotdata_dir, "mesh_nodes.csv"), "w") as handle:
        handle.write("# id, x, y\n")
        for node in instance.nodes:
            c = node.coordinates
            handle.write("%d, %.8f, %.8f\n" % (int(node.label), float(c[0]), float(c[1])))

    with open(os.path.join(plotdata_dir, "mesh_elements.csv"), "w") as handle:
        handle.write("# id, n1, n2, n3\n")
        for elem in instance.elements:
            con = elem.connectivity
            if len(con) >= 3:
                handle.write("%d, %d, %d, %d\n" %
                             (int(elem.label), int(con[0]), int(con[1]), int(con[2])))
    return instance.name

def parse_wall_clock(odb_path):
    base, _ = os.path.splitext(odb_path)
    for path in (base + ".msg", base + ".dat", base + ".sta"):
        if not os.path.exists(path):
            continue
        try:
            lines = open(path, "r").readlines()
        except Exception:
            continue
        for line in reversed(lines):
            up = line.upper()
            if "WALLCLOCK" in up or "WALL CLOCK" in up:
                nums = re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", line)
                if nums:
                    try:
                        return float(nums[-1]), os.path.basename(path)
                    except Exception:
                        pass
    return None, None

def previous_submit_time(results_dir):
    path = os.path.join(results_dir, "abaqus_timing.txt")
    if not os.path.exists(path):
        return None
    try:
        text = open(path, "r").read()
    except Exception:
        return None
    match = re.search(r"Solver wall-clock\s*\(submit->done\):\s*([0-9.]+)\s*s", text)
    return float(match.group(1)) if match else None

def write_timing(results_dir, odb_path, peak_load, peak_cmod, nrows):
    wall, source = parse_wall_clock(odb_path)
    submit = previous_submit_time(results_dir)
    with open(os.path.join(results_dir, "abaqus_timing.txt"), "w") as handle:
        handle.write("Gregoire 3PB job: %s\n" % odb_path)
        handle.write("Peak load:  %.2f N\n" % peak_load)
        handle.write("CMOD@peak:  %.6f mm\n" % peak_cmod)
        handle.write("Rows:       %d\n" % nrows)
        handle.write("Solver wall-clock (submit->done): %s   <-- compare to MATLAB\n" %
                     (("%.2f s" % submit) if submit is not None else "n/a"))
        handle.write("Abaqus WALLCLOCK (.msg/.dat):     %s\n" %
                     (("%.2f s (%s)" % (wall, source)) if wall is not None else "n/a"))

def run_plotter(results_dir):
    plotter = os.path.join(os.path.dirname(results_dir), "plot_results.py")
    if not os.path.exists(plotter):
        print("Plotter not found: %s" % plotter)
        return False

    old_argv = list(sys.argv)
    try:
        sys.argv = [plotter, results_dir]
        runpy.run_path(plotter, run_name="__main__")
        return True
    finally:
        sys.argv = old_argv

def extract_from_odb(odb_path, threshold=DEFAULT_THRESHOLD):
    try:
        from odbAccess import openOdb
    except Exception as exc:
        raise RuntimeError(
            "ODB extraction must be run with Abaqus Python:\n"
            "  abaqus python extract_plots_from_results.py\n"
            "Import error: %s" % exc
        )

    odb_path = find_odb(odb_path)
    case_dir = os.path.dirname(odb_path)
    results_dir = os.path.join(case_dir, "results")
    plotdata_dir = os.path.join(results_dir, "plotdata")
    ensure_dir(results_dir)
    ensure_dir(plotdata_dir)

    print("ODB:     %s" % odb_path)
    print("Results: %s" % results_dir)

    odb = openOdb(path=odb_path, readOnly=True)
    try:
        assembly = odb.rootAssembly
        step = get_by_name_ci(odb.steps, "Loading")
        if step is None:
            step = odb.steps[list(odb.steps.keys())[-1]]

        load_set = get_by_name_ci(assembly.nodeSets, "LOAD_NODES")
        cmod1 = get_by_name_ci(assembly.nodeSets, "CMOD1")
        cmod2 = get_by_name_ci(assembly.nodeSets, "CMOD2")
        if load_set is None or cmod1 is None or cmod2 is None:
            raise RuntimeError("Missing LOAD_NODES, CMOD1, or CMOD2 in the ODB.")

        rf = collect_history(step, labels_from_nodeset(load_set), "RF2", "sum")
        u1 = collect_history(step, labels_from_nodeset(cmod1), "U1", "avg")
        u2 = collect_history(step, labels_from_nodeset(cmod2), "U1", "avg")
        if rf is None or u1 is None or u2 is None:
            raise RuntimeError("History output RF2/U1 not found; cannot build load-CMOD.")

        times, rows = make_load_cmod_rows(rf, u1, u2)
        if not rows:
            raise RuntimeError("No load-CMOD rows could be assembled.")

        load_cmod_csv = os.path.join(results_dir, "abaqus_load_cmod.csv")
        with open(load_cmod_csv, "w") as handle:
            handle.write("# cmod[mm], load[N]\n")
            for cmod, load in rows:
                handle.write("%.6e, %.6e\n" % (cmod, load))

        peak_row = max(range(len(rows)), key=lambda i: rows[i][1])
        peak_cmod, peak_load = rows[peak_row]
        peak_time = times[peak_row]

        postpeak_time = times[-1]
        for idx, (cmod, _) in enumerate(rows):
            if cmod > POSTPEAK_CMOD:
                postpeak_time = times[idx]
                break

        peak_frame = closest_damage_frame_index(step, peak_time)
        postpeak_frame = closest_damage_frame_index(step, postpeak_time, prefer_at_or_after=True)
        final_frame = last_damage_frame_index(step)

        instance_name = dump_mesh(odb, plotdata_dir)
        snapshot_rows = []
        for frame_idx, name in (
            (peak_frame, "omega_peak.csv"),
            (postpeak_frame, "omega_postpeak.csv"),
            (final_frame, "omega_last.csv"),
        ):
            frame = step.frames[frame_idx] if frame_idx is not None else None
            count, max_omega = write_omega_csv(frame, os.path.join(plotdata_dir, name))
            frame_time = float(frame.frameValue) if frame is not None else -1.0
            snapshot_rows.append((name, frame_idx, frame_time, count, max_omega))

        write_timing(results_dir, odb_path, peak_load, peak_cmod, len(rows))

        diagnostic = os.path.join(results_dir, "abaqus_extraction_diagnostic.txt")
        with open(diagnostic, "w") as handle:
            handle.write("ODB: %s\n" % odb_path)
            handle.write("Step: %s\n" % list(odb.steps.keys())[-1])
            handle.write("Instance: %s\n" % instance_name)
            handle.write("Load-CMOD rows: %d\n" % len(rows))
            handle.write("Peak load: %.6f N\n" % peak_load)
            handle.write("CMOD at peak: %.9f mm\n" % peak_cmod)
            handle.write("Peak history time: %.9f\n" % peak_time)
            handle.write("Postpeak target CMOD: %.6f mm\n" % POSTPEAK_CMOD)
            handle.write("Postpeak history time: %.9f\n" % postpeak_time)
            handle.write("Damage threshold used by plotter: %.2f\n" % threshold)
            handle.write("\nSnapshots written from nearest available damage field frame:\n")
            for name, frame_idx, frame_time, count, max_omega in snapshot_rows:
                handle.write("  %-20s frame=%s time=%.9f rows=%d max_omega=%.6f\n" %
                             (name, str(frame_idx), frame_time, count, max_omega))

        print("Wrote %s" % load_cmod_csv)
        print("Wrote plotdata to %s" % plotdata_dir)
        print("Wrote %s" % diagnostic)

    finally:
        odb.close()

    try:
        ok = run_plotter(results_dir)
    except Exception as exc:
        ok = False
        print("Plotting inside this Python failed: %s" % exc)

    if not ok:
        print("Now regenerate plots with normal Python:")
        print('  python "%s" --plot "%s"' % (os.path.abspath(__file__), results_dir))

def plot_only(results_arg=None):
    if results_arg and os.path.isdir(results_arg):
        results_dir = os.path.abspath(results_arg)
    else:
        results_dir = os.path.join(script_dir(), MODEL, "results")
    if not os.path.isdir(results_dir):
        raise RuntimeError("Results folder not found: %s" % results_dir)
    return run_plotter(results_dir)

def main():
    args = sys.argv[1:]
    if "--help" in args or "-h" in args:
        print(__doc__)
        return

    if "--plot" in args:
        idx = args.index("--plot")
        results_arg = args[idx + 1] if idx + 1 < len(args) else None
        plot_only(results_arg)
        return

    odb_arg = args[0] if args else None
    extract_from_odb(odb_arg)

if __name__ == "__main__":
    main()
