# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
run_3pb_abaqus_OLIVER_T3.py
=========================================
ONE-FILE Abaqus workflow for the Gregoire 3PB CDM/UMAT case.
Keeps ONLY what is needed for damage and response figures:
    1) Load vs CMOD curve
    2) Damage plot at peak load: fully damaged elements only
    3) Damage plot at post-peak: fully damaged elements only
    4) Damage plot at last step: fully damaged elements only

Two modes:
  Abaqus mode (default):  abaqus cae noGUI=run_3pb_single_build_run_plot.py
      build -> write Oliver T3 gradN table -> run (cdm_umat_2d_OLIVER_T3.for) -> extract Load-CMOD + omega -> auto-plot
  Plot mode:              python run_3pb_single_build_run_plot.py --plot Gregoire_3PB/results
      draw MATLAB-style figures with matplotlib (also auto-called above)

Memory-lean: field output is SDV only (no U/RF/S/E tensors); the smooth
Load-CMOD curve comes from per-increment HISTORY output, so coarse field
frequency is fine. Damage is plotted on the UNDEFORMED mesh, so no
displacement field is stored or dumped.

Outputs (Gregoire_3PB/results/):
    abaqus_load_cmod.csv
    abaqus_timing.txt
    abaqus_load_cmod_fig.png/.pdf
    abaqus_fig_damage_peak.png/.pdf       fully damaged elements only
    abaqus_fig_damage_postpeak.png/.pdf   fully damaged elements only
    abaqus_fig_damage_last_step.png/.pdf  fully damaged elements only
    plotdata/  (mesh + omega snapshots used by the plotter)
    oliver_t3_gradN.dat  (T3 gradients used by the UMAT for Oliver h(n))

Fast environment options:
    set ABQ_CPUS=4          number of Abaqus CPUs/domains
    set ABQ_FIELD_FREQ=100  write SDV field every 100 increments
    set ABQ_AUTO_PLOT=1     make figures immediately after run; default 0
"""

from __future__ import print_function

# =====================================================================
# VERSION FAST
# Fix: Abaqus MeshElement.connectivity can contain zero-based part.nodes
# indices, not node labels. This file avoids direct node_xy dictionary lookup
# from raw element connectivity.
# =====================================================================

import os
import sys
import re
import time

PLOT_MODE = ('--plot' in sys.argv)

if not PLOT_MODE:
    try:
        from abaqus import mdb
        from abaqusConstants import *
        import mesh as abqMesh
        from odbAccess import openOdb
    except Exception as e:
        print('ERROR: run with Abaqus/CAE Python:')
        print('  abaqus cae noGUI=run_3pb_single_build_run_plot.py')
        print('Import error: %s' % str(e))
        sys.exit(1)


# =====================================================================
# CONFIG
# =====================================================================
MODEL = 'Gregoire_3PB'

# Geometry [mm]
D        = 100.0
S        = 2.5 * D
OVERHANG = 0.5 * D
L        = S + 2.0 * OVERHANG
B_THICK  = 50.0
A0       = 0.2 * D
WN       = D / 40.0
XC       = L / 2.0
XL_NOTCH = XC - WN / 2.0
XR_NOTCH = XC + WN / 2.0

# Mesh (matched to MATLAB solver)
ELEM_SIZE_GLOBAL = D / 16.0
ELEM_SIZE_REFINE = ELEM_SIZE_GLOBAL / 5.0
REFINE_W         = 0.5 * D
REFINE_H         = D

# UMAT material constants: E, nu, ft, GF, fc/ft
E_C, NU_C, FT, GF, FCFT = 37000.0, 0.20, 3.5, 0.090, 10.0

# Loading
U_FINAL = -0.2
N_INC   = int(os.environ.get('ABQ_N_INC', '1000'))  # use 1000 for MATLAB comparison
CPUS    = int(os.environ.get('ABQ_CPUS', '4'))       # set ABQ_CPUS=1 for fair MATLAB timing

# Field output: SDV only, coarse frequency (damage snapshots only).
FIELD_VARS = ('SDV',)
FIELD_FREQ = int(os.environ.get('ABQ_FIELD_FREQ', '100'))  # higher = less ODB output = faster

# Plotting
AUTO_PLOT     = (os.environ.get('ABQ_AUTO_PLOT', '0') == '1')  # 0 fastest; run --plot later
# Full path to a REAL Python (CPython/Anaconda) that has numpy+matplotlib.
# Do NOT use bare 'python' on Windows: it can resolve to Abaqus's ABQcaeK.exe,
# which is a binary, not an interpreter -> "Non-UTF-8 code \x90" SyntaxError.
# Leave as '' to auto-detect, or hardcode e.g. r'C:\Users\you\anaconda3\python.exe'
PYTHON_EXE    = r''         # e.g. r'C:\Users\you\anaconda3\python.exe'
POSTPEAK_CMOD = 0.3         # post-peak snapshot at first CMOD > this
FULL_DAMAGE_THRESHOLD = 0.95   # match MATLAB final figures: show fully damaged elements


# =====================================================================
# UTILITIES
# =====================================================================
def die(msg):
    print('ERROR: ' + str(msg))
    sys.exit(1)


def ensure_dir(path):
    if not os.path.isdir(path):
        os.makedirs(path)


def norm_name(name):
    return str(name).strip().upper().replace('-', '_')


def script_base_dir():
    return os.path.abspath(os.getcwd())


def find_umat_file(base_dir):
    """Prefer the optimized Oliver-table UMAT, but keep fallback names."""
    preferred = (
        'cdm_umat_2d_OLIVER_T3_FAST.for',
        'cdm_umat_2d_OLIVER_T3.for',
        'cdm_umat_2d_oliver_t3.for',
        'cdm_umat_2d.for',
    )
    for nm in preferred:
        pth = os.path.join(base_dir, nm)
        if os.path.exists(pth):
            return pth
    try:
        names = os.listdir(base_dir)
    except Exception:
        names = []
    hits = []
    for nm in names:
        low = nm.lower()
        if low.startswith('cdm_umat_2d') and low.endswith(('.for', '.f', '.f90')):
            hits.append(os.path.join(base_dir, nm))
    # Put Oliver names first if several UMATs are present.
    hits.sort(key=lambda p: (0 if 'oliver' in os.path.basename(p).lower() else 1, os.path.basename(p).lower()))
    return hits[0] if hits else None


def print_last_lines(path, n=60):
    if not os.path.exists(path):
        return
    print('\n--- Last lines of %s ---' % os.path.basename(path))
    try:
        with open(path, 'r') as f:
            lines = f.readlines()
        for ln in lines[-n:]:
            print(ln.rstrip())
    except Exception as e:
        print('Could not read %s: %s' % (path, str(e)))
    print('--- end ---\n')



# =====================================================================
# OLIVER T3 BANDWIDTH TABLE FOR UMAT
# =====================================================================
def _detect_connectivity_mode(part):
    """
    Abaqus/CAE can return MeshElement.connectivity either as node labels
    or as zero-based indices into part.nodes, depending on how the mesh was
    created.  The previous version assumed labels, so it crashed when the
    first connectivity entry was 0.  This detector makes the Oliver table
    writer work for both cases.
    """
    labels = set(int(nd.label) for nd in part.nodes)
    saw_zero = False
    max_index = len(part.nodes) - 1
    checked = 0

    for el in part.elements:
        con = list(el.connectivity)
        if len(con) != 3:
            continue
        checked += 1
        for c in con:
            ic = int(c)
            if ic == 0:
                saw_zero = True
            if ic < 0:
                return 'label'
        if checked >= 100:
            break

    if saw_zero:
        return 'index'

    # If there is no zero, labels are the safest interpretation for orphan
    # meshes and for .inp-style numbering.  Native CAE meshes with indices
    # normally show zero in at least one early element.
    return 'label'


def _node_label_and_xy(part, token, mode):
    """Return (node_label, x, y) from either a label or a zero-based index."""
    itok = int(token)
    if mode == 'index':
        if itok < 0 or itok >= len(part.nodes):
            die('Element connectivity index %d is outside part.nodes range 0..%d' %
                (itok, len(part.nodes)-1))
        nd = part.nodes[itok]
        return int(nd.label), float(nd.coordinates[0]), float(nd.coordinates[1])

    # label mode
    try:
        nd = part.nodes.sequenceFromLabels(labels=(itok,))[0]
    except Exception:
        # Last-resort fallback: if label lookup fails but the value is a valid
        # node-array index, use index mode for this node.  This prevents a hard
        # crash on unusual CAE meshes.
        if 0 <= itok < len(part.nodes):
            nd = part.nodes[itok]
            return int(nd.label), float(nd.coordinates[0]), float(nd.coordinates[1])
        die('Node label %d from element connectivity was not found.' % itok)
    return int(nd.label), float(nd.coordinates[0]), float(nd.coordinates[1])


def write_oliver_t3_gradN_from_part(part, out_path):
    """
    Write element-wise T3 shape-function gradients for the UMAT.

    The UMAT cannot access all element nodal coordinates directly, so this
    table supplies the same gradients used by MATLAB:
        grad(N1) = [g1x, g1y], grad(N2) = [g2x, g2y], grad(N3) = [g3x, g3y]
    Then the UMAT computes the current direction-dependent Oliver bandwidth:
        h(n) = 2 / (|gradN1.n| + |gradN2.n| + |gradN3.n|)

    IMPORTANT FIX:
    Abaqus/CAE part.elements[i].connectivity may be zero-based node indices,
    not node labels.  This function detects that and converts correctly.
    """
    mode = _detect_connectivity_mode(part)
    print('Abaqus element connectivity interpreted as: %s' %
          ('zero-based part.nodes indices' if mode == 'index' else 'node labels'))

    n_written = 0
    with open(out_path, 'w') as f:
        f.write('# NOEL, g1x, g1y, g2x, g2y, g3x, g3y\n')
        for el in part.elements:
            con = list(el.connectivity)
            if len(con) != 3:
                continue

            _, x1, y1 = _node_label_and_xy(part, con[0], mode)
            _, x2, y2 = _node_label_and_xy(part, con[1], mode)
            _, x3, y3 = _node_label_and_xy(part, con[2], mode)

            area_signed = 0.5 * ((x2 - x1) * (y3 - y1) -
                                 (x3 - x1) * (y2 - y1))
            if abs(area_signed) <= 1.0e-14:
                die('Zero/near-zero CPS3 area while writing Oliver table, element %d' % int(el.label))

            b1 = y2 - y3
            b2 = y3 - y1
            b3 = y1 - y2
            c1 = x3 - x2
            c2 = x1 - x3
            c3 = x2 - x1
            inv2A = 1.0 / (2.0 * area_signed)

            g1x = b1 * inv2A
            g1y = c1 * inv2A
            g2x = b2 * inv2A
            g2y = c2 * inv2A
            g3x = b3 * inv2A
            g3y = c3 * inv2A

            f.write('%d, %.16e, %.16e, %.16e, %.16e, %.16e, %.16e\n' %
                    (int(el.label), g1x, g1y, g2x, g2y, g3x, g3y))
            n_written += 1

    print('Wrote Oliver T3 gradN table: %s  (%d CPS3 elements)' % (out_path, n_written))
    if n_written == 0:
        die('Oliver T3 gradN table is empty. Check that the mesh uses CPS3 triangles.')
    return out_path


# =====================================================================
# MODEL BUILDING
# =====================================================================
def pick_single_nearest(part, xt, yt):
    best = None
    bd = 1.0e99
    for nd in part.nodes:
        x, y = nd.coordinates[0], nd.coordinates[1]
        d2 = (x - xt) ** 2 + (y - yt) ** 2
        if d2 < bd:
            bd = d2
            best = nd
    if best is None:
        die('No node near x=%g y=%g' % (xt, yt))
    return part.nodes.sequenceFromLabels(labels=(best.label,))


def pick_n_nearest_top(part, xt, y_top, n_keep):
    tol = 1.0e-3 * abs(y_top) + 1.0e-6
    cands = []
    for nd in part.nodes:
        x, y = nd.coordinates[0], nd.coordinates[1]
        if abs(y - y_top) <= tol:
            cands.append((abs(x - xt), nd.label))
    cands.sort()
    if not cands:
        die('No top-edge nodes for load set.')
    labels = [lbl for _, lbl in cands[:max(1, n_keep)]]
    return part.nodes.sequenceFromLabels(labels=tuple(labels))


def build_model(case_dir, umat_file):
    print('\n==== 1) Build model ====')
    if MODEL in mdb.models.keys():
        del mdb.models[MODEL]
    m = mdb.Model(name=MODEL)

    sk = m.ConstrainedSketch(name='BeamSketch', sheetSize=2.0 * L)
    sk.Line(point1=(0.0, 0.0),       point2=(XL_NOTCH, 0.0))
    sk.Line(point1=(XL_NOTCH, 0.0),  point2=(XL_NOTCH, A0))
    sk.Line(point1=(XL_NOTCH, A0),   point2=(XR_NOTCH, A0))
    sk.Line(point1=(XR_NOTCH, A0),   point2=(XR_NOTCH, 0.0))
    sk.Line(point1=(XR_NOTCH, 0.0),  point2=(L, 0.0))
    sk.Line(point1=(L, 0.0),         point2=(L, D))
    sk.Line(point1=(L, D),           point2=(0.0, D))
    sk.Line(point1=(0.0, D),         point2=(0.0, 0.0))

    part = m.Part(name='Beam2D', dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
    part.BaseShell(sketch=sk)
    del sk

    tform = part.MakeSketchTransform(sketchPlane=part.faces[0], sketchPlaneSide=SIDE1,
                                     origin=(0.0, 0.0, 0.0))
    sk_r = m.ConstrainedSketch(name='RefineZone', sheetSize=L, transform=tform)
    sk_r.rectangle(point1=(XC - REFINE_W / 2.0, 0.0),
                   point2=(XC + REFINE_W / 2.0, REFINE_H))
    part.PartitionFaceBySketch(sketch=sk_r, faces=part.faces)
    del sk_r

    mat = m.Material(name='Concrete')
    mat.UserMaterial(mechanicalConstants=(E_C, NU_C, FT, GF, FCFT))
    mat.Depvar(n=2)   # fastest: kappa, omega only. Use n=4 only when debugging Oliver h/flag
    m.HomogeneousSolidSection(name='BeamSec', material='Concrete', thickness=B_THICK)
    part.SectionAssignment(region=(part.faces,), sectionName='BeamSec')

    elem_t3 = abqMesh.ElemType(elemCode=CPS3, elemLibrary=STANDARD)
    part.setElementType(regions=(part.faces,), elemTypes=(elem_t3,))
    for f in part.faces:
        part.setMeshControls(regions=(f,), technique=FREE, elemShape=TRI,
                             algorithm=ADVANCING_FRONT)
    part.seedPart(size=ELEM_SIZE_GLOBAL, deviationFactor=0.1)
    refine_edges = part.edges.getByBoundingBox(
        xMin=XC - REFINE_W / 2.0 - 1.0e-3, xMax=XC + REFINE_W / 2.0 + 1.0e-3,
        yMin=-1.0e-3, yMax=REFINE_H + 1.0e-3)
    if len(refine_edges):
        part.seedEdgeBySize(edges=refine_edges, size=ELEM_SIZE_REFINE, constraint=FINER)
    part.generateMesh()
    print('Mesh: %d CPS3, %d nodes' % (len(part.elements), len(part.nodes)))

    node_tol = 0.5 * ELEM_SIZE_REFINE
    part.Set(nodes=pick_single_nearest(part, OVERHANG, 0.0),     name='Support_Left')
    part.Set(nodes=pick_single_nearest(part, L - OVERHANG, 0.0), name='Support_Right')
    part.Set(nodes=pick_n_nearest_top(part, XC, D, 3),           name='Load_Nodes')

    cmod_L = part.nodes.getByBoundingBox(xMin=XL_NOTCH - node_tol, xMax=XL_NOTCH + node_tol,
                                         yMin=-node_tol, yMax=node_tol)
    cmod_R = part.nodes.getByBoundingBox(xMin=XR_NOTCH - node_tol, xMax=XR_NOTCH + node_tol,
                                         yMin=-node_tol, yMax=node_tol)
    if len(cmod_L) == 0 or len(cmod_R) == 0:
        die('CMOD notch-lip node sets empty; tolerance too small.')
    part.Set(nodes=cmod_L, name='CMOD1')
    part.Set(nodes=cmod_R, name='CMOD2')

    asm = m.rootAssembly
    asm.DatumCsysByDefault(CARTESIAN)
    inst = asm.Instance(name='BeamInst', part=part, dependent=ON)
    for nm in ('Support_Left', 'Support_Right', 'Load_Nodes', 'CMOD1', 'CMOD2'):
        asm.Set(nodes=inst.sets[nm].nodes, name=nm)

    m.StaticStep(name='Loading', previous='Initial', maxNumInc=20000,
                 initialInc=1.0 / float(N_INC), minInc=1.0e-10,
                 maxInc=1.0 / float(N_INC), nlgeom=OFF)

    m.DisplacementBC(name='BC_SupL', createStepName='Loading',
                     region=asm.sets['Support_Left'], u1=SET, u2=SET)
    m.DisplacementBC(name='BC_SupR', createStepName='Loading',
                     region=asm.sets['Support_Right'], u2=SET)
    m.DisplacementBC(name='BC_Load', createStepName='Loading',
                     region=asm.sets['Load_Nodes'], u2=U_FINAL)

    # Field: SDV only (damage snapshots). History: load + CMOD (per increment).
    m.FieldOutputRequest(name='F-SDV', createStepName='Loading',
                         variables=FIELD_VARS, frequency=FIELD_FREQ)
    m.HistoryOutputRequest(name='H-Load', createStepName='Loading',
                           variables=('RF2',), region=asm.sets['Load_Nodes'])
    m.HistoryOutputRequest(name='H-CMOD1', createStepName='Loading',
                           variables=('U1',), region=asm.sets['CMOD1'])
    m.HistoryOutputRequest(name='H-CMOD2', createStepName='Loading',
                           variables=('U1',), region=asm.sets['CMOD2'])
    print('Build complete.')
    return part


# =====================================================================
# JOB
# =====================================================================
def run_job(case_dir, umat_file):
    print('\n==== 2) Write input + run ====')
    print('UMAT: %s' % umat_file)
    if MODEL in mdb.jobs.keys():
        del mdb.jobs[MODEL]
    job = mdb.Job(name=MODEL, model=MODEL, description='Gregoire 3PB CDM UMAT',
                  userSubroutine=umat_file, numCpus=CPUS, numDomains=CPUS)
    cwd0 = os.getcwd()
    ensure_dir(case_dir)
    os.chdir(case_dir)
    try:
        job.writeInput(consistencyChecking=OFF)
        t0 = time.time()
        job.submit(consistencyChecking=OFF)
        job.waitForCompletion()
        solve_wall = time.time() - t0
        print('Job done in %.2f s.' % solve_wall)
    finally:
        os.chdir(cwd0)
    odb_path = os.path.join(case_dir, MODEL + '.odb')
    if not os.path.exists(odb_path):
        for ext in ('.dat', '.msg', '.sta'):
            print_last_lines(os.path.join(case_dir, MODEL + ext))
        die('ODB not created. Job failed (see lines above).')
    print('ODB: %s' % odb_path)
    return odb_path, solve_wall


# =====================================================================
# ODB EXTRACTION (Load-CMOD via history; omega via SDV field)
# =====================================================================
def get_by_name_ci(container, wanted):
    if wanted in container.keys():
        return container[wanted]
    tgt = norm_name(wanted)
    for k in container.keys():
        if norm_name(k) == tgt:
            return container[k]
    return None


def _flatten_nodes(obj, out):
    if hasattr(obj, 'label'):
        out.add(int(obj.label))
        return
    try:
        for it in obj:
            _flatten_nodes(it, out)
    except TypeError:
        pass


def labels_from_nodeset(ns):
    s = set()
    _flatten_nodes(ns.nodes, s)
    return s


def _last_int(text):
    vals = re.findall(r'\d+', str(text))
    if not vals:
        return None
    try:
        return int(vals[-1])
    except Exception:
        return None


def history_node_label(hkey, hreg):
    lbl = _last_int(hkey)
    if lbl is not None:
        return lbl
    for attr in ('point', 'position'):
        try:
            lbl = _last_int(str(getattr(hreg, attr)))
            if lbl is not None:
                return lbl
        except Exception:
            pass
    return None


def collect_history(step, labels, component, mode):
    labels = set(int(x) for x in labels)
    accum, count, nreg = {}, {}, 0
    for hkey, hreg in step.historyRegions.items():
        lbl = history_node_label(hkey, hreg)
        if lbl is None or lbl not in labels:
            continue
        if component not in hreg.historyOutputs.keys():
            continue
        nreg += 1
        for t, v in hreg.historyOutputs[component].data:
            tt = round(float(t), 12)
            accum[tt] = accum.get(tt, 0.0) + float(v)
            count[tt] = count.get(tt, 0) + 1
    if nreg == 0 or not accum:
        return None
    out = []
    for tt in sorted(accum.keys()):
        val = accum[tt]
        if mode == 'avg':
            val = val / float(max(count.get(tt, 1), 1))
        out.append((tt, val))
    return out


def force_positive_load(rows):
    if not rows:
        return rows
    loads = [r for _, r in rows]
    if max(loads) <= 0.0 and min(loads) < 0.0:
        return [(c, -r) for c, r in rows]
    return rows


def make_rows(rf_hist, u1_hist, u2_hist):
    rf, u1, u2 = dict(rf_hist), dict(u1_hist), dict(u2_hist)
    common = sorted(set(rf).intersection(u1).intersection(u2))
    rows, times = [], []
    if common:
        for t in common:
            rows.append((abs(u2[t] - u1[t]), -rf[t]))
            times.append(t)
        return times, force_positive_load(rows)
    n = min(len(rf_hist), len(u1_hist), len(u2_hist))
    for i in range(n):
        rows.append((abs(u2_hist[i][1] - u1_hist[i][1]), -rf_hist[i][1]))
        times.append(rf_hist[i][0])
    return times, force_positive_load(rows)


def closest_frame_index(step, target_time):
    if len(step.frames) == 0:
        return None
    best_i, best_dt = 0, 1.0e99
    for i, fr in enumerate(step.frames):
        dt = abs(float(fr.frameValue) - float(target_time))
        if dt < best_dt:
            best_dt, best_i = dt, i
    return best_i


def float_component(data, idx):
    try:
        return float(data[idx])
    except Exception:
        return float(data)


def get_damage_field(frame):
    keys = frame.fieldOutputs.keys()
    for k in ('SDV2', 'SDV_2', 'STATEV2', 'STATEV_2'):
        if k in keys:
            return frame.fieldOutputs[k], 0, k
    for k in keys:
        if norm_name(k) in ('SDV2', 'SDV_2', 'STATEV2', 'STATEV_2'):
            return frame.fieldOutputs[k], 0, k
    for k in ('SDV', 'STATEV'):
        if k in keys:
            return frame.fieldOutputs[k], 1, k
    return None, None, None


def write_omega_csv(frame, path):
    if frame is None:
        open(path, 'w').write('# element_id, ip, omega\n')
        return False
    fld, comp, _ = get_damage_field(frame)
    if fld is None:
        open(path, 'w').write('# element_id, ip, omega\n')
        return False
    with open(path, 'w') as f:
        f.write('# element_id, ip, omega\n')
        for v in fld.values:
            try:
                w = float_component(v.data, comp)
            except Exception:
                continue
            ip = getattr(v, 'integrationPoint', 1)
            f.write('%s, %s, %.6e\n' % (str(v.elementLabel), str(ip), w))
    return True


def parse_wall_clock(odb_path):
    base, _ = os.path.splitext(odb_path)
    for path in (base + '.msg', base + '.dat', base + '.sta'):
        if not os.path.exists(path):
            continue
        try:
            with open(path, 'r') as f:
                lines = f.readlines()
        except Exception:
            continue
        for ln in reversed(lines):
            up = ln.upper()
            if 'WALLCLOCK' in up or 'WALL CLOCK' in up:
                nums = re.findall(r'[-+]?\d*\.\d+|[-+]?\d+', ln)
                if nums:
                    try:
                        return float(nums[-1]), os.path.basename(path)
                    except Exception:
                        pass
    return None, None


# =====================================================================
# PLOTDATA DUMP  (mesh + omega snapshots only)
# =====================================================================
def dump_plotdata(odb, step, res_dir, peak_idx, postpeak_idx, last_idx):
    pd_dir = os.path.join(res_dir, 'plotdata')
    ensure_dir(pd_dir)

    asm = odb.rootAssembly
    inst = list(asm.instances.values())[0]
    if len(asm.instances.keys()) > 1:
        for nm in asm.instances.keys():
            if norm_name(nm) != 'ASSEMBLY':
                inst = asm.instances[nm]
                break

    with open(os.path.join(pd_dir, 'mesh_nodes.csv'), 'w') as f:
        f.write('# id, x, y\n')
        for nd in inst.nodes:
            c = nd.coordinates
            f.write('%d, %.8f, %.8f\n' % (int(nd.label), float(c[0]), float(c[1])))

    with open(os.path.join(pd_dir, 'mesh_elements.csv'), 'w') as f:
        f.write('# id, n1, n2, n3\n')
        for el in inst.elements:
            con = el.connectivity
            if len(con) >= 3:
                f.write('%d, %d, %d, %d\n' %
                        (int(el.label), int(con[0]), int(con[1]), int(con[2])))

    def frame_at(idx):
        if idx is None or idx < 0 or idx >= len(step.frames):
            return None
        return step.frames[idx]

    write_omega_csv(frame_at(peak_idx),     os.path.join(pd_dir, 'omega_peak.csv'))
    write_omega_csv(frame_at(postpeak_idx), os.path.join(pd_dir, 'omega_postpeak.csv'))
    write_omega_csv(frame_at(last_idx),     os.path.join(pd_dir, 'omega_last.csv'))
    print('Wrote plotdata to %s' % pd_dir)


def _this_script_path():
    try:
        p = os.path.abspath(__file__)
        if os.path.exists(p):
            return p
    except Exception:
        pass
    for cand in (sys.argv[0] if sys.argv else '',
                 os.path.join(os.getcwd(), 'run_3pb_single_build_run_plot.py')):
        if cand and os.path.exists(cand):
            return os.path.abspath(cand)
    return os.path.abspath('run_3pb_single_build_run_plot.py')


def _find_real_python():
    """Return path to a real CPython/Anaconda interpreter (NOT Abaqus EXE)."""
    import glob
    # 1) user override
    if PYTHON_EXE:
        return PYTHON_EXE
    # 2) common locations / PATH lookups, skipping anything inside SIMULIA/Abaqus
    cands = []
    for ev in ('PYTHON_EXE_OVERRIDE',):
        if os.environ.get(ev):
            cands.append(os.environ[ev])
    # PATH search via 'where' (Windows) / 'which' (posix)
    try:
        import subprocess
        finder = 'where' if os.name == 'nt' else 'which'
        out = subprocess.check_output([finder, 'python'],
                                      stderr=subprocess.STDOUT)
        for ln in out.decode('utf-8', 'ignore').splitlines():
            ln = ln.strip()
            if ln:
                cands.append(ln)
    except Exception:
        pass
    # typical Windows installs
    if os.name == 'nt':
        home = os.environ.get('USERPROFILE', 'C:\\')
        for pat in (os.path.join(home, 'anaconda3', 'python.exe'),
                    os.path.join(home, 'miniconda3', 'python.exe'),
                    os.path.join(home, 'AppData', 'Local', 'Programs',
                                 'Python', 'Python3*', 'python.exe'),
                    'C:\\Python3*\\python.exe',
                    'C:\\ProgramData\\anaconda3\\python.exe'):
            cands.extend(glob.glob(pat))
    # filter out Abaqus/SIMULIA binaries
    for c in cands:
        low = c.lower()
        if 'simulia' in low or 'abaqus' in low or 'abq' in low:
            continue
        if os.path.exists(c):
            return c
    return None


def invoke_plots(res_dir):
    print('\n==== 4) MATLAB-style figures (matplotlib) ====')
    # Preferred: run the plotter IN-PROCESS (no subprocess, no python-path issue).
    try:
        import numpy as _np            # noqa: F401
        import matplotlib as _mpl       # noqa: F401
        print('Plotting in-process (Abaqus Python has matplotlib).')
        plot_main(os.path.abspath(res_dir))
        return
    except Exception as e:
        print('In-process plot unavailable (%s). Trying external Python...'
              % str(e))
    # Fallback: spawn a REAL Python (never the Abaqus EXE).
    py = _find_real_python()
    if not py:
        print('WARNING: no standalone Python with matplotlib found.')
        print('Set PYTHON_EXE at top of script, then run manually:')
        print('  <python> "%s" --plot "%s"'
              % (_this_script_path(), os.path.abspath(res_dir)))
        return
    cmd = '"%s" "%s" --plot "%s"' % (py, _this_script_path(),
                                     os.path.abspath(res_dir))
    print('Running: %s' % cmd)
    rc = os.system(cmd)
    if rc != 0:
        print('WARNING: figure step returned %d. Run manually:' % rc)
        print('  %s' % cmd)


# =====================================================================
# EXTRACT DRIVER
# =====================================================================
def extract_and_plot(odb_path, solve_wall=None):
    print('\n==== 3) Extract Load-CMOD + omega ====')
    res_dir = os.path.join(os.path.dirname(os.path.abspath(odb_path)), 'results')
    ensure_dir(res_dir)

    odb = None
    try:
        odb = openOdb(odb_path, readOnly=True)
        asm = odb.rootAssembly

        step = get_by_name_ci(odb.steps, 'Loading')
        if step is None:
            if not odb.steps.keys():
                die('ODB has no steps.')
            step = odb.steps[list(odb.steps.keys())[-1]]

        load_set = get_by_name_ci(asm.nodeSets, 'LOAD_NODES')
        c1 = get_by_name_ci(asm.nodeSets, 'CMOD1')
        c2 = get_by_name_ci(asm.nodeSets, 'CMOD2')
        if load_set is None or c1 is None or c2 is None:
            die('Missing LOAD_NODES / CMOD1 / CMOD2 in ODB.')

        L_lbl = labels_from_nodeset(load_set)
        c1_lbl = labels_from_nodeset(c1)
        c2_lbl = labels_from_nodeset(c2)

        rf = collect_history(step, L_lbl, 'RF2', 'sum')
        u1 = collect_history(step, c1_lbl, 'U1', 'avg')
        u2 = collect_history(step, c2_lbl, 'U1', 'avg')
        if rf is None or u1 is None or u2 is None:
            die('History output (RF2/U1) not found. Cannot build Load-CMOD.')
        times, rows = make_rows(rf, u1, u2)
        if not rows:
            die('Could not assemble Load-CMOD rows.')

        # Load-CMOD csv
        csv_path = os.path.join(res_dir, 'abaqus_load_cmod.csv')
        with open(csv_path, 'w') as f:
            f.write('# cmod[mm], load[N]\n')
            for c, r in rows:
                f.write('%.6e, %.6e\n' % (c, r))
        print('Wrote ' + csv_path)

        # peak + post-peak frame indices
        peak_idx_row, peak_load = 0, rows[0][1]
        for i, row in enumerate(rows):
            if row[1] > peak_load:
                peak_load, peak_idx_row = row[1], i
        peak_cmod = rows[peak_idx_row][0]
        peak_frame = closest_frame_index(step, times[peak_idx_row])

        pp_time = None
        for i, (c, r) in enumerate(rows):
            if c > POSTPEAK_CMOD:
                pp_time = times[i]
                break
        pp_frame = (len(step.frames) - 1) if pp_time is None \
            else closest_frame_index(step, pp_time)

        # Last available field-output frame: used for the final damage geometry.
        last_frame = len(step.frames) - 1

        dump_plotdata(odb, step, res_dir, peak_frame, pp_frame, last_frame)

        # timing
        wall, src = parse_wall_clock(odb_path)
        with open(os.path.join(res_dir, 'abaqus_timing.txt'), 'w') as f:
            f.write('Gregoire 3PB job: %s\n' % odb_path)
            f.write('Peak load:  %.2f N\n' % peak_load)
            f.write('CMOD@peak:  %.6f mm\n' % peak_cmod)
            f.write('Rows:       %d\n' % len(rows))
            f.write('Solver wall-clock (submit->done): %s   <-- compare to MATLAB\n'
                    % (('%.2f s' % solve_wall) if solve_wall else 'n/a'))
            f.write('Abaqus WALLCLOCK (.msg/.dat):     %s\n'
                    % (('%.2f s (%s)' % (wall, src)) if wall else 'n/a'))

        odb.close()
        odb = None

        print('\n==== Summary ====')
        print('Peak load : %.2f N' % peak_load)
        print('CMOD@peak : %.6f mm' % peak_cmod)
        print('Results   : %s' % res_dir)

        if AUTO_PLOT:
            invoke_plots(res_dir)
    finally:
        try:
            if odb is not None:
                odb.close()
        except Exception:
            pass


# =====================================================================
# =====================================================================
#   MATPLOTLIB SECTION (--plot mode). Two figures, MATLAB style.
# =====================================================================
# =====================================================================
def _mpl():
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import LinearSegmentedColormap
    return np, plt, PolyCollection, LinearSegmentedColormap


def crack_cmap(LinearSegmentedColormap):
    pos = [0.00, 0.10, 0.30, 0.50, 0.70, 0.85, 1.00]
    stops = [(0.84, 0.88, 0.95), (0.18, 0.42, 0.86), (0.05, 0.72, 0.88),
             (0.18, 0.80, 0.32), (0.96, 0.90, 0.08), (0.98, 0.44, 0.04),
             (0.82, 0.04, 0.04)]
    cd = {'red': [], 'green': [], 'blue': []}
    for p, (r, g, b) in zip(pos, stops):
        cd['red'].append((p, r, r))
        cd['green'].append((p, g, g))
        cd['blue'].append((p, b, b))
    return LinearSegmentedColormap('crack', cd, N=256)


def _save_hq(fh, base):
    fh.savefig(base + '.png', dpi=600, facecolor='w')
    fh.savefig(base + '.pdf', dpi=600, facecolor='w')


def _read_csv(np, path, ncol):
    rows = []
    if not os.path.exists(path):
        return np.zeros((0, ncol))
    with open(path, 'r') as f:
        for ln in f:
            s = ln.strip()
            if not s or s.startswith('#'):
                continue
            parts = s.replace(',', ' ').split()
            try:
                vals = [float(x) for x in parts[:ncol]]
            except ValueError:
                continue
            if len(vals) >= ncol:
                rows.append(vals)
    return np.array(rows, dtype=float) if rows else np.zeros((0, ncol))


def _load_plotdata(np, pd_dir):
    nraw = _read_csv(np, os.path.join(pd_dir, 'mesh_nodes.csv'), 3)
    eraw = _read_csv(np, os.path.join(pd_dir, 'mesh_elements.csv'), 4)
    nlbl = nraw[:, 0].astype(int)
    coords = nraw[:, 1:3]
    nmap = {}
    for i in range(len(nlbl)):
        nmap[int(nlbl[i])] = i
    conn = np.array([[nmap[int(c)] for c in row[1:4]] for row in eraw], dtype=int)
    elbl = eraw[:, 0].astype(int)
    emap = {}
    for i in range(len(elbl)):
        emap[int(elbl[i])] = i

    def omega(name):
        d = _read_csv(np, os.path.join(pd_dir, name), 3)
        w = np.zeros(conn.shape[0])
        for row in d:
            lbl = int(row[0])
            if lbl in emap:
                w[emap[lbl]] = row[2]
        return w

    return {'nodes': coords, 'elems': conn,
            'w_peak': omega('omega_peak.csv'),
            'w_pp': omega('omega_postpeak.csv'),
            'w_last': omega('omega_last.csv')}


def _fig_damage(np, plt, PolyCollection, cmap, nodes, elems, omega, label, base):
    """
    Clean damage plot:
      - grey full mesh
      - ONLY fully damaged elements are shown in black
      - no box around the specimen axes
      - no boxed legend or boxed information panel
      - right-side text is kept inside the saved PNG/PDF canvas
    """
    from matplotlib.lines import Line2D

    verts = nodes[elems]  # undeformed mesh, same style as MATLAB
    x_min, x_max = nodes[:, 0].min(), nodes[:, 0].max()
    y_min, y_max = nodes[:, 1].min(), nodes[:, 1].max()

    th = FULL_DAMAGE_THRESHOLD
    idx = np.where(omega >= th)[0]

    fh = plt.figure(figsize=(24.0 / 2.54, 8.0 / 2.54), facecolor='w')

    # Main specimen axis. The right side is reserved for clean text.
    ax = fh.add_axes([0.065, 0.17, 0.715, 0.74])

    mesh_fc = (0.94, 0.94, 0.95)
    mesh_ec = (0.82, 0.84, 0.86)
    crack_fc = (0.00, 0.00, 0.00)

    # Layer 1: complete mesh in light grey
    ax.add_collection(PolyCollection(verts,
                                     facecolors=mesh_fc,
                                     edgecolors=mesh_ec,
                                     linewidths=0.055))

    # Layer 2: ONLY fully damaged elements, plotted in black
    if idx.size:
        ax.add_collection(PolyCollection(verts[idx],
                                         facecolors=crack_fc,
                                         edgecolors=crack_fc,
                                         linewidths=0.16))

    ax.set_aspect('equal')
    ax.set_xlim(x_min - 5, x_max + 5)
    ax.set_ylim(y_min - 5, y_max + 5)
    ax.set_xlabel(r'$x$ [mm]', fontsize=11)
    ax.set_ylabel(r'$y$ [mm]', fontsize=11)
    ax.set_title('Fully damaged elements only', fontsize=13, fontweight='bold')

    # Remove the rectangular plot box. Keep only bottom/left ticks.
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(labelsize=10, direction='out', top=False, right=False)

    # Right-side clean legend/text axis. No frame, no patch box.
    axp = fh.add_axes([0.805, 0.17, 0.18, 0.74])
    axp.set_xlim(0, 1)
    axp.set_ylim(0, 1)
    axp.axis('off')

    axp.text(0.00, 0.98, 'Plot key', fontsize=11, fontweight='bold',
             ha='left', va='top')

    handles = [
        Line2D([0], [0], color=mesh_ec, lw=7, solid_capstyle='butt', label='FE mesh'),
        Line2D([0], [0], color=crack_fc, lw=7, solid_capstyle='butt',
               label=r'$\omega \ge %.2f$' % th)
    ]
    axp.legend(handles=handles,
               loc='upper left',
               bbox_to_anchor=(0.00, 0.86),
               frameon=False,
               fontsize=10.5,
               handlelength=2.0,
               handletextpad=0.8,
               borderaxespad=0.0,
               labelspacing=1.05)

    # No bbox here: this removes the box visible in the previous figures.
    axp.text(0.00, 0.47, label,
             fontsize=10.0, fontweight='bold', ha='left', va='top')
    axp.text(0.00, 0.35, 'Threshold: omega >= %.2f' % th,
             fontsize=9.5, fontweight='bold', ha='left', va='top')
    axp.text(0.00, 0.24, 'Fully damaged: %d elements' % int(idx.size),
             fontsize=9.5, fontweight='bold', ha='left', va='top')

    _save_hq(fh, base)
    plt.close(fh)

def _fig_load_cmod(np, plt, CMOD, F, res_dir):
    ip = int(np.argmax(F))
    pk_load, pk_cmod = F[ip], CMOD[ip]

    fh = plt.figure(figsize=(8.8 / 2.54, 7.0 / 2.54), facecolor='w')
    ax = fh.add_axes([0.14, 0.14, 0.82, 0.78])
    ax.grid(True, linestyle=':', color=(0.80, 0.80, 0.80))
    ax.tick_params(labelsize=9)
    for sp in ax.spines.values():
        sp.set_linewidth(0.7)

    col = (0.08, 0.30, 0.72)
    ax.fill_between(CMOD, F / 1000.0, 0, color=col, alpha=0.10, edgecolor='none')
    ax.plot(CMOD, F / 1000.0, '-', color=col, linewidth=2.0)
    ax.plot(pk_cmod, pk_load / 1000.0, marker='*', markersize=13,
            mfc=(0.98, 0.82, 0.0), mec=(0.40, 0.28, 0.0), linewidth=1.0)

    lbl_x = pk_cmod + 0.35 * 0.03
    lbl_y = pk_load / 1000.0 * 1.04
    ha = 'left'
    if lbl_x > 0.35 * 0.72:
        lbl_x, ha = pk_cmod - 0.35 * 0.03, 'right'
    ax.text(lbl_x, lbl_y,
            '$P_{\\rm peak} = %.2f$ kN\nCMOD $= %.4f$ mm' % (pk_load / 1000.0, pk_cmod),
            fontsize=8.5, fontweight='bold', color=(0.30, 0.20, 0.0), ha=ha, va='bottom',
            bbox=dict(facecolor='w', edgecolor=(0.60, 0.50, 0.20), linewidth=0.5, pad=3))

    ax.set_xlabel('CMOD [mm]', fontsize=10)
    ax.set_ylabel('Load [kN]', fontsize=10)
    ax.set_title('Load--CMOD response', fontsize=10, fontweight='bold')
    ax.set_xlim(0, 0.35)
    ax.set_ylim(0, pk_load / 1000.0 * 1.28)
    _save_hq(fh, os.path.join(res_dir, 'abaqus_load_cmod_fig'))
    plt.close(fh)


def plot_main(res_dir):
    try:
        np, plt, PolyCollection, LSC = _mpl()
    except Exception as e:
        print('ERROR: matplotlib/numpy not available: %s' % str(e))
        print('Set PYTHON_EXE to a Python that has them.')
        return

    lc = _read_csv(np, os.path.join(res_dir, 'abaqus_load_cmod.csv'), 2)
    if lc.shape[0] == 0:
        print('ERROR: abaqus_load_cmod.csv missing/empty.')
        return
    CMOD, F = lc[:, 0], lc[:, 1]
    pk_load = F.max()

    D = _load_plotdata(np, os.path.join(res_dir, 'plotdata'))
    cmap = crack_cmap(LSC)

    _fig_load_cmod(np, plt, CMOD, F, res_dir)
    _fig_damage(np, plt, PolyCollection, cmap, D['nodes'], D['elems'], D['w_peak'],
                'Peak load: %.2f kN' % (pk_load / 1000.0),
                os.path.join(res_dir, 'abaqus_fig_damage_peak'))

    pp_load = F[-1]
    for i in range(len(CMOD)):
        if CMOD[i] > 0.3:
            pp_load = F[i]
            break
    _fig_damage(np, plt, PolyCollection, cmap, D['nodes'], D['elems'], D['w_pp'],
                'Post-peak: %.2f kN' % (pp_load / 1000.0),
                os.path.join(res_dir, 'abaqus_fig_damage_postpeak'))

    _fig_damage(np, plt, PolyCollection, cmap, D['nodes'], D['elems'], D['w_last'],
                'Last step: %.2f kN' % (F[-1] / 1000.0),
                os.path.join(res_dir, 'abaqus_fig_damage_last_step'))

    print('Figures written:')
    for nm in ('abaqus_load_cmod_fig',
               'abaqus_fig_damage_peak',
               'abaqus_fig_damage_postpeak',
               'abaqus_fig_damage_last_step'):
        print('  %s.png / .pdf' % os.path.join(res_dir, nm))


# =====================================================================
# MAIN
# =====================================================================
def main():
    base = script_base_dir()
    case_dir = os.path.join(base, MODEL)
    ensure_dir(case_dir)
    print('=' * 60)
    print('ABAQUS 3PB  Oliver-T3 FAST build + run')
    print('Base: %s' % base)
    print('=' * 60)
    umat = find_umat_file(base)
    if umat is None:
        die('UMAT not found. Put cdm_umat_2d_OLIVER_T3.for in this folder: %s' % base)
    part = build_model(case_dir, os.path.abspath(umat))
    grad_path = os.path.join(case_dir, 'oliver_t3_gradN.dat')
    write_oliver_t3_gradN_from_part(part, grad_path)
    print('Using UMAT: %s' % os.path.abspath(umat))
    print('Oliver table must be in the job folder: %s' % grad_path)
    odb_path, solve_wall = run_job(case_dir, os.path.abspath(umat))
    extract_and_plot(odb_path, solve_wall)
    print('\nDONE. Results: %s' % os.path.join(case_dir, 'results'))


def plot_mode_main():
    res_dir = None
    args = sys.argv[1:]
    for i, a in enumerate(args):
        if a == '--plot' and i + 1 < len(args):
            res_dir = args[i + 1]
            break
    if res_dir is None:
        res_dir = os.path.join(os.getcwd(), MODEL, 'results')
    print('Plot mode. Results: %s' % res_dir)
    plot_main(res_dir)


if __name__ == '__main__':
    if PLOT_MODE:
        plot_mode_main()
    else:
        main()