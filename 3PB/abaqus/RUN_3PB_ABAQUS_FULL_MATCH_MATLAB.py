# -*- coding: utf-8 -*-
# =====================================================================
#  RUN_3PB_ABAQUS_FULL_MATCH_MATLAB.py
#
#  ONE-COMMAND full Abaqus pipeline that reproduces the MATLAB solver
#  result (solver_main_3pb_OLIVER_bandwidth.m) as closely as Abaqus can.
#
#  It does EVERYTHING in one run:
#    PHASE 1  Build the EXACT same Gregoire 3PB CPS3 mesh in Abaqus/CAE
#             (identical seeding/refinement to make_mesh_3pb_MATLAB_ONLY.py,
#              so MATLAB and Abaqus share the same mesh by construction).
#    PHASE 2  Write oliver_bandwidth_data.dat from that exact mesh.
#    PHASE 3  Write a UMAT-ready .inp (UMAT material, MATLAB step size,
#             output requests for Load-CMOD + damage field).
#    PHASE 4  Run the Abaqus/Standard job with the UMAT.
#    PHASE 5  Extract Load-CMOD from the ODB -> CSV.
#    PHASE 6  Extract damage geometry at peak + post-peak.
#    PHASE 7  Plot in EXACT MATLAB style (same colormap, same axes,
#             same shaded fill, same pentagram peak marker).
#
#  REQUIRED in the working folder:
#    scm_umat_2d_OLIVER_MATCH_MATLAB.for   (the fixed UMAT)
#
#  RUN:
#    abaqus cae noGUI=RUN_3PB_ABAQUS_FULL_MATCH_MATLAB.py
#
#  OUTPUTS  ->  abaqus_results/
#    abaqus_load_cmod.csv
#    abaqus_load_cmod.csv            (Load-CMOD raw data only)
#    abaqus_summary.txt
#    damage_peak.dat / damage_postpeak.dat  (raw geometry only)
#
#  NOTE:
#    Plotting is intentionally disabled. This avoids Abaqus/Python
#    exiting or crashing after the simulation when matplotlib/system
#    Python plotting is triggered.
# =====================================================================
from __future__ import print_function

import os
import sys
import csv
import math
import time

# =====================================================================
#  CONFIG  — must match MATLAB solver_main_3pb exactly
# =====================================================================
MODEL      = 'Gregoire_3PB'
JOB_NAME   = 'Gregoire_3PB_MATLAB_MATCH'
STEP_NAME  = 'Loading'
RESULT_DIR = 'abaqus_results'

# --- geometry (mm) — identical to make_mesh_3pb_MATLAB_ONLY.py --------
D        = 100.0
S        = 2.5 * D
overhang = 0.5 * D
L        = S + 2.0 * overhang
b_thick  = 50.0
a0       = 0.2 * D
wn       = D / 40.0
xC       = L / 2.0
xL_notch = xC - wn / 2.0
xR_notch = xC + wn / 2.0

ELEM_SIZE_GLOBAL = D / 16.0
ELEM_SIZE_REFINE = ELEM_SIZE_GLOBAL / 5.0
REFINE_W = 0.5 * D
REFINE_H = D

# --- material (MATLAB p.*) -------------------------------------------
E_C   = 37000.0
NU_C  = 0.20
FT    = 3.50
GF    = 0.090
FC_FT = 10.0          # k = fc/ft = 35/3.5

# --- loading / stepping (MATLAB p.max_disp, p.num_steps) -------------
U_FINAL  = -0.20      # mm total midspan displacement
NUM_STEP = 10000      # MATLAB load steps
INC_SIZE = 1.0 / NUM_STEP   # = 1e-4 pseudo-time per increment
MAX_INC  = 50000

UMAT_CANDIDATES = [
    'scm_umat_2d_OLIVER_MATCH_MATLAB.for',
    'scm_umat_2d_OLIVER_MATCH_MATLAB(1).for',
]

CWD      = os.getcwd()
CASE_DIR = os.path.join(CWD, MODEL)


# =====================================================================
#  PHASE 1 — build the EXACT MATLAB mesh in Abaqus/CAE
# =====================================================================
def phase_build_mesh():
    from abaqus import mdb
    from abaqusConstants import (
        TWO_D_PLANAR, DEFORMABLE_BODY, SIDE1,
        FINER, FREE, TRI, ADVANCING_FRONT,
        CARTESIAN, ON, OFF, CPS3, STANDARD, SET
    )
    import mesh as abqMesh

    if not os.path.exists(CASE_DIR):
        os.makedirs(CASE_DIR)

    if MODEL in mdb.models.keys():
        del mdb.models[MODEL]
    m = mdb.Model(name=MODEL)

    # --- sketch the notched beam outline ---
    sk = m.ConstrainedSketch(name='BeamSketch', sheetSize=2.0 * L)
    sk.Line(point1=(0.0, 0.0),      point2=(xL_notch, 0.0))
    sk.Line(point1=(xL_notch, 0.0), point2=(xL_notch, a0))
    sk.Line(point1=(xL_notch, a0),  point2=(xR_notch, a0))
    sk.Line(point1=(xR_notch, a0),  point2=(xR_notch, 0.0))
    sk.Line(point1=(xR_notch, 0.0), point2=(L, 0.0))
    sk.Line(point1=(L, 0.0),        point2=(L, D))
    sk.Line(point1=(L, D),          point2=(0.0, D))
    sk.Line(point1=(0.0, D),        point2=(0.0, 0.0))

    part = m.Part(name='Beam2D', dimensionality=TWO_D_PLANAR,
                  type=DEFORMABLE_BODY)
    part.BaseShell(sketch=sk)
    del sk

    # --- partition refinement zone around notch ---
    tform = part.MakeSketchTransform(sketchPlane=part.faces[0],
                                     sketchPlaneSide=SIDE1,
                                     origin=(0.0, 0.0, 0.0))
    sk_r = m.ConstrainedSketch(name='RefineZone', sheetSize=L, transform=tform)
    sk_r.rectangle(point1=(xC - REFINE_W / 2.0, 0.0),
                   point2=(xC + REFINE_W / 2.0, REFINE_H))
    part.PartitionFaceBySketch(sketch=sk_r, faces=part.faces)
    del sk_r

    # --- dummy section (UMAT props supplied in .inp) ---
    mat = m.Material(name='Concrete')
    mat.Elastic(table=((E_C, NU_C),))
    m.HomogeneousSolidSection(name='BeamSec', material='Concrete',
                              thickness=b_thick)
    part.SectionAssignment(region=(part.faces,), sectionName='BeamSec')

    # --- CPS3 triangle mesh, identical controls to MATLAB mesh script ---
    elem_t3 = abqMesh.ElemType(elemCode=CPS3, elemLibrary=STANDARD)
    part.setElementType(regions=(part.faces,), elemTypes=(elem_t3,))
    for f in part.faces:
        part.setMeshControls(regions=(f,), technique=FREE,
                             elemShape=TRI, algorithm=ADVANCING_FRONT)

    part.seedPart(size=ELEM_SIZE_GLOBAL, deviationFactor=0.1)

    refine_edges = part.edges.getByBoundingBox(
        xMin=xC - REFINE_W / 2.0 - 1e-3,
        xMax=xC + REFINE_W / 2.0 + 1e-3,
        yMin=-1e-3, yMax=REFINE_H + 1e-3)
    if len(refine_edges):
        part.seedEdgeBySize(edges=refine_edges,
                            size=ELEM_SIZE_REFINE, constraint=FINER)

    part.generateMesh()

    # --- node sets (same picking logic as MATLAB mesh script) ---
    def pick_single_nearest(p, xt, yt):
        best = None; bd = 1.0e30
        for nd in p.nodes:
            x = nd.coordinates[0]; y = nd.coordinates[1]
            d2 = (x - xt) ** 2 + (y - yt) ** 2
            if d2 < bd:
                bd = d2; best = nd
        return p.nodes.sequenceFromLabels(labels=(best.label,))

    def pick_n_nearest_top(p, xt, y_top, n_keep):
        tol = 1.0e-3 * abs(y_top) + 1.0e-6
        cands = []
        for nd in p.nodes:
            x = nd.coordinates[0]; y = nd.coordinates[1]
            if abs(y - y_top) <= tol:
                cands.append((abs(x - xt), nd.label))
        cands.sort()
        keep = [lb for _, lb in cands[:max(1, n_keep)]]
        return p.nodes.sequenceFromLabels(labels=tuple(keep))

    node_tol = 0.5 * ELEM_SIZE_REFINE

    part.Set(nodes=pick_single_nearest(part, overhang, 0.0),
             name='Support_Left')
    part.Set(nodes=pick_single_nearest(part, L - overhang, 0.0),
             name='Support_Right')
    part.Set(nodes=pick_n_nearest_top(part, xC, D, 3),
             name='Load_Nodes')
    part.Set(nodes=part.nodes.getByBoundingBox(
        xMin=xL_notch - node_tol, xMax=xL_notch + node_tol,
        yMin=-node_tol, yMax=node_tol), name='CMOD1')
    part.Set(nodes=part.nodes.getByBoundingBox(
        xMin=xR_notch - node_tol, xMax=xR_notch + node_tol,
        yMin=-node_tol, yMax=node_tol), name='CMOD2')

    # --- assembly + assembly-level sets (needed for history output) ---
    asm = m.rootAssembly
    asm.DatumCsysByDefault(CARTESIAN)
    inst = asm.Instance(name='BeamInst', part=part, dependent=ON)
    for nm in ('Support_Left', 'Support_Right', 'Load_Nodes', 'CMOD1', 'CMOD2'):
        asm.Set(nodes=inst.sets[nm].nodes, name=nm)

    # --- step + BCs (write base .inp, patched later) ---
    m.StaticStep(name=STEP_NAME, previous='Initial',
                 maxNumInc=MAX_INC, initialInc=INC_SIZE,
                 minInc=1e-8, maxInc=INC_SIZE, nlgeom=OFF)

    m.DisplacementBC(name='BC_SupL', createStepName=STEP_NAME,
                     region=asm.sets['Support_Left'], u1=SET, u2=SET)
    m.DisplacementBC(name='BC_SupR', createStepName=STEP_NAME,
                     region=asm.sets['Support_Right'], u2=SET)
    m.DisplacementBC(name='BC_Load', createStepName=STEP_NAME,
                     region=asm.sets['Load_Nodes'], u2=U_FINAL)

    if MODEL in mdb.jobs.keys():
        del mdb.jobs[MODEL]
    job = mdb.Job(name=MODEL, model=MODEL,
                  description='Gregoire 3PB MATLAB-match base mesh')

    cwd0 = os.getcwd()
    os.chdir(CASE_DIR)
    try:
        job.writeInput(consistencyChecking=OFF)
    finally:
        os.chdir(cwd0)

    n_el = len(part.elements); n_nd = len(part.nodes)
    print('  mesh: %d CPS3 elements, %d nodes' % (n_el, n_nd))
    print('  base inp: %s' % os.path.join(CASE_DIR, MODEL + '.inp'))
    return os.path.join(CASE_DIR, MODEL + '.inp')


# =====================================================================
#  PHASE 2 — Oliver bandwidth table from exact mesh
# =====================================================================
def read_inp_nodes_elems(inp_path):
    def is_kw(s): return s.strip().startswith('*')
    def nums(s):
        out = []
        for x in s.replace(',', ' ').split():
            try: out.append(float(x))
            except Exception: pass
        return out
    with open(inp_path) as f:
        lines = f.readlines()
    nodes = []; elems = []
    i = 0
    while i < len(lines):
        low = lines[i].strip().lower()
        if low.startswith('*node') and not low.startswith('*node output'):
            i += 1
            while i < len(lines) and not is_kw(lines[i]):
                v = nums(lines[i])
                if len(v) >= 3: nodes.append((int(v[0]), v[1], v[2]))
                i += 1
            continue
        if low.startswith('*element') and not low.startswith('*element output'):
            i += 1
            while i < len(lines) and not is_kw(lines[i]):
                v = nums(lines[i])
                if len(v) >= 4:
                    elems.append((int(v[0]), int(v[1]), int(v[2]), int(v[3])))
                i += 1
            continue
        i += 1
    return nodes, elems


def write_oliver_bandwidth(inp_path, out_path):
    nodes, elems = read_inp_nodes_elems(inp_path)
    if not nodes or not elems:
        raise RuntimeError('Could not read mesh from %s' % inp_path)
    xy = dict((nid, (x, y)) for nid, x, y in nodes)
    with open(out_path, 'w') as f:
        for eid, n1, n2, n3 in elems:
            x1, y1 = xy[n1]; x2, y2 = xy[n2]; x3, y3 = xy[n3]
            area_s = 0.5 * ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))
            area = abs(area_s)
            if area <= 1e-14:
                raise RuntimeError('zero CPS3 area, elem %d' % eid)
            b1 = y2 - y3; b2 = y3 - y1; b3 = y1 - y2
            c1 = x3 - x2; c2 = x1 - x3; c3 = x2 - x1
            inv2A = 1.0 / (2.0 * area_s)
            g1x = b1 * inv2A; g1y = c1 * inv2A
            g2x = b2 * inv2A; g2y = c2 * inv2A
            g3x = b3 * inv2A; g3y = c3 * inv2A
            f.write('%d %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n' %
                    (eid, g1x, g1y, g2x, g2y, g3x, g3y, area, math.sqrt(area)))
    print('  oliver table: %s  (%d elems)' % (out_path, len(elems)))
    return len(elems)


# =====================================================================
#  PHASE 3 — patch .inp: UMAT material + step + outputs
# =====================================================================
def _skip_block(lines, i, n):
    """Skip a keyword line and all its data lines (until next * keyword)."""
    i += 1
    while i < n and not lines[i].lstrip().startswith('*'):
        i += 1
    return i


def make_run_ready_inp(src_inp, out_inp):
    """Patch the CAE-written base .inp into a UMAT run-ready .inp.

    The CAE base .inp ALREADY contains its own *Restart and *Output
    blocks (from the StaticStep defaults). We must STRIP those first,
    otherwise Abaqus aborts with:
       *RESTART MAY BE SPECIFIED ONLY ONCE PER STEP
    Then we inject exactly one set of MATLAB-match output requests.
    """
    with open(src_inp) as f:
        lines = f.readlines()

    out = []
    i = 0
    n = len(lines)
    while i < n:
        ln = lines[i]
        low = ln.strip().lower()

        # --- fix Solid Section material name -> Concrete ---------------
        if low.startswith('*solid section') or low.startswith('*shell section'):
            import re as _re
            fixed = _re.sub(r'material\s*=\s*[^,\s]+',
                            'material=Concrete', ln, flags=_re.I)
            out.append(fixed)
            i += 1
            continue

        # --- bump *Step increment cap to MAX_INC ----------------------
        if low.startswith('*step'):
            import re as _re
            if _re.search(r'inc\s*=', low):
                fixed = _re.sub(r'inc\s*=\s*\d+', 'inc=%d' % MAX_INC,
                                ln, flags=_re.I)
            else:
                fixed = ln.rstrip('\n') + ', inc=%d\n' % MAX_INC
            out.append(fixed)
            i += 1
            continue

        # --- replace Concrete material block with UMAT material ---------
        if low.startswith('*material') and 'concrete' in low:
            out.append('*Material, name=Concrete\n')
            out.append('*Depvar\n')
            out.append('  3,\n')
            out.append('** 1=kappa  2=omega(damage)  3=H(Oliver bandwidth)\n')
            out.append('*User Material, constants=5\n')
            out.append('** E, nu, ft, GF, fc/ft\n')
            out.append('%.1f, %.2f, %.2f, %.4f, %.1f\n' %
                       (E_C, NU_C, FT, GF, FC_FT))
            i += 1
            # skip old material suboptions until next major block
            while i < n:
                li = lines[i].strip().lower()
                if (li.startswith('*step') or li.startswith('*material')
                        or li.startswith('*boundary')
                        or li.startswith('*assembly')
                        or li.startswith('*end ')):
                    break
                if li.startswith('*') and not li.startswith('**'):
                    i = _skip_block(lines, i, n)
                    continue
                i += 1
            continue

        # --- force MATLAB step size on *Static line --------------------
        if low.startswith('*static'):
            # MATLAB-match stepping: DIRECT=NO STOP gives FIXED increments
            # (no adaptive cut-back), exactly like MATLAB's 10000 equal
            # displacement steps. NO STOP => proceed even if an increment
            # does not converge (mirrors MATLAB's fixed-iteration loop).
            out.append('*Static, direct=NO STOP\n')
            i = _skip_block(lines, i, n)
            out.append('%.1e, 1.0\n' % INC_SIZE)
            # Tighten convergence toward MATLAB's tol (1e-6). With the
            # staggered secant tangent each increment is linear, so this
            # is satisfied in ~1 iteration.
            out.append('*Controls, parameters=field, field=displacement\n')
            out.append('1e-6,\n')
            continue

        # --- STRIP any existing *Restart line (CAE writes one) ---------
        if low.startswith('*restart'):
            i = _skip_block(lines, i, n)
            continue

        # --- STRIP any existing *Output block (field or history) -------
        if low.startswith('*output'):
            i = _skip_block(lines, i, n)
            # also skip the *Node Output / *Element Output suboptions
            # that belong to this output block
            while i < n:
                li = lines[i].strip().lower()
                if (li.startswith('*node output')
                        or li.startswith('*element output')
                        or li.startswith('*contact output')
                        or li.startswith('*energy output')):
                    i = _skip_block(lines, i, n)
                else:
                    break
            continue

        # --- also strip standalone *Node Output / *Element Output ------
        # (in case they appear outside a recognized *Output block)
        if (low.startswith('*node output')
                or low.startswith('*element output')
                or low.startswith('*energy output')
                or low.startswith('*contact output')):
            i = _skip_block(lines, i, n)
            continue

        # --- inject our outputs once, right before *End Step -----------
        if low.startswith('*end step'):
            out.append('** === MATLAB-MATCH OUTPUTS (single set) ===\n')
            out.append('*Restart, write, frequency=0\n')
            out.append('*Output, field, number interval=200\n')
            out.append('*Node Output\n')
            out.append('RF, U\n')
            out.append('*Element Output, directions=YES\n')
            out.append('S, E, SDV\n')
            out.append('*Output, history, frequency=1\n')
            out.append('*Node Output, nset=CMOD1\n')
            out.append('U1,\n')
            out.append('*Node Output, nset=CMOD2\n')
            out.append('U1,\n')
            out.append('*Node Output, nset=Load_Nodes\n')
            out.append('RF2, U2\n')
            out.append('*End Step\n')
            i += 1
            continue

        out.append(ln)
        i += 1

    with open(out_inp, 'w') as f:
        f.writelines(out)
    print('  run-ready inp: %s' % out_inp)
    return out_inp


# =====================================================================
#  PHASE 4 — run Abaqus job with UMAT
# =====================================================================
def choose_existing(cands, what):
    for c in cands:
        if os.path.exists(c):
            return c
    raise RuntimeError('Missing %s. Tried: %s' % (what, ', '.join(cands)))


def cleanup_job(job):
    for ext in ['.odb', '.lck', '.dat', '.msg', '.sta', '.sim', '.prt',
                '.com', '.log', '.res', '.mdl', '.stt', '.pac', '.sel', '.abq']:
        p = job + ext
        if os.path.exists(p):
            try: os.remove(p)
            except Exception: pass


def run_job(inp, umat, job):
    cleanup_job(job)
    print('  job=%s  inp=%s  umat=%s' % (job, inp, umat))
    from abaqus import mdb
    from abaqusConstants import OFF
    if job in mdb.jobs.keys():
        del mdb.jobs[job]
    mdb.JobFromInputFile(name=job, inputFileName=inp, userSubroutine=umat)
    mdb.jobs[job].submit(consistencyChecking=OFF)
    mdb.jobs[job].waitForCompletion()
    odb = job + '.odb'
    if not os.path.exists(odb):
        raise RuntimeError('ODB not found after run: ' + odb)
    print('  job done. odb=%s' % odb)
    return odb


# =====================================================================
#  PHASE 5 — extract Load-CMOD from ODB
# =====================================================================
def find_node_set(odb, wanted):
    sets = odb.rootAssembly.nodeSets
    keys = list(sets.keys())
    uw = wanted.upper()
    for k in keys:
        if k.upper() == uw:
            return sets[k]
    for k in keys:
        ku = k.upper()
        if ku.endswith('.' + uw) or ku.endswith('-' + uw):
            return sets[k]
    for k in keys:
        if uw in k.upper():
            return sets[k]
    raise RuntimeError('node set %s not found. ODB sets: %s' % (wanted, str(keys)))


def flatten_labels(nset):
    out = []
    try:
        for arr in nset.nodes:
            for nd in arr:
                out.append(int(nd.label))
    except Exception:
        try:
            for nd in nset.nodes:
                out.append(int(nd.label))
        except Exception:
            pass
    return sorted(set(out))


def hist_series(step, node_label, varname):
    """Return [(time, value), ...] for one node + variable, or None.

    Crash-proof: skips outputs whose .data is None/empty, and requires
    an EXACT history-output name match (U1, U2, RF2) so we never grab an
    unrelated empty output (the previous cause of the NoneType crash).
    """
    label = str(int(node_label))
    vu = varname.upper().replace(' ', '')
    for rname, region in step.historyRegions.items():
        rn = rname.upper().replace(' ', '')
        # Region must reference this node label
        if not (('.' + label) in rn or rn.endswith('-' + label)
                or rn.endswith(label)):
            continue
        for hname, hobj in region.historyOutputs.items():
            hn = hname.upper().replace(' ', '')
            # EXACT match only (not substring) to avoid empty look-alikes
            if hn != vu:
                continue
            d = getattr(hobj, 'data', None)
            if d is None:
                continue
            try:
                lst = list(d)
            except Exception:
                continue
            if len(lst) == 0:
                continue
            return lst
    return None


def avg_hist(step, labels, varname):
    series = []
    for lb in labels:
        s = hist_series(step, lb, varname)
        if s:
            series.append(s)
    if not series:
        return None
    n = min(len(s) for s in series)
    if n == 0:
        return None
    return [(series[0][i][0], sum(s[i][1] for s in series) / float(len(series)))
            for i in range(n)]


def sum_hist(step, labels, varname):
    series = []
    for lb in labels:
        s = hist_series(step, lb, varname)
        if s:
            series.append(s)
    if not series:
        return None
    n = min(len(s) for s in series)
    if n == 0:
        return None
    return [(series[0][i][0], sum(s[i][1] for s in series))
            for i in range(n)]


def _instance_coord_map(odb):
    """Build {node_label: (x, y)} from the largest instance in the ODB."""
    inst = biggest_instance(odb)
    coord = {}
    for nd in inst.nodes:
        coord[int(nd.label)] = (float(nd.coordinates[0]),
                                float(nd.coordinates[1]))
    return coord


def _nearest_label(coord, xt, yt):
    best = -1; bd = 1.0e30
    for lb, (x, y) in coord.items():
        d = (x - xt) ** 2 + (y - yt) ** 2
        if d < bd:
            bd = d; best = lb
    return best


def _top_labels(coord, xt, yt, n_keep):
    cands = [(abs(x - xt), lb) for lb, (x, y) in coord.items()
             if abs(y - yt) < 0.5]
    cands.sort()
    return [lb for _, lb in cands[:max(1, n_keep)]]


def extract_load_cmod(odb, csv_path):
    """Bulletproof Load-CMOD extraction.

    Strategy (in order):
      1. COORDINATE-BASED field extraction — finds the notch-mouth nodes
         and load nodes by geometry, reads U + RF from field output.
         Immune to node-set naming and to empty history data.
      2. History extraction (finer resolution) if coordinate method gives
         too few points and history is healthy.

    Heavy diagnostics are printed so the ODB contents are always visible
    in the Abaqus message/log if anything is off.
    """
    step = odb.steps[STEP_NAME]

    # ---- diagnostics --------------------------------------------------
    print('  --- ODB diagnostics ---')
    print('  steps        : %s' % str(list(odb.steps.keys())))
    print('  frames       : %d' % len(step.frames))
    print('  hist regions : %d' % len(step.historyRegions))
    try:
        nset_keys = list(odb.rootAssembly.nodeSets.keys())
        print('  node sets    : %s' % str(nset_keys))
    except Exception as e:
        print('  node sets    : <error %s>' % str(e))
    if len(step.frames) > 0:
        fo0 = step.frames[-1].fieldOutputs
        print('  field outputs in last frame: %s' % str(list(fo0.keys())))
    print('  -----------------------')

    if len(step.frames) < 2:
        print('  WARNING: ODB has < 2 frames.')
        print('  This means the Abaqus solve diverged almost immediately.')
        print('  It is a SOLVER problem (not extraction). Check .msg/.sta.')

    # ---- find target nodes by COORDINATE (geometry is known) ----------
    coord = _instance_coord_map(odb)
    print('  instance nodes mapped: %d' % len(coord))

    c1 = _nearest_label(coord, xL_notch, 0.0)   # notch left lip  (173.75, 0)
    c2 = _nearest_label(coord, xR_notch, 0.0)   # notch right lip (176.25, 0)
    ld = _top_labels(coord, xC, D, 3)           # load nodes (175, 100)

    print('  CMOD1 node %d @ %s (target %.3f,0)' %
          (c1, str(coord.get(c1)), xL_notch))
    print('  CMOD2 node %d @ %s (target %.3f,0)' %
          (c2, str(coord.get(c2)), xR_notch))
    print('  Load nodes %s' % str(ld))

    ld_set = set(ld)

    # ---- PRIMARY: coordinate-based field extraction -------------------
    data = []
    for fr in step.frames:
        fo = fr.fieldOutputs
        if 'U' not in fo or 'RF' not in fo:
            continue
        u_c1 = None; u_c2 = None; rf_sum = 0.0; got_rf = False
        try:
            for v in fo['U'].values:
                nl = int(v.nodeLabel)
                if nl == c1: u_c1 = float(v.data[0])
                elif nl == c2: u_c2 = float(v.data[0])
            for v in fo['RF'].values:
                if int(v.nodeLabel) in ld_set:
                    rf_sum += float(v.data[1]); got_rf = True
        except Exception:
            continue
        if (u_c1 is not None) and (u_c2 is not None) and got_rf:
            cmod = u_c2 - u_c1               # CMOD = U1_right - U1_left
            load = -rf_sum                    # reaction -> applied load
            data.append((float(fr.frameValue), cmod, load))

    if data:
        print('  coordinate/field extraction OK: %d points' % len(data))
    else:
        print('  coordinate/field extraction gave 0 points')

    # ---- prefer HISTORY (frequency=1 => one point per increment) ------
    # MATLAB writes all load steps (~10000). History output matches that
    # density, while field output is only `number interval` frames (200),
    # which under-samples the peak. Use history whenever it is denser.
    if True:
        print('  trying history extraction (per-increment resolution)...')
        try:
            u1l = avg_hist(step, [c1], 'U1')
            u1r = avg_hist(step, [c2], 'U1')
            rf2 = sum_hist(step, ld, 'RF2')
            if u1l and u1r and rf2:
                n = min(len(u1l), len(u1r), len(rf2))
                hist_data = []
                for i in range(n):
                    t = u1l[i][0]
                    cmod = u1r[i][1] - u1l[i][1]
                    load = -rf2[i][1]
                    hist_data.append((t, cmod, load))
                if len(hist_data) > len(data):
                    print('  history gave %d points (using history)'
                          % len(hist_data))
                    data = hist_data
        except Exception as e:
            print('  history extraction skipped: %s' % str(e))

    # ---- write csv ----------------------------------------------------
    with open(csv_path, 'w') as f:
        wcsv = csv.writer(f)
        wcsv.writerow(['time', 'cmod_mm', 'load_N'])
        for t, c, p in data:
            wcsv.writerow(['%.12e' % t, '%.12e' % c, '%.12e' % p])
    print('  wrote %s (%d rows)' % (csv_path, len(data)))
    return data


# =====================================================================
#  PHASE 6 — damage geometry at peak + post-peak
# =====================================================================
def biggest_instance(odb):
    best = None; bn = -1
    for nm, inst in odb.rootAssembly.instances.items():
        try:
            k = len(inst.elements)
            if k > bn: bn = k; best = inst
        except Exception:
            pass
    return best


def frame_at_time(step, target):
    bi = 0; bd = 1e100
    for i, fr in enumerate(step.frames):
        d = abs(float(fr.frameValue) - float(target))
        if d < bd: bd = d; bi = i
    return step.frames[bi], bi


def dump_damage_geometry(odb, frame, out_dat):
    inst = biggest_instance(odb)
    coords = {}
    for nd in inst.nodes:
        coords[int(nd.label)] = (float(nd.coordinates[0]),
                                 float(nd.coordinates[1]))
    disp = {}
    if 'U' in frame.fieldOutputs:
        for v in frame.fieldOutputs['U'].values:
            disp[int(v.nodeLabel)] = (float(v.data[0]), float(v.data[1]))
    dmg = {}
    if 'SDV' in frame.fieldOutputs:
        for v in frame.fieldOutputs['SDV'].values:
            try:
                dat = v.data
                if hasattr(dat, '__len__'):
                    om = float(dat[1]) if len(dat) >= 2 else float(dat[0])
                else:
                    om = float(dat)
                eid = int(v.elementLabel)
                dmg[eid] = max(dmg.get(eid, 0.0), om)
            except Exception:
                pass
    elif 'SDV2' in frame.fieldOutputs:
        for v in frame.fieldOutputs['SDV2'].values:
            dmg[int(v.elementLabel)] = float(v.data)

    with open(out_dat, 'w') as f:
        # header: deformed nodal coords; then elements with omega
        f.write('# NODES: id x y ux uy\n')
        for lb in sorted(coords.keys()):
            x, y = coords[lb]
            ux, uy = disp.get(lb, (0.0, 0.0))
            f.write('N %d %.10e %.10e %.10e %.10e\n' % (lb, x, y, ux, uy))
        f.write('# ELEMS: id n1 n2 n3 omega\n')
        for el in inst.elements:
            con = el.connectivity
            if len(con) < 3: continue
            eid = int(el.label)
            f.write('E %d %d %d %d %.10e\n' %
                    (eid, int(con[0]), int(con[1]), int(con[2]),
                     dmg.get(eid, 0.0)))
    print('  damage geometry: %s' % out_dat)


# =====================================================================
#  DRIVER
# =====================================================================
def main():
    print('=' * 70)
    print('Gregoire 3PB — full Abaqus pipeline matching MATLAB solver')
    print('=' * 70)

    umat = choose_existing(UMAT_CANDIDATES, 'UMAT .for file')

    if not os.path.isdir(RESULT_DIR):
        os.makedirs(RESULT_DIR)

    print('\n[1/7] Getting mesh...')
    # CRITICAL for matching MATLAB: reuse the SAME .inp that MATLAB read.
    # MATLAB reads its mesh from Gregoire_3PB/Gregoire_3PB.inp (written by
    # make_mesh_3pb_MATLAB_ONLY.py). If that file exists, REUSE it so both
    # codes run the identical mesh. Only build a fresh mesh if none exists.
    existing_inp = os.path.join(CASE_DIR, MODEL + '.inp')
    if os.path.exists(existing_inp):
        nds, els = read_inp_nodes_elems(existing_inp)
        print('  REUSING existing mesh (same as MATLAB):')
        print('    %s' % existing_inp)
        print('    %d nodes, %d CPS3 elements' % (len(nds), len(els)))
        print('  (delete this file to force a fresh Abaqus/CAE rebuild)')
        base_inp = existing_inp
    else:
        print('  No existing .inp found -> building fresh mesh in CAE')
        print('  WARNING: a fresh build may differ from MATLAB mesh if')
        print('  MATLAB used a different mesh-script version.')
        print('  Best practice: run make_mesh_3pb_MATLAB_ONLY.py FIRST,')
        print('  then run MATLAB and this script on that same mesh.')
        base_inp = phase_build_mesh()

    print('\n[2/7] Writing Oliver bandwidth table...')
    # UMAT reads this from the JOB working dir (cwd), so write it here.
    write_oliver_bandwidth(base_inp, 'oliver_bandwidth_data.dat')

    print('\n[3/7] Patching input for UMAT + outputs...')
    run_inp = JOB_NAME + '_run.inp'
    make_run_ready_inp(base_inp, run_inp)

    print('\n[4/7] Running Abaqus/Standard with UMAT...')
    odb_path = run_job(run_inp, umat, JOB_NAME)

    print('\n[5/7] Extracting Load-CMOD...')
    from odbAccess import openOdb
    odb = openOdb(odb_path, readOnly=True)
    csv_path = os.path.join(RESULT_DIR, 'abaqus_load_cmod.csv')
    data = extract_load_cmod(odb, csv_path)
    if not data:
        odb.close()
        raise RuntimeError('No Load-CMOD data extracted.')
    save_summary(data, os.path.join(RESULT_DIR, 'abaqus_summary.txt'))

    print('\n[6/7] Extracting damage geometry...')
    step = odb.steps[STEP_NAME]
    ip = max(range(len(data)), key=lambda i: data[i][2])
    peak_t = data[ip][0]
    peak_fr, _ = frame_at_time(step, peak_t)
    # post-peak: first frame where cmod > 0.3 (MATLAB snap_pp)
    pp_t = data[-1][0]
    for t, c, p in data:
        if c > 0.3:
            pp_t = t; break
    pp_fr, _ = frame_at_time(step, pp_t)

    dump_damage_geometry(odb, peak_fr,
                         os.path.join(RESULT_DIR, 'damage_peak.dat'))
    dump_damage_geometry(odb, pp_fr,
                         os.path.join(RESULT_DIR, 'damage_postpeak.dat'))
    odb.close()

    print('\n[7/7] Plotting skipped.')
    print('  Load-CMOD PNG plotting is disabled to prevent Abaqus/Python exit/crash.')
    print('  The simulation, Load-CMOD CSV, summary, and raw damage geometry are saved.')

    print('\nDONE. Results in: %s' % os.path.abspath(RESULT_DIR))
    print('  abaqus_load_cmod.csv')
    print('  damage_peak.dat / damage_postpeak.dat')
    print('  abaqus_summary.txt')


if __name__ == '__main__':
    main()
