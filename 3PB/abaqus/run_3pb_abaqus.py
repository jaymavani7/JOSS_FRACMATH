# -*- coding: utf-8 -*-
# =====================================================================
# run_3pb_abaqus.py  –  SINGLE-FILE Abaqus driver for notched 3PB test
#
# USAGE (one command from shell):
#     python3 run_3pb_abaqus.py
#
# OUTPUTS  →  Gregoire_3PB/results/
#     abaqus_load_cmod.csv       CMOD[mm], Load[N] history
#     abaqus_load_cmod.png       Load vs CMOD plot
#
# DISK/MEMORY SAVING VERSION:
#   - No field output is written to the ODB.
#   - No SDV/damage output is extracted.
#   - No timing JSON/TXT or PDF plots are saved.
#   - Abaqus heavy files can be deleted after Load-CMOD extraction.
#
# JOB-SUBMISSION FIXES (this version):
#   FIX-J1  UMAT copied to CASE_DIR, passed as bare filename (relative).
#           Absolute paths with backslashes or spaces are silently mangled
#           by the Abaqus launcher on Windows.
#   FIX-J2  'double=both' added.  The UMAT uses REAL*8 throughout; without
#           this flag Abaqus runs in single precision and the UMAT argument
#           types mismatch, causing an immediate job abort.
#   FIX-J3  'input=<name>.inp' made explicit so Abaqus finds the deck even
#           when cwd is not the directory that owns the job files.
#   FIX-J4  Pre-submission sanity routine checks: .inp exists, UMAT copy is
#           present, and the Fortran compiler is reachable via 'abaqus make'.
#   FIX-J5  subprocess uses shell=True on Windows via _abaqus_call() because
#           the 'abaqus' entry on Windows is a .bat wrapper that needs the
#           Windows shell to execute.
#   FIX-J6  stdout/stderr are NOT suppressed, so Abaqus error lines print
#           immediately (no more silent failures).
#
# PREVIOUS FIXES (carried forward):
#   FIX-P1  SCRIPT_DIR: try __file__, fall back to sys.argv scan so
#           Abaqus CAE noGUI= mode (where __file__ is undefined) works.
#   FIX-P2  UMAT filename: scm_umat_2d.for
#   FIX-P3  Matplotlib runs in plain Python (Phase D), not abaqus python.
#   FIX-P4  Arg detection uses exact match on sys.argv[-1].
# =====================================================================

import os
import sys
import shutil
import subprocess
import json
import time as _time

# =====================================================================
# SCRIPT_DIR  (FIX-P1)
# =====================================================================
try:
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
except NameError:
    # Abaqus CAE noGUI= does not define __file__; find the .py in argv.
    _candidates = [a for a in sys.argv if a.endswith('.py')]
    if _candidates:
        SCRIPT_DIR = os.path.dirname(os.path.abspath(_candidates[0]))
    else:
        SCRIPT_DIR = os.getcwd()

MODEL      = 'Gregoire_3PB'
CASE_DIR   = os.path.join(SCRIPT_DIR, MODEL)
RES_DIR    = os.path.join(CASE_DIR, 'results')
UMAT_SRC   = os.path.join(SCRIPT_DIR, 'scm_umat_2d.for')   # FIX-P2
UMAT_LOCAL = os.path.join(CASE_DIR,   'scm_umat_2d.for')   # FIX-J1 destination
TIMING_JSON = os.path.join(RES_DIR, 'abaqus_timing.json')  # not written in this lite version

# Save disk space after extraction. Keep False if you need to inspect the ODB later.
DELETE_HEAVY_ABAQUS_FILES = True

# =====================================================================
# GEOMETRY  [mm]
# =====================================================================
D         = 100.0
S         = 2.5 * D
overhang  = 0.5 * D
L         = S + 2.0 * overhang
b_thick   = 50.0
a0        = 0.2 * D
wn        = D / 40.0
xC        = L / 2.0
xL_notch  = xC - wn / 2.0
xR_notch  = xC + wn / 2.0

ELEM_SIZE_GLOBAL = D / 10.0
ELEM_SIZE_REFINE = ELEM_SIZE_GLOBAL / 5.0
REFINE_W = 0.5 * D
REFINE_H = D

# =====================================================================
# MATERIAL
# =====================================================================
E_C, NU_C, FT, GF, FCFT = 37000.0, 0.20, 3.5, 0.090, 10.0

# =====================================================================
# LOADING
# =====================================================================
U_FINAL = -0.2

# User-requested target loading resolution.
# maxInc=1/N_INC means one normal converged increment = 1/10000 of the step.
N_INC   = 10000

# Abaqus may cut back increments during nonlinear damage/softening.
# Therefore maxNumInc must be larger than N_INC, otherwise the job can stop with:
# "Too many increments needed to complete the step".
MAX_NUM_INC = 5 * N_INC

# =====================================================================
# PLATFORM HELPER  (FIX-J5)
# =====================================================================
_ON_WINDOWS = sys.platform.startswith('win')

def _abaqus_call(args, cwd=None):
    """
    Run an Abaqus command robustly on Windows and Linux/macOS.

    On Windows 'abaqus' is a .bat wrapper; subprocess needs shell=True
    to find and invoke it.  On Linux the script is a real executable so
    shell=False (safer) is used.
    stdout/stderr are inherited (FIX-J6) — Abaqus messages appear live
    in the terminal instead of being swallowed silently.
    """
    if _ON_WINDOWS:
        cmd_str = ' '.join(args)
        return subprocess.call(cmd_str, cwd=cwd, shell=True)
    else:
        return subprocess.call(args, cwd=cwd)

# =====================================================================
# TIMING HELPERS
# =====================================================================
def _load_timing():
    # Lite version: do not write/read timing files.
    return {}

def _save_timing(d):
    # Lite version: intentionally do nothing.
    return

# =====================================================================
# PRE-SUBMISSION SANITY CHECK  (FIX-J4)
# =====================================================================
def _check_before_solve():
    """
    Returns a list of error strings.  Empty list = all clear.
    Checks: .inp written, UMAT copy present, Fortran compiler reachable.
    """
    errors = []

    inp_path = os.path.join(CASE_DIR, MODEL + '.inp')
    if not os.path.exists(inp_path):
        errors.append('  [MISSING .inp]  ' + inp_path)
        errors.append('                  Phase A (build) did not complete.')

    if not os.path.exists(UMAT_LOCAL):
        errors.append('  [MISSING UMAT]  ' + UMAT_LOCAL)
        errors.append('                  UMAT was not copied to CASE_DIR.')

    # Quick compiler probe — compile a trivial subroutine.
    probe_src = os.path.join(CASE_DIR, '_probe.for')
    try:
        with open(probe_src, 'w') as fh:
            fh.write('      SUBROUTINE PROBE\n      RETURN\n      END\n')
        if _ON_WINDOWS:
            rc = subprocess.call(
                'abaqus make library=' + probe_src,
                cwd=CASE_DIR, shell=True,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        else:
            rc = subprocess.call(
                ['abaqus', 'make', 'library=' + probe_src],
                cwd=CASE_DIR,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if rc != 0:
            errors.append('  [NO COMPILER]  abaqus make probe returned rc=%d' % rc)
            errors.append('  Fix on Windows : install Intel oneAPI Fortran and')
            errors.append('                   launch from the Abaqus CAE shortcut.')
            errors.append('  Fix on Linux   : set compile_fortran in abaqus_v6.env')
            errors.append('                   (usually ifort or gfortran).')
    except Exception as ex:
        errors.append('  [COMPILER PROBE FAILED]  ' + str(ex))
    finally:
        for p in (probe_src,
                  os.path.join(CASE_DIR, '_probe.o'),
                  os.path.join(CASE_DIR, '_probe.obj'),
                  os.path.join(CASE_DIR, 'StandardU.lib'),
                  os.path.join(CASE_DIR, 'StandardU.dll')):
            try:
                if os.path.exists(p):
                    os.remove(p)
            except Exception:
                pass

    return errors

# =====================================================================
# PHASE A  — build model + write .inp        (abaqus cae noGUI=…)
# =====================================================================
def phase_build():
    from abaqus import mdb
    from abaqusConstants import (TWO_D_PLANAR, DEFORMABLE_BODY, SIDE1,
                                  FINER, FREE, TRI, ADVANCING_FRONT,
                                  CARTESIAN, ON, OFF,
                                  CPS3, STANDARD, SET)
    import mesh as abqMesh

    t0 = _time.time()

    for d in (CASE_DIR, RES_DIR):
        if not os.path.isdir(d):
            os.makedirs(d)

    if MODEL in mdb.models.keys():
        del mdb.models[MODEL]
    m = mdb.Model(name=MODEL)

    # ---- sketch ----
    sk = m.ConstrainedSketch(name='BeamSketch', sheetSize=2.0 * L)
    for p1, p2 in [
        ((0.0, 0.0),        (xL_notch, 0.0)),
        ((xL_notch, 0.0),   (xL_notch, a0)),
        ((xL_notch, a0),    (xR_notch, a0)),
        ((xR_notch, a0),    (xR_notch, 0.0)),
        ((xR_notch, 0.0),   (L, 0.0)),
        ((L, 0.0),          (L, D)),
        ((L, D),            (0.0, D)),
        ((0.0, D),          (0.0, 0.0)),
    ]:
        sk.Line(point1=p1, point2=p2)

    part = m.Part(name='Beam2D', dimensionality=TWO_D_PLANAR,
                   type=DEFORMABLE_BODY)
    part.BaseShell(sketch=sk)
    del sk

    # ---- refine partition ----
    tform = part.MakeSketchTransform(sketchPlane=part.faces[0],
                                      sketchPlaneSide=SIDE1,
                                      origin=(0.0, 0.0, 0.0))
    sk_r = m.ConstrainedSketch(name='RefineZone', sheetSize=L, transform=tform)
    sk_r.rectangle(point1=(xC - REFINE_W / 2.0, 0.0),
                   point2=(xC + REFINE_W / 2.0, REFINE_H))
    part.PartitionFaceBySketch(sketch=sk_r, faces=part.faces)
    del sk_r

    # ---- material + section ----
    mat = m.Material(name='Concrete')
    mat.UserMaterial(mechanicalConstants=(E_C, NU_C, FT, GF, FCFT))
    mat.Depvar(n=2)
    m.HomogeneousSolidSection(name='BeamSec', material='Concrete',
                              thickness=b_thick)
    part.SectionAssignment(region=(part.faces,), sectionName='BeamSec')

    # ---- mesh ----
    elem_t3 = abqMesh.ElemType(elemCode=CPS3, elemLibrary=STANDARD)
    part.setElementType(regions=(part.faces,), elemTypes=(elem_t3,))
    for f in part.faces:
        part.setMeshControls(regions=(f,), technique=FREE,
                              elemShape=TRI, algorithm=ADVANCING_FRONT)
    part.seedPart(size=ELEM_SIZE_GLOBAL, deviationFactor=0.1)
    refine_edges = part.edges.getByBoundingBox(
        xMin=xC - REFINE_W / 2.0 - 1e-3, xMax=xC + REFINE_W / 2.0 + 1e-3,
        yMin=-1e-3, yMax=REFINE_H + 1e-3)
    if len(refine_edges):
        part.seedEdgeBySize(edges=refine_edges, size=ELEM_SIZE_REFINE,
                             constraint=FINER)
    part.generateMesh()
    n_elem = len(part.elements)
    n_node = len(part.nodes)
    print('  mesh: %d CPS3 elements, %d nodes' % (n_elem, n_node))

    # ---- node sets ----
    def pick_nearest(p, xt, yt):
        best, bd = None, 1e30
        for nd in p.nodes:
            dd = (nd.coordinates[0]-xt)**2 + (nd.coordinates[1]-yt)**2
            if dd < bd:
                bd, best = dd, nd
        return p.nodes.sequenceFromLabels((best.label,))

    def pick_top_nearest(p, xt, y_top, n):
        tol = 1e-3 * abs(y_top) + 1e-6
        cands = [(abs(nd.coordinates[0]-xt), nd.label)
                 for nd in p.nodes
                 if abs(nd.coordinates[1]-y_top) <= tol]
        cands.sort()
        return p.nodes.sequenceFromLabels(
            tuple(lb for _, lb in cands[:max(1, n)]))

    node_tol = 0.5 * ELEM_SIZE_REFINE

    part.Set(nodes=pick_nearest(part, overhang, 0.0),     name='Support_Left')
    part.Set(nodes=pick_nearest(part, L - overhang, 0.0), name='Support_Right')
    part.Set(nodes=pick_top_nearest(part, xC, D, 3),      name='Load_Nodes')
    part.Set(nodes=part.nodes.getByBoundingBox(
        xMin=xL_notch-node_tol, xMax=xL_notch+node_tol,
        yMin=-node_tol, yMax=node_tol),                   name='CMOD1')
    part.Set(nodes=part.nodes.getByBoundingBox(
        xMin=xR_notch-node_tol, xMax=xR_notch+node_tol,
        yMin=-node_tol, yMax=node_tol),                   name='CMOD2')

    # ---- assembly ----
    asm = m.rootAssembly
    asm.DatumCsysByDefault(CARTESIAN)
    inst = asm.Instance(name='BeamInst', part=part, dependent=ON)
    for nm in ('Support_Left', 'Support_Right', 'Load_Nodes', 'CMOD1', 'CMOD2'):
        asm.Set(nodes=inst.sets[nm].nodes, name=nm)

    # ---- step + BCs ----
    m.StaticStep(name='Loading', previous='Initial',
                  maxNumInc=MAX_NUM_INC, initialInc=1.0/N_INC,
                  minInc=1e-12, maxInc=1.0/N_INC, nlgeom=OFF)
    m.DisplacementBC(name='BC_SupL', createStepName='Loading',
                      region=asm.sets['Support_Left'], u1=SET, u2=SET)
    m.DisplacementBC(name='BC_SupR', createStepName='Loading',
                      region=asm.sets['Support_Right'], u2=SET)
    m.DisplacementBC(name='BC_Load', createStepName='Loading',
                      region=asm.sets['Load_Nodes'], u2=U_FINAL)

    # ---- output requests ----
    # IMPORTANT: delete all field outputs to keep the ODB very small.
    # We only need history output for Load-CMOD.
    for _key in list(m.fieldOutputRequests.keys()):
        del m.fieldOutputRequests[_key]

    # Keep only nodal history needed for load vs CMOD.
    m.HistoryOutputRequest(name='H-Load',  createStepName='Loading',
                            variables=('U2', 'RF2'),
                            region=asm.sets['Load_Nodes'])
    m.HistoryOutputRequest(name='H-CMOD1', createStepName='Loading',
                            variables=('U1',), region=asm.sets['CMOD1'])
    m.HistoryOutputRequest(name='H-CMOD2', createStepName='Loading',
                            variables=('U1',), region=asm.sets['CMOD2'])

    # ---- write .inp ----
    if MODEL in mdb.jobs.keys():
        del mdb.jobs[MODEL]
    mdb.Job(name=MODEL, model=MODEL,
            description='Gregoire D=100 a/D=0.2 3PB, SCM CDM UMAT')

    cwd0 = os.getcwd()
    os.chdir(CASE_DIR)
    try:
        mdb.jobs[MODEL].writeInput(consistencyChecking=OFF)
    finally:
        os.chdir(cwd0)

    inp_path = os.path.join(CASE_DIR, MODEL + '.inp')
    t1 = _time.time()
    print('  .inp written: ' + inp_path)

    # FIX-J1: copy UMAT next to the .inp immediately after build
    if os.path.exists(UMAT_SRC):
        shutil.copy2(UMAT_SRC, UMAT_LOCAL)
        print('  UMAT copied : ' + UMAT_LOCAL)
    else:
        print('  WARNING: UMAT source not found at ' + UMAT_SRC)

    td = _load_timing()
    td['build_wall_s'] = round(t1 - t0, 2)
    td['n_elements']   = n_elem
    td['n_nodes']      = n_node
    _save_timing(td)

# =====================================================================
# PHASE B  — submit solver job               (plain Python subprocess)
#
# All FIX-J fixes live here.
# =====================================================================
def phase_solve():
    t0 = _time.time()

    # FIX-J4 — sanity checks before touching the licence token
    print('  Pre-submission checks ...')
    errs = _check_before_solve()
    if errs:
        print('\nERROR: Cannot submit job. Reasons:')
        for e in errs:
            print(e)
        print('\nFix the issues above and re-run.')
        return 1
    print('  All checks passed.')

    umat_name = os.path.basename(UMAT_LOCAL)  # 'scm_umat_2d.for'

    # FIX-J1 (relative path) + FIX-J2 (double=both) + FIX-J3 (input=)
    cmd = [
        'abaqus',
        'job='   + MODEL,
        'input=' + MODEL + '.inp',  # FIX-J3: explicit deck name
        'user='  + umat_name,       # FIX-J1: relative — no spaces/backslash
        'double=both',              # FIX-J2: REAL*8 UMAT needs full double
        'interactive',
        'cpus=4',
        'ask_delete=OFF',
    ]

    print('  Solver command : ' + ' '.join(cmd))
    print('  Working dir    : ' + CASE_DIR)

    # FIX-J5: _abaqus_call uses shell=True on Windows
    rc = _abaqus_call(cmd, cwd=CASE_DIR)

    solver_wall = round(_time.time() - t0, 2)
    td = _load_timing()
    td['solver_subprocess_wall_s'] = solver_wall
    _save_timing(td)
    print('  done in %.1f s' % solver_wall)

    if rc != 0:
        # Show the tail of the log so the user sees the real error
        log_path = os.path.join(CASE_DIR, MODEL + '.log')
        if os.path.exists(log_path):
            print('\n--- last 30 lines of %s ---' % log_path)
            with open(log_path) as fh:
                lines = fh.readlines()
            for ln in lines[-30:]:
                print(ln, end='')
            print('--- end of log ---\n')
        print('WARNING: solver rc=%d. Check .log and .dat for details.' % rc)

    return rc

# =====================================================================
# PHASE C  — extract only Load-CMOD history  (abaqus python …)
# =====================================================================
def _cleanup_heavy_abaqus_files():
    """
    Delete heavy Abaqus files after extracting Load-CMOD.
    The results folder keeps only abaqus_load_cmod.csv and abaqus_load_cmod.png.
    """
    if not DELETE_HEAVY_ABAQUS_FILES:
        return

    heavy_ext = (
        '.odb', '.sim', '.stt', '.mdl', '.pac', '.res', '.sel',
        '.dat', '.msg', '.sta', '.prt', '.com', '.log', '.lck'
    )
    for ext in heavy_ext:
        p = os.path.join(CASE_DIR, MODEL + ext)
        try:
            if os.path.exists(p):
                os.remove(p)
                print('  deleted heavy file: ' + p)
        except Exception as ex:
            print('  could not delete ' + p + ': ' + str(ex))

def phase_extract():
    from odbAccess import openOdb

    odb_path = os.path.join(CASE_DIR, MODEL + '.odb')

    if not os.path.exists(odb_path):
        print('ERROR: ODB not found: ' + odb_path)
        return 1

    if not os.path.isdir(RES_DIR):
        os.makedirs(RES_DIR)

    odb  = openOdb(odb_path, readOnly=True)
    step = odb.steps['Loading']
    asm  = odb.rootAssembly

    def gather_labels(name):
        out = set()
        for nlist in asm.nodeSets[name].nodes:
            for nd in nlist:
                out.add(nd.label)
        return out

    load_lbls  = gather_labels('LOAD_NODES')
    cmod1_lbls = gather_labels('CMOD1')
    cmod2_lbls = gather_labels('CMOD2')

    def sum_series(labels, key):
        series = None
        for hkey, hreg in step.historyRegions.items():
            if not hkey.startswith('Node '):
                continue
            try:
                lbl = int(hkey.split('.')[-1])
            except Exception:
                continue
            if lbl not in labels or key not in hreg.historyOutputs:
                continue
            data = hreg.historyOutputs[key].data
            if series is None:
                series = list(data)
            else:
                n = min(len(series), len(data))
                for i in range(n):
                    series[i] = (series[i][0], series[i][1] + data[i][1])
        return series

    def avg_series(labels, key):
        series = None
        n_nodes = 0
        for hkey, hreg in step.historyRegions.items():
            if not hkey.startswith('Node '):
                continue
            try:
                lbl = int(hkey.split('.')[-1])
            except Exception:
                continue
            if lbl not in labels or key not in hreg.historyOutputs:
                continue
            data = hreg.historyOutputs[key].data
            if series is None:
                series = list(data)
            else:
                n = min(len(series), len(data))
                for i in range(n):
                    series[i] = (series[i][0], series[i][1] + data[i][1])
            n_nodes += 1

        if series is None or n_nodes == 0:
            return None
        return [(t, v / float(n_nodes)) for t, v in series]

    rf = sum_series(load_lbls,  'RF2')
    u1 = avg_series(cmod1_lbls, 'U1')
    u2 = avg_series(cmod2_lbls, 'U1')

    if rf is None or u1 is None or u2 is None:
        print('ERROR: missing history outputs in ODB.')
        odb.close()
        return 1

    n = min(len(rf), len(u1), len(u2))
    rows = []
    for i in range(n):
        load_N  = -rf[i][1]
        cmod_mm = abs(u2[i][1] - u1[i][1])
        rows.append((i + 1, load_N, cmod_mm))

    csv_path = os.path.join(RES_DIR, 'abaqus_load_cmod.csv')
    with open(csv_path, 'w') as fh:
        fh.write('# inc, load[N], cmod[mm]\n')
        for inc, load, cmod in rows:
            fh.write('%d, %.6e, %.6e\n' % (inc, load, cmod))
    print('  wrote ' + csv_path)

    pk_idx, pk = 0, -1.0e99
    for i, (inc, load, cmod) in enumerate(rows):
        if load > pk:
            pk, pk_idx = load, i
    pk_cmod = rows[pk_idx][2]

    print('  Peak load: %.2f N  (%.4f kN)' % (pk, pk / 1000.0))
    print('  CMOD@peak: %.6f mm' % pk_cmod)

    odb.close()

    # Delete heavy Abaqus files only after the CSV has been safely written.
    _cleanup_heavy_abaqus_files()

    return 0

# =====================================================================
# PHASE D  — save only Load-CMOD plot PNG     (plain Python / matplotlib)
# =====================================================================
def phase_plot():
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('ERROR: matplotlib not found. Install: pip install matplotlib')
        return 1

    csv_path = os.path.join(RES_DIR, 'abaqus_load_cmod.csv')

    rows = []
    try:
        with open(csv_path) as fh:
            for ln in fh:
                if ln.startswith('#') or not ln.strip():
                    continue
                parts = ln.split(',')
                rows.append((int(parts[0]), float(parts[1]), float(parts[2])))
    except Exception as ex:
        print('ERROR reading ' + csv_path + ': ' + str(ex))
        return 1

    if not rows:
        print('ERROR: no Load-CMOD rows found in ' + csv_path)
        return 1

    load_v  = [r[1] / 1000.0 for r in rows]
    cmod_v  = [r[2] for r in rows]
    pk      = max(load_v)
    pk_idx  = load_v.index(pk)
    pk_cmod = cmod_v[pk_idx]

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(cmod_v, load_v, '-', lw=2.0, label='Abaqus + SCM UMAT')
    ax.plot(pk_cmod, pk, marker='*', markersize=16, linestyle='None',
            label='Peak: %.3f kN @ %.4f mm' % (pk, pk_cmod))
    ax.set_xlabel('CMOD [mm]', fontsize=12)
    ax.set_ylabel('Load [kN]', fontsize=12)
    ax.set_xlim(left=0, right=max(cmod_v) * 1.05 if max(cmod_v) > 0 else 1.0)
    ax.set_ylim(bottom=0, top=max(load_v) * 1.20 if max(load_v) > 0 else 1.0)
    ax.grid(True, alpha=0.3, linestyle=':')
    ax.legend(loc='upper right', frameon=True, framealpha=0.95, fontsize=10)
    ax.set_title('Abaqus 3PB – Load vs CMOD', fontsize=12)
    plt.tight_layout()

    png_path = os.path.join(RES_DIR, 'abaqus_load_cmod.png')
    plt.savefig(png_path, dpi=200)
    plt.close()
    print('  wrote ' + png_path)

    return 0

# =====================================================================
# DRIVER  — plain Python orchestrator
# =====================================================================
def driver():
    t_start = _time.time()

    print('=' * 60)
    print('Abaqus 3PB driver  –  single-file (SCM CDM UMAT)')
    print('  Script dir : ' + SCRIPT_DIR)
    print('  UMAT source: ' + UMAT_SRC)
    print('  UMAT local : ' + UMAT_LOCAL)
    print('  Output     : ' + RES_DIR)
    print('  Platform   : ' + ('Windows' if _ON_WINDOWS else sys.platform))
    print('=' * 60)

    if not os.path.exists(UMAT_SRC):
        sys.exit(
            '\nERROR: UMAT source not found: ' + UMAT_SRC +
            '\nMake sure scm_umat_2d.for is in the same folder as this script.')

    for d in (CASE_DIR, RES_DIR):
        if not os.path.isdir(d):
            os.makedirs(d)

    # ---- Phase A: build model + .inp ----
    print('\n[1/4] Building model + writing .inp ...')
    t0 = _time.time()
    # driver() is plain Python so __file__ is defined here.
    rc = _abaqus_call(
        ['abaqus', 'cae',
         'noGUI=' + os.path.abspath(__file__),
         '--', 'build'])
    if rc != 0:
        sys.exit('Abaqus CAE build failed (rc=%d).' % rc)
    print('  done in %.1f s' % (_time.time() - t0))

    # ---- Phase B: solve with UMAT ----
    print('\n[2/4] Submitting Abaqus solver job ...')
    rc = phase_solve()
    if rc != 0:
        print('WARNING: solver rc=%d — attempting extraction anyway.' % rc)

    # ---- Phase C: extract ----
    print('\n[3/4] Extracting results from .odb ...')
    t0 = _time.time()
    rc = _abaqus_call(
        ['abaqus', 'python',
         os.path.abspath(__file__),
         '--', 'extract'])
    if rc != 0:
        sys.exit('Extraction failed (rc=%d).' % rc)
    print('  done in %.1f s' % (_time.time() - t0))

    # ---- Phase D: plot (plain Python) ----
    print('\n[4/4] Generating Load-CMOD plot only ...')
    t0 = _time.time()
    rc = phase_plot()
    plot_wall = round(_time.time() - t0, 2)
    if rc != 0:
        print('WARNING: plot step returned error.')
    print('  done in %.1f s' % plot_wall)

    t_total = round(_time.time() - t_start, 2)

    print('\n' + '=' * 60)
    print('DONE.  Wall-clock: %.1f s  (%.1f min)' % (t_total, t_total / 60.0))
    print('Outputs in: ' + RES_DIR)
    for f in ('abaqus_load_cmod.csv', 'abaqus_load_cmod.png'):
        tick = '[OK]' if os.path.exists(os.path.join(RES_DIR, f)) else '[--]'
        print('  %s %s' % (tick, f))
    print('=' * 60)

# =====================================================================
# MODE DETECTION
# =====================================================================
def _mode():
    last = sys.argv[-1] if len(sys.argv) > 1 else ''
    if last == 'build':
        return 'build'
    if last == 'extract':
        return 'extract'
    try:
        import abaqus  # noqa
        return 'build'
    except ImportError:
        return 'driver'

if __name__ == '__main__':
    mode = _mode()
    if   mode == 'driver':  driver()
    elif mode == 'build':   phase_build()
    elif mode == 'extract': sys.exit(phase_extract())