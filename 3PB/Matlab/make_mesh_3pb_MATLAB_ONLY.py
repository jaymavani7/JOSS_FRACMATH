# -*- coding: utf-8 -*-
# =====================================================================
# make_mesh_3pb_MATLAB_ONLY.py
#
# Purpose:
#   This script is ONLY for making the Abaqus/CAE mesh and exporting the
#   mesh/sets needed by the MATLAB solver.
#
# It does NOT:
#   - run Abaqus analysis
#   - open an ODB
#   - extract Load-CMOD from Abaqus
#
# It DOES:
#   1) Build the Gregoire 3PB notched beam mesh in Abaqus/CAE.
#   2) Write Gregoire_3PB/Gregoire_3PB.inp.
#   3) Parse the .inp into MATLAB-readable files:
#        nodes.txt
#        elements.txt
#        top_nodes.txt
#        left_nodes.txt
#        right_nodes.txt
#        cmod1.txt
#        cmod2.txt
#        boundary_conditions.txt
#        oliver_bandwidth_data.dat
#
# Run option A, from normal Windows/Linux terminal with Python available:
#     python make_mesh_3pb_MATLAB_ONLY.py
#
# Run option B, using only Abaqus commands:
#     abaqus cae noGUI=make_mesh_3pb_MATLAB_ONLY.py -- build
#     abaqus python make_mesh_3pb_MATLAB_ONLY.py parse
#
# Then in MATLAB:
#     Run solver_main_3pb.m / solver_main_3pb_OLIVER_bandwidth()
#
# =====================================================================

from __future__ import print_function

import os
import re
import sys
import math
import subprocess


# =====================================================================
# CONFIG: must match MATLAB solver
# =====================================================================
MODEL = 'Gregoire_3PB'

# Geometry in mm
D = 100.0
S = 2.5 * D
overhang = 0.5 * D
L = S + 2.0 * overhang

b_thick = 50.0
a0 = 0.2 * D
wn = D / 40.0

xC = L / 2.0
xL_notch = xC - wn / 2.0
xR_notch = xC + wn / 2.0

# Mesh
ELEM_SIZE_GLOBAL = D / 16.0
ELEM_SIZE_REFINE = ELEM_SIZE_GLOBAL / 5.0
REFINE_W = 0.5 * D
REFINE_H = D

# Dummy Abaqus material/step values.
# MATLAB solver uses its own material constants; these are only written
# so the .inp is complete and the BC sets are clear.
E_C = 37000.0
NU_C = 0.20
U_FINAL = -0.2
N_INC = 200

CWD = os.getcwd()
CASE_DIR = os.path.join(CWD, MODEL)


# =====================================================================
# PHASE A: build Abaqus/CAE mesh and write .inp
# =====================================================================
def phase_build():
    from abaqus import mdb
    from abaqusConstants import (
        TWO_D_PLANAR, DEFORMABLE_BODY, SIDE1,
        FINER, FREE, TRI, ADVANCING_FRONT,
        CARTESIAN, ON, OFF,
        CPS3, STANDARD, SET
    )
    import mesh as abqMesh

    if not os.path.exists(CASE_DIR):
        os.makedirs(CASE_DIR)

    if MODEL in mdb.models.keys():
        del mdb.models[MODEL]

    m = mdb.Model(name=MODEL)

    # ---------------------------------------------------------------
    # Geometry: 3PB notched beam
    # ---------------------------------------------------------------
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

    # Local partition around notch for refined mesh
    tform = part.MakeSketchTransform(sketchPlane=part.faces[0],
                                     sketchPlaneSide=SIDE1,
                                     origin=(0.0, 0.0, 0.0))
    sk_r = m.ConstrainedSketch(name='RefineZone', sheetSize=L, transform=tform)
    sk_r.rectangle(point1=(xC - REFINE_W / 2.0, 0.0),
                   point2=(xC + REFINE_W / 2.0, REFINE_H))
    part.PartitionFaceBySketch(sketch=sk_r, faces=part.faces)
    del sk_r

    # Dummy Abaqus section/material.
    # MATLAB ignores this and supplies material properties in solver_main.
    mat = m.Material(name='Concrete_dummy_for_mesh_only')
    mat.Elastic(table=((E_C, NU_C),))

    m.HomogeneousSolidSection(name='BeamSec',
                              material='Concrete_dummy_for_mesh_only',
                              thickness=b_thick)
    part.SectionAssignment(region=(part.faces,), sectionName='BeamSec')

    # CPS3 triangle mesh, same type expected by MATLAB T3 solver
    elem_t3 = abqMesh.ElemType(elemCode=CPS3, elemLibrary=STANDARD)
    part.setElementType(regions=(part.faces,), elemTypes=(elem_t3,))

    for f in part.faces:
        part.setMeshControls(regions=(f,), technique=FREE,
                             elemShape=TRI, algorithm=ADVANCING_FRONT)

    part.seedPart(size=ELEM_SIZE_GLOBAL, deviationFactor=0.1)

    refine_edges = part.edges.getByBoundingBox(
        xMin=xC - REFINE_W / 2.0 - 1e-3,
        xMax=xC + REFINE_W / 2.0 + 1e-3,
        yMin=-1e-3,
        yMax=REFINE_H + 1e-3
    )

    if len(refine_edges):
        part.seedEdgeBySize(edges=refine_edges,
                            size=ELEM_SIZE_REFINE,
                            constraint=FINER)

    part.generateMesh()

    # ---------------------------------------------------------------
    # Node set picking
    # ---------------------------------------------------------------
    def pick_single_nearest(p, xt, yt):
        best = None
        bd = 1.0e30
        for nd in p.nodes:
            x = nd.coordinates[0]
            y = nd.coordinates[1]
            d2 = (x - xt) ** 2 + (y - yt) ** 2
            if d2 < bd:
                bd = d2
                best = nd
        return p.nodes.sequenceFromLabels(labels=(best.label,))

    def pick_n_nearest_top(p, xt, y_top, n_keep):
        tol = 1.0e-3 * abs(y_top) + 1.0e-6
        cands = []
        for nd in p.nodes:
            x = nd.coordinates[0]
            y = nd.coordinates[1]
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
        xMin=xL_notch - node_tol,
        xMax=xL_notch + node_tol,
        yMin=-node_tol,
        yMax=node_tol),
        name='CMOD1'
    )

    part.Set(nodes=part.nodes.getByBoundingBox(
        xMin=xR_notch - node_tol,
        xMax=xR_notch + node_tol,
        yMin=-node_tol,
        yMax=node_tol),
        name='CMOD2'
    )

    # Assembly sets, so the .inp has normal Abaqus NSET names
    asm = m.rootAssembly
    asm.DatumCsysByDefault(CARTESIAN)
    inst = asm.Instance(name='BeamInst', part=part, dependent=ON)

    for nm in ('Support_Left', 'Support_Right', 'Load_Nodes', 'CMOD1', 'CMOD2'):
        asm.Set(nodes=inst.sets[nm].nodes, name=nm)

    # Dummy step/BCs only to write boundary_conditions.txt consistently.
    # MATLAB uses the exported node-set files for actual solve.
    m.StaticStep(name='Loading',
                 previous='Initial',
                 maxNumInc=4000,
                 initialInc=1.0 / N_INC,
                 minInc=1e-10,
                 maxInc=1.0 / N_INC,
                 nlgeom=OFF)

    m.DisplacementBC(name='BC_SupL',
                     createStepName='Loading',
                     region=asm.sets['Support_Left'],
                     u1=SET, u2=SET)

    m.DisplacementBC(name='BC_SupR',
                     createStepName='Loading',
                     region=asm.sets['Support_Right'],
                     u2=SET)

    m.DisplacementBC(name='BC_Load',
                     createStepName='Loading',
                     region=asm.sets['Load_Nodes'],
                     u2=U_FINAL)

    if MODEL in mdb.jobs.keys():
        del mdb.jobs[MODEL]

    job = mdb.Job(name=MODEL,
                  model=MODEL,
                  description='Mesh-only Gregoire 3PB export for MATLAB solver')

    cwd0 = os.getcwd()
    os.chdir(CASE_DIR)
    try:
        job.writeInput(consistencyChecking=OFF)
    finally:
        os.chdir(cwd0)

    print('  mesh created: %d CPS3 elements, %d nodes' %
          (len(part.elements), len(part.nodes)))
    print('  inp written : %s' % os.path.join(CASE_DIR, MODEL + '.inp'))
    print('  no Abaqus analysis was run; no ODB is needed')


# =====================================================================
# PHASE B: parse .inp into MATLAB files
# =====================================================================
def phase_parse():
    inp_path = os.path.join(CASE_DIR, MODEL + '.inp')

    if not os.path.exists(inp_path):
        sys.exit('ERROR: .inp not found: %s' % inp_path)

    nodes, elems, sets, bcs = read_inp_mesh_sets_bcs(inp_path)

    def W(name, lines_):
        path = os.path.join(CASE_DIR, name)
        with open(path, 'w') as f:
            if lines_:
                f.write('\n'.join(lines_) + '\n')
        return path

    W('nodes.txt', ['%d %.16e %.16e' % (nid, x, y)
                    for (nid, x, y) in nodes])

    W('elements.txt', ['%d %d %d %d' % e for e in elems])

    W('top_nodes.txt', ['%d' % i for i in sets['Load_Nodes']])
    W('left_nodes.txt', ['%d' % i for i in sets['Support_Left']])
    W('right_nodes.txt', ['%d' % i for i in sets['Support_Right']])
    W('cmod1.txt', ['%d' % i for i in sets['CMOD1']])
    W('cmod2.txt', ['%d' % i for i in sets['CMOD2']])

    W('boundary_conditions.txt',
      ['%s %d %.16e' % (nm, dof, val) for (nm, dof, val) in bcs])

    write_oliver_bandwidth_from_mesh(
        nodes, elems,
        os.path.join(CASE_DIR, 'oliver_bandwidth_data.dat')
    )

    print('  parsed nodes    : %d' % len(nodes))
    print('  parsed elements : %d' % len(elems))
    print('  Load_Nodes      : %s' % sets['Load_Nodes'])
    print('  Support_Left    : %s' % sets['Support_Left'])
    print('  Support_Right   : %s' % sets['Support_Right'])
    print('  CMOD1           : %s' % sets['CMOD1'])
    print('  CMOD2           : %s' % sets['CMOD2'])
    print('  MATLAB files written in: %s' % CASE_DIR)


def read_inp_mesh_sets_bcs(inp_path):
    with open(inp_path) as f:
        lines = f.readlines()

    def is_kw(s):
        return s.strip().startswith('*')

    def to_nums(s):
        out = []
        for x in s.replace(',', ' ').split():
            if x.strip():
                out.append(float(x))
        return out

    nodes = []
    elems = []
    sets = {
        'Support_Left': [],
        'Support_Right': [],
        'Load_Nodes': [],
        'CMOD1': [],
        'CMOD2': []
    }
    bcs = []

    i = 0
    n = len(lines)

    while i < n:
        ln = lines[i].strip()
        low = ln.lower()

        if low.startswith('*node') and not low.startswith('*node output'):
            i += 1
            while i < n and not is_kw(lines[i]):
                v = to_nums(lines[i])
                if len(v) >= 3:
                    nodes.append((int(v[0]), v[1], v[2]))
                i += 1
            continue

        if low.startswith('*element') and not low.startswith('*element output'):
            i += 1
            while i < n and not is_kw(lines[i]):
                v = to_nums(lines[i])
                if len(v) >= 4:
                    elems.append((int(v[0]), int(v[1]), int(v[2]), int(v[3])))
                i += 1
            continue

        if low.startswith('*nset'):
            tok = re.search(r'nset\s*=\s*([^,\s]+)', ln, re.I)
            name = tok.group(1).strip() if tok else None
            is_gen = ('generate' in low)

            i += 1
            acc = []
            while i < n and not is_kw(lines[i]):
                nv = to_nums(lines[i])
                if is_gen:
                    for j in range(0, len(nv), 3):
                        if j + 2 < len(nv):
                            start = int(nv[j])
                            stop = int(nv[j + 1])
                            step = int(nv[j + 2])
                            acc.extend(range(start, stop + 1, step))
                else:
                    acc.extend(int(x) for x in nv)
                i += 1

            if name in sets:
                sets[name] = sorted(set(sets[name] + acc))
            continue

        if low.startswith('*boundary') and not low.startswith('*boundary output'):
            i += 1
            while i < n and not is_kw(lines[i]):
                pt = [p.strip() for p in lines[i].split(',')]
                if len(pt) >= 3:
                    nm = pt[0]
                    dof1 = int(pt[1])
                    dof2 = int(pt[2])
                    val = 0.0
                    if len(pt) >= 4 and pt[3]:
                        try:
                            val = float(pt[3])
                        except Exception:
                            val = 0.0
                    for dof in range(dof1, dof2 + 1):
                        bcs.append((nm, dof, val))
                i += 1
            continue

        i += 1

    return nodes, elems, sets, bcs


# =====================================================================
# Oliver crack-band gradient table for MATLAB / UMAT comparison
# =====================================================================
def write_oliver_bandwidth_from_mesh(nodes, elems, out_path):
    """Write Oliver gradient table from the exact mesh.

    Columns:
      elem_id g1x g1y g2x g2y g3x g3y area sqrt_area

    These gradients are the same T3 shape-function derivatives used in
    MATLAB precompute_T3(). The direction-dependent bandwidth itself is:

      h(n) = 2 / sum_a |grad(N_a) dot n|

    where n is the crack normal / max principal strain direction.
    """
    node_xy = dict((nid, (x, y)) for nid, x, y in nodes)

    with open(out_path, 'w') as f:
        for eid, n1, n2, n3 in elems:
            x1, y1 = node_xy[n1]
            x2, y2 = node_xy[n2]
            x3, y3 = node_xy[n3]

            area_signed = 0.5 * ((x2 - x1) * (y3 - y1) -
                                 (x3 - x1) * (y2 - y1))
            area = abs(area_signed)

            if area <= 1.0e-14:
                raise RuntimeError('Zero/near-zero CPS3 area at element %d' % eid)

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

            f.write('%d %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n' %
                    (eid, g1x, g1y, g2x, g2y, g3x, g3y, area, math.sqrt(area)))

    print('  Oliver bandwidth table written: %s' % out_path)


# =====================================================================
# DRIVER
# =====================================================================
def driver():
    print('=' * 70)
    print('Gregoire 3PB mesh-only exporter for MATLAB')
    print('Geometry: D=100 mm, a/D=0.2, span=250 mm, CPS3 triangles')
    print('=' * 70)

    print('\n[1/2] Building mesh in Abaqus/CAE...')
    rc = subprocess.call([
        'abaqus', 'cae',
        'noGUI=' + os.path.abspath(__file__),
        '--', 'build'
    ])

    if rc != 0:
        sys.exit('ERROR: Abaqus CAE mesh build failed.')

    print('\n[2/2] Parsing .inp into MATLAB mesh files...')
    phase_parse()

    print('\nDONE.')
    print('Use this MATLAB folder:')
    print('  %s' % CASE_DIR)
    print('\nDo NOT run Abaqus extractor for this mesh-only workflow.')
    print('No ODB file is required.')


def _mode():
    # Abaqus passes script arguments after "--"; keep this simple and robust.
    args = [a.lower() for a in sys.argv]
    if 'build' in args:
        return 'build'
    if 'parse' in args:
        return 'parse'

    # If imported by Abaqus/CAE without arguments, build only.
    try:
        import abaqus  # noqa
        return 'build'
    except Exception:
        return 'driver'


if __name__ == '__main__':
    mode = _mode()

    if mode == 'build':
        phase_build()
    elif mode == 'parse':
        phase_parse()
    else:
        driver()
