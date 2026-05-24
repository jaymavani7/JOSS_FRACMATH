# -*- coding: utf-8 -*-
# =====================================================================
# make_mesh_3pb.py
#
# SINGLE file: build the 3PB geometry in Abaqus/CAE, write the .inp,
# THEN parse that .inp into plain-text files for the MATLAB solver to
# read (nodes.txt, elements.txt, boundary conditions, CMOD nodes).
#
# Two-phase script:
#   Phase A (inside abaqus cae) -- build geometry, mesh, write .inp
#   Phase B (plain python)      -- parse .inp -> .txt files
#
# Run from a normal shell:
#     python3 make_mesh_3pb.py
#
# Output (in folder Gregoire_3PB/):
#     Gregoire_3PB.inp        Abaqus deck
#     nodes.txt               id  x  y
#     elements.txt            id  n1  n2  n3
#     top_nodes.txt           IDs of load nodes
#     left_nodes.txt          IDs of left-support nodes
#     right_nodes.txt         IDs of right-support nodes
#     cmod1.txt               IDs of left-notch-lip nodes
#     cmod2.txt               IDs of right-notch-lip nodes
#     boundary_conditions.txt nset_name, dof, value
#
# Same data as the MATLAB solver expects:
#     Gregoire D=100, a/D=0.2 (FifthNotched), Medium mesh
#
# Authors: [Name to be added]
# =====================================================================

import os
import re
import sys
import subprocess


# =====================================================================
# CONFIG
# =====================================================================
MODEL    = 'Gregoire_3PB'

# Geometry (mm)
D         = 100.0
S         = 2.5 * D
overhang  = 0.5 * D
L         = S + 2.0 * overhang
b_thick   = 50.0
a0        = 0.2 * D
wn        = D / 40.0
xC        = L / 2.0
xL_notch  = xC - wn/2.0
xR_notch  = xC + wn/2.0

ELEM_SIZE_GLOBAL = D / 16
ELEM_SIZE_REFINE = ELEM_SIZE_GLOBAL / 5.0
REFINE_W = 0.5 * D
REFINE_H = D

E_C, NU_C, FT, GF, FCFT = 37000.0, 0.20, 3.5, 0.090, 10.0
U_FINAL = -0.2
N_INC   = 200

CWD       = os.getcwd()
CASE_DIR  = os.path.join(CWD, MODEL)


# =====================================================================
# PHASE A: build the model in CAE, write .inp
# =====================================================================
def phase_build():
    from abaqus import mdb
    from abaqusConstants import (TWO_D_PLANAR, DEFORMABLE_BODY, SIDE1,
                                  FINER, FREE, TRI, ADVANCING_FRONT,
                                  CARTESIAN, ON, OFF,
                                  CPS3, STANDARD, SET)
    import mesh as abqMesh

    if not os.path.exists(CASE_DIR):
        os.makedirs(CASE_DIR)

    if MODEL in mdb.models.keys():
        del mdb.models[MODEL]
    m = mdb.Model(name=MODEL)

    sk = m.ConstrainedSketch(name='BeamSketch', sheetSize=2.0*L)
    sk.Line(point1=(0.0,0.0),         point2=(xL_notch,0.0))
    sk.Line(point1=(xL_notch,0.0),    point2=(xL_notch,a0))
    sk.Line(point1=(xL_notch,a0),     point2=(xR_notch,a0))
    sk.Line(point1=(xR_notch,a0),     point2=(xR_notch,0.0))
    sk.Line(point1=(xR_notch,0.0),    point2=(L,0.0))
    sk.Line(point1=(L,0.0), point2=(L,D))
    sk.Line(point1=(L,D),   point2=(0.0,D))
    sk.Line(point1=(0.0,D), point2=(0.0,0.0))

    part = m.Part(name='Beam2D', dimensionality=TWO_D_PLANAR,
                   type=DEFORMABLE_BODY)
    part.BaseShell(sketch=sk)
    del sk

    tform = part.MakeSketchTransform(sketchPlane=part.faces[0],
                                      sketchPlaneSide=SIDE1,
                                      origin=(0.0,0.0,0.0))
    sk_r = m.ConstrainedSketch(name='RefineZone', sheetSize=L, transform=tform)
    sk_r.rectangle(point1=(xC - REFINE_W/2.0, 0.0),
                   point2=(xC + REFINE_W/2.0, REFINE_H))
    part.PartitionFaceBySketch(sketch=sk_r, faces=part.faces)
    del sk_r

    # Placeholder elastic material (the MATLAB solver supplies its own)
    mat = m.Material(name='Concrete')
    mat.Elastic(table=((E_C, NU_C),))
    m.HomogeneousSolidSection(name='BeamSec', material='Concrete',
                              thickness=b_thick)
    part.SectionAssignment(region=(part.faces,), sectionName='BeamSec')

    elem_t3 = abqMesh.ElemType(elemCode=CPS3, elemLibrary=STANDARD)
    part.setElementType(regions=(part.faces,), elemTypes=(elem_t3,))
    for f in part.faces:
        part.setMeshControls(regions=(f,), technique=FREE,
                              elemShape=TRI, algorithm=ADVANCING_FRONT)
    part.seedPart(size=ELEM_SIZE_GLOBAL, deviationFactor=0.1)
    refine_edges = part.edges.getByBoundingBox(
        xMin=xC - REFINE_W/2.0 - 1e-3, xMax=xC + REFINE_W/2.0 + 1e-3,
        yMin=-1e-3, yMax=REFINE_H + 1e-3)
    if len(refine_edges):
        part.seedEdgeBySize(edges=refine_edges, size=ELEM_SIZE_REFINE,
                             constraint=FINER)
    part.generateMesh()

    def pick_single_nearest(p, xt, yt):
        best, bd = None, 1e30
        for nd in p.nodes:
            x,y = nd.coordinates[0], nd.coordinates[1]
            d = (x-xt)**2 + (y-yt)**2
            if d < bd: bd, best = d, nd
        return p.nodes.sequenceFromLabels(labels=(best.label,))

    def pick_n_nearest_top(p, xt, y_top, n):
        tol = 1e-3*abs(y_top) + 1e-6
        cands = []
        for nd in p.nodes:
            x,y = nd.coordinates[0], nd.coordinates[1]
            if abs(y-y_top) > tol: continue
            cands.append((abs(x-xt), nd.label))
        cands.sort()
        keep = [lb for _,lb in cands[:max(1,n)]]
        return p.nodes.sequenceFromLabels(labels=tuple(keep))

    node_tol = 0.5 * ELEM_SIZE_REFINE

    part.Set(nodes=pick_single_nearest(part, overhang,    0.0), name='Support_Left')
    part.Set(nodes=pick_single_nearest(part, L-overhang,  0.0), name='Support_Right')
    part.Set(nodes=pick_n_nearest_top(part, xC, D, 3),          name='Load_Nodes')
    part.Set(nodes=part.nodes.getByBoundingBox(
        xMin=xL_notch-node_tol, xMax=xL_notch+node_tol,
        yMin=-node_tol, yMax=node_tol), name='CMOD1')
    part.Set(nodes=part.nodes.getByBoundingBox(
        xMin=xR_notch-node_tol, xMax=xR_notch+node_tol,
        yMin=-node_tol, yMax=node_tol), name='CMOD2')

    asm = m.rootAssembly
    asm.DatumCsysByDefault(CARTESIAN)
    inst = asm.Instance(name='BeamInst', part=part, dependent=ON)
    for nm in ('Support_Left','Support_Right','Load_Nodes','CMOD1','CMOD2'):
        asm.Set(nodes=inst.sets[nm].nodes, name=nm)

    m.StaticStep(name='Loading', previous='Initial',
                  maxNumInc=4000, initialInc=1.0/N_INC,
                  minInc=1e-10, maxInc=1.0/N_INC, nlgeom=OFF)
    m.DisplacementBC(name='BC_SupL', createStepName='Loading',
                      region=asm.sets['Support_Left'], u1=SET, u2=SET)
    m.DisplacementBC(name='BC_SupR', createStepName='Loading',
                      region=asm.sets['Support_Right'], u2=SET)
    m.DisplacementBC(name='BC_Load', createStepName='Loading',
                      region=asm.sets['Load_Nodes'], u2=U_FINAL)

    if MODEL in mdb.jobs.keys():
        del mdb.jobs[MODEL]
    job = mdb.Job(name=MODEL, model=MODEL,
                   description='Gregoire D=100 a/D=0.2 3PB (mesh only)')

    cwd0 = os.getcwd()
    os.chdir(CASE_DIR)
    try:
        job.writeInput(consistencyChecking=OFF)
    finally:
        os.chdir(cwd0)

    print('  mesh: %d CPS3 elements, %d nodes'
          % (len(part.elements), len(part.nodes)))
    print('  .inp: ' + os.path.join(CASE_DIR, MODEL + '.inp'))


# =====================================================================
# PHASE B: parse the .inp into plain-text files
# =====================================================================
def phase_parse():
    inp_path = os.path.join(CASE_DIR, MODEL + '.inp')
    if not os.path.exists(inp_path):
        sys.exit('ERROR: .inp not found at ' + inp_path)

    with open(inp_path) as f:
        lines = f.readlines()

    is_kw   = lambda s: s.strip().startswith('*')
    to_nums = lambda s: [float(x) for x in s.replace(',', ' ').split()
                          if x.strip()]

    nodes = []   # (id, x, y)
    elems = []   # (id, n1, n2, n3)
    sets  = {
        'Support_Left':  [], 'Support_Right': [],
        'Load_Nodes':    [], 'CMOD1':         [], 'CMOD2': [],
    }
    bcs   = []   # (nset, dof, value)

    i, n = 0, len(lines)
    while i < n:
        ln  = lines[i].strip()
        low = ln.lower()

        if low.startswith('*node') and not low.startswith('*node o'):
            i += 1
            while i < n and not is_kw(lines[i]):
                v = to_nums(lines[i])
                if len(v) >= 3:
                    nodes.append((int(v[0]), v[1], v[2]))
                i += 1
            continue

        if low.startswith('*element') and not low.startswith('*element o'):
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
            is_gen = 'generate' in low
            i += 1
            acc = []
            while i < n and not is_kw(lines[i]):
                nv = to_nums(lines[i])
                if is_gen:
                    # 'start, stop, step' triplets
                    for j in range(0, len(nv), 3):
                        if j+2 < len(nv):
                            acc.extend(range(int(nv[j]), int(nv[j+1])+1,
                                              int(nv[j+2])))
                else:
                    acc.extend(int(x) for x in nv)
                i += 1
            if name and name in sets:
                # set is unique
                sets[name] = sorted(set(sets[name] + acc))
            continue

        if low.startswith('*boundary') and not low.startswith('*boundary o'):
            i += 1
            while i < n and not is_kw(lines[i]):
                pt = [p.strip() for p in lines[i].split(',')]
                if len(pt) >= 3:
                    v = 0.0
                    if len(pt) >= 4 and pt[3]:
                        try: v = float(pt[3])
                        except ValueError: pass
                    bcs.append((pt[0], int(pt[1]), v))
                i += 1
            continue

        i += 1

    # --- write the .txt files ----------------------------------------
    def W(name, lines_):
        path = os.path.join(CASE_DIR, name)
        with open(path, 'w') as f:
            f.write('\n'.join(lines_) + ('\n' if lines_ else ''))
        return path

    W('nodes.txt',    ['%d %.6f %.6f' % (i,x,y) for (i,x,y) in nodes])
    W('elements.txt', ['%d %d %d %d' % e for e in elems])
    W('top_nodes.txt',   ['%d' % i for i in sets['Load_Nodes']])
    W('left_nodes.txt',  ['%d' % i for i in sets['Support_Left']])
    W('right_nodes.txt', ['%d' % i for i in sets['Support_Right']])
    W('cmod1.txt',       ['%d' % i for i in sets['CMOD1']])
    W('cmod2.txt',       ['%d' % i for i in sets['CMOD2']])
    W('boundary_conditions.txt',
      ['%s %d %.6g' % (nm, dof, v) for (nm, dof, v) in bcs])

    print('  parsed: %d nodes, %d CPS3 elements' % (len(nodes), len(elems)))
    print('  load nodes: %s' % sets['Load_Nodes'])
    print('  left sup  : %s' % sets['Support_Left'])
    print('  right sup : %s' % sets['Support_Right'])
    print('  CMOD1     : %s' % sets['CMOD1'])
    print('  CMOD2     : %s' % sets['CMOD2'])
    print('  .txt files written in ' + CASE_DIR)


# =====================================================================
# DRIVER
# =====================================================================
def driver():
    print('=' * 60)
    print('MATLAB-side mesh generator')
    print('  Geometry: D=100 mm, a/D=0.2, span=250 mm, CPS3')
    print('=' * 60)

    # --- Phase A: build inside Abaqus CAE ----------------------------
    print('\n[1/2] Building geometry in Abaqus/CAE...')
    rc = subprocess.call(['abaqus', 'cae', 'noGUI=' + os.path.abspath(__file__),
                          '--', 'build'])
    if rc != 0:
        sys.exit('Abaqus CAE build failed.')

    # --- Phase B: parse .inp -> .txt (plain python) ------------------
    print('\n[2/2] Parsing .inp into .txt files for MATLAB...')
    phase_parse()

    print('\nDone. MATLAB-ready files in ' + CASE_DIR + '/')


def _mode():
    if 'build' in sys.argv: return 'build'
    if 'parse' in sys.argv: return 'parse'
    try:
        import abaqus  # noqa: F401
        return 'build'
    except ImportError:
        return 'driver'


if __name__ == '__main__':
    mode = _mode()
    if   mode == 'driver': driver()
    elif mode == 'build':  phase_build()
    elif mode == 'parse':  phase_parse()