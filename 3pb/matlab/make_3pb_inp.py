
from __future__ import print_function
import os
import sys

try:
    from abaqus import mdb
    from abaqusConstants import *
    import mesh as abqMesh
except Exception as e:
    print('ERROR: run with Abaqus/CAE Python:')
    print('  abaqus cae noGUI=make_3pb_inp.py')
    print('Import error: %s' % str(e))
    sys.exit(1)

MODEL = 'Gregoire_3PB'

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

ELEM_SIZE_GLOBAL = D / 16.0
ELEM_SIZE_REFINE = ELEM_SIZE_GLOBAL / 5.0
REFINE_W         = 0.5 * D
REFINE_H         = D

E_C, NU_C, FT, GF, FCFT = 37000.0, 0.20, 3.5, 0.090, 10.0

U_FINAL = -0.2
N_INC   = 1000

FIELD_VARS = ('SDV',)
FIELD_FREQ = 10

def die(msg):
    print('ERROR: ' + str(msg))
    sys.exit(1)

def ensure_dir(path):
    if not os.path.isdir(path):
        os.makedirs(path)

def pick_single_nearest(part, xt, yt):
    best, bd = None, 1.0e99
    for nd in part.nodes:
        x, y = nd.coordinates[0], nd.coordinates[1]
        d2 = (x - xt) ** 2 + (y - yt) ** 2
        if d2 < bd:
            bd, best = d2, nd
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

def build_model():
    print('\n==== Build model ====')
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
    mat.Depvar(n=2)
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

    m.FieldOutputRequest(name='F-SDV', createStepName='Loading',
                         variables=FIELD_VARS, frequency=FIELD_FREQ)
    m.HistoryOutputRequest(name='H-Load', createStepName='Loading',
                           variables=('RF2',), region=asm.sets['Load_Nodes'])
    m.HistoryOutputRequest(name='H-CMOD1', createStepName='Loading',
                           variables=('U1',), region=asm.sets['CMOD1'])
    m.HistoryOutputRequest(name='H-CMOD2', createStepName='Loading',
                           variables=('U1',), region=asm.sets['CMOD2'])
    print('Build complete.')
    return m

def write_inp(case_dir):
    print('\n==== Write .inp (no run) ====')
    if MODEL in mdb.jobs.keys():
        del mdb.jobs[MODEL]
    job = mdb.Job(name=MODEL, model=MODEL, description='Gregoire 3PB CDM UMAT (inp only)')
    cwd0 = os.getcwd()
    ensure_dir(case_dir)
    os.chdir(case_dir)
    try:
        job.writeInput(consistencyChecking=OFF)
    finally:
        os.chdir(cwd0)
    inp_path = os.path.join(case_dir, MODEL + '.inp')
    if not os.path.exists(inp_path):
        die('.inp was not written.')
    print('Wrote: %s' % inp_path)
    return inp_path

def main():
    base = os.path.abspath(os.getcwd())
    case_dir = os.path.join(base, MODEL)
    ensure_dir(case_dir)
    print('=' * 60)
    print('Gregoire 3PB  ->  WRITE .inp ONLY')
    print('Base: %s' % base)
    print('=' * 60)
    build_model()
    inp_path = write_inp(case_dir)
    print('\nDONE.')
    print('Run later with:')
    print('  abaqus job=%s input=%s user=cdm_umat_2d.for cpus=4'
          % (MODEL, MODEL))

if __name__ == '__main__':
    main()
