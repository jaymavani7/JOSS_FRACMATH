
from __future__ import print_function
import os
import sys
import glob

try:
    from odbAccess import openOdb
except Exception as e:
    print('ERROR: run this INSIDE Abaqus python:')
    print('    abaqus python extract_damage.py [odb]')
    print('Import error: %s' % str(e))
    sys.exit(1)

def norm_name(s):
    return ''.join(ch for ch in str(s).upper() if ch.isalnum())

def float_component(data, comp):

    try:
        return float(data)
    except (TypeError, ValueError):
        pass
    try:
        return float(data[comp])
    except Exception:
        return float(data[0])

def find_odb(arg):
    if arg and os.path.exists(arg):
        return os.path.abspath(arg)
    cands = glob.glob('*.odb') + glob.glob(os.path.join('Gregoire_3PB', '*.odb'))
    cands = [c for c in cands if os.path.exists(c)]
    if not cands:
        print('ERROR: no .odb found. Pass the path explicitly.')
        sys.exit(1)
    cands.sort(key=os.path.getmtime)
    return os.path.abspath(cands[-1])

def ensure_dir(p):
    if not os.path.isdir(p):
        os.makedirs(p)

def get_damage_field(frame):

    keys = frame.fieldOutputs.keys()
    for k in ('SDV2', 'SDV_2', 'STATEV2', 'STATEV_2'):
        if k in keys:
            return frame.fieldOutputs[k], 0, k
    for k in keys:
        if norm_name(k) in ('SDV2', 'SDV2', 'STATEV2'):
            return frame.fieldOutputs[k], 0, k
    for k in ('SDV', 'STATEV'):
        if k in keys:
            return frame.fieldOutputs[k], 1, k
    return None, None, None

def frame_max_omega(frame):

    fld, comp, key = get_damage_field(frame)
    if fld is None:
        return None, None
    mx = -1.0e99
    for v in fld.values:
        try:
            w = float_component(v.data, comp)
        except Exception:
            continue
        if w > mx:
            mx = w
    if mx <= -1.0e98:
        return None, key
    return mx, key

def write_omega_csv(frame, path):
    if frame is None:
        open(path, 'w').write('# element_id, ip, omega\n')
        return 0, -1.0
    fld, comp, _ = get_damage_field(frame)
    if fld is None:
        open(path, 'w').write('# element_id, ip, omega\n')
        return 0, -1.0
    n, mx = 0, -1.0e99
    f = open(path, 'w')
    f.write('# element_id, ip, omega\n')
    for v in fld.values:
        try:
            w = float_component(v.data, comp)
        except Exception:
            continue
        ip = getattr(v, 'integrationPoint', 1)
        f.write('%s, %s, %.6e\n' % (str(v.elementLabel), str(ip), w))
        n += 1
        if w > mx:
            mx = w
    f.close()
    return n, (mx if mx > -1.0e98 else -1.0)

def dump_mesh(odb, pd_dir):
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
    return inst.name

def pick_indices(step):

    nf = len(step.frames)
    if nf == 0:
        return None, None, None
    last = nf - 1

    peak = last
    try:
        best = None
        for ro in step.historyRegions.values():
            for kk, out in ro.historyOutputs.items():
                if kk.upper().startswith('RF2'):

                    series = out.data

                    times = [fr.frameValue for fr in step.frames]
                    vmax, tmax = -1.0e99, None
                    for (t, val) in series:
                        av = abs(val)
                        if av > vmax:
                            vmax, tmax = av, t
                    if tmax is not None:

                        di = min(range(nf), key=lambda i: abs(times[i] - tmax))
                        best = di
        if best is not None:
            peak = best
    except Exception:
        peak = last
    postpeak = min(last, int(round(0.5 * (peak + last))))
    if postpeak <= peak and last > peak:
        postpeak = peak + 1
    return peak, postpeak, last

def try_plot(pd_dir, out_dir, thresh):
    try:
        import numpy as np
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.collections import PolyCollection
    except Exception as e:
        print('matplotlib not available in Abaqus python (%s).' % str(e))
        print('CSVs were written. Make the figures with normal Python:')
        print('    python plot_damage_from_csv.py "%s"' % out_dir)
        return False

    nodes = {}
    for ln in open(os.path.join(pd_dir, 'mesh_nodes.csv')):
        if ln.startswith('#') or not ln.strip():
            continue
        p = ln.split(',')
        nodes[int(p[0])] = (float(p[1]), float(p[2]))
    elems = []
    for ln in open(os.path.join(pd_dir, 'mesh_elements.csv')):
        if ln.startswith('#') or not ln.strip():
            continue
        p = ln.split(',')
        elems.append((int(p[0]), int(p[1]), int(p[2]), int(p[3])))
    eid_order = [e[0] for e in elems]
    polys_all = [[nodes[e[1]], nodes[e[2]], nodes[e[3]]] for e in elems]

    def load_omega(name):
        d = {}
        path = os.path.join(pd_dir, name)
        if not os.path.exists(path):
            return d
        for ln in open(path):
            if ln.startswith('#') or not ln.strip():
                continue
            p = ln.split(',')
            try:
                eid, w = int(p[0]), float(p[2])
            except Exception:
                continue
            d[eid] = max(d.get(eid, -1.0), w)
        return d

    def one_fig(csv_name, label, base):
        wmap = load_omega(csv_name)
        omega = np.array([wmap.get(eid, 0.0) for eid in eid_order])
        idx = np.where(omega >= thresh)[0]
        fh = plt.figure(figsize=(16.0 / 2.54, 7.0 / 2.54), facecolor='w')
        ax = fh.add_axes([0.10, 0.16, 0.87, 0.74])
        grey = PolyCollection(polys_all, facecolors=(0.86, 0.86, 0.86),
                              edgecolors=(0.62, 0.62, 0.62), linewidths=0.25)
        ax.add_collection(grey)
        if idx.size:
            blk = PolyCollection([polys_all[i] for i in idx],
                                 facecolors=(0.0, 0.0, 0.0),
                                 edgecolors=(0.0, 0.0, 0.0), linewidths=0.25)
            ax.add_collection(blk)
        xs = [p[0] for p in nodes.values()]
        ys = [p[1] for p in nodes.values()]
        ax.set_xlim(min(xs), max(xs))
        ax.set_ylim(min(ys), max(ys))
        ax.set_aspect('equal')
        ax.set_xlabel('x [mm]', fontsize=11)
        ax.set_ylabel('y [mm]', fontsize=11)
        ax.set_title('Fully damaged elements only', fontsize=12, fontweight='bold')
        ax.tick_params(labelsize=9)
        ax.text(0.015, 0.94,
                '%s\ncracked elems: %d / %d\nmax omega: %.4f\nthreshold: %.2f'
                % (label, idx.size, len(eid_order),
                   (omega.max() if omega.size else 0.0), thresh),
                transform=ax.transAxes, fontsize=8.5, va='top',
                bbox=dict(facecolor='w', edgecolor=(0.6, 0.6, 0.6), lw=0.6, pad=3))
        import matplotlib.patches as mp
        ax.legend(handles=[mp.Patch(facecolor='k', edgecolor='k',
                                    label='omega >= %.2f' % thresh)],
                  loc='upper right', frameon=False, fontsize=9)
        fh.savefig(base + '.png', dpi=600)
        fh.savefig(base + '.pdf', dpi=600)
        plt.close(fh)
        print('  wrote %s.png  (cracked elems: %d, max omega %.4f)'
              % (os.path.basename(base), idx.size,
                 (omega.max() if omega.size else 0.0)))

    one_fig('omega_peak.csv',     'peak load',  os.path.join(out_dir, 'abaqus_fig_damage_peak'))
    one_fig('omega_postpeak.csv', 'post-peak',  os.path.join(out_dir, 'abaqus_fig_damage_postpeak'))
    one_fig('omega_last.csv',     'last step',  os.path.join(out_dir, 'abaqus_fig_damage_last_step'))
    return True

def main():
    odb_path = find_odb(sys.argv[1] if len(sys.argv) > 1 else None)
    thresh = float(sys.argv[2]) if len(sys.argv) > 2 else 0.99
    out_dir = os.path.abspath(sys.argv[3]) if len(sys.argv) > 3 \
        else os.path.join(os.path.dirname(odb_path), 'results')
    pd_dir = os.path.join(out_dir, 'plotdata')
    ensure_dir(out_dir)
    ensure_dir(pd_dir)

    print('ODB:       %s' % odb_path)
    print('Threshold: %.2f' % thresh)
    print('Output:    %s' % out_dir)

    odb = openOdb(path=odb_path, readOnly=True)
    diag = []
    diag.append('ODB: %s' % odb_path)
    diag.append('Threshold (crack = omega >= this): %.2f' % thresh)

    stepname = list(odb.steps.keys())[-1]
    step = odb.steps[stepname]
    nf = len(step.frames)
    diag.append('Step: %s   frames: %d' % (stepname, nf))
    if nf == 0:
        diag.append('FATAL: step has 0 frames -> job wrote no field output.')
        odb.close()
        open(os.path.join(out_dir, 'damage_diagnostic.txt'), 'w').write('\n'.join(diag) + '\n')
        print('\n'.join(diag))
        return

    last = step.frames[-1]
    fkeys = list(last.fieldOutputs.keys())
    diag.append('Field outputs on last frame: %s' % ', '.join(fkeys) if fkeys
                else 'Field outputs on last frame: (NONE)')
    fld, comp, key = get_damage_field(last)
    if fld is None:
        diag.append('FATAL: no SDV/SDV2 field found -> omega cannot be read.')
        diag.append('  Fix: in the run script ensure FIELD_VARS=("SDV",) and the')
        diag.append('  material has Depvar(n>=2); re-run the job.')
    else:
        diag.append('Damage field used: %s (component %d)' % (key, comp))

    gmax, gmax_i = -1.0e99, -1
    for i, fr in enumerate(step.frames):
        mx, _ = frame_max_omega(fr)
        if mx is not None and mx > gmax:
            gmax, gmax_i = mx, i
    if gmax > -1.0e98:
        diag.append('MAX omega anywhere in ODB: %.6f  (frame %d of %d)'
                    % (gmax, gmax_i, nf - 1))
        if gmax < thresh:
            diag.append('==> THIS is why the crack is blank: max omega %.4f < %.2f.'
                        % (gmax, thresh))
            diag.append('    Either lower the threshold, or run further into softening')
            diag.append('    (larger |U_FINAL| / more increments) so damage localizes.')
        else:
            diag.append('==> Damage reaches the threshold: crack SHOULD render.')
            diag.append('    If you still saw no figure, matplotlib was missing in')
            diag.append('    Abaqus python (CSVs are written; use plot_damage_from_csv.py).')
    else:
        diag.append('MAX omega: could not be determined (no readable SDV values).')

    inst_name = dump_mesh(odb, pd_dir)
    diag.append('Instance: %s' % inst_name)
    peak_i, pp_i, last_i = pick_indices(step)
    diag.append('Frame indices -> peak: %s  postpeak: %s  last: %s'
                % (str(peak_i), str(pp_i), str(last_i)))

    for idx, nm in ((peak_i, 'omega_peak.csv'),
                    (pp_i, 'omega_postpeak.csv'),
                    (last_i, 'omega_last.csv')):
        fr = step.frames[idx] if (idx is not None and 0 <= idx < nf) else None
        n, mx = write_omega_csv(fr, os.path.join(pd_dir, nm))
        diag.append('  %-20s rows: %5d   max omega: %.4f' % (nm, n, mx))

    odb.close()

    open(os.path.join(out_dir, 'damage_diagnostic.txt'), 'w').write('\n'.join(diag) + '\n')
    print('\n----- diagnostic -----')
    print('\n'.join(diag))
    print('----------------------\n')

    try_plot(pd_dir, out_dir, thresh)

    print('\nDONE. Read: %s' % os.path.join(out_dir, 'damage_diagnostic.txt'))

if __name__ == '__main__':
    main()
