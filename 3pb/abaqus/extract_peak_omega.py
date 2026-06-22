from __future__ import print_function
import os, sys, glob

try:
    from odbAccess import openOdb
except Exception as e:
    print('ERROR: run inside Abaqus python:')
    print('  abaqus python extract_peak_omega.py')
    sys.exit(1)

script_dir = os.path.abspath(os.path.dirname(sys.argv[0]) if sys.argv[0] else os.getcwd())
odb_path   = None
for cand in glob.glob(os.path.join(script_dir, '*.odb')):
    odb_path = cand
if odb_path is None:
    print('ERROR: no .odb found in %s' % script_dir)
    sys.exit(1)
print('ODB: %s' % odb_path)

res_dir = os.path.join(script_dir, 'results')
pd_dir  = os.path.join(res_dir, 'plotdata')
for d in (res_dir, pd_dir):
    if not os.path.isdir(d):
        os.makedirs(d)

lc_csv = os.path.join(res_dir, 'abaqus_load_cmod.csv')
peak_load, peak_time_approx = -1e99, None
if os.path.exists(lc_csv):
    rows = []
    with open(lc_csv, 'r') as f:
        for ln in f:
            s = ln.strip()
            if not s or s.startswith('#'):
                continue
            p = s.replace(',', ' ').split()
            try:
                rows.append((float(p[0]), float(p[1])))
            except Exception:
                continue
    if rows:

        best = max(rows, key=lambda r: r[1])
        peak_cmod, peak_load = best
        peak_idx  = rows.index(best)
        print('Peak load from CSV: %.4f N  at CMOD=%.6f mm  (row %d)' %
              (peak_load, peak_cmod, peak_idx))

        pp_cmod = 0.3
        pp_idx  = peak_idx
        for i, (c, r) in enumerate(rows):
            if c > pp_cmod:
                pp_idx = i
                break
        print('Post-peak row %d: CMOD=%.4f mm' % (pp_idx, rows[pp_idx][0]))
else:
    print('WARNING: %s not found; will use last frame for peak.' % lc_csv)

odb  = openOdb(path=odb_path, readOnly=True)
step = odb.steps[list(odb.steps.keys())[-1]]
nf   = len(step.frames)
print('Step: %s   frames: %d' % (list(odb.steps.keys())[-1], nf))

frame_times = [float(step.frames[i].frameValue) for i in range(nf)]

def nearest_frame(target_frac):

    t_max = frame_times[-1] if frame_times else 1.0
    target_t = target_frac * t_max
    best_i, best_d = 0, 1e99
    for i, t in enumerate(frame_times):
        d = abs(t - target_t)
        if d < best_d:
            best_d, best_i = d, i
    return best_i

if rows:
    n_rows   = len(rows)
    pk_frac  = float(peak_idx) / max(n_rows - 1, 1)
    pp_frac  = float(pp_idx)   / max(n_rows - 1, 1)
    peak_fi  = nearest_frame(pk_frac)
    pp_fi    = nearest_frame(pp_frac)
else:
    peak_fi  = nf - 1
    pp_fi    = nf - 1
last_fi = nf - 1

print('Field-frame indices → peak: %d (t=%.4f)  postpeak: %d (t=%.4f)  last: %d (t=%.4f)'
      % (peak_fi, frame_times[peak_fi],
         pp_fi,   frame_times[pp_fi],
         last_fi, frame_times[last_fi]))

def get_damage_field(frame):
    keys = list(frame.fieldOutputs.keys())
    for k in ('SDV2', 'SDV_2', 'STATEV2'):
        if k in keys:
            return frame.fieldOutputs[k], 0
    for k in ('SDV', 'STATEV'):
        if k in keys:
            return frame.fieldOutputs[k], 1
    return None, None

def write_omega(frame, path):
    fld, comp = get_damage_field(frame)
    n, mx = 0, -1.0
    with open(path, 'w') as f:
        f.write('# element_id, ip, omega\n')
        if fld is None:
            print('  WARNING: no SDV field in frame %d' % frame.incrementNumber)
            return 0, 0.0
        for v in fld.values:
            try:
                data = v.data
                try:
                    w = float(data[comp])
                except Exception:
                    w = float(data)
            except Exception:
                continue
            ip = getattr(v, 'integrationPoint', 1)
            f.write('%s, %s, %.6e\n' % (str(v.elementLabel), str(ip), w))
            n += 1
            if w > mx:
                mx = w
    return n, mx

for fi, name in ((peak_fi, 'omega_peak.csv'),
                 (pp_fi,   'omega_postpeak.csv'),
                 (last_fi, 'omega_last.csv')):
    fr = step.frames[fi]
    n, mx = write_omega(fr, os.path.join(pd_dir, name))
    print('  %-22s  rows: %5d   max omega: %.6f   (frame %d, t=%.4f)'
          % (name, n, mx, fi, frame_times[fi]))

odb.close()
print('\nDone. Now run:')
print('  python plot_results.py')
