#!/usr/bin/env python3
"""
DSMC / AMR monitor (single file).

Drop into an OpenFOAM dsmcFoam+ case dir and run:  python3 monitor.py

Buttons: Run (./Allrun) | Stop (SIGINT) | Continue (./Allcont) | Clean (./Allclean)
         | Reconstruct (./Allreconstruct) | .foam | del procs

Right hand side is a notebook with three tabs:

  mass flow   live mass-flow rate + net residual from postProcessing/massFlow.csv
  live view   Mach-number contour map of the newest written time directory,
              redrawn automatically at every new write
  mesh        the mesh itself (faces coloured by AMR refinement level) plus
              cell / level / parcel counters, redrawn at every new write

Case parameters (number density, nEquivalentParticles, cells, deltaT, gas) come
from system/dsmcInitialiseDict + constant/dsmcProperties; wall models from
system/boundariesDict.

Rendering uses matplotlib with the Agg backend in a background thread, so the
GUI never freezes.  If matplotlib is missing the graphs still work and the
image tabs show a message.
"""

import glob
import os
import queue
import re
import signal
import subprocess
import sys
import threading
import time
import tkinter as tk
from tkinter import ttk, messagebox

CASEDIR = os.path.dirname(os.path.abspath(__file__))
POLL_MS = 500                 # GUI refresh
PLOT_STEP_FRAC = 0.001        # redraw mass-flow graph every 0.1 % of the run
CONFIG_FILE = os.path.join(CASEDIR, 'monitor_config.json')
SOLVE_SCRIPTS = ('Allrun', 'Allcont')
IMG_DIR = os.path.join(CASEDIR, '.monitor')
FIELD_IMG = os.path.join(IMG_DIR, 'field.png')
MESH_IMG = os.path.join(IMG_DIR, 'mesh.png')


# ------------------------------------------------------------- readers -----
def read_end_time():
    try:
        s = open(os.path.join(CASEDIR, 'system', 'controlDict')).read()
        m = re.search(r'endTime\s+([-+0-9.eE]+)', s)
        return float(m.group(1)) if m else None
    except Exception:
        return None


def solver_application():
    """Application name from controlDict (dsmcFoam, dsmcFoamPlus, ...)."""
    try:
        s = open(os.path.join(CASEDIR, 'system', 'controlDict')).read()
        m = re.search(r'application\s+(\w+)\s*;', s)
        return m.group(1) if m else 'dsmcFoam'
    except Exception:
        return 'dsmcFoam'


def solver_log_path():
    """log.<application>, falling back to any log.dsmcFoam* (parallel runs write
    the same log through mpirun)."""
    p = os.path.join(CASEDIR, 'log.' + solver_application())
    if os.path.exists(p):
        return p
    cand = sorted(glob.glob(os.path.join(CASEDIR, 'log.dsmcFoam*')))
    return cand[-1] if cand else p


def log_time():
    try:
        s = open(solver_log_path(), errors='replace').read()
        ts = re.findall(r'(?<![\w])Time\s*=\s*([-+0-9.eE]+)', s)
        return float(ts[-1]) if ts else None
    except Exception:
        return None


def time_dirs():
    """[(t, 'name')] for every written time directory, ascending."""
    out = []
    try:
        for d in os.listdir(CASEDIR):
            if not re.match(r'^[0-9]+(\.[0-9]+)?$', d):
                continue
            if os.path.isdir(os.path.join(CASEDIR, d)):
                try:
                    out.append((float(d), d))
                except ValueError:
                    pass
    except OSError:
        pass
    return sorted(out)


def proc_time_dirs():
    """[(t, dir)] of the time directories present in processor0 (parallel run)."""
    p0 = os.path.join(CASEDIR, 'processor0')
    out = []
    if not os.path.isdir(p0):
        return out
    for d in os.listdir(p0):
        if re.match(r'^[0-9]+(\.[0-9]+)?$', d):
            try:
                out.append((float(d), os.path.join(p0, d)))
            except ValueError:
                pass
    return sorted(out)


def current_solution_time():
    """Solver progress: the log is authoritative while a run is writing it
    (a stale, newer time directory from an earlier run must not be shown);
    without a log, fall back to the newest written time directory."""
    t = log_time()
    if t is not None:
        return t
    td = time_dirs()
    return td[-1][0] if td else None


def count_cells(tdir):
    """(nCells, {level: count}) from <time>/cellLevel."""
    f = os.path.join(CASEDIR, tdir, 'cellLevel')
    if not os.path.isfile(f):
        return None, {}
    try:
        s = open(f, errors='replace').read()
    except OSError:
        return None, {}
    m = re.search(r'(?:internalField|value)\s+nonuniform\s+List<[\w:]+>\s*\n(\d+)\s*\n\(', s)
    if not m:
        # no AMR yet: count the base mesh cells (all at level 0)
        n = _owner_cells()
        return (n, {0: n}) if n else (None, {})
    vals = []
    for line in s[m.end():].split('\n'):
        t = line.strip()
        if t == ')':
            break
        try:
            vals.append(int(float(t)))
        except ValueError:
            pass
    if not vals:
        return None, {}
    counts = {}
    for v in vals:
        counts[v] = counts.get(v, 0) + 1
    return len(vals), counts


def _owner_cells():
    """Cell count of the base mesh (constant/polyMesh/owner)."""
    try:
        owner = open(os.path.join(CASEDIR, 'constant', 'polyMesh', 'owner')).read()
        bm = re.search(r'\d+\s*\((.*?)\)\s*//', owner, re.S)
        if bm:
            return max(int(x) for x in bm.group(1).split()) + 1
    except Exception:
        pass
    return None


def log_int(pattern, cast=int):
    f = solver_log_path()
    if not os.path.isfile(f):
        return None
    try:
        s = open(f, errors='replace').read()
    except OSError:
        return None
    m = re.findall(pattern, s)
    return cast(m[-1]) if m else None


def parcels():
    return log_int(r'Number of DSMC particles\s*=\s*(\d+)')


def refinements():
    f = solver_log_path()
    if not os.path.isfile(f):
        return 0
    try:
        s = open(f, errors='replace').read()
    except OSError:
        return 0
    return len(re.findall(r'Selected\s+\d+\s+cells for refinement', s))


def last_selection():
    f = solver_log_path()
    if not os.path.isfile(f):
        return None
    try:
        s = open(f, errors='replace').read()
    except OSError:
        return None
    m = re.findall(r'Selected\s+(\d+)\s+cells for refinement out of (\d+)', s)
    return m[-1] if m else None


def last_fatal():
    """First FOAM FATAL / ERROR line from the newest log, if the run died."""
    for name in ('log.Allrun', 'log.dsmcFoamPlus', 'log.' + solver_application()):
        p = os.path.join(CASEDIR, name)
        if not os.path.isfile(p):
            continue
        try:
            s = open(p, errors='replace').read()
        except OSError:
            continue
        for line in s.split('\n'):
            if 'FOAM FATAL' in line or re.match(r'\s*--> FOAM FATAL', line):
                return line.strip()[:160]
            if 'prun:exe-not-found' in line or 'exe-not-found' in line:
                return 'MPI launcher could not find the solver executable'
    return None


def post_dir(name):
    p = os.path.join(CASEDIR, 'postProcessing', name)
    return p if os.path.isdir(p) else None


def read_config():
    try:
        return json_load(CONFIG_FILE)
    except Exception:
        return {}


def section_planes():
    """x-planes used for the mass-flow integration (monitor_config.json)."""
    c = read_config().get('sections') or [1.0, 2.0]
    try:
        return [float(v) for v in c]
    except Exception:
        return [1.0, 2.0]


def _vec_field_vals(tname, name):
    """vector volField values [(x, y, z), ...] from <time>/name."""
    p = field_path(tname, name)
    if not p:
        return None
    txt = open(p, errors='replace').read()
    out = []
    m = re.search(r'(?:internalField|value)\s+nonuniform\s+List<[\w:]+>\s*\n(\d+)\s*\n\(', txt)
    if m:
        for line in txt[m.end():].split('\n'):
            t = line.strip()
            if t == ')':
                break
            mm = re.match(r'^\(([^)]*)\)$', t)
            if mm:
                parts = mm.group(1).split()
                if len(parts) == 3:
                    try:
                        out.append(tuple(float(x) for x in parts))
                    except ValueError:
                        pass
        return out
    m = re.search(r'(?:internalField|value)\s+uniform\s+\(([^)]*)\)', txt)
    if m:
        parts = m.group(1).split()
        if len(parts) == 3:
            return [tuple(float(x) for x in parts)]
    return None


def _face_sf(pts, fp):
    """Face area vector of one face (3-D polygon, pointer winding)."""
    sx = sy = sz = 0.0
    n = len(fp)
    for i in range(n):
        a = pts[fp[i]]
        b = pts[fp[(i + 1) % n]]
        sx += a[1] * b[2] - a[2] * b[1]
        sy += a[2] * b[0] - a[0] * b[2]
        sz += a[0] * b[1] - a[1] * b[0]
    return (0.5 * sx, 0.5 * sy, 0.5 * sz)


def _cell_geometry(pts, faces, own, nb):
    """Cell centres (x, y), x-extent and volume, approximated from the faces."""
    n = max(own) + 1
    nf = len(faces)
    cx = [0.0] * n
    cy = [0.0] * n
    cnt = [0] * n
    xmin = [1e300] * n
    xmax = [-1e300] * n
    vol = [0.0] * n
    for f in range(nf):
        fp = faces[f]
        sf = _face_sf(pts, fp)
        fx = sum(pts[i][0] for i in fp) / len(fp)
        fy = sum(pts[i][1] for i in fp) / len(fp)
        fz = sum(pts[i][2] for i in fp) / len(fp)
        # owner: outward normal = +Sf; neighbour: outward = -Sf
        c = own[f]
        cx[c] += fx
        cy[c] += fy
        cnt[c] += 1
        for i in fp:
            if pts[i][0] < xmin[c]:
                xmin[c] = pts[i][0]
            if pts[i][0] > xmax[c]:
                xmax[c] = pts[i][0]
        vol[c] += (fx * sf[0] + fy * sf[1] + fz * sf[2]) / 3.0
        c2 = nb[f] if f < len(nb) else -1
        if c2 >= 0:
            cx[c2] += fx
            cy[c2] += fy
            cnt[c2] += 1
            for i in fp:
                if pts[i][0] < xmin[c2]:
                    xmin[c2] = pts[i][0]
                if pts[i][0] > xmax[c2]:
                    xmax[c2] = pts[i][0]
            vol[c2] -= (fx * sf[0] + fy * sf[1] + fz * sf[2]) / 3.0
    for c in range(n):
        if cnt[c]:
            cx[c] /= cnt[c]
            cy[c] /= cnt[c]
        if not vol[c]:
            vol[c] = 1.0
        if xmin[c] > xmax[c]:
            xmin[c] = xmax[c] = 0.0
    return cx, cy, xmin, xmax, vol


def compute_mass_flow(tname):
    """(sections, mdots, net): mass flow through interior x-planes of one time
    directory.

    mdot = integral of rhoM*U_x over the plane.  When the mesh has faces on the
    plane the face fluxes are used (same definition as dsmcMassFlowRate); when
    it has none (a graded mesh need not put a face exactly on x0) the cell
    column containing the plane is integrated instead.  The mesh and fields of
    *that* time directory are used (root or processorN), so AMR changes are
    handled."""
    dirs, src = state_dirs(tname)
    if not dirs:
        return None
    rho = _field_vals_multi(dirs, 'rhoM_mix')
    U = _vec_vals_multi(dirs, 'UMean_mix')
    if not rho or not U:
        return None
    pts, faces, own, nb = _poly_mesh_multi(dirs)
    n = max(own) + 1
    if len(rho) < n or len(U) < n:
        return None
    secs = section_planes()
    cx, cy, xmin, xmax, vol = _cell_geometry(pts, faces, own, nb)
    face_cx = [sum(pts[i][0] for i in fp) / len(fp) for fp in faces]
    mdots = []
    for x0 in secs:
        # Integrate the cell column that contains the plane:
        #   mdot = sum_cells rho * U_x * (V/dx)   over |x_c - x0| <= dx/2
        # The graded/refined meshes do not put a face exactly on x0 (and after
        # AMR only part of the plane may have faces), so a column integral is
        # used for every plane; it is consistent between times.
        m = 0.0
        for c in range(n):
            dxc = xmax[c] - xmin[c]
            if dxc <= 0.0:
                continue
            # exactly the cells the plane passes through (half-open interval:
            # a face lying on x0 is not counted twice)
            if xmin[c] <= x0 < xmax[c]:
                m += rho[c] * U[c][0] * vol[c] / dxc
        mdots.append(m)
    net = (mdots[-1] - mdots[0]) if len(mdots) >= 2 else sum(mdots)
    return secs, mdots, net


def write_mass_flow_csv(history):
    """history: {float(time): (mdots, net)} -> postProcessing/massFlow.csv."""
    if not history:
        return None
    secs = section_planes()
    d = os.path.join(CASEDIR, 'postProcessing')
    try:
        os.makedirs(d, exist_ok=True)
    except OSError:
        return None
    path = os.path.join(d, 'massFlow.csv')
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        f.write('#time ' + ' '.join('x=%g' % s for s in secs) + ' net\n')
        for t in sorted(history):
            mdots, net = history[t]
            f.write('%.6g %s %.6g\n' % (t, ' '.join('%.6g' % v for v in mdots), net))
    os.replace(tmp, path)
    return path


class MassFlowWorker(threading.Thread):
    """Computes the section mass flows of written time directories off the GUI
    thread; results arrive in .done as (time, secs, mdots, net)."""

    def __init__(self):
        super().__init__(daemon=True)
        self.q = queue.Queue()
        self.done = queue.Queue()

    def submit(self, tname):
        self.q.put(tname)

    def run(self):
        while True:
            tname = self.q.get()
            try:
                res = compute_mass_flow(tname)
            except Exception:
                res = None
            self.done.put((tname, res))


def read_mass_flow():
    """(names, rows); names = column labels, rows = (t, v1, v2, net)."""
    names = ('sec1', 'sec2')
    rows = []
    p = os.path.join(CASEDIR, 'postProcessing', 'massFlow.csv')
    if not os.path.exists(p):
        return names, rows
    for line in open(p):
        if line.lstrip().startswith('#'):
            parts = line.split()
            if len(parts) >= 4:
                names = (parts[1], parts[2])
            continue
        parts = line.split()
        if len(parts) < 4:
            continue
        try:
            rows.append(tuple(float(x) for x in parts[:4]))
        except ValueError:
            pass
    rows.sort(key=lambda r: r[0])
    return names, rows


def parse_case_params():
    params = {'numberDensity': '—', 'nEquivalent': '—', 'nCells': '—',
              'deltaT': '—', 'gas': '—'}

    def num(s):
        return '%.5e' % float(s)

    try:
        dp = open(os.path.join(CASEDIR, 'constant', 'dsmcProperties')).read()
    except Exception:
        dp = ''
    try:
        m = re.search(r'nEquivalentParticles\s+([-+0-9.eE]+)', dp)
        if m:
            params['nEquivalent'] = num(m.group(1))
    except Exception:
        pass
    try:
        m = re.search(r'typeIdList\s*\(([^)]*)\)', dp)
        if m:
            params['gas'] = m.group(1).strip()
    except Exception:
        pass
    try:
        m = re.search(r'numberDensities\s*\{([^}]*)\}', dp)
        if m:
            toks = re.findall(r'([A-Za-z0-9_.+-]+)\s+([-+0-9.eE]+)', m.group(1))
            if toks:
                params['numberDensity'] = num(toks[0][1])
    except Exception:
        pass
    try:
        init = open(os.path.join(CASEDIR, 'system', 'dsmcInitialiseDict')).read()
        m = re.search(r'numberDensities\s*\{([^}]*)\}', init)
        if m and params['numberDensity'] == '—':
            toks = re.findall(r'([A-Za-z0-9_.+-]+)\s+([-+0-9.eE]+)', m.group(1))
            if toks:
                params['numberDensity'] = num(toks[0][1])
    except Exception:
        pass
    try:
        cd = open(os.path.join(CASEDIR, 'system', 'controlDict')).read()
        m = re.search(r'deltaT\s+([-+0-9.eE]+)', cd)
        if m:
            params['deltaT'] = m.group(1)
    except Exception:
        pass
    # cells: newest written time dir (AMR) first, then constant/polyMesh
    n = None
    td = time_dirs()
    for _t, name in reversed(td):
        n, _lvl = count_cells(name)
        if n:
            break
    if not n:
        try:
            owner = open(os.path.join(CASEDIR, 'constant', 'polyMesh', 'owner')).read()
            bm = re.search(r'\d+\s*\((.*)\)\s*//', owner, re.S)
            if bm:
                n = max(int(x) for x in bm.group(1).split()) + 1
        except Exception:
            pass
    if n:
        params['nCells'] = str(n)
    return params


def wall_models():
    """(patchName, model-description) from system/boundariesDict."""
    out = []
    try:
        txt = open(os.path.join(CASEDIR, 'system', 'boundariesDict')).read()
    except Exception:
        txt = ''
    if txt:
        for blk in re.finditer(r'patchBoundaryProperties\s*\{(.*?)\}', txt, re.S):
            tail = txt[blk.end():]
            mn = re.search(r'boundaryModel\s+(\w+)\s*;', tail)
            pn = re.search(r"patchName\s+([^\s;]+)\s*;", blk.group(1))
            if not (mn and pn):
                continue
            extra = ''
            if mn.group(1).endswith('MixedDiffuseSpecularWallPatch'):
                dm = re.search(r'diffuseFraction\s+([-+0-9.eE]+)', tail)
                extra = '  diffuseFraction %s' % dm.group(1) if dm else '  diffuseFraction ?'
            out.append('%s: %s%s' % (pn.group(1), mn.group(1), extra))
        for blk in re.finditer(r'generalBoundaryProperties\s*\{(.*?)\}', txt, re.S):
            tail = txt[blk.end():]
            mn = re.search(r'boundaryModel\s+(\w+)\s*;', tail)
            pn = re.search(r"patchName\s+([^\s;]+)\s*;", blk.group(1))
            if mn and pn:
                out.append('%s: %s (inflow)' % (pn.group(1), mn.group(1)))
        if out:
            return out
    try:
        dp = open(os.path.join(CASEDIR, 'constant', 'dsmcProperties')).read()
        m = re.search(r'WallInteractionModel\s+(\S+)\s*;', dp)
        if m:
            model = m.group(1)
            if 'MixedDiffuseSpecular' in model:
                dm = re.search(r'diffuseFraction\s+([-+0-9.eE]+)', dp)
                out.append('all walls: MixedDiffuseSpecular (diffuseFraction %s)'
                           % (dm.group(1) if dm else '?'))
            else:
                out.append('all walls: %s' % model)
    except Exception:
        pass
    return out or ['—']


def load_threshold():
    try:
        return float(json_load(CONFIG_FILE).get('threshold', 0.01))
    except Exception:
        return 0.01


def json_load(path):
    import json
    return json.load(open(path))


def read_stat_summary():
    p = os.path.join(CASEDIR, 'postProcessing', 'statisticalUncertainty', 'summary.dat')
    try:
        return open(p).read()
    except Exception:
        return ''


def relative_residual(row):
    v1, v2, net = abs(row[1]), abs(row[2]), abs(row[3])
    mean_flow = 0.5 * (v1 + v2)
    return net / mean_flow if mean_flow > 1e-300 else 1e9


def is_converged(rows, threshold):
    n = len(rows)
    if n < 2:
        return False
    window = max(20, int(0.05 * n))
    if n < window:
        return False
    return all(relative_residual(r) < threshold for r in rows[-window:])


def newest_mesh_dir(tname=None):
    """Directory that holds a polyMesh for time tname (AMR writes per-time meshes)."""
    cands = []
    if tname:
        cands.append(os.path.join(CASEDIR, tname, 'polyMesh'))
    for _t, name in reversed(time_dirs()):
        cands.append(os.path.join(CASEDIR, name, 'polyMesh'))
    cands.append(os.path.join(CASEDIR, 'constant', 'polyMesh'))
    for c in cands:
        if os.path.isfile(os.path.join(c, 'points')) and os.path.isfile(os.path.join(c, 'owner')):
            return c
    return None


def field_path(tname, name):
    p = os.path.join(CASEDIR, tname, name)
    return p if os.path.isfile(p) else None


# --------------------------------------------------------- rendering -------
class Renderer(threading.Thread):
    """Runs the matplotlib renderer in a *separate process* (monitor.py --render),
    one job at a time.  The GUI thread only loads the finished PNG, so a slow or
    crashing render can never freeze or kill the monitor."""

    def __init__(self):
        super().__init__(daemon=True)
        self.q = queue.Queue()
        self.done = queue.Queue()      # (kind, tdir, path, error)

    def submit(self, kind, tdir):
        self.q.put((kind, tdir))

    def run(self):
        os.makedirs(IMG_DIR, exist_ok=True)
        while True:
            kind, tdir = self.q.get()
            out = FIELD_IMG if kind == 'field' else MESH_IMG
            try:
                res = subprocess.run([sys.executable, os.path.abspath(__file__),
                                      '--render', kind, tdir, out],
                                     cwd=CASEDIR, capture_output=True, text=True,
                                     timeout=600)
                if res.returncode != 0 or not os.path.isfile(out):
                    msg = (res.stderr or res.stdout or 'render failed').strip().split('\n')[-1]
                    self.done.put((kind, tdir, None, msg[:160]))
                else:
                    self.done.put((kind, tdir, out, None))
            except Exception as err:
                self.done.put((kind, tdir, None, '%s: %s' % (type(err).__name__, err)))


def n_processors():
    """Number of processorN directories present in the case."""
    n = 0
    for d in os.listdir(CASEDIR):
        m = re.match(r'^processor(\d+)$', d)
        if m:
            n = max(n, int(m.group(1)) + 1)
    return n


def state_dirs(tname):
    """Directories holding the mesh and fields of time tname.

    After a write the data may only exist in the parallel time directories
    (processorN/<time>), because AMR writes a new mesh and reconstructPar
    cannot merge those meshes (the decomposed meshes have no
    *ProcAddressing).  So fall back to the processor directories."""
    root = os.path.join(CASEDIR, tname)
    if os.path.isdir(root):
        # a case-root time directory: the mesh may still be in constant/polyMesh
        # (static case) and the fields may not exist yet (fresh case)
        return [root], 'root'
    dirs = []
    for n in range(n_processors()):
        p = os.path.join(CASEDIR, 'processor%d' % n, tname)
        if os.path.isdir(p):
            dirs.append(p)
    return (dirs, 'proc') if dirs else ([], None)


def state_times():
    """Time directory names available in the root or in processor0."""
    out = set()
    for _t, name in time_dirs():
        out.add(name)
    p0 = os.path.join(CASEDIR, 'processor0')
    if os.path.isdir(p0):
        for d in os.listdir(p0):
            if re.match(r'^[0-9]+(\.[0-9]+)?$', d) and d != '0':
                out.add(d)
    return sorted(out, key=lambda x: float(x))


def mesh_dir_of(d):
    """polyMesh of a case/processor time directory.

    A write that did not change the mesh has no polyMesh of its own, so fall
    back to the newest earlier time directory of the same processor that has
    one, and finally to the constant mesh."""
    own = os.path.join(d, 'polyMesh')
    if os.path.isfile(os.path.join(own, 'points')):
        return own
    parent = os.path.dirname(d)
    try:
        t_now = float(os.path.basename(d))
    except ValueError:
        t_now = None
    if t_now is not None:
        best, best_t = None, -1.0
        try:
            names = os.listdir(parent)
        except OSError:
            names = []
        for n in names:
            try:
                tv = float(n)
            except ValueError:
                continue
            if tv > t_now:
                continue
            cand = os.path.join(parent, n, 'polyMesh')
            if os.path.isfile(os.path.join(cand, 'points')) and tv > best_t:
                best, best_t = cand, tv
        if best is not None:
            return best
    const = os.path.join(parent, 'constant', 'polyMesh')
    return const if os.path.isfile(os.path.join(const, 'points')) else None


def _poly_mesh_multi(dirs):
    """Concatenate the meshes of several (processor) directories."""
    pts, faces, own, nb = [], [], [], []
    p_off = c_off = 0
    for d in dirs:
        mdir = mesh_dir_of(d)
        if mdir is None:
            raise RuntimeError('no polyMesh for %s' % d)
        p, f, o, n = _poly_mesh(mdir)
        pts.extend(p)
        for fp in f:
            faces.append([i + p_off for i in fp])
        own.extend([c + c_off for c in o])
        nb.extend([c + c_off if c >= 0 else -1 for c in n])
        p_off = len(pts)
        c_off += (max(o) + 1) if o else 0
    return pts, faces, own, nb


def _field_vals_multi(dirs, name):
    """Concatenated values of scalar field <name> over all source directories."""
    out = []
    for d in dirs:
        v = _scalar_field_in_dir(d, name)
        if not v:
            return None
        out.extend(v)
    return out


def _scalar_field_in_dir(d, name):
    return _field_vals_dir(d, name)


def _vec_vals_multi(dirs, name):
    if len(dirs) == 1:
        return _vec_field_vals_dir(dirs[0], name)
    out = []
    for d in dirs:
        v = _vec_field_vals_dir(d, name)
        if not v:
            return None
        out.extend(v)
    return out


def _vec_field_vals_dir(d, name):
    p = os.path.join(d, name)
    if not os.path.isfile(p):
        return None
    txt = open(p, errors='replace').read()
    out = []
    m = re.search(r'(?:internalField|value)\s+nonuniform\s+List<[\w:]+>\s*\n(\d+)\s*\n\(', txt)
    if m:
        for line in txt[m.end():].split('\n'):
            t = line.strip()
            if t == ')':
                break
            mm = re.match(r'^\(([^)]*)\)$', t)
            if mm:
                parts = mm.group(1).split()
                if len(parts) == 3:
                    try:
                        out.append(tuple(float(x) for x in parts))
                    except ValueError:
                        pass
        return out
    m = re.search(r'(?:internalField|value)\s+uniform\s+\(([^)]*)\)', txt)
    if m:
        parts = m.group(1).split()
        if len(parts) == 3:
            return [tuple(float(x) for x in parts)]
    return None


def count_cells_state(tname):
    """(nCells, {level: count}) over the root or processor directories."""
    dirs, src = state_dirs(tname)
    if not dirs:
        return None, {}
    total = {}
    n = 0
    for d in dirs:
        _c, lvl = count_cells_in_dir(d)
        if _c:
            n += _c
            for k, v in lvl.items():
                total[k] = total.get(k, 0) + v
    return (n, total) if n else (None, {})


def count_cells_in_dir(d):
    f = os.path.join(d, 'cellLevel')
    if not os.path.isfile(f):
        return None, {}
    try:
        s = open(f, errors='replace').read()
    except OSError:
        return None, {}
    m = re.search(r'(?:internalField|value)\s+nonuniform\s+List<[\w:]+>\s*\n(\d+)\s*\n\(', s)
    vals = []
    if m:
        for line in s[m.end():].split('\n'):
            t = line.strip()
            if t == ')':
                break
            try:
                vals.append(int(float(t)))
            except ValueError:
                pass
    if not vals:
        mdir = mesh_dir_of(d)
        if mdir:
            try:
                pts, faces, own, _nb = _poly_mesh(mdir)
                n = max(own) + 1
                return n, {0: n}
            except Exception:
                pass
        return None, {}
    counts = {}
    for v in vals:
        counts[v] = counts.get(v, 0) + 1
    return len(vals), counts


def _foam_list(path, as_int=False):
    txt = open(path, errors='replace').read()
    m = re.search(r'\n(\d+)\s*\n\(', txt)
    if not m:
        raise RuntimeError('cannot parse ' + path)
    body = txt[m.end():]
    out = []
    for line in body.split('\n'):
        t = line.strip()
        if t == ')':
            break
        if not t or t.startswith('(') is False and not as_int:
            continue
        if as_int:
            for tok in t.split():
                try:
                    out.append(int(float(tok)))
                except ValueError:
                    pass
        else:
            if t.startswith('('):
                out.append([float(x) for x in t.strip('()').split()])
    return out


def _poly_mesh(mdir):
    """points (N,3), faces (list of lists), owner, neighbour."""
    pts = _foam_list(os.path.join(mdir, 'points'))
    own = _foam_list(os.path.join(mdir, 'owner'), as_int=True)
    txt = open(os.path.join(mdir, 'faces'), errors='replace').read()
    m = re.search(r'\n(\d+)\s*\n\(', txt)
    faces = []
    for line in txt[m.end():].split('\n'):
        s = line.strip()
        mm = re.match(r'^(\d+)\((.*)\)$', s)
        if mm:
            faces.append([int(x) for x in mm.group(2).split()])
        elif s == ')':
            break
    nb = []
    p = os.path.join(mdir, 'neighbour')
    if os.path.isfile(p):
        nb = _foam_list(p, as_int=True)
    if len(nb) < len(faces):
        nb = nb + [-1] * (len(faces) - len(nb))     # boundary faces
    return pts, faces, own, nb


def _field_vals_dir(d, name):
    """scalar field <name> straight from a case or processor directory."""
    p = os.path.join(d, name)
    if not os.path.isfile(p):
        return None
    txt = open(p, errors='replace').read()
    m = re.search(r'(?:internalField|value)\s+nonuniform\s+List<[\w:]+>\s*\n(\d+)\s*\n\(', txt)
    vals = []
    if m:
        for line in txt[m.end():].split('\n'):
            t = line.strip()
            if t == ')':
                break
            try:
                vals.append(float(t))
            except ValueError:
                pass
        return vals
    m = re.search(r'(?:internalField|value)\s+uniform\s+([-+0-9.eE]+)', txt)
    return [float(m.group(1))] if m else None


def _field_vals(tname, name):
    p = field_path(tname, name)
    if not p:
        return None
    txt = open(p, errors='replace').read()
    m = re.search(r'(?:internalField|value)\s+nonuniform\s+List<[\w:]+>\s*\n(\d+)\s*\n\(', txt)
    vals = []
    if m:
        for line in txt[m.end():].split('\n'):
            s = line.strip()
            if s == ')':
                break
            try:
                vals.append(float(s))
            except ValueError:
                pass
        return vals
    m = re.search(r'(?:internalField|value)\s+uniform\s+([-+0-9.eE]+)', txt)
    return [float(m.group(1))] if m else None


def _cell_centres(pts, faces, own, nb):
    n = max(own) + 1
    cx = [0.0] * n
    cy = [0.0] * n
    cnt = [0] * n
    for f, fp in enumerate(faces):
        mp = [sum(pts[i][0] for i in fp) / len(fp),
              sum(pts[i][1] for i in fp) / len(fp)]
        c = own[f]
        cx[c] += mp[0]
        cy[c] += mp[1]
        cnt[c] += 1
    return [cx[i] / cnt[i] if cnt[i] else 0.0 for i in range(n)], \
           [cy[i] / cnt[i] if cnt[i] else 0.0 for i in range(n)]


def load_state(tname):
    """(pts, faces, own, nb, dirs) for one time directory."""
    dirs, _src = state_dirs(tname)
    if not dirs:
        raise RuntimeError('no mesh for t = %s yet (run Allrun first)' % tname)
    pts, faces, own, nb = _poly_mesh_multi(dirs)
    return pts, faces, own, nb, dirs


def render_field(plt, np, tname, out):
    """Mach-number contour of one written time directory (horizontal colorbar)."""
    pts, faces, own, nb, dirs = load_state(tname)
    x, y = _cell_centres(pts, faces, own, nb)
    vals = _field_vals_multi(dirs, 'Ma_mix')
    if not vals:
        raise RuntimeError('no Ma_mix field in %s' % tname)
    n = min(len(x), len(vals))
    xs, ys, vs = np.array(x[:n]), np.array(y[:n]), np.array(vals[:n])
    fig, ax = plt.subplots(figsize=(7.2, 4.4), dpi=100)
    try:
        art = ax.tricontourf(xs, ys, vs, levels=40, cmap='turbo')
    except Exception:
        art = ax.scatter(xs, ys, c=vs, s=4, cmap='turbo', marker='s', linewidths=0)
    ax.set_aspect('equal')
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_title('Ma at t = %s s   (%d cells, max %.2f)' % (tname, n, float(vs.max())))
    cb = fig.colorbar(art, ax=ax, orientation='horizontal',
                      fraction=0.05, pad=0.14, aspect=40)
    cb.set_label('Ma')
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def render_mesh(plt, np, LineCollection, tname, out):
    """Black face wireframe of one written time directory."""
    pts, faces, own, nb, dirs = load_state(tname)
    n = max(own) + 1
    maxlev = 0
    lv = _field_vals_multi(dirs, 'cellLevel')
    if lv:
        try:
            maxlev = max(int(v) for v in lv)
        except Exception:
            maxlev = 0
    segs = []
    for fp in faces:
        ps = [pts[i] for i in fp]
        for a in range(len(ps)):
            b = (a + 1) % len(ps)
            if abs(ps[a][2] - ps[b][2]) > 1e-12:
                continue                    # skip through-thickness (z) edges
            segs.append([(ps[a][0], ps[a][1]), (ps[b][0], ps[b][1])])
    fig, ax = plt.subplots(figsize=(7.2, 4.4), dpi=100)
    ax.add_collection(LineCollection(segs, linewidths=0.25, colors='black'))
    ax.autoscale()
    ax.set_aspect('equal')
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_title('mesh at t = %s s   (%d cells, %d faces, AMR level 0..%d)'
                 % (tname, n, len(faces), maxlev))
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


# ---------------------------------------------------------------- plot -----
def fmt(v):
    return '0' if v == 0 else '%.3g' % v


class MiniPlot:
    def __init__(self, parent, height=220):
        self.cv = tk.Canvas(parent, height=height, bg='white',
                            highlightthickness=1, highlightbackground='#bbb')
        self.cv.pack(fill='both', expand=True, padx=8, pady=2)

    def plot(self, title, series, xmax_hint):
        cv = self.cv
        cv.delete('all')
        W = max(cv.winfo_width(), 300)
        H = max(cv.winfo_height(), 150)
        pad_l, pad_r, pad_t, pad_b = 64, 14, 26, 30
        pw, ph = W - pad_l - pad_r, H - pad_t - pad_b
        xs = [p[0] for _n, _c, pts in series for p in pts]
        ys = [p[1] for _n, _c, pts in series for p in pts]
        if not xs:
            cv.create_text(W // 2, H // 2, text='waiting for postProcessing/massFlow.csv',
                           fill='#666', font=('Helvetica', 11))
            return
        xmax = max(max(xs), xmax_hint or 0.0) * 1.02 or 1e-9
        ymax, ymin = max(ys), min(ys)
        span = (ymax - ymin) or 1e-9
        ymax += 0.08 * span
        ymin -= 0.08 * span

        def X(x):
            return pad_l + x / xmax * pw

        def Y(y):
            return pad_t + ph - (y - ymin) / (ymax - ymin) * ph

        for i in range(5):
            yv = ymin + (ymax - ymin) * i / 4.0
            cv.create_line(pad_l, Y(yv), W - pad_r, Y(yv), fill='#eee')
            cv.create_text(pad_l - 6, Y(yv), text=fmt(yv), anchor='e',
                           font=('Helvetica', 8), fill='#444')
        for i in range(5):
            xv = xmax * i / 4.0
            cv.create_text(X(xv), H - pad_b + 12, text=fmt(xv), anchor='n',
                           font=('Helvetica', 8), fill='#444')
        for name, color, pts in series:
            if len(pts) < 2:
                continue
            coords = []
            for x, y in pts:
                coords += [X(x), Y(y)]
            cv.create_line(*coords, fill=color, width=2)
        cv.create_text(pad_l + 8, pad_t - 12, text=title, anchor='w',
                       font=('Helvetica', 10, 'bold'), fill='#111')
        ly = pad_t + 14
        for name, color, _pts in series:
            cv.create_line(W - pad_r - 40, ly, W - pad_r - 20, ly, fill=color, width=2)
            cv.create_text(W - pad_r - 16, ly, text=name, anchor='w',
                           font=('Helvetica', 8), fill='#111')
            ly += 16


# ------------------------------------------------------------------ app ----
class ImageTab(ttk.Frame):
    """One image tab: caption + scaled picture + its own Fit / 100 % / Refresh
    controls (so every graph can be refreshed and re-fitted to the window)."""

    def __init__(self, parent, app, kind, caption):
        super().__init__(parent)
        self.app = app
        self.kind = kind
        self.caption_var = caption
        self.pil = None          # PIL image of the rendered PNG
        self.path = None
        self.fit = True
        self.tkimg = None
        self._after = None

        bar = ttk.Frame(self)
        bar.pack(fill='x', padx=6, pady=(4, 0))
        ttk.Button(bar, text='Fit to window',
                   command=lambda: self.set_fit(True)).pack(side='left', padx=2)
        ttk.Button(bar, text='100 %',
                   command=lambda: self.set_fit(False)).pack(side='left', padx=2)
        ttk.Button(bar, text='Refresh', command=self.refresh).pack(side='left', padx=2)
        self.zoom_var = tk.StringVar(value='fit')
        ttk.Label(bar, textvariable=self.zoom_var).pack(side='left', padx=8)

        ttk.Label(self, textvariable=self.caption_var, anchor='w').pack(fill='x', padx=6)
        self.label = tk.Label(self, bg='white', anchor='center', width=1, height=1)
        self.label.pack(fill='both', expand=True, padx=6, pady=4)
        self.label.bind('<Configure>', self._on_configure)

    # ---- sizing ----
    def set_fit(self, fit):
        self.fit = fit
        self.zoom_var.set('fit' if fit else '100 %')
        self._show_scaled()

    def refresh(self):
        self.app.rerender(self.kind)

    def _on_configure(self, _event=None):
        if not self.fit or self.path is None:
            return
        if self._after is not None:
            try:
                self.after_cancel(self._after)
            except Exception:
                pass
        self._after = self.after(250, self._show_scaled)

    def load(self, path):
        self.path = path
        try:
            from PIL import Image
            self.pil = Image.open(path)
        except Exception:
            self.pil = None
        self._show_scaled()

    def _target_size(self, iw, ih):
        if not self.fit:
            return iw, ih
        W = max(self.label.winfo_width(), 120)
        H = max(self.label.winfo_height(), 80)
        if W <= 1 or H <= 1:
            return None
        scale = min(W / iw, H / ih)
        return max(1, int(iw * scale)), max(1, int(ih * scale))

    def _show_scaled(self):
        if self.path is None:
            return
        if self.pil is None:
            # no PIL: show the rendered image as it is
            try:
                self.tkimg = tk.PhotoImage(file=self.path)
                self.label.configure(image=self.tkimg)
                self.label.image = self.tkimg
            except Exception:
                pass
            return
        try:
            iw, ih = self.pil.size
            size = self._target_size(iw, ih)
            tmp = os.path.join(IMG_DIR, '%s_view.png' % self.kind)
            if size is None or size == (iw, ih):
                self.tkimg = tk.PhotoImage(file=self.path)
            else:
                img = self.pil.resize(size)
                img.save(tmp)
                self.tkimg = tk.PhotoImage(file=tmp)
            self.label.configure(image=self.tkimg)
            self.label.image = self.tkimg
            self.zoom_var.set('%dx%d (%s)' % (size[0], size[1], 'fit' if self.fit else '100 %')
                              if size else '')
        except Exception as err:
            self.caption_var.set('display error: %s' % err)



class MonitorApp:
    def __init__(self, root):
        self.root = root
        self.end_time = read_end_time()
        self.proc = None
        self.script = None
        self._stopped_by_user = False
        self._was_proc_running = None
        self._last_plot_t = None
        self._force_plot = False
        self._logfile = None
        self.threshold = load_threshold()

        # rendering state
        self.renderer = Renderer()
        self.renderer.start()
        self.mf_worker = MassFlowWorker()
        self.mf_worker.start()
        self.mf = {}            # {time: (mdots, net)}
        self._mf_done = set()   # times already handed to the worker
        self._mf_tries = {}     # attempts per time (retries incomplete fields)
        self._shown = {'field': None, 'mesh': None}   # time dir last displayed
        self._want = {'field': None, 'mesh': None}    # time dir last submitted

        root.title('DSMC monitor — ' + os.path.basename(os.path.normpath(CASEDIR)))
        root.geometry('1180x760')

        # ---- buttons ----
        bar = ttk.Frame(root)
        bar.pack(fill='x', padx=8, pady=6)
        self.btn_run = ttk.Button(bar, text='Run', command=self.run)
        self.btn_stop = ttk.Button(bar, text='Stop', command=self.stop, state='disabled')
        self.btn_cont = ttk.Button(bar, text='Continue', command=self.cont)
        self.btn_clean = ttk.Button(bar, text='Clean', command=self.clean)
        self.btn_recon = ttk.Button(bar, text='Reconstruct', command=self.reconstruct)
        self.btn_foam = ttk.Button(bar, text='.foam', command=self.make_foam)
        self.btn_proc = ttk.Button(bar, text='del procs', command=self.clean_procs)
        for b in (self.btn_run, self.btn_stop, self.btn_cont, self.btn_clean,
                  self.btn_recon, self.btn_foam, self.btn_proc):
            b.pack(side='left', padx=3)
        self.status = tk.StringVar(value='idle')
        self.status_lbl = ttk.Label(bar, textvariable=self.status)
        self.status_lbl.pack(side='right', padx=6)

        # ---- progress ----
        pframe = ttk.Frame(root)
        pframe.pack(fill='x', padx=8, pady=(0, 4))
        self.progress = ttk.Progressbar(pframe, maximum=100.0)
        self.progress.pack(fill='x')
        self.progress_label = tk.StringVar(value='progress: —')
        ttk.Label(pframe, textvariable=self.progress_label).pack(anchor='w')

        # ---- left panel ----
        left = ttk.LabelFrame(root, text='Case parameters')
        left.pack(side='left', fill='y', padx=8, pady=4)
        self.param_vars = {}
        rows = [('number density (1/m3)', 'numberDensity'),
                ('nEquivalentParticles', 'nEquivalent'),
                ('number of cells', 'nCells'),
                ('time step (s)', 'deltaT'),
                ('gas', 'gas')]
        for label, key in rows:
            r = ttk.Frame(left)
            r.pack(fill='x', padx=6, pady=2)
            ttk.Label(r, text=label + ':', width=22, anchor='w').pack(side='left')
            v = tk.StringVar(value='—')
            self.param_vars[key] = v
            ttk.Label(r, textvariable=v, width=16, anchor='e').pack(side='right')
        wf = ttk.LabelFrame(left, text='Boundaries (system/boundariesDict)')
        wf.pack(fill='x', padx=6, pady=(8, 2))
        self.walls_var = tk.StringVar(value='—')
        ttk.Label(wf, textvariable=self.walls_var, anchor='w', justify='left',
                  font=('Menlo', 7)).pack(fill='x', padx=4, pady=2)
        self.update_params()

        fr = ttk.Frame(left)
        fr.pack(fill='x', padx=6, pady=(12, 2))
        ttk.Label(fr, text='flow time (s):', width=22, anchor='w').pack(side='left')
        self.flow_time_var = tk.StringVar(value='—')
        ttk.Label(fr, textvariable=self.flow_time_var, width=16, anchor='e').pack(side='right')

        th = ttk.Frame(left)
        th.pack(fill='x', padx=6, pady=(14, 2))
        ttk.Label(th, text='Convergence threshold (frac):').pack(side='left')
        self.thresh_var = tk.StringVar(value=str(self.threshold))
        ttk.Entry(th, textvariable=self.thresh_var, width=8).pack(side='left', padx=4)
        ttk.Button(th, text='Set', command=self.set_threshold).pack(side='left')
        rfb = ttk.Frame(left)
        rfb.pack(anchor='w', padx=6, pady=(2, 0))
        ttk.Button(rfb, text='Refresh parameters', command=self.update_params).pack(side='left')
        ttk.Button(rfb, text='Refresh all', command=self.refresh_all).pack(side='left', padx=4)
        ttk.Button(rfb, text='Reset axes', command=self.reset_axes).pack(side='left')
        self.residual_label = tk.StringVar(value='residual: —')
        ttk.Label(left, textvariable=self.residual_label).pack(anchor='w', padx=6, pady=(4, 0))

        statframe = ttk.LabelFrame(left, text='Statistical convergence (dsmcStatistics)')
        statframe.pack(fill='both', expand=True, padx=6, pady=(10, 4))
        self.statbox = tk.Text(statframe, height=9, font=('Menlo', 7), state='disabled',
                               bg='#f7f7f7', wrap='none')
        self.statbox.pack(fill='both', expand=True, padx=3, pady=3)

        # ---- right: notebook with mass flow / live view / mesh ----
        right = ttk.Frame(root)
        right.pack(side='left', fill='both', expand=True, padx=4, pady=4)
        nb = ttk.Notebook(right)
        nb.pack(fill='both', expand=True)

        tab_plot = ttk.Frame(nb)
        nb.add(tab_plot, text='mass flow')
        self.plot_m = MiniPlot(tab_plot)
        self.plot_net = MiniPlot(tab_plot)

        tab_field = ttk.Frame(nb)
        nb.add(tab_field, text='live Ma')
        self.field_caption = tk.StringVar(value='waiting for the first write with Ma_mix')
        self.field_tab = ImageTab(tab_field, self, 'field', self.field_caption)
        self.field_tab.pack(fill='both', expand=True)
        self.field_label = self.field_tab.label

        tab_mesh = ttk.Frame(nb)
        nb.add(tab_mesh, text='mesh')
        self.mesh_caption = tk.StringVar(value='waiting for the first written time directory')
        self.mesh_tab = ImageTab(tab_mesh, self, 'mesh', self.mesh_caption)
        self.mesh_tab.pack(fill='both', expand=True)
        self.mesh_label = self.mesh_tab.label

        self.nb = nb
        self.refresh()

    # --------------------------------------------------------- actions -----
    def _launch(self, script):
        if self.proc and self.proc.poll() is None:
            return
        path = os.path.join(CASEDIR, script)
        if not os.path.isfile(path):
            self.status.set('%s not found in %s' % (script, CASEDIR))
            return
        # dsmcFoam+ cases have no 0/ directory until dsmcInitialisePlus creates
        # it (fields, parcels, parcel dictionaries); Allrun does that itself, so
        # a missing 0/ is not a reason to refuse to start.

        # the case scripts delete log.* themselves, so keep the monitor's own
        # log of the launched script under .monitor/
        try:
            os.makedirs(IMG_DIR, exist_ok=True)
        except OSError:
            pass
        logpath = os.path.join(IMG_DIR, script + '.log')
        try:
            logf = open(logpath, 'w')
        except OSError as err:
            self.status.set('cannot open %s: %s' % (logpath, err))
            return

        env = dict(os.environ)
        env.setdefault('HOME', os.path.expanduser('~') or '/home/ubuntu')
        env.setdefault('USER', os.environ.get('USER', 'ubuntu'))
        try:
            self.proc = subprocess.Popen(['bash', path], cwd=CASEDIR, env=env,
                                         start_new_session=True,
                                         stdout=logf, stderr=subprocess.STDOUT)
        except Exception as err:
            self.status.set('launch failed: %s' % err)
            logf.close()
            return
        finally:
            self._logfile = logf

        self.script = script
        self._stopped_by_user = False
        self.status.set('running %s (pid %d) -> .monitor/%s.log'
                        % (script, self.proc.pid, script))

    def reconstruct(self):
        self._launch('Allreconstruct')

    def make_foam(self):
        name = os.path.basename(os.path.normpath(CASEDIR)) + '.foam'
        try:
            open(os.path.join(CASEDIR, name), 'w').close()
            self.status.set('created %s' % name)
        except Exception as err:
            self.status.set('touch failed: %s' % err)

    def clean_procs(self):
        import shutil
        n = 0
        for d in glob.glob(os.path.join(CASEDIR, 'processor*')):
            if os.path.isdir(d):
                shutil.rmtree(d, ignore_errors=True)
                n += 1
        self.status.set('removed %d processor dir(s)' % n)

    def run(self):
        self._launch('Allrun')

    def cont(self):
        self._launch('Allcont')

    def clean(self):
        self._launch('Allclean')

    def stop(self):
        if self.proc and self.proc.poll() is None:
            try:
                os.killpg(os.getpgid(self.proc.pid), signal.SIGINT)
            except Exception:
                self.proc.terminate()
            self._stopped_by_user = True
            self.status.set('sent Ctrl+C (SIGINT)')

    def set_threshold(self):
        try:
            v = float(self.thresh_var.get())
            if v <= 0:
                raise ValueError
            self.threshold = v
            cfg = read_config()
            cfg['threshold'] = v
            import json
            json.dump(cfg, open(CONFIG_FILE, 'w'), indent=2)
            self.thresh_var.set(str(v))
            self.status.set('threshold set: %g' % v)
        except Exception:
            messagebox.showerror('Set threshold', 'Enter a positive number, e.g. 0.01')

    def refresh_plot(self):
        self._force_plot = True
        self.refresh()
        self.status.set('plots refreshed')

    def update_params(self):
        self.end_time = read_end_time()
        p = parse_case_params()
        for k, v in p.items():
            self.param_vars[k].set(v)
        self.walls_var.set('\n'.join(wall_models()))

    # -------------------------------------------------------- rendering ----
    def _kick_renders(self, latest):
        """Submit field/mesh renders for a newly written time directory.

        Times that only exist in the parallel directories are rendered too
        (they cannot be reconstructed once AMR has changed the mesh)."""
        field_dir = None
        for name in reversed(state_times()):
            dirs, _src = state_dirs(name)
            if not dirs:
                continue
            if os.path.isfile(os.path.join(dirs[0], 'Ma_mix')):
                field_dir = name
                break
        jobs = (('field', field_dir), ('mesh', latest))
        for kind, tdir in jobs:
            if tdir is None:
                if kind == 'field':
                    self.field_caption.set('waiting for the first write with Ma_mix')
                continue
            if self._want[kind] != tdir:
                self._want[kind] = tdir
                self.renderer.submit(kind, tdir)

    def _collect_images(self):
        while True:
            try:
                kind, tdir, path, err = self.renderer.done.get_nowait()
            except queue.Empty:
                break
            tab = self.field_tab if kind == 'field' else self.mesh_tab
            if path is None:
                # never keep a failed render "done": allow a retry next poll
                self._want[kind] = None
                tab.caption_var.set('t = %s s: %s' % (tdir, err))
                continue
            if kind == 'field':
                tab.load(path)
                tab.caption_var.set('Mach number at t = %s s' % tdir)
            else:
                tab.load(path)
                n, lvl = count_cells_state(tdir)
                txt = 'mesh at t = %s s' % tdir
                if n:
                    txt += '   %d cells (%s)' % (n, ', '.join(
                        'L%d: %d' % (k, lvl[k]) for k in sorted(lvl)))
                p = parcels()
                if p:
                    txt += '   %d parcels' % p
                sel = last_selection()
                if sel:
                    txt += '   last AMR: %s of %s cells' % sel
                tab.caption_var.set(txt)
            self._shown[kind] = tdir

    # ---- per-tab / global refresh ----
    def rerender(self, kind):
        """Force a fresh render of one image tab at the newest time."""
        self._want[kind] = None
        self.refresh()

    def refresh_all(self):
        """Refresh every graph: mass flow, live Ma and mesh."""
        self._want['field'] = None
        self._want['mesh'] = None
        self._force_plot = True
        self.refresh()
        self.status.set('refreshed all graphs')

    def reset_axes(self):
        """Rescale the mass-flow axes and re-fit both image tabs to the window."""
        self._force_plot = True
        self.field_tab.set_fit(True)
        self.mesh_tab.set_fit(True)
        self.refresh()

    # ---------------------------------------------------------- refresh ----
    def refresh(self):
        self.end_time = read_end_time()
        proc_running = self.proc is not None and self.proc.poll() is None
        t = current_solution_time()

        self.btn_run.config(state='disabled' if proc_running else 'normal')
        self.btn_stop.config(state='normal' if proc_running else 'disabled')
        self.btn_cont.config(state='disabled' if proc_running else 'normal')
        self.btn_clean.config(state='disabled' if proc_running else 'normal')

        if t is not None:
            self.flow_time_var.set('%.6g' % t)
        if t is not None and self.end_time:
            frac = min(1.0, max(0.0, t / self.end_time))
            self.progress['value'] = 100.0 * frac
            self.progress_label.set('t = %.6g / %.6g s  (%.2f%%)' % (t, self.end_time, 100.0 * frac))
        else:
            self.progress['value'] = 0.0
            self.progress_label.set('progress: waiting for log.%s' % solver_application())

        # ---- mass flow through the section planes -------------------------
        # Computed here from the written time directories (rhoM_mix, UMean_mix
        # and the mesh of *that* time), because the dsmcMassFlowRate function
        # object caches its section faces and cannot follow an AMR mesh.
        # The result is written to postProcessing/massFlow.csv in the usual
        # format, so the same file the other cases produce exists here too.
        mf_write = read_config().get('massflow_write', True)
        if mf_write:
            for tname in state_times():
                if tname != '0' and tname not in self._mf_done:
                    self._mf_done.add(tname)
                    self.mf_worker.submit(tname)
        changed = False
        while True:
            try:
                tname, res = self.mf_worker.done.get_nowait()
            except queue.Empty:
                break
            if res is None:
                # fields not complete yet (reconstruct in progress): retry a few
                n = self._mf_tries.get(tname, 0) + 1
                self._mf_tries[tname] = n
                if n < 5:
                    self._mf_done.discard(tname)
                continue
            _secs, mdots, net = res
            self.mf[float(tname)] = (mdots, net)
            changed = True
        if changed:
            write_mass_flow_csv(self.mf)
            self._force_plot = True

        names, rows = read_mass_flow()
        if self.param_vars['nCells'].get() == '—' and \
                os.path.exists(os.path.join(CASEDIR, '0')):
            self.update_params()

        s = read_stat_summary()
        self.statbox.config(state='normal')
        self.statbox.delete('1.0', 'end')
        self.statbox.insert('1.0', s if s else 'waiting for dsmcStatistics summary.dat')
        self.statbox.config(state='disabled')

        conv = is_converged(rows, self.threshold)
        if rows:
            self.residual_label.set('residual: %.4g  (threshold %.3g)'
                                    % (relative_residual(rows[-1]), self.threshold))

        redraw = getattr(self, '_force_plot', False)
        self._force_plot = False
        if self.end_time:
            step = PLOT_STEP_FRAC * self.end_time
            if self._last_plot_t is None:
                redraw = True
            elif t is not None and (t - self._last_plot_t) >= step - 1e-12:
                redraw = True
        if redraw:
            if t is not None:
                self._last_plot_t = t
            inlet = [(r[0], r[1]) for r in rows]
            outlet = [(r[0], r[2]) for r in rows]
            net = [(r[0], r[3]) for r in rows]
            self.plot_m.plot('mass flow rate (kg/s)',
                             [(names[0], 'royalblue', inlet),
                              (names[1], 'forestgreen', outlet)], self.end_time)
            self.plot_net.plot('net residual (kg/s)', [('net', 'crimson', net)], self.end_time)

        # live view / mesh: render whenever a new time directory appears
        td_now = state_times()
        latest = td_now[-1] if td_now else None
        self._kick_renders(latest)
        self._collect_images()

        was_proc = self._was_proc_running
        finished = (was_proc and not proc_running and
                    self.script in SOLVE_SCRIPTS and not self._stopped_by_user)
        if finished:
            msg = 'simulation finished'
            fatal = last_fatal()
            if fatal:
                msg = 'run failed: %s' % fatal
            self.status.set(msg)
        elif was_proc and not proc_running and self.script in SOLVE_SCRIPTS:
            self.status.set('stopped by user')
        self._was_proc_running = proc_running

        if conv and proc_running:
            self.status.set('CONVERGED (residual < %g)' % self.threshold)

        color = 'green' if (conv or finished) else \
                ('red' if self._stopped_by_user else
                 ('dark orange' if proc_running else 'black'))
        self.status_lbl.config(foreground=color)

        self.root.after(POLL_MS, self.refresh)


def main():
    if len(sys.argv) > 1 and sys.argv[1] == '--render':
        # internal helper for the background renderer:
        #   monitor.py --render field|mesh <timeDir> <out.png>
        kind, tname, out = sys.argv[2], sys.argv[3], sys.argv[4]
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import numpy as np
        from matplotlib.collections import LineCollection
        if kind == 'field':
            render_field(plt, np, tname, out)
        else:
            render_mesh(plt, np, LineCollection, tname, out)
        return

    root = tk.Tk()
    MonitorApp(root)
    root.mainloop()


if __name__ == '__main__':
    main()
