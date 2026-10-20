#!/usr/bin/env python3
"""
DSMC convergence monitor (portable, single file, stdlib only).

Drop into an OpenFOAM dsmcFoam case dir and run:  python3 monitor.py

Buttons: Run (./Allrun) | Stop (SIGINT) | Continue (./Allcont) | Clean (./Allclean)
         | Reconstruct (./Allreconstruct, every written time) | .foam (touch case.foam)
         | del procs (rm -rf processor*)

Works as-is for the dsmcFoam+ cases: case parameters are read from
system/dsmcInitialiseDict and the boundary models are listed patch-by-patch
from system/boundariesDict (diffuse / specular / mixed + diffuseFraction /
deletion / inflow).

Features:
  * Live mass-flow plot from postProcessing/massFlow.csv (written every step by
    the dsmcMassFlowRate FO), redrawn every 0.1% of solution progress.
  * Progress bar (log.<application> "Time =" vs controlDict endTime); the
    application comes from controlDict, so it works for dsmcFoam (v2606) and
    dsmcFoamPlus (v2406) alike, serial or mpirun-parallel.
  * Convergence alert: set a threshold (% residual) in the box + "Set"; when the
    sustained |sec2-sec1|/mean < threshold, an alert shows (simulation keeps runng).
  * Finish alert when the solve subprocess stops for any reason.
  * Case parameters panel (number density, nEq, cells, deltaT, gas, diffusion frac).
  * Color status (amber running / green converged / red stopped).
"""

import os, re, json, signal, subprocess
import tkinter as tk
from tkinter import ttk, messagebox

CASEDIR = os.path.dirname(os.path.abspath(__file__))
POLL_MS = 200
PLOT_STEP_FRAC = 0.001
CONFIG_FILE = os.path.join(CASEDIR, 'monitor_config.json')
SOLVE_SCRIPTS = ('Allrun', 'Allcont')


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
    """log.<application>, falling back to any log.dsmcFoam* (also handles
    parallel runs, where the same log file is written by mpirun)."""
    app = solver_application()
    p = os.path.join(CASEDIR, 'log.' + app)
    if os.path.exists(p):
        return p
    import glob
    cand = sorted(glob.glob(os.path.join(CASEDIR, 'log.dsmcFoam*')))
    return cand[-1] if cand else p


def current_solution_time():
    try:
        s = open(solver_log_path(), errors='replace').read()
        ts = re.findall(r'(?<![\w])Time\s*=\s*([-+0-9.eE]+)', s)
        return float(ts[-1]) if ts else None
    except Exception:
        return None


def read_mass_flow():
    """(names, rows); names = the two measured column labels, rows = (t,v1,v2,net)."""
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
    rows.sort(key=lambda r: r[0])      # rows may interleave after a restart
    return names, rows


def parse_case_params():
    """Extract the case parameters shown in the side panel (best-effort)."""
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
        if m: params['nEquivalent'] = num(m.group(1))
    except Exception: pass
    try:
        m = re.search(r'typeIdList\s*\(([^)]*)\)', dp)
        if m: params['gas'] = m.group(1).strip()
    except Exception: pass
    try:
        m = re.search(r'numberDensities\s*\{([^}]*)\}', dp)
        if m:
            toks = re.findall(r'([A-Za-z0-9_.+-]+)\s+([-+0-9.eE]+)', m.group(1))
            if toks: params['numberDensity'] = num(toks[0][1])
    except Exception: pass
    # dsmcFoam+ keeps the initial state in system/dsmcInitialiseDict
    try:
        init = open(os.path.join(CASEDIR, 'system', 'dsmcInitialiseDict')).read()
        m = re.search(r'numberDensities\s*\{([^}]*)\}', init)
        if m and params['numberDensity'] == '—':
            toks = re.findall(r'([A-Za-z0-9_.+-]+)\s+([-+0-9.eE]+)', m.group(1))
            if toks: params['numberDensity'] = num(toks[0][1])
    except Exception:
        pass
    try:
        cd = open(os.path.join(CASEDIR, 'system', 'controlDict')).read()
        m = re.search(r'deltaT\s+([-+0-9.eE]+)', cd)
        if m: params['deltaT'] = m.group(1)
    except Exception: pass
    try:
        owner = open(os.path.join(CASEDIR, 'constant', 'polyMesh', 'owner')).read()
        bm = re.search(r'\d+\s*\((.*)\)\s*//', owner, re.S)
        if bm:
            vals = [int(x) for x in bm.group(1).split()]
            params['nCells'] = str(max(vals) + 1)
    except Exception: pass
    return params


def wall_models():
    """List (patchName, model-description) from the dsmcFoam+ boundariesDict.

    Reads the dsmcPatchBoundaries entries, so what is shown is exactly the
    case setup: diffuse / specular / mixed (with its diffuseFraction) /
    deletion.  Falls back to the v2606 WallInteractionModel when there is no
    boundariesDict.
    """
    out = []
    path = os.path.join(CASEDIR, 'system', 'boundariesDict')
    try:
        txt = open(path).read()
    except Exception:
        txt = ''

    if txt:
        # every entry: patchBoundaryProperties { patchName X; } ... boundaryModel Y;
        for blk in re.finditer(r'patchBoundaryProperties\s*\{(.*?)\}', txt, re.S):
            tail = txt[blk.end():]
            mn = re.search(r'boundaryModel\s+(\w+)\s*;', tail)
            pn = re.search(r"patchName\s+([^\s;]+)\s*;", blk.group(1))
            if not (mn and pn):
                continue
            patch, model = pn.group(1), mn.group(1)
            extra = ''
            if model.endswith('MixedDiffuseSpecularWallPatch'):
                dm = re.search(r'diffuseFraction\s+([-+0-9.eE]+)', tail)
                extra = '  diffuseFraction %s' % dm.group(1) if dm else '  diffuseFraction ?'
            out.append(('%s: %s%s' % (patch, model, extra)))
        for blk in re.finditer(r'generalBoundaryProperties\s*\{(.*?)\}', txt, re.S):
            tail = txt[blk.end():]
            mn = re.search(r'boundaryModel\s+(\w+)\s*;', tail)
            pn = re.search(r"patchName\s+([^\s;]+)\s*;", blk.group(1))
            if mn and pn:
                out.append('%s: %s (inflow)' % (pn.group(1), mn.group(1)))
        if out:
            return out

    # v2606 dsmcFoam fallback
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
        return float(json.load(open(CONFIG_FILE)).get('threshold', 0.01))
    except Exception:
        return 0.01


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
    """True if the relative residual stayed below threshold for a sustained run."""
    n = len(rows)
    if n < 2:
        return False
    window = max(20, int(0.05 * n))
    if n < window:
        return False
    return all(relative_residual(r) < threshold for r in rows[-window:])


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

        def X(x): return pad_l + (x) / xmax * pw
        def Y(y): return pad_t + ph - (y - ymin) / (ymax - ymin) * ph

        for i in range(5):
            yv = ymin + (ymax - ymin) * i / 4.0
            cv.create_line(pad_l, Y(yv), W - pad_r, Y(yv), fill='#eee')
            cv.create_text(pad_l - 6, Y(yv), text=fmt(yv), anchor='e', font=('Helvetica', 8), fill='#444')
        for i in range(5):
            xv = xmax * i / 4.0
            cv.create_text(X(xv), H - pad_b + 12, text=fmt(xv), anchor='n', font=('Helvetica', 8), fill='#444')
        for name, color, pts in series:
            if len(pts) < 2:
                continue
            coords = []
            for x, y in pts:
                coords += [X(x), Y(y)]
            cv.create_line(*coords, fill=color, width=2)
        cv.create_text(pad_l + 8, pad_t - 12, text=title, anchor='w', font=('Helvetica', 10, 'bold'), fill='#111')
        ly = pad_t + 14
        for name, color, _pts in series:
            cv.create_line(W - pad_r - 40, ly, W - pad_r - 20, ly, fill=color, width=2)
            cv.create_text(W - pad_r - 16, ly, text=name, anchor='w', font=('Helvetica', 8), fill='#111')
            ly += 16


# ------------------------------------------------------------------ app ----
class MonitorApp:
    def __init__(self, root):
        self.root = root
        self.end_time = read_end_time()
        self.proc = None
        self.script = None
        self._last_running = None
        self._last_plot_t = None
        self._last_running = None
        self._stopped_by_user = False
        self._end_alerted = False
        self._was_proc_running = None
        self.threshold = load_threshold()

        root.title('DSMC convergence monitor — mass flow rate')
        root.geometry('980x720')

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

        # ---- left panel: parameters + threshold ----
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
        # walls, straight from system/boundariesDict
        wf = ttk.LabelFrame(left, text='Boundaries (system/boundariesDict)')
        wf.pack(fill='x', padx=6, pady=(8, 2))
        self.walls_var = tk.StringVar(value='—')
        ttk.Label(wf, textvariable=self.walls_var, anchor='w', justify='left',
                  font=('Menlo', 7)).pack(fill='x', padx=4, pady=2)

        self.update_params()

        # live flow time (under case parameters)
        fr = ttk.Frame(left)
        fr.pack(fill='x', padx=6, pady=(12, 2))
        ttk.Label(fr, text='flow time (s):', width=22, anchor='w').pack(side='left')
        self.flow_time_var = tk.StringVar(value='—')
        ttk.Label(fr, textvariable=self.flow_time_var, width=16, anchor='e').pack(side='right')

        # threshold + refresh
        th = ttk.Frame(left)
        th.pack(fill='x', padx=6, pady=(14, 2))
        ttk.Label(th, text='Convergence threshold (frac):').pack(side='left')
        self.thresh_var = tk.StringVar(value=str(self.threshold))
        ttk.Entry(th, textvariable=self.thresh_var, width=8).pack(side='left', padx=4)
        ttk.Button(th, text='Set', command=self.set_threshold).pack(side='left')
        ttk.Button(left, text='Refresh parameters', command=self.update_params)\
            .pack(anchor='w', padx=6, pady=(2, 0))
        self.residual_label = tk.StringVar(value='residual: —')
        ttk.Label(left, textvariable=self.residual_label).pack(anchor='w', padx=6, pady=(4, 0))

        # statistical-uncertainty readout (from dsmcStatistics FO summary.dat)
        statframe = ttk.LabelFrame(left, text='Statistical convergence (dsmcStatistics)')
        statframe.pack(fill='both', expand=True, padx=6, pady=(10, 4))
        self.statbox = tk.Text(statframe, height=9, font=('Menlo', 7), state='disabled',
                               bg='#f7f7f7', wrap='none')
        self.statbox.pack(fill='both', expand=True, padx=3, pady=3)

        # ---- plots ----
        right = ttk.Frame(root)
        right.pack(side='left', fill='both', expand=True, padx=4, pady=4)
        self.plot_m = MiniPlot(right)
        self.plot_net = MiniPlot(right)

        self.refresh()

    # --------------------------------------------------------- actions -----
    def _launch(self, script):
        if self.proc and self.proc.poll() is None:
            return
        path = os.path.join(CASEDIR, script)
        if not os.path.isfile(path):
            self.status.set('%s not found in %s' % (script, CASEDIR))
            return

        # The case scripts source the OpenFOAM environment themselves, so no
        # shell setup is needed here. Output goes to log.<script> so that a
        # failing step is visible instead of vanishing into /dev/null.
        logpath = os.path.join(CASEDIR, 'log.' + script)
        try:
            logf = open(logpath, 'w')
        except OSError as err:
            self.status.set('cannot open %s: %s' % (logpath, err))
            return

        try:
            self.proc = subprocess.Popen(['bash', path], cwd=CASEDIR,
                                         start_new_session=True,
                                         stdout=logf, stderr=subprocess.STDOUT)
        except Exception as err:
            self.status.set('launch failed: %s' % err)
            logf.close()
            return
        finally:
            self._logfile = logf

        self.script = script
        self.last_script = script
        self._stopped_by_user = False
        self._end_alerted = False
        self.status.set('running %s (pid %d) -> log.%s' % (script, self.proc.pid, script))

    def reconstruct(self):
        # reconstruct every written time (not just the latest)
        self._launch('Allreconstruct')

    def make_foam(self):
        name = os.path.basename(os.path.normpath(CASEDIR)) + '.foam'
        try:
            open(os.path.join(CASEDIR, name), 'w').close()
            self.status.set('created %s' % name)
        except Exception as err:
            self.status.set('touch failed: %s' % err)

    def clean_procs(self):
        import glob as _glob, shutil
        n = 0
        for d in _glob.glob(os.path.join(CASEDIR, 'processor*')):
            if os.path.isdir(d):
                shutil.rmtree(d, ignore_errors=True)
                n += 1
        self.status.set('removed %d processor dir(s)' % n)

    def run(self): self._launch('Allrun')
    def cont(self): self._launch('Allcont')
    def clean(self): self._launch('Allclean')

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
            json.dump({'threshold': v}, open(CONFIG_FILE, 'w'))
            self.thresh_var.set(str(v))
            self.status.set('threshold set: %g' % v)
        except Exception:
            messagebox.showerror('Set threshold', 'Enter a positive number, e.g. 0.01')

    def update_params(self):
        self.end_time = read_end_time()   # refresh final flow time too
        p = parse_case_params()
        for k, v in p.items():
            self.param_vars[k].set(v)
        self.walls_var.set('\n'.join(wall_models()))

    # ---------------------------------------------------------- refresh ----
    def refresh(self):
        # always read the final flow time fresh from the case dictionaries
        # (controlDict endTime) so a mid-run change is picked up.
        self.end_time = read_end_time()
        proc_running = self.proc is not None and self.proc.poll() is None
        t = current_solution_time()
        running = proc_running   # simulation mode ends when the subprocess exits

        # buttons auto-refresh each poll (so Stop disables once it finishes)
        self.btn_run.config(state='disabled' if running else 'normal')
        self.btn_stop.config(state='normal' if running else 'disabled')
        self.btn_cont.config(state='disabled' if running else 'normal')
        self.btn_clean.config(state='disabled' if running else 'normal')

        # live flow time (panel) + progress bar
        if t is not None:
            self.flow_time_var.set('%.6g' % t)
            self.param_vars['nCells'].set(self.param_vars['nCells'].get())
        if t is not None and self.end_time:
            frac = min(1.0, max(0.0, t / self.end_time))
            self.progress['value'] = 100.0 * frac
            self.progress_label.set('t = %.6g / %.6g s  (%.2f%%)' % (t, self.end_time, 100.0 * frac))
        else:
            self.progress['value'] = 0.0
            self.progress_label.set('progress: no log.dsmcFoam yet (t = %s)' % (t if t is not None else '—'))

        names, rows = read_mass_flow()

        # re-parse nCells once the mesh is generated (blockMesh creates polyMesh)
        if self.param_vars['nCells'].get() == '—' and \
           os.path.exists(os.path.join(CASEDIR, 'constant', 'polyMesh', 'owner')):
            self.update_params()

        # statistical-uncertainty readout
        s = read_stat_summary()
        self.statbox.config(state='normal')
        self.statbox.delete('1.0', 'end')
        self.statbox.insert('1.0', s if s else 'waiting for dsmcStatistics summary.dat')
        self.statbox.config(state='disabled')

        # convergence: green status text only (no popup)
        conv = is_converged(rows, self.threshold)
        if rows:
            self.residual_label.set('residual: %.4g  (threshold %.3g)' % (relative_residual(rows[-1]), self.threshold))

        # redraw plot every 0.1% progress / on state changes / on first data
        redraw = False
        if self.end_time:
            step = PLOT_STEP_FRAC * self.end_time
            if self._last_plot_t is None:
                redraw = True
            elif t is not None and (t - self._last_plot_t) >= step - 1e-12:
                redraw = True
        if running != self._last_running:
            redraw = True
        if redraw:
            if t is not None:
                self._last_plot_t = t
            inlet = [(r[0], r[1]) for r in rows]
            outlet = [(r[0], r[2]) for r in rows]
            net = [(r[0], r[3]) for r in rows]
            name1, name2 = names[0], names[1]
            self.plot_m.plot('mass flow rate (kg/s)',
                             [(name1, 'royalblue', inlet), (name2, 'forestgreen', outlet)], self.end_time)
            self.plot_net.plot('net residual (kg/s)', [('net', 'crimson', net)], self.end_time)

        # simulation finished (subprocess exited) -> text-only, no alert.
        # The buttons refresh automatically from `running` (Stop disables).
        was_proc = self._was_proc_running
        finished = (was_proc and not proc_running and
                    self.script in SOLVE_SCRIPTS and not self._stopped_by_user)
        if finished:
            self.status.set('simulation finished')
        elif was_proc and not proc_running and self.script in SOLVE_SCRIPTS:
            self.status.set('stopped by user')

        self._was_proc_running = proc_running
        self._last_running = running

        if conv and running:
            self.status.set('CONVERGED (residual < %g)' % self.threshold)

        color = 'green' if (conv or finished) else \
                ('red' if self._stopped_by_user else
                 ('dark orange' if running else 'black'))
        self.status_lbl.config(foreground=color)

        self.root.after(POLL_MS, self.refresh)


def main():
    root = tk.Tk()
    MonitorApp(root)
    root.mainloop()


if __name__ == '__main__':
    main()
