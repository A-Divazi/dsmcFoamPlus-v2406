#!/usr/bin/env python3
"""Prepare the fine-mesh DSMC cases (2, 1, 0.5, 0.1, 0.05, 0.01 Pa).

    python3 utils-make-fine-cases.py            # all
    python3 utils-make-fine-cases.py num20-2pa  # one

The mean free path, the mean collision time and the time step are computed
with the functions of the DSMC Setup Calculator
(/home/ubuntu/Desktop/DSMC-Calculator/dsmc_calc.py):

    sigma = pi dref^2 (Tref/T)^(omega-1/2)                 (VHS cross-section)
    lambda = 1 / (sqrt(2) n sigma)                         (mean free path)
    cbar   = sqrt(8 k T / (pi m))                          (mean thermal speed)
    tau    = lambda / cbar                                 (mean collision time)
    dt     = min(f1*tau, f2*dx/(U+cbar))   f1 = 0.2, f2 = 0.3 (calculator defaults)
    dx     = (Vdomain / Ncells)^(1/3)

Mesh: the 116 000-cell ramp mesh (3.5 x 1 x 0.001 m, V = 2.6165e-3 m3) refined
by the factor f in x and y so that dx <= lambda; cells = 116000 f^2.
Parcels: 20 per cell when that fits the RAM budget (PARCEL_MAX), otherwise
capped at PARCEL_MAX (the resulting parcels/cell is reported and flagged).
"""
import math
import os
import re
import shutil
import sys

CALC = '/home/ubuntu/Desktop/DSMC-Calculator'
sys.path.insert(0, CALC)
try:
    import dsmc_calc as CALC_MOD           # the user's DSMC Setup Calculator
except Exception as err:                   # fallback: identical formulas
    CALC_MOD = None
    print('warning: DSMC calculator not importable (%s) - using the built-in '
          'formulas' % err)

B = '/home/ubuntu/OpenFOAM/dsmcFoamPlus-v2306'
V = 2.6165e-3             # domain volume [m3]
K_B = 1.380649e-23
T = 150.0                 # gas temperature [K]
U = 736.5                 # freestream velocity [m/s] (Mach 3.23)
MASS = 66.3e-27           # argon [kg]
DREF = 4.17e-10
OMEGA = 0.81
TREF = 298.15
PARCEL_MAX = 9.0e6        # RAM budget of this machine
F1, F2 = 0.2, 0.3         # calculator safety factors
RUNNERS = ('Allrun', 'Allcont', 'Allclean', 'Allreconstruct')

RUNNER_FILES = {
    'Allrun': r"""#!/bin/bash
cd "$(dirname -- "$0")" || exit 1
#------------------------------------------------------------------------------
# Parallel AMR run.  Ranks from system/decomposeParDict (NP=n overrides).
#   ./Allrun            # full endTime from controlDict
#   END=0.005 ./Allrun  # shorter
#   NP=4 ./Allrun
#------------------------------------------------------------------------------
: "${HOME:=/home/ubuntu}"; export HOME
: "${USER:=ubuntu}"; export USER
# always use the v2406 environment this case was built for, even if the
# shell that started us has v2606 (or nothing) sourced
if [ -z "${WM_PROJECT_DIR:-}" ] || [ "${WM_PROJECT_VERSION:-}" != "v2406" ]; then
    . /home/ubuntu/opt/setup-openfoam2406.sh
fi
. "${WM_PROJECT_DIR:?}/bin/tools/RunFunctions"

[ -n "${END:-}" ] && foamDictionary -entry endTime -set "$END" system/controlDict
[ -n "${NP:-}" ]  && foamDictionary -entry numberOfSubdomains -set "$NP" system/decomposeParDict
nProcs=$(foamDictionary -entry numberOfSubdomains -value system/decomposeParDict)
app=$(getApplication)

rm -f log.*
rm -rf processor*

foamDictionary -entry startFrom -set startTime system/controlDict >/dev/null
# A mesh produced by the adaptation workflow (adapt_mesh.py / one-shot AMR
# splitting) already carries cellLevel/refinementHistory: do not overwrite it
# with blockMesh.
if [ -f constant/polyMesh/cellLevel ]; then
    echo "adapted mesh found in constant/polyMesh (cellLevel) - skipping blockMesh"
else
    runApplication blockMesh
fi

# The initial parcels and per-cell fields in 0/ belong to one particular mesh.
# After a change of the blockMesh resolution they must be regenerated, otherwise
# decomposePar stops with "Size N is not equal to the expected length M".
meshCells=$(python3 - "$PWD/constant/polyMesh/owner" <<'PYEOF'
import re, sys
s = open(sys.argv[1]).read()
m = re.search(r'\n(\d+)\s*\n\((.*?)\)\s*\n', s, re.S)
print(max(int(x) for x in m.group(2).split()) + 1 if m else 0)
PYEOF
)
v0=0
[ -f 0/V0 ] && v0=$(python3 - "$PWD/0/V0" <<'PYEOF'
import re, sys
m = re.search(r'internalField\s+nonuniform\s+List<[\w:]+>\s*\n(\d+)', open(sys.argv[1]).read())
print(m.group(1) if m else 0)
PYEOF
)
if [ -d 0 ] && [ "$v0" != "0" ] && [ "$v0" != "$meshCells" ]; then
    echo "0/ was written for ${v0} cells but the mesh has ${meshCells}: re-initialising" >&2
    rm -rf 0
fi
if [ ! -d 0 ]; then
    # fresh case: dsmcInitialisePlus creates 0/ (fields, parcels, parcel dicts)
    echo "no 0/ directory - creating it with a serial dsmcInitialisePlus" >&2
    runApplication dsmcInitialisePlus
    rm -f log.dsmcInitialisePlus
fi

runApplication decomposePar -force
runParallel -np "$nProcs" dsmcInitialisePlus
runParallel -np "$nProcs" "$app"
# AMR writes a new mesh per write, and reconstructed meshes have no
# *ProcAddressing, so reconstructPar can only merge times before the first
# refinement.  Do not let that failure look like a failed run.
runApplication reconstructPar -no-lagrangian -newTimes || \
    echo "note: reconstructPar could not merge the AMR times (data is complete in processor*/)"
echo "=== $PWD done (endTime $(foamDictionary -entry endTime -value system/controlDict)) ==="
""",
    'Allcont': r"""#!/bin/bash
cd "$(dirname -- "$0")" || exit 1
# Continue from the latest written time (parallel, ranks from decomposeParDict).
: "${HOME:=/home/ubuntu}"; export HOME
: "${USER:=ubuntu}"; export USER
# always use the v2406 environment this case was built for, even if the
# shell that started us has v2606 (or nothing) sourced
if [ -z "${WM_PROJECT_DIR:-}" ] || [ "${WM_PROJECT_VERSION:-}" != "v2406" ]; then
    . /home/ubuntu/opt/setup-openfoam2406.sh
fi
. "${WM_PROJECT_DIR:?}/bin/tools/RunFunctions"
[ -n "${END:-}" ] && foamDictionary -entry endTime -set "$END" system/controlDict
[ -n "${NP:-}" ]  && foamDictionary -entry numberOfSubdomains -set "$NP" system/decomposeParDict
nProcs=$(foamDictionary -entry numberOfSubdomains -value system/decomposeParDict)
app=$(getApplication)
foamDictionary -entry startFrom -set latestTime system/controlDict
rm -f "log.$(basename "$app")"
runParallel -np "$nProcs" "$app"
runApplication reconstructPar -no-lagrangian
echo "=== $PWD continued ==="
""",
    'Allclean': r"""#!/bin/bash
cd "$(dirname -- "$0")" || exit 1
# Remove solver output but KEEP 0/ (the initial fields) and the mesh in constant/.
: "${HOME:=/home/ubuntu}"; export HOME
: "${USER:=ubuntu}"; export USER
# always use the v2406 environment this case was built for, even if the
# shell that started us has v2606 (or nothing) sourced
if [ -z "${WM_PROJECT_DIR:-}" ] || [ "${WM_PROJECT_VERSION:-}" != "v2406" ]; then
    . /home/ubuntu/opt/setup-openfoam2406.sh
fi

# every numeric time directory except 0 (and not e.g. 0.014.old)
for d in [0-9]*; do
    [ -d "$d" ] || continue
    [[ "$d" =~ ^[0-9]+(\.[0-9]+)?$ ]] || continue
    [ "$d" = "0" ] && continue
    rm -rf "$d"
done
rm -rf processor* postProcessing dynamicCode log.* .monitor
echo "cleaned $PWD (0/ and constant/polyMesh kept)"
""",
    'Allreconstruct': r"""#!/bin/bash
cd "$(dirname -- "$0")" || exit 1
#------------------------------------------------------------------------------
# Reconstruct EVERY written time (not only the latest) — use after stopping a
# parallel run by hand.  Lagrangian/parcel data is skipped (field data only).
#------------------------------------------------------------------------------

: "${HOME:=/home/ubuntu}"; export HOME
: "${USER:=ubuntu}"; export USER

if [ -z "${WM_PROJECT_DIR:-}" ] || [ "${WM_PROJECT_VERSION:-}" != "v2406" ]
then
    . /home/ubuntu/opt/setup-openfoam2406.sh
fi
. "${WM_PROJECT_DIR:?}/bin/tools/RunFunctions"

rm -f log.reconstructPar
# only the times written before the first AMR refinement can be merged (their
# decomposed meshes still match the decomposed base mesh)
runApplication reconstructPar -no-lagrangian || \
    echo "note: reconstructPar could not merge the AMR times (see log.reconstructPar)"
""",
}


# case directory, pressure [Pa], mesh refinement factor
SPEC = [
    ('num20-2pa',                 2.0,    6),
    ('num20-20particle-1pa',      1.0,    4),
    ('num19.5-20particle-0.5pa',  0.5,    2),
    ('num19',                     0.1,    1),
    ('num18.5-20particle-0.05pa', 0.05,   1),
    ('num18-20particle-0.01pa',   0.01,   1),
]


def sigma(Tg):
    if CALC_MOD:
        return CALC_MOD.vhs_cross_section(DREF, TREF, OMEGA, Tg)
    return math.pi * DREF ** 2 * (TREF / Tg) ** (OMEGA - 0.5)


def lam(n, Tg):
    if CALC_MOD:
        return CALC_MOD.mean_free_path(n, sigma(Tg))
    return 1.0 / (math.sqrt(2.0) * n * sigma(Tg))


def cbar(Tg):
    if CALC_MOD:
        return CALC_MOD.mean_thermal_speed(MASS, Tg)
    return math.sqrt(8.0 * K_B * Tg / (math.pi * MASS))


def patch_blockmesh(path, f):
    if f == 1:
        return
    s = open(path).read()
    s = re.sub(r'\((\d+)\s+(\d+)\s+(\d+)\)(\s+simpleGrading)',
               lambda m: '(%d %d %d)%s' % (int(m.group(1)) * f,
                                           int(m.group(2)) * f,
                                           int(m.group(3)), m.group(4)), s)
    open(path, 'w').write(s)


def set_entry(path, pattern, value):
    s = open(path).read()
    if not re.search(pattern, s):
        print('   ! %s: pattern %r not found' % (path, pattern))
        return
    open(path, 'w').write(re.sub(pattern, value, s))


def build(src_name, P, f, only=None):
    src = os.path.join(B, src_name)
    dst = src + '-fine'
    if only and only not in (src_name, os.path.basename(dst)):
        return None
    if not os.path.isdir(src):
        print('missing source case', src)
        return None
    if os.path.isdir(dst):
        shutil.rmtree(dst)
    os.makedirs(os.path.join(dst, 'system'))
    os.makedirs(os.path.join(dst, 'constant'))
    for x in os.listdir(os.path.join(src, 'system')):
        shutil.copy2(os.path.join(src, 'system', x), os.path.join(dst, 'system', x))
    for x in os.listdir(os.path.join(src, 'constant')):
        p = os.path.join(src, 'constant', x)
        if os.path.isfile(p):
            shutil.copy2(p, os.path.join(dst, 'constant', x))
    for x in RUNNERS:
        with open(os.path.join(dst, x), 'w') as fh:
            fh.write(RUNNER_FILES[x])
        os.chmod(os.path.join(dst, x), 0o755)
    sp = os.path.join(B, 'monitor.py.new')
    if os.path.isfile(sp):
        shutil.copy2(sp, os.path.join(dst, 'monitor.py'))
    open(os.path.join(dst, 'monitor_config.json'), 'w').write(
        '{\n  "threshold": 0.01,\n  "sections": [1.0, 2.0],\n'
        '  "massflow_write": false\n}\n')

    # --- mesh ---------------------------------------------------------------
    patch_blockmesh(os.path.join(dst, 'system', 'blockMeshDict'), f)
    cells = 116000 * f * f
    dx = (V / cells) ** (1.0 / 3.0)

    # --- DSMC parameters (calculator) ---------------------------------------
    n = P / (K_B * T)
    s_T = sigma(T)
    lbd = lam(n, T)
    cb = cbar(T)
    tau = lbd / cb
    dt_coll = F1 * tau
    dt_trans = F2 * dx / (U + cb)
    dt = min(dt_coll, dt_trans)

    want = 20 * cells
    parcels = min(want, PARCEL_MAX)
    capped = parcels < want
    nEq = n * V / parcels

    cd = os.path.join(dst, 'system', 'controlDict')
    set_entry(cd, r'(?m)^deltaT\s+.*;', 'deltaT          %.6g;' % dt)
    set_entry(cd, r'(?m)^endTime\s+.*;', 'endTime         0.02;')
    # the ASCII parcel files of a 9e6-parcel case are ~1.5 GB per write, so the
    # heavy meshes write less often and keep only the last few times
    heavy = parcels >= 4.0e6
    write_int = '0.005' if heavy else '0.0005'
    set_entry(cd, r'(?m)^writeInterval\s+.*;', 'writeInterval   %s;' % write_int)
    set_entry(cd, r'(?m)^purgeWrite\s+.*;', 'purgeWrite      %d;' % (6 if heavy else 0))
    set_entry(os.path.join(dst, 'constant', 'dsmcProperties'),
              r'nEquivalentParticles\s+[0-9.eE+-]+;',
              'nEquivalentParticles            %.6e;' % nEq)
    init = os.path.join(dst, 'system', 'dsmcInitialiseDict')
    s = open(init).read()
    s = re.sub(r'(numberDensities\s*\{\s*\w+\s+)[0-9.eE+-]+;',
               lambda m: m.group(1) + '%.6e;' % n, s)
    open(init, 'w').write(s)

    note = ('**Warning:** 20 parcels/cell needs %d parcels; only %d fit in the '
            'RAM budget, so this case runs with %.1f parcels/cell. Lower the '
            'mesh factor or use a bigger machine for better statistics.'
            % (want, parcels, parcels / cells)) if capped else \
           'Statistics: 20 parcels per cell.'

    open(os.path.join(dst, 'README.md'), 'w').write('''# %s - %g Pa on the fine mesh (factor %d)

Prepared with the DSMC Setup Calculator
(`/home/ubuntu/Desktop/DSMC-Calculator/dsmc_calc.py`, functions
`vhs_cross_section`, `mean_free_path`, `mean_thermal_speed`).

| quantity | value |
|---|---|
| pressure | %g Pa |
| number density n | %.4e 1/m^3 (T = %g K) |
| VHS cross-section sigma | %.4e m^2 |
| mean free path lambda | %.4e m |
| mean thermal speed cbar | %.3f m/s |
| mean collision time tau | %.4e s |
| mesh | 116 000 x %d^2 = **%d cells** |
| cell size dx = (V/N)^(1/3) | %.4e m |
| dx / lambda | %.2f |
| dt collision limit f1*tau (f1 = %g) | %.3e s |
| dt transit limit f2*dx/(U+cbar) (f2 = %g) | %.3e s |
| **dt used** | **%.3e s** |
| parcels | %d (%s) |
| nEquivalentParticles | %.4e |

%s

Run:

    ./Allrun            # blockMesh -> dsmcInitialisePlus -> dsmcFoamPlus (mpirun)
    END=2e-3 ./Allrun   # short
    python3 monitor.py  # mass flow | live Ma | mesh

The mesh is graded in y (simpleGrading 10/-10): dx above is the bulk
equivalent cell size used by the calculator; the near-wall cells are smaller.
''' % (os.path.basename(dst), P, f, P, n, T, s_T, lbd, cb, tau, f, cells, dx,
       dx / lbd, F1, dt_coll, F2, dt_trans, dt, parcels,
       'capped at the RAM budget PARCEL_MAX' if capped else '20 per cell',
       nEq, note))

    return dict(case=os.path.basename(dst), P=P, f=f, cells=cells, dx=dx,
                lbd=lbd, tau=tau, dt=dt, parcels=parcels, nEq=nEq,
                capped=capped, ppc=parcels / cells, dxlam=dx / lbd)


def main():
    only = sys.argv[1] if len(sys.argv) > 1 else None
    rows = []
    for src, P, f in SPEC:
        print('building %s-fine (%g Pa, factor %d)' % (src, P, f))
        r = build(src, P, f, only)
        if r:
            rows.append(r)
    if not rows:
        return
    print()
    print('%-28s %7s %3s %11s %10s %9s %10s %10s %6s %7s'
          % ('case', 'P[Pa]', 'f', 'cells', 'dx[m]', 'dx/lam', 'dt[s]',
             'parcels', 'ppc', 'flag'))
    for r in rows:
        print('%-28s %7g %3d %11d %10.3e %9.2f %10.2e %10d %6.1f %7s'
              % (r['case'], r['P'], r['f'], r['cells'], r['dx'], r['dxlam'],
                 r['dt'], r['parcels'], r['ppc'],
                 'CAPPED' if r['capped'] else ''))


if __name__ == '__main__':
    main()
