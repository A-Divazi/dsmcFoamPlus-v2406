# dsmcFoam+ for OpenFOAM v2406

`dsmcFoam+` is the extended DSMC solver of the IAP Kühlungsborn group
(Staszak et al., see the original repository
[`igit.iap-kborn.de/staszak/dsmcfoamplus-v2306`](https://igit.iap-kborn.de/staszak/dsmcfoamplus-v2306)).
This repository is a **ready-to-build backup for OpenFOAM v2406 (ESI)**:
the upstream `v2306` sources plus the fixes needed to compile and run them on
v2406, install scripts, the upstream tutorials, the v2406 `dsmcFoam`
tutorials for reference, two fully converted working cases, and the custom
function objects developed here (mass flow rate, wall shear stress, density
gradient, Mach number, statistics).

Everything below was verified on Ubuntu 26.04 / arm64 with a user-space
OpenFOAM v2406 install.

## What is in here

```
src/                      dsmcFoam+ sources (fixed for v2406)
applications/             dsmcInitialisePlus + dsmcFoamPlus
utils/functionObjects/    libUserDsmcFoam: dsmcMassFlowRate, wallShearStressDsmc,
                          dsmcGradRho, dsmcMachNumber, dsmcStatistics
tutorials/                upstream dsmcFoam+ tutorials
                          (simpleSphere, simpleDiffusiveMesh_AMR, statisticalMeshDemos)
v2406-dsmcFoam-tutorials/ the stock OpenFOAM v2406 dsmcFoam tutorials (reference)
cases/Mach-3.23/          Ma 3.23 ramp case, converted to dsmcFoam+ (system, constant,
                          Allrun/Allcont/Allclean/Allreconstruct, monitor.py)
cases/wedge15Ma5/         Ma 5 15-degree wedge, converted from the v2406 tutorial
cases/*-fine/             the ramp case at 2, 1, 0.5, 0.1, 0.05 and 0.01 Pa on the
                          fine mesh, with the time step from the DSMC Setup
                          Calculator (see "4. Fine-mesh cases" below)
scripts/                  setup / build / scotch-compat helpers, make-fine-cases.py
```

## 1. Get OpenFOAM v2406

Any v2406 installation works (ESI precompiled debs for Ubuntu 22.04/24.04,
or a source build). With the debs unpacked into a directory of your choice:

```bash
export OF2406_ROOT=/path/to/openfoam2406root     # contains usr/lib/openfoam/openfoam2406
source scripts/setup-openfoam2406.sh             # exports the OpenFOAM v2406 environment
```

`scripts/setup-openfoam2406.sh` also adds `scripts/compat-lib` to
`LD_LIBRARY_PATH`. Create that directory once with

```bash
./scripts/make-compat-lib.sh
```

It only makes symlinks such as `libscotch-7.0.so -> libscotch.so.7.0`: the ESI
binaries are linked against ESI's sonames, while the distribution ships the
Debian ones under different names. Without it `decomposePar` fails with
`libscotch-7.0.so: cannot open shared object file`.

## 2. Build

```bash
source scripts/setup-openfoam2406.sh
./scripts/build-dsmcfoamplus.sh          # JOBS=1 by default (low-RAM friendly)
```

Result: `libgeneralMoleculeDynamics.so`, `libdsmcFoamPlus.so`,
`libUserDsmcFoam.so` in `$FOAM_USER_LIBBIN`; `dsmcInitialisePlus` and
`dsmcFoamPlus` in `$FOAM_USER_APPBIN`.

`FOAM_USER_SRC` must point at this repository's `src/` (the build script sets
it). If you build by hand, build in this order:

```bash
export FOAM_USER_SRC=$PWD/src
wmake -j1 src/lagrangian/molecularDynamics/general
wmake -j1 src/lagrangian/dsmc
wmake -j1 applications/utilities/preProcessing/dsmc/dsmcInitialisePlus
wmake -j1 applications/solvers/discreteMethods/dsmc/dsmcFoamPlus
wmake -j1 utils/functionObjects
```

## 3. Run a case

```bash
cd cases/Mach-3.23
./Allrun            # always parallel; ranks from system/decomposeParDict
NP=4 ./Allrun       # one-off override (writes numberOfSubdomains into the dict)
./Allcont           # continue from the latest written time
./Allreconstruct    # reconstruct ALL written times (after Ctrl-C)
./Allclean          # remove results, logs, processor dirs, mesh
python3 monitor.py  # GUI monitor: Run/Stop/Continue/Clean/Reconstruct/.foam/del procs
```

The case scripts source the v2406 environment themselves, so they also work
when started from the monitor GUI or a plain shell.

## 4. What differs from the upstream `v2306` sources

Compilation on v2406 (each of these is a real error you would otherwise hit):

1. `src/lagrangian/dsmc/Make/options` — added
   `-I$(FOAM_USER_SRC)/lagrangian/dsmc/boundaries/derived/cyclicBoundaries`.
   The three `dsmcCyclic*StatisticalSpecular/Diffusive/IsotropicPatch` headers
   include `regularGridInterpolator.H`, which lives one directory above them,
   and the committed `lnInclude` symlink for it is broken.
2. `src/lagrangian/dsmc/Make/files` — `dsmcSpecularWallPatch` and
   `dsmcMixedDiffuseSpecularWallPatch` were commented out; they are enabled
   and now compile (both were also missing the two pure virtuals
   `autoBoundaryModelMap()` / `writeField()` required by the current
   `dsmcPatchBoundary` interface — added as empty implementations).
3. `dsmcBoundaries.C` — general (inflow) boundary output path: the absolute
   case path was appended inside each time directory
   (`.../generalBoundaryModels/home/user/.../boundaries/...`); now uses the
   model name like the patch models do.
4. `dsmcVolFields.C` —
   * removed an uncommented debug `Pout << rhoNInstantaneous_ << ...` that
     dumped the whole per-cell field on every rank every step (a 5-minute
     8-rank run produced a 459 MB log);
   * clamped a negative translational temperature (statistical noise in cells
     with almost no parcels) which otherwise raised `FE_INVALID` inside a
     later `sqrt()` and aborted with `FOAM_SIGFPE`.

Known limitation: `measureMeanFreePath true` in `system/fieldPropertiesDict`
still trips that floating-point trap when a cell has (almost) no parcels —
keep it `false` (the shipped cases do).

## 5. Custom function objects (`utils/functionObjects`, `libUserDsmcFoam.so`)

| object | output |
|---|---|
| `dsmcMassFlowRate` | mass flow through interior section planes or patches → `postProcessing/massFlow.csv`; builds `momentum = rhoM_mix * UMean_mix` because dsmcFoam+ writes no `momentum`, and processes each *written* time directory once (`executeControl timeStep; writeControl none;`) |
| `wallShearStressDsmc` | `wallShearStress`, `wallShearStressMag` vol fields from `fD_mix` |
| `dsmcGradRho` | `gradRho`, `gradRhoMag` from `rhoN_mix` |
| `dsmcMachNumber` | `Ma` from a velocity and temperature field |
| `dsmcStatistics` | statistical uncertainty (tau_int, N_eff, Student-t CI, block averaging) |

`cases/Mach-3.23/system/controlDict` shows how they are enabled.

## 6. Conventions of dsmcFoam+ cases

* `0/` is **generated** by `dsmcInitialisePlus` — never ship it.
* boundaries live in `system/boundariesDict` (three lists:
  `dsmcPatchBoundaries`, `dsmcGeneralBoundaries`, `dsmcCyclicBoundaries`);
  every non-`empty`/`symmetry`/`wedge`/`cyclic` patch needs exactly one model.
* `system/fieldPropertiesDict` replaces `fieldAverage`/`dsmcFields` function
  objects; `system/controllersDict`, `system/chemReactDict`,
  `system/loadBalanceDict` must exist (they may be empty).
* a `symmetryPlane` patch is **not** recognised as symmetry by this solver —
  use `symmetry` in `blockMeshDict`.

A full model reference (wall interaction models and their implementation,
collision models, controllers, initialisation configurations, field output
keys, and how to add a new model) is in
[`cases/Mach-3.23/README.md`](cases/Mach-3.23/README.md) and in the
`MODELS.md` section of this repository's documentation.

## License

The dsmcFoam+ sources are GPL v3 (OpenFOAM / OpenCFD Ltd. and the dsmcFoam+
authors: C. White et al. 2018, T. Staszak 2024). The converted cases,
scripts and function objects in this repository are provided under the same
terms.

## 4. Fine-mesh cases (2 Pa ... 0.01 Pa)

`cases/*-fine/` contain the ramp case prepared with the DSMC Setup Calculator
(`dsmc_calc.py`: `vhs_cross_section`, `mean_free_path`, `mean_thermal_speed`):

    sigma  = pi dref^2 (Tref/T)^(omega-1/2)
    lambda = 1/(sqrt(2) n sigma)          cbar = sqrt(8 k T / (pi m))
    tau    = lambda / cbar
    dt     = min(0.2 tau, 0.3 dx/(U+cbar))        (calculator defaults f1, f2)
    dx     = (Vdomain/Ncells)^(1/3)

Argon, T = 150 K, U = 736.5 m/s (Mach 3.23), V = 2.6165e-3 m3 (3.5 x 1 x 0.001 m),
116 000-cell base mesh refined by f in x and y so that dx <= lambda.

| case | P [Pa] | f | cells | dx [m] | dx/lambda | lambda [m] | tau [s] | dt [s] | parcels |
|---|---|---|---|---|---|---|---|---|---|
| num20-2pa-fine | 2 | 6 | 4 176 000 | 8.56e-4 | 0.79 | 1.08e-3 | 3.84e-6 | 2.52e-7 | 9.0e6 (2.2/cell, RAM cap) |
| num20-20particle-1pa-fine | 1 | 4 | 1 856 000 | 1.12e-3 | 0.52 | 2.17e-3 | 7.68e-6 | 3.30e-7 | 9.0e6 (4.8/cell, RAM cap) |
| num19.5-20particle-0.5pa-fine | 0.5 | 2 | 464 000 | 1.78e-3 | 0.41 | 4.33e-3 | 1.54e-5 | 5.24e-7 | 9.0e6 (19.4/cell, RAM cap) |
| num19-fine | 0.1 | 1 | 116 000 | 2.83e-3 | 0.13 | 2.17e-2 | 7.68e-5 | 8.32e-7 | 2.32e6 (20/cell) |
| num18.5-20particle-0.05pa-fine | 0.05 | 1 | 116 000 | 2.83e-3 | 0.07 | 4.33e-2 | 1.54e-4 | 8.32e-7 | 2.32e6 (20/cell) |
| num18-20particle-0.01pa-fine | 0.01 | 1 | 116 000 | 2.83e-3 | 0.01 | 2.17e-1 | 7.68e-4 | 8.32e-7 | 2.32e6 (20/cell) |

`dt` is the smaller of the collision limit (0.2 tau) and the transit limit
(0.3 dx/(U+cbar)); on the fine meshes the transit limit dominates. The mesh is
graded in y, so `dx` is the bulk equivalent cell size used by the calculator
while the near-wall cells are smaller.

Parcels are 20 per cell where that fits the RAM budget of the very small
development VM (9e6); for the three refined meshes 20/cell would need 37 M / 84 M
parcels, so they are capped and the statistics are 2.2 - 19.4 parcels/cell.
Rebuild the family with

    python3 scripts/make-fine-cases.py            # all
    python3 scripts/make-fine-cases.py num20-2pa  # one

(this writes system/, constant/dsmcProperties, the runner scripts, monitor.py
and a README with the calculated values into `<case>-fine`).
