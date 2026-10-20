# wedge15Ma5 — dsmcFoam+ conversion of the v2406 `dsmcFoam` tutorial

Source case: `/home/ubuntu/opt/openfoam2406/tutorials/discreteMethods/dsmcFoam/wedge15Ma5`
(dsmcFoam, OpenFOAM v2406). Converted to the dsmcFoam+ (v2306 branch) input
format following the layout of the dsmcFoam+ tutorials in `../tutorials/`.

## Physics (unchanged)

| Quantity | Value |
|---|---|
| Gas | air, N2 `0.777e20` + O2 `0.223e20` 1/m^3 |
| Free stream | U = `(1736 0 0)` m/s, T = 300 K → Ma ≈ 5 |
| Wedge | 15° ramp (`obstacle`, wall), diffuse wall at **550 K** |
| VHS | N2: 46.5e-27 kg, 4.17e-10 m, ω=0.74; O2: 53.12e-27 kg, 4.07e-10 m, ω=0.77 |
| Internal DOF | rotational only (2), no vibration/electronic |
| `nEquivalentParticles` | 5e12 → 20 949 parcels on 2400 cells |
| Time | `deltaT 2e-6`, `endTime 0.02` (10 000 steps), write every `1e-3` |

## Conversion map

| v2406 `dsmcFoam` | here (dsmcFoamPlus) |
|---|---|
| `0.orig/` (rhoN, momentum, boundaryT, …) | **none** — `dsmcInitialisePlus` creates `0/` |
| `InflowBoundaryModel FreeStream` | `dsmcGeneralBoundaries` → `dsmcFreeStreamInflowPatch` on `flow` |
| `WallInteractionModel MaxwellianThermal` | `dsmcPatchBoundaries` → `dsmcDiffuseWallPatch` on `obstacle`, `temperature 550` |
| outflow through `flow` | `dsmcPatchBoundaries` → `dsmcDeletionPatch` on `flow` |
| `fieldAverage` + `dsmcFields` function objects | `system/fieldPropertiesDict` (`dsmcVolFields`, `resetAtOutput true`) |
| `forces1` (fD on `obstacle`) | not used — wall fluxes are in `uniform/boundaries/` and `surfaceShearStress_mix` / `surfaceHeatTransfer_mix` |
| — | new required dicts: `boundariesDict`, `fieldPropertiesDict`, `controllersDict`, `chemReactDict`, `loadBalanceDict` |

`constant/dsmcProperties` was rewritten in the dsmcFoam+ layout (full
`moleculeProperties` keys, `collisionPartnerSelectionModel`,
`axisymmetricSimulation`, `rotationalRelaxationCollisionNumber`).

## Run — always parallel

```bash
./Allrun            # uses numberOfSubdomains from system/decomposeParDict
NP=4 ./Allrun       # one-off override, written into the dictionary
./Allcont
./Allclean
```

Self-contained (sources the v2406 environment) and always parallel:
`blockMesh → decomposePar → dsmcInitialisePlus → dsmcFoamPlus →
reconstructPar -latestTime -no-lagrangian`.

`Allrun` (parallel branch) does
`blockMesh -> decomposePar -> dsmcInitialisePlus -> dsmcFoamPlus -> reconstructPar -latestTime -no-lagrangian`
and first

* removes `log.*` (OpenFOAM's `runApplication`/`runParallel` **silently skip**
  a step when its log file already exists),
* sets `numberOfSubdomains = $NP` and `method scotch` in
  `system/decomposeParDict` (the solver aborts if the dict and the MPI rank
  count disagree; `foamDictionary` is used, so the file is rewritten),

## Timings measured here (2400 cells, ~21 k parcels, 10 000 steps)

| run | clock | steps/s | notes |
|---|---|---|---|
| serial | 9 min | ~62 | reference fields |
| `NP=4` | 7.3 min (439 s solver) | ~23 | `reconstructPar` ok, no OOM |
| `NP=8` | — | ~25 | **OOM-killed** (`/proc/vmstat oom_kill`), only ~0.022 reached |

Parallel is *not* faster for this small case (only ~600 cells/rank); the
overhead is MPI + per-rank memory. Use it for larger cases.

Or step by step:

```bash
blockMesh && dsmcInitialisePlus && dsmcFoamPlus
```

## Result

Serial run, single core, exit 0 in ≈ 9 min for the full `endTime 0.02`
(~62 steps/s on 2400 cells / ~21 k parcels). Time directories `0.001` … `0.02`
each hold `Ma_mix`, `rhoN_mix`, `rhoM_mix`, `p_mix`, `UMean_mix`,
`translationalT_mix`, `dsmcRhoN_mix`, mean-free-path and error fields, plus the
`lagrangian` parcel data.
