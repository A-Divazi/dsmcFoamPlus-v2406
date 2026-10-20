# Mach-3.23 ramp case — dsmcFoam+ (v2406)

Source: `/home/ubuntu/Desktop/New-Runs/3.23` (dsmcFoam, OpenFOAM v2606),
converted to dsmcFoam+ following the layout of the working `wedge15Ma5` case.

## Setup (unchanged physics)

| Item | Value |
|---|---|
| Gas | Ar, `mass 66.3e-27`, `d 4.17e-10`, VHS `omega 0.81` |
| Free stream | n = 4.8286e19 1/m^3, U = (736.5 0 0) m/s, T = 150 K (Ma ≈ 3.23) |
| Walls | diffuse (`dsmcDiffuseWallPatch`), 150 K |
| `nEquivalentParticles` | 5.5214e11 → 229 723 parcels on 116 000 cells |
| Time | `deltaT 5e-6`, `endTime 0.05`, `writeInterval 1e-3` |

## dsmcFoam+ extra function objects (`system/controlDict`)

`libs ("libUserDsmcFoam.so")` (built against v2406 from
`/home/ubuntu/OpenFOAM-notes/dsmc/src/wallShearStressDsmc`) adds:

| entry | type | output |
|---|---|---|
| `massFlow1` | `dsmcMassFlowRate` | `postProcessing/massFlow.csv` — mass flow through the interior planes x = 1.0 and x = 2.0 (m·U summed over section faces), plus the `net` residual |
| `wallShear1` | `wallShearStressDsmc` | `wallShearStress`, `wallShearStressMag` (vol fields from the DSMC `fD_mix` force density, visible in ParaView) |
| `gradRho1` | `dsmcGradRho` | `gradRho`, `gradRhoMag` from `rhoN_mix` |

dsmcFoam+ writes no `momentum` field, so `dsmcMassFlowRate` builds it as
`rhoM_mix × UMean_mix` (configurable with `rhoMField`/`UField`).

Because dsmcFoam+ creates the time directory and *then* writes the fields, the
object is configured with `executeControl timeStep; writeControl none;`: it
polls for a time directory that already contains the fields, processes it once
(`process(foundTime)`), and skips directories that are still incomplete. Each
row of `massFlow.csv` therefore corresponds exactly to its time stamp (verified:
2.025e-6 → 1.961e-6 → 1.967e-6 kg/s over 0.0002…0.0006).

dsmcFoam+ itself already provides `Ma_mix`, `p_mix`, `surfaceShearStress_mix`,
`surfaceHeatTransfer_mix`, `cellMeanFreePathRatio_mix` and the statistical
error fields `densityError_mix`, `pressureError_mix`, `temperatureError_mix`,
`velocityError_mix` — the separate v2606 `dsmcMachNumber`/`dsmcStatistics`
objects are therefore not enabled here.

## Run — always parallel

```bash
source /home/ubuntu/opt/setup-openfoam2406.sh   # optional: the scripts do it
./Allrun            # uses numberOfSubdomains from system/decomposeParDict
NP=8 ./Allrun       # one-off override, written into the dictionary
./Allcont           # resume from the latest written time, same rank count
./Allclean
```

`Allrun` is self-contained (it sources the v2406 environment, setting `HOME`/
`USER` itself so the GUI or a plain shell can start it) and always
parallel:

`rm -f log.*` → `rm -rf processor*` → `blockMesh` → `foamDictionary
numberOfSubdomains=$NP, method=scotch` → `decomposePar -force` →
`dsmcInitialisePlus` → `dsmcFoamPlus` → `reconstructPar -latestTime
-no-lagrangian`.

`Allrun` removes stale `log.*` (OpenFOAM skips silently otherwise), forces
`numberOfSubdomains = $NP` + `method scotch` in `system/decomposeParDict`,
runs `blockMesh → decomposePar → dsmcInitialisePlus → dsmcFoamPlus →
reconstructPar -latestTime -no-lagrangian`.

## Monitoring

```bash
python3 monitor.py      # inside the case directory
```

`monitor.py` was updated to detect the solver from `controlDict`
(`application dsmcFoamPlus` → `log.dsmcFoamPlus`), to read the initial state
from `system/dsmcInitialiseDict`, the wall model from `system/boundariesDict`,
and to plot `postProcessing/massFlow.csv`.

## Known dsmcFoam+ issues hit with this case

* `measureMeanFreePath true` in `system/fieldPropertiesDict` aborts with a
  floating point exception in `dsmcVolFields::calculateField()` when a cell has
  (almost) no parcels — the pressure goes slightly negative and a `sqrt()`
  traps. It is therefore `false` here.
* The same negative-temperature path was clamped in the solver source
  (`if (translationalT_[cell] < SMALL) translationalT_[cell] = SMALL;`).
