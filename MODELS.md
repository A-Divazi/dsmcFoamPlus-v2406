MODEL REFERENCE — every selectable model in dsmcFoam+ (v2306 sources, built on OpenFOAM v2406)

Everything in dsmcFoam+ that is *selectable* is a run-time-selection model:
each derived class registers itself with
`defineTypeNameAndDebug(<name>, 0)` + `addToRunTimeSelectionTable(...)`, and the
keyword you write in a dictionary is exactly the class name. Model-specific
settings are read from the sub-dictionary named `<typeName>Properties`
(e.g. `dsmcDiffuseWallPatch` → `dsmcDiffuseWallPatchProperties`).

Source paths below are relative to
`/home/ubuntu/opt/dsmcFoamPlus-src/src/lagrangian/dsmc/`.

## 1. Wall interaction models (`system/boundariesDict` → `dsmcPatchBoundaries`)

Applies to every **wall** patch (and to `patch`-type patches such as an
outlet). Each non-`empty`/`symmetry`/`wedge`/`cyclic`/`processor` mesh patch
must have **exactly one** entry — `dsmcBoundaries::checkPatchBoundaryModels()`
aborts otherwise ("Number of poly-patches = N … not equal to the number of
patch models = M").

| Model (dict keyword) | Source | Properties | Physics |
|---|---|---|---|
| `dsmcDiffuseWallPatch` | `boundaries/derived/patchBoundaries/dsmcDiffuseWallPatch` | `temperature` (K, required), `velocity` (wall velocity, default `(0 0 0)`) | fully thermal-diffuse (Maxwellian) reflection at `temperature`; velocities resampled from a half-Maxwellian plus wall drift, rotational/vibrational/electronic energy equipartitioned at `temperature` |
| `dsmcSpecularWallPatch` | `…/dsmcSpecularWallPatch` | none | mirror reflection, `U -= 2 (U·n) n`; energy unchanged (no accommodation) |
| `dsmcMixedDiffuseSpecularWallPatch` | `…/dsmcMixedDiffuseSpecularWallPatch` | `diffuseFraction` (0–1, required), `temperature`, `velocity` | per collision: diffuse if `rnd < diffuseFraction`, else specular — the dsmcFoam+ equivalent of `MixedDiffuseSpecular` |
| `dsmcDeletionPatch` | `…/dsmcDeletionPatch` | `allSpecies yes/no` | particle is removed (`td.keepParticle = false`) — used for outlets/absorbing surfaces |

Skeleton (`dsmcMassFlowRate`-style example):

```cpp
boundary
{
    patchBoundaryProperties { patchName obstacle; }
    boundaryModel   dsmcMixedDiffuseSpecularWallPatch;
    dsmcMixedDiffuseSpecularWallPatchProperties
    {
        diffuseFraction 0.5;
        temperature     300;
        velocity        (0 0 0);
    }
}
```

Implementation notes:

* `controlParticle(dsmcParcel&, trackingData&)` is called from
  `dsmcParcel::hitWallPatch` for each molecule hitting the patch; it modifies
  `p.U()`, `p.ERot()`, `p.vibLevel()`, `p.ELevel()` in place.
* `measurePropertiesBeforeControl()` / `measurePropertiesAfterControl()`
  accumulate the momentum/energy fluxes used for wall shear stress and heat
  transfer (`surfaceShearStress_<field>`, `surfaceHeatTransfer_<field>`, and
  `fD_<field>`, the force density on walls).
* Wall-normal at the hit point comes from `p.normal()`; diffuse reflection
  needs a tangential direction, and if the molecule arrives exactly along the
  normal the code perturbs `U` by `(0.8 + 0.2·rnd)` per component.
* Time-series per patch are written to `<time>/uniform/boundaries/dsmc/patchBoundaryModels/<patchName>`.
* Species-specific energies use `cloud_.equipartitionRotationalEnergy`,
  `equipartitionVibrationalEnergyLevel`, `equipartitionElectronicLevel`
  (so `vibrationalModes`/`rotationalDegreesOfFreedom`/electronic level lists in
  `constant/dsmcProperties` matter here).

## 2. Non-wall (open) boundaries (`dsmcGeneralBoundaries`)

| Model | Source | Properties |
|---|---|---|
| `dsmcFreeStreamInflowPatch` | `boundaries/derived/generalBoundaries/dsmcFreeStreamInflowPatch` | `typeIds (…); numberDensities { … }; velocity; translationalTemperature; rotationalTemperature; vibrationalTemperature; electronicTemperature` |
| `dsmcFreeStreamInflowFieldPatch` | `…/dsmcFreeStreamInflowFieldPatch` | only `typeIds (… )`; the inflow state is read **per face** from the boundary fields `boundaryT`, `boundaryU`, `boundaryNumberDensity` (so a spatially varying inflow can be prescribed) |

They *inject* molecules every step (`controlParcelsBeforeMove()` computes the
Bird eq. 4.22 flux, `exp(-s²) + √π·s·(1+erf s)`); outgoing molecules are **not**
removed by these models, so an outlet on the same patch needs
`dsmcDeletionPatch` as its `dsmcPatchBoundaries` entry (this is how the
Mach-3.23 and wedge cases handle a combined inlet/outlet `flow` patch).

## 3. Cyclic / statistical-mesh boundaries (`dsmcCyclicBoundaries`)

Require a `cyclic` patch pair (e.g. made with `createBaffles`).

| Model | Source | Properties | Behaviour |
|---|---|---|---|
| `dsmcCyclicStatisticalSpecularPatch` | `boundaries/derived/cyclicBoundaries/…` | `Tr` (transmission probability), `velocity` | particle passing the mesh is reflected specularly with probability `1-Tr` |
| `dsmcCyclicStatisticalDiffusivePatch` | 〃 | `Tr`, `velocity`, `temperature` | reflected diffusely |
| `dsmcCyclicForwardScatterStatisticalSpecularPatch` | 〃 | `Tr`, `velocity` | forward-scattering statistics (uses `grid.h` tables + `regularGridInterpolator`) |
| `dsmcCyclicForwardScatterStatisticalDiffusivePatch` | 〃 | `Tr`, `velocity`, `temperature` | as above, diffuse scatter |
| `dsmcCyclicForwardScatterStatisticalIsotropicPatch` | 〃 | `Tr`, `velocity` | isotropic scatter |

## 4. Binary collision models (`constant/dsmcProperties`)

| Keyword | Source | Coeffs sub-dict | Keys |
|---|---|---|---|
| `NoBinaryCollision` | `collisions/derived/NoBinaryCollision` | — | free-molecular flow |
| `VariableHardSphere` | `collisions/derived/VariableHardSphere` | `VariableHardSphereCoeffs` | `Tref` |
| `VariableSoftSphere` | `collisions/derived/VariableSoftSphere` | `VariableSoftSphereCoeffs` | `Tref` (+ `alpha` from the molecule properties) |
| `LarsenBorgnakkeVariableHardSphere` | `collisions/derived/LarsenBorgnakkeVariableHardSphere` | `LarsenBorgnakkeVariableHardSphereCoeffs` | `Tref`, `rotationalRelaxationCollisionNumber`, `electronicRelaxationCollisionNumber` |
| `LarsenBorgnakkeVariableSoftSphere` | `collisions/derived/LarsenBorgnakkeVariableSoftSphere` | `LarsenBorgnakkeVariableSoftSphereCoeffs` | same as above |

Selection key in `dsmcProperties`: `BinaryCollisionModel <keyword>;`.
`dsmcParcel::constantProperties` supplies `mass`, `diameter`, `omega`, `alpha`,
the DOF counts and the electronic/vibrational level data used by
`pow(2kT/(mR cR²), ω)` etc.

Collision partner selection: `collisionPartnerSelectionModel noTimeCounter;`
(only model implemented — `collisionPartnerSelection/derived/noTimeCounter`,
Bird's NTC, driven by the `dsmcSigmaTcRMax` field).

## 5. Initial configurations (`system/dsmcInitialiseDict`)

```cpp
configurations
(
    configuration
    {
        type            dsmcMeshFill;          // dsmcConfiguration type
        numberDensities { Ar 4.8286e+19; }     // per species
        translationalTemperature 150;
        rotationalTemperature    150;
        vibrationalTemperature   150;
        electronicTemperature    150;
        velocity                 (736.5 0 0);
    }
);
```

Only `dsmcMeshFill` is compiled in (`initialiseDsmcParcels/derived/dsmcMeshFill`);
`dsmcConfiguration` is the abstract base (`initialiseDsmcParcels/basic`), and
`dsmcAllConfigurations` holds the list. `dsmcInitialisePlus` writes `0/`
(fields + parcels), so **never** create `0/` by hand.

## 6. Macroscopic output (`system/fieldPropertiesDict`)

`dsmcFields ( field { fieldModel dsmcVolFields; timeProperties {...}
dsmcVolFieldsProperties {...} } )` — the only compiled field model
(`macroscopicProperties/derived/combined/dsmcVolFields`). Keys:

| Key | Meaning |
|---|---|
| `fieldName` | suffix of the written fields, e.g. `mix` → `rhoN_mix`, `Ma_mix`, `UMean_mix`, `p_mix`, `translationalT_mix`, `dsmcRhoN_mix`, `fD_mix`, `surfaceShearStress_mix`, `surfaceHeatTransfer_mix`, … |
| `typeIds` | species list |
| `densityOnly` | skip the momentum/energy fields |
| `measureErrors` | writes `densityError_*`, `pressureError_*`, `temperatureError_*`, `velocityError_*` (from `1/sqrt(N)`) |
| `measureHeatFluxShearStress` | enables the pressure tensor/heat-flux machinery |
| `measureClassifications` | collision classification fractions |
| `measureMeanFreePath` + `mfpReferenceTemperature` | `cellMeanFreePathRatio_*`, `variableHardSphereMeanFreePath_*` — **aborts (FOAM_SIGFPE) when a cell has almost no parcels** on v2406 (see the bug note above); keep it `false` |
| `sampleInterval` | sampling stride |
| `timeProperties { timeOption write\|decoupledFromWriteInterval\|numberWrite; resetAtOutput true; }` | averaging window; `resetAtOutput true` restarts the averages after every write |

## 7. Controllers (`system/controllersDict`)

Two lists: `dsmcStateControllers ( … )` and `dsmcFluxControllers ( … )`.
Each entry has `controllerProperties { … }` + `timeProperties { timeOption… }`.

| Model | Source (`controllers/derived/…`) | Main keys |
|---|---|---|
| `densityController`, `densityZoneController`, `densityControlZone` | `density/…` | `controlSwitch`, `typeId`, `numberDensity`/`massDensity`, `temperature`, `velocity`, `zoneName` |
| `temperatureController`, `temperatureZoneController`, `temperatureBinsController`, `temperatureBerendsenBins` | `temperature/…` | `controlSwitch`, `typeIds`, `controlTemperature`, `tauT`, `componentControl`, `zoneName` |
| `crbsController`, `permeabilityController`, `gravitationalController`, `gravitationalAccelerationController`, `forcingFunctionController`, `timeVaryingOnlyForcingFunctionController` | `others/…` | model specific (`controlSwitch`, `zoneName`, …) |

Both controller families share `timeProperties`
(`timeOption`, `nSamples`, `nAverages`, `resetAtOutput`, `readStateFromFile`,
`initialTimePeriod`) handled by `timeDataMeas`
(`src/lagrangian/molecularDynamics/general/time`).

## 8. Other selectable settings

`constant/dsmcProperties` (read by `dsmcCloud`):

| Key | Meaning |
|---|---|
| `nEquivalentParticles` | real molecules per simulated particle F |
| `typeIdList (… )` | species present |
| `moleculeProperties { <species> { mass, diameter, rotationalDegreesOfFreedom, vibrationalModes, omega, alpha, characteristicVibrationalTemperature, dissociationTemperature, Zref, referenceTempForZref, charDissQuantumLevel, ionisationTemperature, numberOfElectronicLevels, electronicEnergyList, degeneracyList, charge } }` | all read by `dsmcParcel::constantProperties`; the vibrational/electronic lists must be `()` when the corresponding DOF count is 0 |
| `axisymmetricSimulation`, `radialExtentOfDomain`, `maxRadialWeightingFactor` | axisymmetric DSMC (radial weighting) |
| `chargedParticles`, `adsorption` | enable electron/ion and adsorption features |
| `collisionPartnerSelectionModel` | see §4 |
| `BinaryCollisionModel` + coeffs | see §4 |

Other dictionaries: `system/loadBalanceDict`
(`enableBalancing`, `maximumAllowableImbalance`, plus
`system/balanceParDict`/`topoSetDict` for the balancing run),
`system/chemReactDict` (`reactions ( … )`; no reaction models are compiled in
this build — the cloud warns “there are no reactions”),
`system/controlDict` (`nTerminalOutputs` is required).

## 9. Adding a new model (pattern)

1. Copy a sibling directory under the relevant `derived/` folder.
2. Rename the class + `defineTypeNameAndDebug` + `addToRunTimeSelectionTable`
   (`dsmcPatchBoundary`, `dsmcGeneralBoundary`, `dsmcCyclicBoundary`,
   `BinaryCollisionModel`, `collisionPartnerSelection`, `dsmcField`,
   `dsmcConfiguration`, `dsmcStateController`, `dsmcFluxController`).
3. Add the `.C` to `src/lagrangian/dsmc/Make/files`, rebuild with
   `wmake -j1 src/lagrangian/dsmc/`.
4. Set the coefficient sub-dictionary name to `<typeName>Properties` and read
   it with `dict.subDict(typeName + "Properties")`, then implement
   `controlParticle()` / `controlParcelsBeforeMove()` / `calculateField()` as
   appropriate.
