# num20-2pa-fine - 2 Pa on the fine mesh (factor 6)

Prepared with the DSMC Setup Calculator
(`/home/ubuntu/Desktop/DSMC-Calculator/dsmc_calc.py`, functions
`vhs_cross_section`, `mean_free_path`, `mean_thermal_speed`).

| quantity | value |
|---|---|
| pressure | 2 Pa |
| number density n | 9.6573e+20 1/m^3 (T = 150 K) |
| VHS cross-section sigma | 6.7594e-19 m^2 |
| mean free path lambda | 1.0832e-03 m |
| mean thermal speed cbar | 282.033 m/s |
| mean collision time tau | 3.8408e-06 s |
| mesh | 116 000 x 6^2 = **4176000 cells** |
| cell size dx = (V/N)^(1/3) | 8.5570e-04 m |
| dx / lambda | 0.79 |
| dt collision limit f1*tau (f1 = 0.2) | 7.682e-07 s |
| dt transit limit f2*dx/(U+cbar) (f2 = 0.3) | 2.520e-07 s |
| **dt used** | **2.520e-07 s** |
| parcels | 9000000 (capped at the RAM budget PARCEL_MAX) |
| nEquivalentParticles | 2.8076e+11 |

**Warning:** 20 parcels/cell needs 83520000 parcels; only 9000000 fit in the RAM budget, so this case runs with 2.2 parcels/cell. Lower the mesh factor or use a bigger machine for better statistics.

Run:

    ./Allrun            # blockMesh -> dsmcInitialisePlus -> dsmcFoamPlus (mpirun)
    END=2e-3 ./Allrun   # short
    python3 monitor.py  # mass flow | live Ma | mesh

The mesh is graded in y (simpleGrading 10/-10): dx above is the bulk
equivalent cell size used by the calculator; the near-wall cells are smaller.
