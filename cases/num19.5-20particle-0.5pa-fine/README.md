# num19.5-20particle-0.5pa-fine - 0.5 Pa on the fine mesh (factor 2)

Prepared with the DSMC Setup Calculator
(`/home/ubuntu/Desktop/DSMC-Calculator/dsmc_calc.py`, functions
`vhs_cross_section`, `mean_free_path`, `mean_thermal_speed`).

| quantity | value |
|---|---|
| pressure | 0.5 Pa |
| number density n | 2.4143e+20 1/m^3 (T = 150 K) |
| VHS cross-section sigma | 6.7594e-19 m^2 |
| mean free path lambda | 4.3329e-03 m |
| mean thermal speed cbar | 282.033 m/s |
| mean collision time tau | 1.5363e-05 s |
| mesh | 116 000 x 2^2 = **464000 cells** |
| cell size dx = (V/N)^(1/3) | 1.7799e-03 m |
| dx / lambda | 0.41 |
| dt collision limit f1*tau (f1 = 0.2) | 3.073e-06 s |
| dt transit limit f2*dx/(U+cbar) (f2 = 0.3) | 5.243e-07 s |
| **dt used** | **5.243e-07 s** |
| parcels | 9000000 (capped at the RAM budget PARCEL_MAX) |
| nEquivalentParticles | 7.0190e+10 |

**Warning:** 20 parcels/cell needs 9280000 parcels; only 9000000 fit in the RAM budget, so this case runs with 19.4 parcels/cell. Lower the mesh factor or use a bigger machine for better statistics.

Run:

    ./Allrun            # blockMesh -> dsmcInitialisePlus -> dsmcFoamPlus (mpirun)
    END=2e-3 ./Allrun   # short
    python3 monitor.py  # mass flow | live Ma | mesh

The mesh is graded in y (simpleGrading 10/-10): dx above is the bulk
equivalent cell size used by the calculator; the near-wall cells are smaller.
