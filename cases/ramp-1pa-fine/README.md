# ramp-1pa-fine - 1 Pa on the fine mesh (factor 4)

Prepared with the DSMC Setup Calculator
(`/home/ubuntu/Desktop/DSMC-Calculator/dsmc_calc.py`, functions
`vhs_cross_section`, `mean_free_path`, `mean_thermal_speed`).

| quantity | value |
|---|---|
| pressure | 1 Pa |
| number density n | 4.8286e+20 1/m^3 (T = 150 K) |
| VHS cross-section sigma | 6.7594e-19 m^2 |
| mean free path lambda | 2.1665e-03 m |
| mean thermal speed cbar | 282.033 m/s |
| mean collision time tau | 7.6816e-06 s |
| mesh | 116 000 x 4^2 = **1856000 cells** |
| cell size dx = (V/N)^(1/3) | 1.1213e-03 m |
| dx / lambda | 0.52 |
| dt collision limit f1*tau (f1 = 0.2) | 1.536e-06 s |
| dt transit limit f2*dx/(U+cbar) (f2 = 0.3) | 3.303e-07 s |
| **dt used** | **3.303e-07 s** |
| parcels | 37120000 (20 per cell) |
| nEquivalentParticles | 3.4036e+10 |

**Statistics:** 20 parcels per cell (DSMC calculator default).

Estimated memory for the initial parcel set: **~22 GB** (0.55 GB per million parcels + ~1 GB per million cells) - run it on a machine with that much RAM, or rebuild with a smaller mesh factor.

Run:

    ./Allrun            # blockMesh -> dsmcInitialisePlus -> dsmcFoamPlus (mpirun)
    END=2e-3 ./Allrun   # short
    python3 monitor.py  # mass flow | live Ma | mesh

The mesh is graded in y (simpleGrading 10/-10): dx above is the bulk
equivalent cell size used by the calculator; the near-wall cells are smaller.
