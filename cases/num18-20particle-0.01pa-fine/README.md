# num18-20particle-0.01pa-fine - 0.01 Pa on the fine mesh (factor 1)

Prepared with the DSMC Setup Calculator
(`/home/ubuntu/Desktop/DSMC-Calculator/dsmc_calc.py`, functions
`vhs_cross_section`, `mean_free_path`, `mean_thermal_speed`).

| quantity | value |
|---|---|
| pressure | 0.01 Pa |
| number density n | 4.8286e+18 1/m^3 (T = 150 K) |
| VHS cross-section sigma | 6.7594e-19 m^2 |
| mean free path lambda | 2.1665e-01 m |
| mean thermal speed cbar | 282.033 m/s |
| mean collision time tau | 7.6816e-04 s |
| mesh | 116 000 x 1^2 = **116000 cells** |
| cell size dx = (V/N)^(1/3) | 2.8254e-03 m |
| dx / lambda | 0.01 |
| dt collision limit f1*tau (f1 = 0.2) | 1.536e-04 s |
| dt transit limit f2*dx/(U+cbar) (f2 = 0.3) | 8.322e-07 s |
| **dt used** | **8.322e-07 s** |
| parcels | 2320000 (20 per cell) |
| nEquivalentParticles | 5.4458e+09 |

Statistics: 20 parcels per cell.

Run:

    ./Allrun            # blockMesh -> dsmcInitialisePlus -> dsmcFoamPlus (mpirun)
    END=2e-3 ./Allrun   # short
    python3 monitor.py  # mass flow | live Ma | mesh

The mesh is graded in y (simpleGrading 10/-10): dx above is the bulk
equivalent cell size used by the calculator; the near-wall cells are smaller.
