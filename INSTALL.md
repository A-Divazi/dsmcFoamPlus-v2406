# Installing dsmcFoam+ (v2406 build) on another computer

This installs the **solver** (`dsmcFoamPlus`, `dsmcInitialisePlus`, the
`libUserDsmcFoam` function objects) and then runs the cases of `cases/`.

Tested on Ubuntu 24.04 / 26.04 (arm64 and amd64) with OpenFOAM v2406
(ESI, patch 241212). Everything is CPU-only; no GPU is needed.

Requirements: Linux, ~3 GB disk for OpenFOAM + repo, 8 GB RAM minimum
(the prepared 2 Pa / 1 Pa cases with 9 M parcels want 12–16 GB), and the
packages

```bash
sudo apt install -y build-essential git flex bison zlib1g-dev libreadline-dev \
                    libopenmpi-dev openmpi-bin python3 python3-tk \
                    python3-numpy python3-matplotlib
```

(`python3-tk`/`numpy`/`matplotlib` are only for `monitor.py`.)

---

## 1. Get OpenFOAM v2406

dsmcFoam+ v2406 **must** be built against OpenFOAM v2406 (the ESI version,
`WM_PROJECT_VERSION=v2406`). Three ways, pick one.

### 1a. ESI packages (easiest on Ubuntu)

Download and install the **v2406** package set from the official page
<https://www.openfoam.com/download> (or clone their apt repository and
`apt install openfoam2406-*`). The installation ends up either in
`/opt/openfoam2406` (packaged install) or in a directory containing
`usr/lib/openfoam/openfoam2406` (user-space unpack).

Check it:

```bash
source /opt/openfoam2406/etc/bashrc       # or the etc/bashrc of your unpack
echo $WM_PROJECT_VERSION                  # must print v2406
```

### 1b. Build OpenFOAM v2406 from source (any distribution, any CPU)

```bash
git clone https://develop.openfoam.com/Development/openfoam.git openfoam-v2406
cd openfoam-v2406
git checkout OpenFOAM-v2406
source etc/bashrc
./Allwmake -j 4 -s -l        # 30-60 min; use -j 2 on a small machine
```

Result: `$WM_PROJECT_DIR` points at that directory, `WM_PROJECT_VERSION=v2406`.

### 1c. Copy the install from the machine it already runs on (same CPU arch)

The working install on the development VM is 548 MB:

```
/home/ubuntu/opt/openfoam2406root/        (contains usr/lib/openfoam/openfoam2406)
/home/ubuntu/opt/setup-openfoam2406.sh
```

Copy the first directory over (e.g. `rsync -a`) on the new machine and source
the setup script with `OF2406_ROOT` pointing at it. **Only works if the new
machine has the same CPU architecture** (the binaries are
`linuxARM64GccDPInt32Opt` — the development VM is arm64/Apple-M4).

---

## 2. Clone and build dsmcFoam+

```bash
git clone https://github.com/A-Divazi/dsmcFoamPlus-v2406.git
cd dsmcFoamPlus-v2406
```

Make OpenFOAM v2406 available (one of):

```bash
source /opt/openfoam2406/etc/bashrc                  # 1a / 1b, then:
export LD_LIBRARY_PATH=$PWD/scripts/compat-lib:$LD_LIBRARY_PATH

# or, with the user-space layout of 1c:
export OF2406_ROOT=$HOME/opt/openfoam2406root
source scripts/setup-openfoam2406.sh
```

Create the scotch/boost compatibility sonames (Ubuntu/Debian only — needed by
`decomposePar` and the parallel runs):

```bash
./scripts/make-compat-lib.sh
```

Build everything (installs into `$FOAM_USER_APPBIN` / `$FOAM_USER_LIBBIN`,
i.e. `~/OpenFOAM/<user>-v2406/platforms/...`):

```bash
JOBS=1 ./scripts/build-dsmcfoamplus.sh        # JOBS=1 for low-RAM machines
```

Expected output at the end:

```
libdsmcFoamPlus.so  libgeneralMoleculeDynamics.so  libUserDsmcFoam.so
dsmcFoamPlus  dsmcInitialisePlus
```

### If the build fails

| symptom | fix |
|---|---|
| `wmake: command not found` | OpenFOAM not sourced — `source …/etc/bashrc` |
| `cannot find -lPstream` / MPI headers | `sudo apt install libopenmpi-dev openmpi-bin` |
| `regularGridInterpolator.H: No such file` | build order — run `./scripts/build-dsmcfoamplus.sh`, not hand-made `wmake` calls |
| `invalid new-expression of abstract class` | you are building against the wrong OpenFOAM (not v2406) |
| killed during compile | `JOBS=1` (default), or add swap |

---

## 3. Verify

```bash
dsmcFoamPlus -help | head
dsmcInitialisePlus -help | head
ls $FOAM_USER_LIBBIN
```

Quick physics test (shipped v2406 tutorial, not the ramp case):

```bash
cp -r v2406-dsmcFoam-tutorials/freeSpacePeriodic $HOME/test-freeSpace
cd $HOME/test-freeSpace
blockMesh > log.blockMesh && dsmcInitialisePlus && dsmcFoamPlus
```

---

## 4. Run the cases

The prepared cases are the two converted ones (`cases/Mach-3.23`,
`cases/wedge15Ma5`) and the six fine-mesh ramp cases
`cases/num20-2pa-fine`, `num20-20particle-1pa-fine`,
`num19.5-20particle-0.5pa-fine`, `num19-fine`,
`num18.5-20particle-0.05pa-fine`, `num18-20particle-0.01pa-fine`
(parameters in §4 of `README.md`).

Each case directory contains `system/`, `constant/dsmcProperties`, the runner
scripts and `monitor.py` — no mesh, no results.

```bash
cd cases/num19-fine            # 0.1 Pa, 116 000 cells, 2.3M parcels
./Allrun                       # blockMesh -> dsmcInitialisePlus -> decomposePar
                               #   -> dsmcFoamPlus (mpirun) -> reconstructPar
END=2e-3 ./Allrun              # short test run
NP=4 ./Allrun                  # fewer MPI ranks (default from decomposeParDict)
./Allcont                      # continue from the latest written time
./Allclean                     # remove time dirs/logs, keep 0/ and the mesh
python3 monitor.py             # GUI: mass flow | live Ma | mesh
```

Parallel runs need `mpirun` from the system OpenMPI (installed above); the
number of ranks comes from `system/decomposeParDict` (`numberOfSubdomains`).

Results: time directories in the case root (after `reconstructPar`), mass flow
in `postProcessing/massFlow.csv`, parcels in `<time>/lagrangian/dsmc`. For
ParaView: `touch case.foam` and open that file.

Cost hints (8 ranks, 14 GB machine): 464 000 cells with 9M parcels ≈ 3–4 s per
step; a 0.02 s run is 24 000–80 000 steps on the fine meshes (see the Δt column
of the README table) — use `END=` for short tests first. Each 9M-parcel write is
~1.5 GB (ASCII); the heavy cases write every 0.005 s with `purgeWrite 6`.

The GUI needs a display. On a remote machine use `ssh -X`, or run the script
head-less for a single picture: `python3 monitor.py --render mesh <time> out.png`.

---

## 5. Everyday use

```bash
# in every new shell:
source /opt/openfoam2406/etc/bashrc        # or: source scripts/setup-openfoam2406.sh
export LD_LIBRARY_PATH=$PWD/scripts/compat-lib:$LD_LIBRARY_PATH   # after the build

# always start a case run from its own directory:
cd cases/<case> && ./Allrun
```

Common pitfalls

* **`prun:exe-not-found` / “cannot open case directory processor0”** — the
  shell has a *different* OpenFOAM sourced (e.g. v2606) while `dsmcFoamPlus`
  lives in the v2406 tree. Make sure `WM_PROJECT_VERSION=v2406`; the case
  scripts re-source v2406 automatically only when the path
  `/home/ubuntu/opt/setup-openfoam2406.sh` exists, so on another machine put
  `source …/etc/bashrc` into your `~/.bashrc`.
* **`libscotch-7.0.so: cannot open shared object file`** — run
  `./scripts/make-compat-lib.sh` (or install the ESI packages, which ship
  their own scotch).
* **`Size N is not equal to the expected length M` in `decomposePar`** — `0/`
  was generated for a different mesh resolution; `rm -rf 0` (or `./Allclean`)
  and re-run `./Allrun`, which regenerates it.
* **Out of memory** — lower the parcel count: `nEquivalentParticles` in
  `constant/dsmcProperties` (parcels = n·V/nEq; ×Y = fewer parcels by Y), or
  use `NP=` to change the rank count.
