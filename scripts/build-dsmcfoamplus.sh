#!/bin/bash
# Build dsmcFoam+ against OpenFOAM v2406 and install the user function objects.
#
#   source scripts/setup-openfoam2406.sh
#   ./scripts/build-dsmcfoamplus.sh
#
# Use -j1 on machines with little RAM (parallel wmake is easily OOM-killed).
set -o pipefail
here="$(cd "$(dirname "$0")" && pwd)"
root="$(dirname "$here")"

[ -n "${WM_PROJECT_DIR:-}" ] || { echo "source scripts/setup-openfoam2406.sh first" >&2; exit 1; }
echo "OpenFOAM ${WM_PROJECT_VERSION} from ${WM_PROJECT_DIR}"

export FOAM_USER_SRC="$root/src"
JOBS=${JOBS:-1}

build() {
    echo "=== $1"
    wmake -j"$JOBS" "$root/$1" || { echo "FAILED: $1" >&2; exit 1; }
}

# 1. general molecule dynamics library (creates the lnInclude used by the rest)
build src/lagrangian/molecularDynamics/general
# 2. the dsmc library
build src/lagrangian/dsmc
# 3. convert / initialise utility
build applications/utilities/preProcessing/dsmc/dsmcInitialisePlus
# 4. the solver
build applications/solvers/discreteMethods/dsmc/dsmcFoamPlus
# 5. optional user function objects (mass flow, wall shear, grad rho, ...)
build utils/functionObjects

echo
echo "installed:"
ls -l "$FOAM_USER_LIBBIN"/libdsmcFoamPlus.so "$FOAM_USER_LIBBIN"/libgeneralMoleculeDynamics.so \
      "$FOAM_USER_LIBBIN"/libUserDsmcFoam.so 2>/dev/null
ls -l "$FOAM_USER_APPBIN"/dsmcFoamPlus "$FOAM_USER_APPBIN"/dsmcInitialisePlus 2>/dev/null
