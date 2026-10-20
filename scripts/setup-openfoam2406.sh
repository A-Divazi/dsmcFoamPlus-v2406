#!/bin/bash
# Source this to use a user-space OpenFOAM v2406 (ESI) installation.
# Override the install root with OF2406_ROOT if it lives elsewhere.
#
#   export OF2406_ROOT=/path/to/openfoam2406root
#   source scripts/setup-openfoam2406.sh
#
: "${OF2406_ROOT:=/home/ubuntu/opt/openfoam2406root}"
_foam="${OF2406_ROOT}/usr/lib/openfoam/openfoam2406/etc/bashrc"
[ -f "$_foam" ] || { echo "OpenFOAM v2406 not found under $OF2406_ROOT" >&2; return 1 2>/dev/null || exit 1; }
source "$_foam"

# The ESI binaries need scotch/boost sonames that the debs do not ship; point
# them at the distribution libraries under the ESI names (7.0.x series).
_compat="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/compat-lib"
if [ -d "$_compat" ]; then
    case ":$LD_LIBRARY_PATH:" in
    *":$_compat:"*) ;;
    *) export LD_LIBRARY_PATH="$_compat:$LD_LIBRARY_PATH" ;;
    esac
fi
unset _foam _compat
