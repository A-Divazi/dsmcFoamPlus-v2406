#!/bin/bash
# Create ESI-named sonames (libscotch-7.0.so, ...) from the distribution
# libraries so that decomposePar and friends can load. Ubuntu/Debian only.
set -e
here="$(cd "$(dirname "$0")" && pwd)"
out="$here/compat-lib"
mkdir -p "$out"
for pair in \
    libscotch-7.0.so:libscotch.so.7.0 \
    libscotcherrexit-7.0.so:libscotcherrexit.so.7.0 \
    libptscotch-7.0.so:libptscotch.so.7.0 \
    libptscotcherrexit-7.0.so:libptscotcherrexit.so.7.0
do
    esi=${pair%%:*}; sys=${pair##*:}
    for d in /usr/lib/aarch64-linux-gnu /usr/lib/x86_64-linux-gnu /usr/lib; do
        if [ -e "$d/$sys" ]; then ln -sf "$d/$sys" "$out/$esi"; echo "linked $esi"; break; fi
    done
done
