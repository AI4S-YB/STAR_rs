#!/usr/bin/env bash
# Diff STAR outputs between the C++ reference and the Rust port.
# Allowed differences (per plan §0, §7):
#   - compilation-time/server/dir header line
#   - version/run-time lines in Log.out / Log.progress.out
#   - absolute paths inside log files
set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "usage: $0 <orig-dir-or-file> <rs-dir-or-file>" >&2
    exit 64
fi

orig="$1"; rs="$2"

filter() {
    sed \
        -E -e '/^STAR compilation time,server,dir=/d' \
           -e '/^[A-Z][a-z]{2}[[:space:]]+[0-9]{1,2}[[:space:]]+[0-9]{2}:[0-9]{2}:[0-9]{2}/d'
}

if [[ -d "$orig" && -d "$rs" ]]; then
    # Directory mode: for every file in orig, look for counterpart in rs.
    rc=0
    for f in $(cd "$orig" && find . -type f | sort); do
        o="$orig/$f"; r="$rs/$f"
        if [[ ! -e "$r" ]]; then echo "MISSING: $f"; rc=1; continue; fi
        if ! diff -u <(filter <"$o") <(filter <"$r") > "/tmp/diff.$$"; then
            echo "DIFF: $f"
            head -n 20 "/tmp/diff.$$"
            rc=1
        fi
        rm -f "/tmp/diff.$$"
    done
    exit $rc
else
    diff -u <(filter <"$orig") <(filter <"$rs")
fi
