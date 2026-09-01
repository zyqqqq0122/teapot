#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

JAR=${1:-}
if [ -z "$JAR" ]; then
    echo "usage: assets/build_masslist_decoys.sh <path to encyclopedia jar>" >&2
    echo "   or: apptainer exec <image> cat /opt/encyclopedia/encyclopedia.jar > enc.jar" >&2
    exit 1
fi

tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
javac --release 17 -nowarn -cp "$JAR" -d "$tmp" assets/MassListDecoyGenerator.java
( cd "$tmp" && jar cf masslist-decoys.jar ./*.class )
mv "$tmp/masslist-decoys.jar" assets/masslist-decoys.jar

echo "built assets/masslist-decoys.jar"
echo "verify with:  apptainer exec <image> java -cp <enc.jar>:assets/masslist-decoys.jar MassListDecoyGenerator"
