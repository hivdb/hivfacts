#! /bin/bash

BASEDIR=`realpath $(dirname $0)/..`
DATADIR="$BASEDIR/data/aapcnt"

set -e

for rx in all naive art; do
    for subtype in all A B C D F G CRF01_AE CRF02_AG other; do
        hivdbql export-aapcnt \
            --subtype $subtype \
            --rx-type $rx \
            --format "json" \
            "$DATADIR/rx-${rx}_subtype-${subtype}.json"
        ls -lh "$DATADIR/rx-${rx}_subtype-${subtype}.json"
    done
done
