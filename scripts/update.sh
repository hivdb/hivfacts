#! /bin/bash

BASEDIR=`realpath $(dirname $0)/..`
DATADIR="$BASEDIR/data"

set -e

for rx in all naive art; do
    for subtype in All A B C D F G CRF01_AE CRF02_AG Other; do
        hivdbql export_aapcnt \
            --subtype $subtype \
            --rx-type $rx \
            --format "json" \
            "$DATADIR/${rx}${subtype}.json"
        ls -lh "$DATADIR/${rx}${subtype}.json"
    done
done
