#! /bin/bash

BASEDIR=`realpath $(dirname $0)/..`
DATADIR="$BASEDIR/data/codonpcnt"

set -e

for rx in all naive art; do
    for subtype in all A B C D F G CRF01_AE CRF02_AG; do
        hivdbql export-codonpcnt \
            --species HIV1 \
            --subtype $subtype \
            --rx-type $rx \
            --format json \
            --filter PUBLISHED_ONLY \
            --filter NO_DNACHIP \
            --filter NO_QA_ISSUES \
            "$DATADIR/rx-${rx}_subtype-${subtype}.json"
        ls -lh "$DATADIR/rx-${rx}_subtype-${subtype}.json"
    done
done
