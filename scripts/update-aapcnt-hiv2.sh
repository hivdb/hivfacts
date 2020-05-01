#! /bin/bash

BASEDIR=`realpath $(dirname $0)/..`
DATADIR="$BASEDIR/data/aapcnt-hiv2"

set -e

cd $BASEDIR
mkdir -p $DATADIR

for rx in all naive art; do
    for subtype in all; do
        pipenv run hivdbql export-aapcnt2 \
            --species HIV2 \
            --subtype $subtype \
            --rx-type $rx \
            --filter NO_QA_ISSUES \
            --format "json" \
            --unusual-threshold 0.005 \
            --drm-list $BASEDIR/data/drms_hiv2.json \
            "$DATADIR/rx-${rx}_subtype-${subtype}.json"
        ls -lh "$DATADIR/rx-${rx}_subtype-${subtype}.json"
    done
done
