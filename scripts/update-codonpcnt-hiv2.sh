#! /bin/bash

BASEDIR=`realpath $(dirname $0)/..`
DATADIR="$BASEDIR/data/codonpcnt-hiv2"

set -e
cd $BASEDIR

mkdir -p $DATADIR

for rx in all naive art; do
    pipenv run hivdbql export-codonpcnt \
        --species HIV2 \
        --subtype all \
        --rx-type $rx \
        --format json \
        --filter NO_QA_ISSUES \
        "$DATADIR/rx-${rx}_subtype-all.json"
    ls -lh "$DATADIR/rx-${rx}_subtype-all.json"
done
