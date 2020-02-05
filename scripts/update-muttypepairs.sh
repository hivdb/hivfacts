#! /bin/bash

BASEDIR=`realpath $(dirname $0)/..`
DATADIR="$BASEDIR/data"

set -e
cd $BASEDIR

pipenv run python3 scripts/export_muttypepairs.py HIV1 $DATADIR/mutation-type-pairs_hiv1.json
pipenv run python3 scripts/export_muttypepairs.py HIV2 $DATADIR/mutation-type-pairs_hiv2.json
