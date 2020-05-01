#! /bin/bash

BASEDIR=`realpath $(dirname $0)/..`
DATADIR="$BASEDIR/data"

set -e
cd $BASEDIR

pipenv run python3 scripts/export_drmlist.py HIV1 $DATADIR/drms_hiv1.json
pipenv run python3 scripts/export_drmlist.py HIV2 $DATADIR/drms_hiv2.json
