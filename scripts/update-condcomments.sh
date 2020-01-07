#! /bin/bash

BASEDIR=`realpath $(dirname $0)/..`
DATADIR="$BASEDIR/data"

set -e
cd $BASEDIR

pipenv run python3 scripts/export_condcomments.py $DATADIR/conditional-comments.json
