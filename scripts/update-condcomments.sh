#! /bin/bash

BASEDIR=`realpath $(dirname $0)/..`
DATADIR="$BASEDIR/data"

set -e
cd $BASEDIR

pipenv run python3 scripts/export_condcomments.py HIV1 $DATADIR/conditional-comments_hiv1.json
pipenv run python3 scripts/export_condcomments.py HIV2 $DATADIR/conditional-comments_hiv2.json
