#! /bin/bash

BASEDIR=`realpath $(dirname $0)/..`
DATADIR="$BASEDIR/data"

set -e
cd $BASEDIR

pipenv run python3 scripts/export_mut_scores_indiv.py $DATADIR/mutation-scores-indiv.json
pipenv run python3 scripts/export_mut_scores_combi.py $DATADIR/mutation-scores-combi.json