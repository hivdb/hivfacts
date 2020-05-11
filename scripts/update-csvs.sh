#! /bin/bash

BASEDIR=`realpath $(dirname $0)/..`
DATADIR="$BASEDIR/data/"

set -e

cd $BASEDIR

pipenv run python scripts/json2csv.py $DATADIR/aapcnt/rx-all_subtype-all.json
pipenv run python scripts/json2csv.py $DATADIR/aapcnt-hiv2/rx-all_subtype-all.json
pipenv run python scripts/json2csv.py $DATADIR/drms_hiv1.json
pipenv run python scripts/json2csv.py $DATADIR/drms_hiv2.json
pipenv run python scripts/json2csv.py $DATADIR/sdrms_hiv1.json
pipenv run python scripts/json2csv.py $DATADIR/sdrms_hiv2.json
pipenv run python scripts/json2csv.py $DATADIR/mutation-type-pairs_hiv1.json
pipenv run python scripts/json2csv.py $DATADIR/mutation-type-pairs_hiv2.json
