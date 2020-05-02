#! /bin/sh

BASEDIR=`realpath $(dirname $0)/..`
mkdir -p $BASEDIR/local/apobecs-hiv2
AAPCNT_PATH="$BASEDIR/data/aapcnt-hiv2/rx-all_subtype-all.json"
CODONPCNT_PATH="$BASEDIR/local/apobecs-hiv2/codonpcnt_for_apobecs.json"
STOPS_PLASMA_PATH="$BASEDIR/local/apobecs-hiv2/cooccur_stops_plasma.json"
STOPS_PBMC_PATH="$BASEDIR/local/apobecs-hiv2/cooccur_stops_pbmc.json"
STOPS_ALLSRC_PATH="$BASEDIR/local/apobecs-hiv2/cooccur_stops_all_sources.json"
APOBEC_JSON_PATH="$BASEDIR/local/apobecs-hiv2/apobecs_all.json"
APOBEC_CSV_PATH="$BASEDIR/local/apobecs-hiv2/apobecs_all.csv"

set -e
cd $BASEDIR

pipenv run hivdbql export-codonpcnt \
    --species HIV2 \
    --subtype all \
    --rx-type all \
    --format json \
    --include-by-subtypes \
    --filter NO_QA_ISSUES \
    $CODONPCNT_PATH

pipenv run hivdbql export-cooccur-stops \
    --species HIV2 \
    --subtype all \
    --rx-type all \
    --format json \
    --filter PBMC_ONLY \
    $STOPS_PBMC_PATH

pipenv run hivdbql export-cooccur-stops \
    --species HIV2 \
    --subtype all \
    --rx-type all \
    --format json \
    --filter PLASMA_ONLY \
    $STOPS_PLASMA_PATH

pipenv run hivdbql export-cooccur-stops \
    --species HIV2 \
    --subtype all \
    --rx-type all \
    --format json \
    --no-filter \
    $STOPS_ALLSRC_PATH

pipenv run hivdbql discover-apobecs \
    --species HIV2 \
    --aapcnt-json $AAPCNT_PATH \
    --codonpcnt-json $CODONPCNT_PATH \
    --cooccur-stops-plasma-json $STOPS_PLASMA_PATH \
    --cooccur-stops-pbmc-json $STOPS_PBMC_PATH \
    --cooccur-stops-allsrc-json $STOPS_ALLSRC_PATH \
    --threshold 0.95 \
    --subtype-threshold 0.9 \
    --format json \
    $APOBEC_JSON_PATH

pipenv run hivdbql discover-apobecs \
    --species HIV2 \
    --aapcnt-json $AAPCNT_PATH \
    --codonpcnt-json $CODONPCNT_PATH \
    --cooccur-stops-plasma-json $STOPS_PLASMA_PATH \
    --cooccur-stops-pbmc-json $STOPS_PBMC_PATH \
    --cooccur-stops-allsrc-json $STOPS_ALLSRC_PATH \
    --threshold 0.95 \
    --subtype-threshold 0.9 \
    --format csv \
    $APOBEC_CSV_PATH
