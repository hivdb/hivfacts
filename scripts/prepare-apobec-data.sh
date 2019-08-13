#! /bin/sh

BASEDIR=`realpath $(dirname $0)/..`
mkdir -p $BASEDIR/local/apobecs
CODONPCNT_PATH="$BASEDIR/local/apobecs/codonpcnt_for_apobecs.json"
STOPS_PLASMA_PATH="$BASEDIR/local/apobecs/cooccur_stops_plasma.json"
STOPS_PBMC_PATH="$BASEDIR/local/apobecs/cooccur_stops_pbmc.json"
STOPS_ALLSRC_PATH="$BASEDIR/local/apobecs/cooccur_stops_all_sources.json"
APOBEC_JSON_PATH="$BASEDIR/local/apobecs/apobecs_all.json"
APOBEC_CSV_PATH="$BASEDIR/local/apobecs/apobecs_all.csv"

set -e

hivdbql export-codonpcnt \
    --species HIV1 \
    --subtype all \
    --rx-type all \
    --format json \
    --include-by-subtypes \
    --filter PLASMA_ONLY \
    --filter NO_DNACHIP \
    --filter NO_QA_ISSUES \
    --filter NO_CLONES \
    $CODONPCNT_PATH

hivdbql export-cooccur-stops \
    --species HIV1 \
    --subtype all \
    --rx-type all \
    --format json \
    --filter PBMC_ONLY \
    --filter NO_CLONES \
    $STOPS_PBMC_PATH

hivdbql export-cooccur-stops \
    --species HIV1 \
    --subtype all \
    --rx-type all \
    --format json \
    --filter PLASMA_ONLY \
    --filter NO_CLONES \
    $STOPS_PLASMA_PATH

hivdbql export-cooccur-stops \
    --species HIV1 \
    --subtype all \
    --rx-type all \
    --format json \
    --filter NO_CLONES \
    $STOPS_ALLSRC_PATH

hivdbql discover-apobecs \
    --codonpcnt-json $CODONPCNT_PATH \
    --cooccur-stops-plasma-json $STOPS_PLASMA_PATH \
    --cooccur-stops-pbmc-json $STOPS_PBMC_PATH \
    --cooccur-stops-allsrc-json $STOPS_ALLSRC_PATH \
    --threshold 0.98 \
    --subtype-threshold 0.9 \
    --format json \
    $APOBEC_JSON_PATH

hivdbql discover-apobecs \
    --codonpcnt-json $CODONPCNT_PATH \
    --cooccur-stops-plasma-json $STOPS_PLASMA_PATH \
    --cooccur-stops-pbmc-json $STOPS_PBMC_PATH \
    --cooccur-stops-allsrc-json $STOPS_ALLSRC_PATH \
    --threshold 0.98 \
    --subtype-threshold 0.9 \
    --format csv \
    $APOBEC_CSV_PATH
