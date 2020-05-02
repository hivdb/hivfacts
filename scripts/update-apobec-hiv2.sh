#! /bin/sh

BASEDIR=`realpath $(dirname $0)/..`
mkdir -p $BASEDIR/data/apobecs-hiv2
APOBEC_ALL_JSON_PATH="$BASEDIR/local/apobecs-hiv2/apobecs_all.json"
APOBEC_RESULT_JSON="$BASEDIR/data/apobecs-hiv2/apobecs.json"
APOBEC_DRM_RESULT_JSON="$BASEDIR/data/apobecs-hiv2/apobec_drms.json"
APOBEC_RESULT_CSV="$BASEDIR/data/apobecs-hiv2/apobecs.csv"
APOBEC_DRM_RESULT_CSV="$BASEDIR/data/apobecs-hiv2/apobec_drms.csv"

set -e

cd $BASEDIR

pipenv run python3 $BASEDIR/scripts/apply_apobec_filter.py \
    HIV2 $APOBEC_ALL_JSON_PATH \
    $APOBEC_RESULT_JSON $APOBEC_DRM_RESULT_JSON \
    $APOBEC_RESULT_CSV $APOBEC_DRM_RESULT_CSV
