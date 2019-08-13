#! /bin/sh

BASEDIR=`realpath $(dirname $0)/..`
mkdir -p $BASEDIR/data/apobecs
APOBEC_ALL_JSON_PATH="$BASEDIR/local/apobecs/apobecs_all.json"
APOBEC_RESULT_JSON="$BASEDIR/data/apobecs/apobecs.json"
APOBEC_DRM_RESULT_JSON="$BASEDIR/data/apobecs/apobec_drms.json"
APOBEC_RESULT_CSV="$BASEDIR/data/apobecs/apobecs.csv"
APOBEC_DRM_RESULT_CSV="$BASEDIR/data/apobecs/apobec_drms.csv"

set -e


python3 $BASEDIR/scripts/apply_apobec_filter.py \
    $APOBEC_ALL_JSON_PATH \
    $APOBEC_RESULT_JSON $APOBEC_DRM_RESULT_JSON \
    $APOBEC_RESULT_CSV $APOBEC_DRM_RESULT_CSV
