#! /bin/bash

set -e

cd `dirname $0`/..
if [ -z "$1" ]; then
    echo "Usage: $0 <VERSION>" >&2
    exit 1
fi

echo $1 > VERSION
sed -i "" "s/^version = '.*'$/version = '$1'/g" hiv-aapcnt-java/build.gradle
sed -i "" "s/^version = '.*'$/version = '$1'/g" hiv-aapcnt-python/setup.py
sed -i "" "s/^version = '.*'$/version = '$1'/g" hiv-aapcnt-python/hivaapcnt/__init__.py
cat VERSION
