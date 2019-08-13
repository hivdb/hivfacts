#! /bin/bash

set -e

cd `dirname $0`/..
if [ -z "$1" ]; then
    echo "Usage: $0 <VERSION>" >&2
    exit 1
fi

echo $1 > VERSION
sed -i "" "s/^version = '.*'$/version = '$1'/g" hivfacts-java/build.gradle
sed -i "" "s/^version = '.*'$/version = '$1'/g" hivfacts-python/setup.py
sed -i "" "s/^version = '.*'$/version = '$1'/g" hivfacts-python/hivfacts/__init__.py
cat VERSION
