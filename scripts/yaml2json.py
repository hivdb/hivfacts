#! /usr/bin/env python3
import os
import sys
import json
import yaml  # type: ignore
from yaml import Loader

BASEDIR = os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))
)
DATADIR = os.path.join(BASEDIR, 'data')


def yaml2json(yamlfile):
    jsonfile = os.path.splitext(yamlfile)[0] + '.json'
    with open(yamlfile) as infile, open(jsonfile, 'w') as outfile:
        data = yaml.load(infile, Loader=Loader)
        json.dump(data, outfile, indent=2)
        print('Write to', outfile.name)


def main():
    if len(sys.argv) > 1:
        for filename in sys.argv[1:]:
            yaml2json(filename)
    else:
        for dirpath, _, filenames in os.walk(DATADIR):
            for filename in filenames:
                fnl = filename.lower()
                if fnl.endswith('.yml') or fnl.endswith('.yaml'):
                    yaml2json(os.path.join(dirpath, filename))


if __name__ == '__main__':
    main()
