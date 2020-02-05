#! /usr/bin/env python3
import os
import sys
import json
import yaml

BASEDIR = os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))
)
DATADIR = os.path.join(BASEDIR, 'data')


def json2yaml(jsonfile):
    yamlfile = os.path.splitext(jsonfile)[0] + '.yml'
    with open(jsonfile) as infile, open(yamlfile, 'w') as outfile:
        data = json.load(infile)
        yaml.dump(data, outfile)


def main():
    filename = sys.argv[1]
    json2yaml(filename)


if __name__ == '__main__':
    main()
