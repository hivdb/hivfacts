#! /usr/bin/env python3
import os
import json
import yaml
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


def main():
    for dirpath, _, filenames in os.walk(DATADIR):
        for filename in filenames:
            fnl = filename.lower()
            if fnl.endswith('.yml') or fnl.endswith('.yaml'):
                yaml2json(os.path.join(dirpath, filename))


if __name__ == '__main__':
    main()
