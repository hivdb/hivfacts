#! /usr/bin/env python3
import os
import sys
import json
import csv

BASEDIR = os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))
)
DATADIR = os.path.join(BASEDIR, 'data')


def json2csv(jsonfile):
    csvfile = os.path.splitext(jsonfile)[0] + '.csv'
    with open(jsonfile) as infile, open(
        csvfile, 'w', encoding='utf-8-sig') as outfile:

        data = json.load(infile)
        fieldnames = data[0].keys()

        csvwriter = csv.DictWriter(outfile, fieldnames=fieldnames)
        csvwriter.writeheader()

        for record in data:
            csvwriter.writerow(record)


def main():
    filename = sys.argv[1]
    json2csv(filename)


if __name__ == '__main__':
    main()
