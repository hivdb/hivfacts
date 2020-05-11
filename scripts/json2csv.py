#! /usr/bin/env python3
import os
import sys
import json
import csv

BASEDIR = os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))
)
DATADIR = os.path.join(BASEDIR, 'data')


def guess_fieldname(values):
    fieldname = set()
    for value in values:
        if value in ('PI', 'NRTI', 'NNRTI', 'INSTI'):
            fieldname.add('drug_class')
        elif value in ('PR', 'RT', 'IN'):
            fieldname.add('gene')
    if len(fieldname) == 1:
        return fieldname.pop()
    raise RuntimeError('Unable to guess field name for {!r}'.format(value))


def get_fieldnames(data):
    if isinstance(data, list) and data:
        return list(data[0].keys())
    elif isinstance(data, dict) and data:
        fieldnames = [guess_fieldname(data.keys())]
        for subset in data.values():
            if isinstance(subset, list) and subset:
                fieldnames.extend(subset[0].keys())
                return fieldnames
    raise RuntimeError('Unable to convert JSON data to CSV')


def dict2list(data, fieldnames):
    firstfield = fieldnames[0]
    results = []
    for key, subset in data.items():
        for value in subset:
            results.append({firstfield: key, **value})
    return results


def json2csv(jsonfile):
    csvfile = os.path.splitext(jsonfile)[0] + '.csv'
    with open(jsonfile) as infile, open(
        csvfile, 'w', encoding='utf-8-sig'
    ) as outfile:
        data = json.load(infile)
        fieldnames = get_fieldnames(data)
        if isinstance(data, dict):
            data = dict2list(data, fieldnames)

        csvwriter = csv.DictWriter(outfile, fieldnames=fieldnames)
        csvwriter.writeheader()

        for record in data:
            csvwriter.writerow(record)


def main():
    filename = sys.argv[1]
    json2csv(filename)


if __name__ == '__main__':
    main()
