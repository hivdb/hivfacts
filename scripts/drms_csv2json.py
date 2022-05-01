#! /usr/bin/env python3

import os
import sys
import csv
import json
from itertools import groupby


def main():
    csvfile = sys.argv[1]
    payload = {}
    with open(csvfile, encoding='utf-8-sig') as fp, \
            open(
                os.path.splitext(csvfile)[0] + '.json',
                'w'
            ) as out:
        rows = csv.DictReader(fp)
        for dc, dcrows in groupby(rows, lambda r: r['drug_class']):
            payload[dc] = list(dcrows)
            for row in payload[dc]:
                row.pop('drug_class')
                row['position'] = int(row['position'])
        json.dump(payload, out, indent=2)
        print('Write to', out.name)


if __name__ == '__main__':
    main()
