#! /usr/bin/env python3

import os
import sys
import csv
import json


def main():
    csvfile = sys.argv[1]
    with open(csvfile, encoding='utf-8-sig') as fp, \
            open(
                os.path.splitext(csvfile)[0] + '.json',
                'w'
            ) as out:
        rows = list(csv.DictReader(fp))
        for row in rows:
            if 'position' in row:
                row['position'] = int(row['position'])
            if 'percent' in row:
                row['percent'] = float(row['percent'])
            if 'count' in row:
                row['count'] = int(row['count'])
            if 'total' in row:
                row['total'] = int(row['total'])
            if 'isUnusual' in row:
                row['isUnusual'] = row['isUnusual'].lower() == 'true'
        json.dump(rows, out, indent=2)
        print('Write to', out.name)


if __name__ == '__main__':
    main()
