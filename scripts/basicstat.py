#! /usr/bin/env python
import sys
import csv
from collections import defaultdict

from hivfacts import hivaapcnt

GENES = ['PR', 'RT', 'IN']


def main():
    if len(sys.argv) != 2:
        print('Usage: {} <OUTPUT_CSV>'.format(sys.argv[0]), file=sys.stderr)
        exit(1)
    stat = defaultdict(dict)
    aapcnt_all_all = hivaapcnt.HIVAAPcnt('all', 'all')
    aapcnt_all_b = hivaapcnt.HIVAAPcnt('all', 'B')
    aapcnt_art_all = hivaapcnt.HIVAAPcnt('art', 'all')
    for gene in GENES:
        stat['# Total'][gene] = max(one['total']
                                    for one in aapcnt_all_all.get(gene))
        stat['# B'][gene] = max(one['total'] for one in aapcnt_all_b.get(gene))
        stat['# Non-B'][gene] = stat['# Total'][gene] - stat['# B'][gene]
        stat['% Treated'][gene] = (
            max(one['total'] for one in aapcnt_art_all.get(gene)) /
            stat['# Total'][gene])
    with open(sys.argv[1], 'w') as outfp:
        writer = csv.writer(outfp)
        writer.writerow(['', 'PR', 'RT', 'IN'])
        for key, val in stat.items():
            writer.writerow([key, *[val[g] for g in GENES]])


if __name__ == '__main__':
    main()
