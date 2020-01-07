#! /usr/bin/env python3
import os
import json
import click
from collections import defaultdict
from sqlalchemy import create_engine

BASEDIR = os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))
)
VERSIONS_PATH = os.path.join(
    BASEDIR, 'data', 'algorithms', 'versions.yml')

DATABASE_URI = os.environ.get(
    'DATABASE_URI_HIVDBRULES',
    'mysql+pymysql://rshafer:rshafer@10.77.6.244/HIVDB_Results'
)
QUERY_TBL_TSMS = (
    'SELECT DrugClass, Gene, Pos, AA '
    'FROM tblTSMs '
    "WHERE Include = 'Y' "
    'ORDER BY Gene, Pos, AA'
)
GENE_ORDER = {
    'PR': 0,
    'RT': 1,
    'IN': 2
}


def mutation_sortkey(mut):
    return GENE_ORDER[mut[1]], mut[0], mut[2], mut[3]


@click.command()
@click.argument('output_json', type=click.File('w'))
def main(output_json):
    engine = create_engine(DATABASE_URI)
    engine.connect()
    results = engine.execute(QUERY_TBL_TSMS)
    tsmlists = defaultdict(list)
    for row in sorted(results, key=mutation_sortkey):
        tsmlists[row['DrugClass']].append({
            'gene': row['Gene'],
            'position': row['Pos'],
            'aa': row['AA'].replace('#', '_').replace('~', '-')
        })
    json.dump(tsmlists, output_json, indent=2)


if __name__ == '__main__':
    main()
