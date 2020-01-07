#! /usr/bin/env python3
import os
import json
import click
from itertools import groupby
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
QUERY_TBL_SDRMS = (
    'SELECT Class AS DrugClass, Gene, Pos, AAs '
    'FROM tblSDRMs '
    'ORDER BY Gene, Pos, AAs'
)
GENE_ORDER = {
    'PR': 0,
    'RT': 1,
    'IN': 2
}


def mutation_sortkey(mut):
    return GENE_ORDER[mut['Gene']], mut['DrugClass'], mut['Pos'], mut['AAs']


@click.command()
@click.argument('output_json', type=click.File('w'))
def main(output_json):
    engine = create_engine(DATABASE_URI)
    engine.connect()
    results = engine.execute(QUERY_TBL_SDRMS)
    sdrmlist = sorted(results, key=mutation_sortkey)
    sdrmdict = {
        dc: [{
            'gene': gene,
            'position': pos,
            'aa': aa.replace('#', '_').replace('~', '-')
        } for _, gene, pos, aa in sdrms]
        for dc, sdrms in groupby(sdrmlist, key=lambda m: m['DrugClass'])
    }
    json.dump(sdrmdict, output_json, indent=2)


if __name__ == '__main__':
    main()
