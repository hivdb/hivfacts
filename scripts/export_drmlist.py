#! /usr/bin/env python3
import os
import re
import json
import yaml
import click
from yaml import Loader
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
QUERY_TBL_SCORES = (
    'SELECT DrugClass, Gene, Pos, AA '
    'FROM tblScoresWithVersions WHERE Version=%s '
    'GROUP BY DrugClass, Gene, Pos, AA ORDER BY Gene, DrugClass, Pos, AA'
)
QUERY_TBL_COMPSCORES = (
    'SELECT DrugClass, Gene, Rule '
    'FROM tblCombinationScoresWithVersions WHERE Version=%s '
    'Group BY Gene, DrugClass, Rule ORDER BY Gene, DrugClass, Rule'
)
QUERY_HIV2_DRMS = (
    "SELECT DrugClass, Gene, Pos, AAs "
    "FROM tblMutationTypesWithVersions WHERE "
    "Version='V9_0a1' AND Type IN ('Major', 'Accessory')"
    'GROUP BY DrugClass, Gene, Pos, AAs ORDER BY Gene, DrugClass, Pos, AAs'
)
GENE_ORDER = {
    'PR': 0,
    'RT': 1,
    'IN': 2
}


def get_latest_db_version():
    with open(VERSIONS_PATH) as fp:
        data = yaml.load(fp, Loader=Loader)
        ver, *_ = data['HIVDB'][-1]
        return 'V' + ver.split('-', 1)[0].replace('.', '_')


def mutation_sortkey(mut):
    return GENE_ORDER[mut[1]], mut[0], mut[2], mut[3]


@click.command()
@click.argument('species', type=click.Choice(['HIV1', 'HIV2']))
@click.argument('output_json', type=click.File('w'))
def main(species, output_json):
    engine = create_engine(DATABASE_URI)
    engine.connect()
    dbver = get_latest_db_version()
    if species == 'HIV1':
        result_scores = engine.execute(QUERY_TBL_SCORES, dbver)
        drmset = {(
            row['DrugClass'],
            row['Gene'],
            row['Pos'],
            row['AA'].replace('#', '_').replace('~', '-')
        ) for row in result_scores}
        result_compscores = engine.execute(QUERY_TBL_COMPSCORES, dbver)
        for row in result_compscores:
            dc = row['DrugClass']
            gene = row['Gene']
            muts = [re.match(r'^(\d+)([A-Z]+)$', r.strip())
                    for r in row['Rule'].split('+')]
            for m in muts:
                pos = int(m.group(1))
                drmset |= {(
                    dc,
                    gene,
                    pos,
                    aa.replace('#', '_').replace('~', '-')
                ) for aa in m.group(2)}
    elif species == 'HIV2':
        result_muts = engine.execute(QUERY_HIV2_DRMS)
        drmset = set()
        for row in result_muts:
            for aa in row['AAs']:
                drmset.add((
                    row['DrugClass'],
                    row['Gene'],
                    row['Pos'],
                    aa.replace('#', '_').replace('~', '-')
                ))
    drmlist = sorted(drmset, key=mutation_sortkey)
    drmdict = {
        dc: [{'gene': gene,
              'position': pos,
              'aa': aa} for _, gene, pos, aa in drms]
        for dc, drms in groupby(drmlist, key=lambda drm: drm[0])
    }
    json.dump(drmdict, output_json, indent=2)


if __name__ == '__main__':
    main()
