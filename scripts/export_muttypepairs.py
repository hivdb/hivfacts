#! /usr/bin/env python3
import os
import json
import yaml
import click
from yaml import Loader
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
QUERY_HIV1_MUTTYPES = (
    "SELECT Strain, Gene, DrugClass, Pos, AAs, Type, IsUnusual "
    "FROM tblMutationTypesWithVersions WHERE "
    "Strain='HIV1' AND Version=%s "
    "ORDER BY Gene, DrugClass, Pos, "
    "(CASE Type WHEN 'Major' THEN 0"
    " WHEN 'Accessory' THEN 1"
    " WHEN 'NRTI' THEN 0"
    " WHEN 'NNRTI' THEN 1"
    " ELSE 3 END), AAs"
)
QUERY_HIV2_MUTTYPES = (
    "SELECT Strain, Gene, DrugClass, Pos, AAs, Type, IsUnusual "
    "FROM tblMutationTypesWithVersions WHERE "
    "Version='V9_0a3' "
    "ORDER BY Gene, DrugClass, Pos, "
    "(CASE Type WHEN 'Major' THEN 0"
    " WHEN 'Accessory' THEN 1"
    " WHEN 'NRTI' THEN 0"
    " WHEN 'NNRTI' THEN 1"
    " ELSE 3 END), AAs"
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
    return GENE_ORDER[mut[0]], mut[1], mut[2]


@click.command()
@click.argument('species', type=click.Choice(['HIV1', 'HIV2']))
@click.argument('output_json', type=click.File('w'))
def main(species, output_json):
    engine = create_engine(DATABASE_URI)
    engine.connect()
    dbver = get_latest_db_version()
    if species == 'HIV1':
        results = engine.execute(QUERY_HIV1_MUTTYPES, dbver)
        muttypepairs = [{
            'strain': row['Strain'],
            'gene': row['Gene'],
            'drugClass': row['DrugClass'],
            'position': row['Pos'],
            'aas': row['AAs'].replace('#', '_').replace('~', '-'),
            'mutationType': row['Type'],
            'isUnusual': row['IsUnusual'] > 0
        } for row in results]
    else:  # species == 'HIV2'
        results = engine.execute(QUERY_HIV2_MUTTYPES)
        muttypepairs = [{
            'strain': row['Strain'],
            'gene': row['Gene'],
            'drugClass': row['DrugClass'],
            'position': row['Pos'],
            'aas': row['AAs'].replace('#', '_').replace('~', '-'),
            'mutationType': row['Type'],
            'isUnusual': row['IsUnusual'] > 0
        } for row in results]
    json.dump(muttypepairs, output_json, indent=2)


if __name__ == '__main__':
    main()
