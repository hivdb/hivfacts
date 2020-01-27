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
QUERY_TBL_SCORE = (
    "SELECT Gene, DrugClass, Pos, AA, Drug, Score "
    "FROM tblScoresWithVersions WHERE Version=%s "
    "ORDER BY Gene, Pos, AA, Drug;"
)


def get_latest_db_version():
    with open(VERSIONS_PATH) as fp:
        data = yaml.load(fp, Loader=Loader)
        ver, *_ = data['HIVDB'][-1]
        return 'V' + ver.split('-', 1)[0].replace('.', '_')


def to_internalformat(aa):
    aa = aa.replace(
        'i', '_').replace(
        '#', '_').replace(
        'd', '-').replace(
        '~', '-')

    return aa


@click.command()
@click.argument('output_json', type=click.File('w'))
def main(output_json):
    engine = create_engine(DATABASE_URI)
    engine.connect()
    dbver = get_latest_db_version()
    results = engine.execute(QUERY_TBL_SCORE, dbver)
    cmtlist = []
    for row in results:
        cmtlist.append({
            'gene': row['Gene'],
            'drugClass': row['DrugClass'],
            'pos': row['Pos'],
            'aa': to_internalformat(row['AA']),
            'drug': row['Drug'],
            'score': row['Score'],
        })
    json.dump(cmtlist, output_json, indent=2)


if __name__ == '__main__':
    main()
