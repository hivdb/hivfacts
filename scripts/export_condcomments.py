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
QUERY_TBL_CMTS_HIV1 = (
    'SELECT CommentName, DrugClass, ConditionType, ConditionValue, Comment '
    'FROM tblConditionalCommentsWithVersions WHERE Version=%s '
    'ORDER BY CommentName'
)

QUERY_TBL_CMTS_HIV2 = (
    'SELECT CommentName, DrugClass, ConditionType, ConditionValue, Comment '
    "FROM tblConditionalCommentsWithVersions WHERE Version='V9_0a3' "
    'ORDER BY CommentName'
)


def get_latest_db_version():
    with open(VERSIONS_PATH) as fp:
        data = yaml.load(fp, Loader=Loader)
        ver, *_ = data['HIVDB'][-1]
        return 'V' + ver.split('-', 1)[0].replace('.', '_')


@click.command()
@click.argument('species', type=click.Choice(['HIV1', 'HIV2']))
@click.argument('output_json', type=click.File('w'))
def main(species, output_json):
    engine = create_engine(DATABASE_URI)
    engine.connect()
    dbver = get_latest_db_version()
    if species == 'HIV1':
        results = engine.execute(QUERY_TBL_CMTS_HIV1, dbver)
    else:  # species == 'HIV2'
        results = engine.execute(QUERY_TBL_CMTS_HIV2)
    cmtlist = []
    for row in results:
        condval = json.loads(row['ConditionValue'])
        cmtlist.append({
            'strain': condval.get('strain', 'HIV1'),
            'commentName': row['CommentName'],
            'drugClass': row['DrugClass'],
            'conditionType': row['ConditionType'],
            'conditionValue': condval,
            'comment': row['Comment']
        })
    json.dump(cmtlist, output_json, indent=2)


if __name__ == '__main__':
    main()
