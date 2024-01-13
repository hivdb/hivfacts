#! /usr/bin/env python3
import re
import csv
import json
import click  # type: ignore

from typing import TextIO, Dict, List, Any, Optional

DRUG_NAMES: Dict[str, str] = {
    '3TC': 'LMV',
    'ATV/r': 'ATV',
    'DRV/r': 'DRV',
    'FPV/r': 'FPV',
    'IDV/r': 'IDV',
    'LPV/r': 'LPV',
    'SQV/r': 'SQV',
    'TPV/r': 'TPV'
}


def parseCondition(row: Dict[str, str]) -> Dict[str, Any]:
    cond_type = row['conditionType']
    cond_text = row['condition']
    if cond_type == 'MUTATION':
        conds: List[Dict[str, Any]] = []
        for mut in cond_text.split(' OR '):
            m: Optional[re.Match] = re.match(r'^(\d+)([^\d]+)$', mut)
            if not m:
                raise RuntimeError('invalid mutation {}'.format(mut))
            conds.append({
                'gene': row['gene'],
                'pos': int(m.group(1)),
                'aas': m.group(2).replace('i', '_').replace('d', '-')
            })
        if len(conds) == 1:
            return conds[0]
        else:
            return {'or': conds}
    else:  # if cond_type == 'DRUGLEVEL':
        conds = []
        for druglevel in cond_text.split(' AND '):
            drug, level = druglevel.split('=', 1)
            conds.append({
                'drug': DRUG_NAMES.get(drug, drug),
                'levels': [level]
            })
        if len(conds) == 1:
            return conds[0]
        else:
            return {'and': conds}


@click.command()
@click.argument('species', type=click.Choice(['HIV1', 'HIV2']))
@click.argument('input_csv', type=click.File('r'))
@click.argument('output_json', type=click.File('w'))
def main(species: str, input_csv: TextIO, output_json: TextIO) -> None:
    reader: csv.DictReader = csv.DictReader(input_csv)
    cmtlist = []
    for row in reader:
        cmtid = row['id']
        if species == 'HIV1':
            strain = 'HIV1'
        elif cmtid.startswith('HIV2A'):
            strain = 'HIV2A'
        elif cmtid.startswith('HIV2B'):
            strain = 'HIV2B'
        else:
            raise RuntimeError('unknown strain for {}'.format(cmtid))
        cmtlist.append({
            'strain': strain,
            'commentName': row['id'],
            'drugClass': row['drugClass'],
            'conditionType': row['conditionType'],
            'conditionValue': parseCondition(row),
            'comment': row['text']
        })
    json.dump(cmtlist, output_json, indent=2)
    click.echo(output_json.name)


if __name__ == '__main__':
    main()
