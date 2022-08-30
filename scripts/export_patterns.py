#! /usr/bin/env python3
import os
import csv
# import json
from itertools import groupby
from collections import Counter
import click  # type: ignore
from sqlalchemy import create_engine  # type: ignore

BASEDIR = os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))
)
VERSIONS_PATH = os.path.join(
    BASEDIR, 'data', 'algorithms', 'versions.yml')

DATABASE_URI = os.environ.get(
    'DATABASE_URI_HIVDBRULES',
    'mysql+pymysql://rshafer:rshafer@10.77.6.244/HIVDB2'
)
QUERY_TBL_MUTATIONS = """
    SELECT
      M.SequenceID AS sequence_id,
      Gene AS gene,
      CodonPos AS position,
      CASE
        WHEN Insertion = 'Yes' THEN '_'
        WHEN MutAA = '~' THEN '-'
        WHEN MutAA = '*' THEN '*'
        ELSE MutAA
      END AS amino_acid,
      CONCAT(
        Consensus,
        CodonPos,
        CASE
          WHEN Insertion = 'Yes' THEN '_'
          WHEN MutAA = '~' THEN '-'
          WHEN MutAA = '*' THEN '*'
          ELSE MutAA
        END
      ) AS mutation,
      CONCAT(
        Gene,
        ':',
        CodonPos,
        CASE
          WHEN Insertion = 'Yes' THEN '_'
          WHEN MutAA = '~' THEN '-'
          WHEN MutAA = '*' THEN '*'
          ELSE MutAA
        END
      ) AS mutation_key
    FROM _Mutations M, tblSequences S, tblIsolates I, tblSpecies SP
    WHERE
      M.SequenceID = S.SequenceID AND
      S.IsolateID = I.IsolateID AND
      I.IsolateID = SP.IsolateID AND
      SP.Species = 'HIV1' AND
      ({})
"""
DRUG_CLASS_ORDER = {
    'PI': 0,
    'NRTI': 1,
    'NNRTI': 2,
    'INSTI': 3
}


@click.command()
@click.argument('output_csv', type=click.File('w', encoding='utf-8-sig'))
def main(output_csv):
    drm_positions = {}
    drm_lookup = {}
    with open(os.path.join(BASEDIR, 'data', 'drms_hiv1.csv'),
              encoding='utf-8-sig') as fp:
        for row in csv.DictReader(fp):
            drm_positions.setdefault(row['gene'], set())
            drm_positions[row['gene']].add(row['position'])
            drm_lookup[
                '{gene}:{position}{aa}'.format(**row)
            ] = row['drug_class']
    engine = create_engine(DATABASE_URI)
    conn = engine.connect()
    results = (
        conn
        .execution_options(stream_results=True)
        .engine.execute(
            QUERY_TBL_MUTATIONS
            .format(' OR '.join([
                '(Gene = {!r} AND CodonPos IN ({}))'
                .format(gene, ', '.join(pos))
                for gene, pos in drm_positions.items()
            ]))
        )
    )
    drm_results = sorted(
        ({
            'drug_class': drm_lookup[row.mutation_key],
            **row
        } for row in results if row.mutation_key in drm_lookup),
        key=lambda row: (
            row['sequence_id'],
            DRUG_CLASS_ORDER[row['drug_class']],
            row['position'],
            row['amino_acid']
        )
    )
    patterns = Counter(
        (gene, drug_class, '+'.join([mut['mutation'] for mut in muts]))
        for (_, gene, drug_class), muts in groupby(
            drm_results,
            key=lambda row: (row['sequence_id'],
                             row['gene'],
                             row['drug_class'])
        )
    )

    writer = csv.writer(output_csv)
    writer.writerow(['gene', 'drug_class', 'pattern', 'count'])
    for row, count in sorted(
        patterns.items(),
        key=lambda item: (DRUG_CLASS_ORDER[item[0][1]], -item[1])
    ):
        writer.writerow([*row, count])
    click.echo('Write to {}'.format(output_csv.name))


if __name__ == '__main__':
    main()
