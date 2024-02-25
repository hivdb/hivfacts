from typing import Literal

import os
import re
import csv
import json
import click
from lxml import etree

MUTATION_PATTERN = re.compile(r'[A-Z]?(\d+)([A-Zid]+)')
ASI_XML = {
    'HIV1': 'HIVDB_latest.xml',
    'HIV2': 'HIVDB-HIV2_9.0.xml'
}
DRUGS_JSON = 'drugs.json'
DRUGCLASSES_JSON = {
    'HIV1': 'drug-classes_hiv1.json',
    'HIV2': 'drug-classes_hiv2.json'
}
GENES_JSON = {
    'HIV1': 'genes_hiv1.json',
    'HIV2': 'genes_hiv2.json'
}

DRMS_CSV = {
    'HIV1': 'drms_hiv1.csv',
    'HIV2': 'drms_hiv2.csv'
}


def read_drug_class_lookup(data_dir: str) -> dict[str, str]:
    with open(os.path.join(data_dir, DRUGS_JSON), 'r') as drugs_json:
        payload = json.load(drugs_json)
    return {drug['displayAbbr']: drug['drugClass'] for drug in payload}


def read_gene_lookup(
    virus: Literal['HIV1', 'HIV2'],
    data_dir: str
) -> dict[str, str]:
    with open(os.path.join(
        data_dir, DRUGCLASSES_JSON[virus]
    ), 'r') as drugclasses_json:
        payload = json.load(drugclasses_json)
    return {dc['name']: dc['abstractGene'] for dc in payload}


def sort_drms(
    drms: set[tuple[str, str, int, str]],
    virus: Literal['HIV1', 'HIV2'],
    data_dir: str
) -> list[tuple[str, str, int, str]]:
    ordered_drug_classes = list(read_gene_lookup(virus, data_dir).keys())
    return sorted(
        drms,
        key=lambda x: (ordered_drug_classes.index(x[0]), x[2], x[3])
    )


def write_drms_csv(
    virus: Literal['HIV1', 'HIV2'],
    sorted_drms: list[tuple[str, str, int, str]],
    data_dir: str
) -> None:
    with open(os.path.join(data_dir, DRMS_CSV[virus]), 'w') as drms_csv:
        writer = csv.writer(drms_csv)
        writer.writerow(['drug_class', 'gene', 'position', 'aa'])
        writer.writerows(sorted_drms)
        print(f'wrtie to {drms_csv.name}')


def extract_drms(
    asi: etree._Element,
    virus: Literal['HIV1', 'HIV2'],
    data_dir: str
) -> set[tuple[str, str, int, str]]:
    drug_lookup = read_drug_class_lookup(data_dir)
    gene_lookup = read_gene_lookup(virus, data_dir)
    unique_drms = set()
    drugs = asi.xpath('//DRUG')
    for drug in drugs:
        drug_name = drug.find('NAME').text
        drug_class = drug_lookup[drug_name]
        gene = gene_lookup[drug_class]
        conditions = drug.xpath('.//CONDITION')
        for cond in conditions:
            cond_text = cond.text
            mutations = MUTATION_PATTERN.findall(cond_text)
            for pos, aas in mutations:
                for aa in aas:
                    if aa == 'i':
                        aa = '_'
                    elif aa == 'd':
                        aa = '-'
                    unique_drms.add((drug_class, gene, int(pos), aa))
    return unique_drms


@click.command()
@click.argument('virus', type=click.Choice(['HIV1', 'HIV2']))
@click.argument('data_dir', type=click.Path(exists=True, file_okay=False))
def main(virus: Literal['HIV1', 'HIV2'], data_dir: str) -> None:
    with open(os.path.join(
        data_dir, 'algorithms',
        ASI_XML[virus]
    ), 'rb') as asi_xml:
        asi = etree.parse(asi_xml)
        unique_drms = extract_drms(asi, virus, data_dir)
    sorted_drms = sort_drms(unique_drms, virus, data_dir)
    write_drms_csv(virus, sorted_drms, data_dir)


if __name__ == '__main__':
    main()
