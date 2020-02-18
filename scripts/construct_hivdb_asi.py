#! /usr/bin/env python3
from pathlib import Path
from collections import defaultdict
from collections import OrderedDict
from functools import lru_cache
import yaml
from yaml import Loader
import json
import re

from lxml import etree
from lxml.builder import E

import click


# File configuration & ASI file version
BASEDIR = Path(__file__).parent.parent.absolute()
DATADIR = BASEDIR / 'data'
VERSIONS_PATHS = {
    'HIV1': BASEDIR / 'data' / 'algorithms' / 'versions.yml',
    'HIV2': BASEDIR / 'data' / 'algorithms' / 'versions_hiv2.yml'
}


@lru_cache()
def get_versions(species):
    with VERSIONS_PATHS[species].open() as fp:
        versions = yaml.load(fp, Loader=Loader)

        return versions


_VERSION = None


def validate_version(ctx, param, value):
    species = ctx.params['species']
    all_versions = get_versions(species)['HIVDB']
    if value == 'latest':
        version = all_versions[-1]
    else:
        match = re.match(r'^((?:\d+[.\-_]){1,2}\d+(?:[-ap]\d+)?)$', value)
        if not match:
            raise click.ClickException(
                'Malformed version text {!r}'.format(value))
        search_version = match.group(1)

        for item in all_versions:
            if item[0] == search_version:
                version = item
                break
        else:
            raise click.ClickException(
                ('Version {} for {} is not defined. '
                 'Available Choice(s): {}')
                .format(
                    value,
                    species,
                    ', '.join(v for v, _, _ in reversed(all_versions))
                )
            )

    return dict(zip(
        ['version', 'date', 'species'],
        version
        ))


# Load data
class DataSource(type):

    def __init__(self, name, bases, attrs):
        filenames = attrs.get('FILENAMES')
        self.DATA = {}
        if not filenames:
            return
        else:
            for species, filename in filenames.items():
                filepath = DATADIR / filename

                with open(filepath) as fd:
                    self.DATA[species] = json.load(fd)


class DrugClassData(metaclass=DataSource):

    FILENAMES = {
        'HIV1': 'drug-classes_hiv1.json',
        'HIV2': 'drug-classes_hiv2.json'
    }

    @classmethod
    def gene_definition(cls, species):
        _gene_definition = defaultdict(list)
        for drugclass in cls.DATA[species]:
            gene = drugclass['abstractGene']
            _gene_definition[gene].append(drugclass['name'])

        return _gene_definition


class DrugData(metaclass=DataSource):

    FILENAMES = {
        '*': 'drugs.json'
    }

    DRUGCLASS_ORDER = {
        'HIV1': [
            'NRTI',
            'NNRTI',
            'PI',
            'INSTI'
        ],
        'HIV2': [
            'NRTI',
            'PI',
            'INSTI'
        ]
    }

    @classmethod
    def drugclass(cls, species):
        _drugclass = defaultdict(list)
        all_drugclasses = cls.DRUGCLASS_ORDER[species]

        for drug in cls.DATA['*']:
            drugclass = drug['drugClass']
            if drugclass not in all_drugclasses:
                continue
            _drugclass[drugclass].append(drug)

        return OrderedDict(sorted(
                _drugclass.items(),
                key=lambda i: all_drugclasses.index(i[0]))
            )

    @classmethod
    def get_display_abbr(cls, species, drugname):
        for item in cls.DATA['*']:
            if item['name'] == drugname:
                return item['displayAbbr']


class ConditionalCommentData(metaclass=DataSource):

    FILENAMES = {
        'HIV1': 'conditional-comments_hiv1.json',
        'HIV2': 'conditional-comments_hiv2.json'
    }

    @classmethod
    def get_comments_by(cls, species, condtype):
        _data = [
            item for item in cls.DATA[species]
            if item['conditionType'] == condtype
        ]
        if condtype == 'DRUGLEVEL':
            comments = _data
        else:
            comments = defaultdict(list)
            for comment in _data:
                gene = comment['conditionValue']['gene']
                comments[gene].append(comment)

        return comments


class MutationComments(metaclass=DataSource):

    @classmethod
    def get_by_drugname(cls, species, drugname):
        rules = cls._filtered_rules = []

        for rule in cls.DATA[species]:
            if rule['drug'] == drugname:
                rules.append(rule)

        return cls

    @classmethod
    def group_by_position(cls):
        rules = defaultdict(list)

        for rule in cls._filtered_rules:
            pos = rule['pos']
            rules[pos].append(rule)

        return rules


class MutationIndivComments(MutationComments):

    FILENAMES = {
        'HIV1': 'mutation-scores-indiv.json',
        'HIV2': 'mutation-scores-indiv_hiv2.json'
    }


class MutationCombiComments(MutationComments):

    FILENAMES = {
        'HIV1': 'mutation-scores-combi.json',
        'HIV2': 'mutation-scores-combi_hiv2.json'
    }
    MUTATION_PATTERN = r'(\d+)([A-Z_-]+)'

    @classmethod
    def group_and_order_by_position(cls):
        rules = defaultdict(list)

        for item in cls._filtered_rules:
            rule_list = item['rule'].split('+')
            positions = []
            for mut in rule_list:
                position = re.match(cls.MUTATION_PATTERN, mut).group(1)
                positions.append(position)

            rules[','.join(positions)].append(item)

        single_rule = []
        multi_rule = []

        for key, item in rules.items():
            if len(item) == 1:
                single_rule.append([key, item])
            else:
                multi_rule.append([key, item])

        return single_rule + multi_rule


# Algorithm constants
ALG_NAME = "HIVDB"


LEVEL_DEFINITION = [
    ["Susceptible", "S", '1'],
    ["Potential Low-Level Resistance", "S", '2'],
    ["Low-Level Resistance", "I", '3'],
    ["Intermediate Resistance", "I", '4'],
    ["High-Level Resistance", "R", '5'],
]

GLOBALRANGE_CONTENTS = (
    "(-INF TO 9 => 1,  "
    "10 TO 14 => 2,  "
    "15 TO 29 => 3,  "
    "30 TO 59 => 4,  "
    "60 TO INF => 5)"
)
DRUG_RULE_INDENTATION = ' ' * 25


# Util function
def to_asi_format(aa):
    aa = aa.replace(
        '#', 'i').replace(
        '_', 'i').replace(
        '~', 'd').replace(
        '-', 'd')

    return aa


# ASI file parts
def doc_type(name, pubid, system):
    if pubid:
        return '<!DOCTYPE {name} {type} "{pubid}" "{system}">'.format(
            name=name,
            type='PUBLIC',
            pubid=pubid,
            system=system
        )
    else:
        return '<!DOCTYPE {name} {type} "{system}">'.format(
            name=name,
            type='SYSTEM',
            system=system
        )


# Definitions
def gene_definition(species):
    nodes = []
    for gene, drugs in DrugClassData.gene_definition(species).items():
        nodes.append(
            E.GENE_DEFINITION(
                E.NAME(gene),
                E.DRUGCLASSLIST(','.join(drugs))
            )
        )
    return nodes


def level_definition(species):
    nodes = []
    for item in LEVEL_DEFINITION:
        nodes.append(
            E.LEVEL_DEFINITION(
                E.ORDER(item[2]),
                E.ORIGINAL(item[0]),
                E.SIR(item[1])
            )
        )
    return nodes


def drugclass_definition(species):
    nodes = []
    for drugclass, drugs in DrugData.drugclass(species).items():
        drugs = [i['displayAbbr'] for i in drugs]
        nodes.append(
            E.DRUGCLASS(
                E.NAME(drugclass),
                E.DRUGLIST(
                    ','.join(drugs)
                )
            )
        )
    return nodes


def global_range_definition():
    return [
        E.GLOBALRANGE(
            etree.CDATA(GLOBALRANGE_CONTENTS)
        )
    ]


def comment_definitions(species):
    tree = E.COMMENT_DEFINITIONS()

    for comment in ConditionalCommentData.DATA[species]:
        tree.append(
            E.COMMENT_STRING(
                E.TEXT(etree.CDATA(comment['comment'])),
                E.SORT_TAG('1'),
                id=comment['commentName'],
            )
        )

    return [tree]


def definition(species):
    tree = E('DEFINITIONS')

    tree.extend(gene_definition(species))
    tree.extend(level_definition(species))
    tree.extend(drugclass_definition(species))

    tree.extend(global_range_definition())
    tree.extend(comment_definitions(species))

    return tree


# Drug
def drug_rules(species, drugname):
    all_rules = []

    # individual rule
    rules = MutationIndivComments.get_by_drugname(
        species, drugname).group_by_position()

    for pos, _ in rules.items():
        _rules = []
        for rule in _:
            aa = to_asi_format(rule['aa'])
            score = rule['score']
            _rules.append(
                '{}{} => {}'.format(pos, aa, score)
            )
        if len(_rules) == 1:
            all_rules.append(''.join(_rules))
        else:
            all_rules.append('MAX ( ' + ', '.join(_rules) + ' )')

    # Combinational rule
    rules = MutationCombiComments.get_by_drugname(
        species, drugname).group_and_order_by_position()

    for _, rule_list in rules:
        _rules = []
        for _ in rule_list:
            rule = _['rule'].replace('+', ' AND ')
            rule = to_asi_format(rule)
            score = _['score']

            rule = '({}) => {}'.format(rule, score)
            _rules.append(rule)

        if len(_rules) == 1:
            all_rules.append(''.join(_rules))
        else:
            all_rules.append('MAX (' + ', '.join(_rules) + ')')

    ident = ',\n' + DRUG_RULE_INDENTATION

    if all_rules:
        return (
            "SCORE FROM(" +
            ident.join(all_rules) +
            ")"
        )
    else:
        return None


def drug(species):
    drug_nodes = []

    for drugclass, drugs in DrugData.drugclass(species).items():
        for drug in drugs:
            rulestext = drug_rules(species, drug['name'])
            rule_node = []
            if rulestext:
                rule_node = [E.RULE(
                    E.CONDITION(
                        etree.CDATA(rulestext)
                    ),
                    E(
                        'ACTIONS',
                        E(
                            'SCORERANGE',
                            E('USE_GLOBALRANGE')
                        )
                    )
                )]
            drug_node = E(
                'DRUG',
                E.NAME(drug['displayAbbr']),
                E.FULLNAME(drug['fullName']),
                *rule_node
            )

            drug_nodes.append(drug_node)

    return drug_nodes


def mutation_comment(species):
    comments_node = E('MUTATION_COMMENTS')

    mutation_comments = ConditionalCommentData.get_comments_by(
        species, 'MUTATION')

    for gene, comments in mutation_comments.items():
        gene_node = E(
            'GENE',
            E('NAME', gene)
        )

        for comment in comments:
            rule_node = E(
                'RULE',
                E(
                    'CONDITION',
                    '{}{}'.format(
                        comment['conditionValue']['pos'],
                        to_asi_format(comment['conditionValue']['aas'])
                    )
                ),
                E(
                    'ACTIONS',
                    E('COMMENT', ref=comment['commentName'])
                )
            )

            gene_node.append(rule_node)

        comments_node.append(gene_node)

    return [comments_node]


def result_comments(species):
    tree = E('RESULT_COMMENTS')

    for comment in ConditionalCommentData.get_comments_by(
            species, 'DRUGLEVEL'):
        rule_node = E('RESULT_COMMENT_RULE')
        condition = E('DRUG_LEVEL_CONDITIONS')

        cond_value = comment['conditionValue']
        if 'and' not in cond_value:
            cond_value = [cond_value]
        else:
            cond_value = cond_value['and']

        for value in cond_value:
            for level in value['levels']:
                drugname = DrugData.get_display_abbr(species, value['drug'])
                condition.append(E.DRUG_LEVEL_CONDITION(
                    E.DRUG_NAME(drugname),
                    E.EQ(str(level))
                ))

        rule_node.append(condition)
        rule_node.append(
            E.LEVEL_ACTION(
                E.COMMENT(ref=comment['commentName'])
            )
        )
        tree.append(rule_node)

    return [tree]


def algorithm(species, version):
    """Assemble content"""
    species_suffix = '-HIV2' if species == 'HIV2' else ''

    tree = E.ALGORITHM(
        E.ALGNAME(ALG_NAME + species_suffix),
        E.ALGVERSION(version['version']),
        E.ALGDATE(version['date']),
        definition(species)
    )

    tree.extend(drug(species))
    tree.extend(mutation_comment(species))
    tree.extend(result_comments(species))

    return tree


DEFAULT_OUTPUT_FOLDER = BASEDIR / 'data' / 'algorithms'


def output(content, species, output_path, version):
    # First line, change from single quote to double quote
    lines = content.split('\n')
    for _, line in enumerate(lines):
        if _ == 0:
            line = line.replace('\'', '"')
            lines[_] = line
    species_suffix = '-HIV2' if species == 'HIV2' else ''

    if not output_path:
        output_path = 'HIVDB{suffix}_{version}.xml'.format(
            suffix=species_suffix, **version)
    else:
        output_path = output_path[0]

    output_path = DEFAULT_OUTPUT_FOLDER / output_path

    click.echo('Write to {}'.format(output_path))
    with output_path.open('w') as fd:
        fd.write('\n'.join(lines))


@click.command()
@click.argument('species', is_eager=True,
                type=click.Choice(['HIV1', 'HIV2']))
@click.option('--output-path',
              type=click.Path(dir_okay=False, writable=True))
@click.option('--target-version',
              type=str, required=True,
              callback=validate_version)
def work(species, output_path, target_version):
    result = etree.tostring(
        algorithm(species, target_version),
        pretty_print=True,
        xml_declaration=True,
        encoding='UTF-8',
        standalone=False,
        doctype=doc_type(
            'ALGORITHM',
            '',
            'https://cms.hivdb.org/prod/downloads/asi/ASI2.2.dtd'
        )
    )

    output(result.decode('utf-8'), species, output_path, target_version)


if __name__ == "__main__":
    work()
