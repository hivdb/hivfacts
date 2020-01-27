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
VERSIONS_PATH = BASEDIR / 'data' / 'algorithms' / 'versions.yml'


@lru_cache()
def get_versions():
    with VERSIONS_PATH.open() as fp:
        versions = yaml.load(fp, Loader=Loader)

        return versions


_VERSION = None


def get_cur_version():
    if not _VERSION:
        version = get_versions()['HIVDB'][-1]
    else:
        match = re.match(r'(\d+)[.\-_](\d+)[.\-_]?(\d+)?', _VERSION)
        if not match:
            raise Exception(
                'Version {} is not a well-formed version.'.format(version))
        major, middle, minor = match.groups()
        search_version = major + '.' + middle
        if minor:
            search_version += '-' + minor

        for item in get_versions()['HIVDB']:
            if item[0] == search_version:
                version = item
                break
        else:
            raise Exception('Version {} not found.'.format(version))

    return dict(zip(
        ['version', 'date', 'type'],
        version
        ))


# Load data
class DataSource(type):

    def __init__(self, name, bases, attrs):
        file_name = attrs.get('FILENAME')
        if not file_name:
            return
        else:
            file_path = DATADIR / file_name

            with open(file_path) as fd:
                self.DATA = json.load(fd)


class DrugClassData(metaclass=DataSource):

    FILENAME = 'drug-classes.json'

    @classmethod
    def gene_definition(cls):
        _gene_definition = defaultdict(list)
        for drugclass in cls.DATA:
            gene = drugclass['abstractGene']
            _gene_definition[gene].append(drugclass['name'])

        return _gene_definition


class DrugData(metaclass=DataSource):

    FILENAME = 'drugs.json'
    DRUGCLASS_ORDER = [
        'NRTI',
        'NNRTI',
        'PI',
        'INSTI'
    ]

    @classmethod
    def drugclass(cls):
        _drugclass = defaultdict(list)

        for drug in cls.DATA:
            drugclass = drug['drugClass']

            _drugclass[drugclass].append(drug)

        return OrderedDict(sorted(
                _drugclass.items(),
                key=lambda i: cls.DRUGCLASS_ORDER.index(i[0]))
            )

    @classmethod
    def get_display_abbr(cls, drugname):
        for item in cls.DATA:
            if item['name'] == drugname:
                return item['displayAbbr']


class ConditionalCommentData(metaclass=DataSource):

    FILENAME = 'conditional-comments.json'

    @classmethod
    def get_comments_by(cls, condtype):
        _data = [
            item for item in cls.DATA
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


def filter_out_HIV2(cls):
    _DATA = []
    for comment in cls.DATA:
        if comment['strain'].startswith('HIV2'):
            continue
        else:
            _DATA.append(comment)

    return _DATA


ConditionalCommentData.DATA = filter_out_HIV2(ConditionalCommentData)


class MutationComments(metaclass=DataSource):

    @classmethod
    def get_by_drugname(cls, drugname):
        rules = cls._filtered_rules = []

        for rule in cls.DATA:
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

    FILENAME = 'mutation-scores-indiv.json'


class MutationCombiComments(MutationComments):

    FILENAME = 'mutation-scores-combi.json'
    MUTATION_PATTERN = r'(\d+)([A-Z]+)'

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
def gene_definition(tree):
    for gene, drugs in DrugClassData.gene_definition().items():
        tree.append(
            E.GENE_DEFINITION(
                E.NAME(gene),
                E.DRUGCLASSLIST(','.join(drugs))
            )
        )


def level_definition(tree):
    for item in LEVEL_DEFINITION:
        _ = etree.SubElement(tree, 'LEVEL_DEFINITION')

        etree.SubElement(_, 'ORDER').text = item[-1]
        etree.SubElement(_, 'ORIGINAL').text = item[0]
        etree.SubElement(_, 'SIR').text = item[1]


def drugclass_definition(tree):
    for drugclass, drugs in DrugData.drugclass().items():
        drugs = [i['displayAbbr'] for i in drugs]
        tree.append(
            E.DRUGCLASS(
                E.NAME(drugclass),
                E.DRUGLIST(
                    ','.join(drugs)
                )
            )

        )


def global_range_definition():
    return E.GLOBALRANGE(
            etree.CDATA(GLOBALRANGE_CONTENTS)
        )


def comment_definitions():
    tree = E.COMMENT_DEFINITIONS()

    for comment in ConditionalCommentData.DATA:
        tree.append(
            E.COMMENT_STRING(
                E.TEXT(etree.CDATA(comment['comment'])),
                E.SORT_TAG('1'),
                id=comment['commentName'],
            )
        )

    return tree


def definition():
    tree = E('DEFINITIONS')

    gene_definition(tree)
    level_definition(tree)
    drugclass_definition(tree)

    tree.append(global_range_definition())
    tree.append(comment_definitions())

    return tree


# Drug
def drug_rules(drugname):
    all_rules = []

    # individual rule
    rules = MutationIndivComments.get_by_drugname(
        drugname).group_by_position()

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
        drugname).group_and_order_by_position()

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

    cc = (
        "SCORE FROM(" +
        ident.join(all_rules) +
        ")"
    )
    return cc


def drug(tree):

    for drugclass, drugs in DrugData.drugclass().items():
        for drug in drugs:
            drug_node = E(
                'DRUG',
                E.NAME(drug['displayAbbr']),
                E.FULLNAME(drug['fullName']),
                E.RULE(
                    E.CONDITION(
                        etree.CDATA(drug_rules(drug['name']))
                    ),
                    E(
                        'ACTIONS',
                        E(
                            'SCORERANGE',
                            E('USE_GLOBALRANGE')
                        )
                    )
                )
            )

            tree.append(drug_node)

    return None


def mutation_comment(tree):
    comments_node = E('MUTATION_COMMENTS')

    mutation_comments = ConditionalCommentData.get_comments_by('MUTATION')

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

    tree.append(comments_node)


def result_comments():
    tree = E('RESULT_COMMENTS')

    for comment in ConditionalCommentData.get_comments_by('DRUGLEVEL'):
        rule_node = E('RESULT_COMMENT_RULE')
        condition = E('DRUG_LEVEL_CONDITIONS')

        cond_value = comment['conditionValue']
        if 'and' not in cond_value:
            cond_value = [cond_value]
        else:
            cond_value = cond_value['and']

        for value in cond_value:
            for level in value['levels']:
                drugname = DrugData.get_display_abbr(value['drug'])
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

    return tree


# Assemble content
def algorithm():

    tree = E.ALGORITHM(
        E.ALGNAME(ALG_NAME),
        E.ALGVERSION(get_cur_version()['version']),
        E.ALGDATE(get_cur_version()['date']),
        definition()
    )

    drug(tree)
    mutation_comment(tree)
    tree.append(result_comments())

    return tree


DEFAULT_OUTPUT_FOLDER = BASEDIR / 'data' / 'algorithms'
DEFAULT_OUTPUT_PATH = 'HIVDB_{}.xml'.format(get_cur_version()['version'])


def output(content, output_path=None):
    # First line, change from single quote to double quote
    lines = content.split('\n')
    for _, line in enumerate(lines):
        if _ == 0:
            line = line.replace('\'', '"')
            lines[_] = line

    if not output_path:
        output_path = DEFAULT_OUTPUT_PATH
    else:
        output_path = output_path[0]

    output_path = DEFAULT_OUTPUT_FOLDER / output_path

    with output_path.open('w') as fd:
        fd.write('\n'.join(lines))


@click.command()
@click.argument('output_path', nargs=-1)
@click.option('--version')
def work(output_path, version):
    global _VERSION
    _VERSION = version

    result = etree.tostring(
        algorithm(),
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

    output(result.decode('utf-8'), output_path)


if __name__ == "__main__":
    work()
