# -*- coding: utf-8 -*-

import os
import json
from copy import deepcopy

from . import data

RES_DIRECTORY = os.path.dirname(data.__file__)


class HIVCodonPcnt:

    singletons = {}

    def __new__(cls, treatment, subtype):
        lrx = treatment.lower()
        resource_name = 'rx-{}_subtype-{}.json'.format(lrx, subtype)
        if resource_name in cls.singletons:
            return cls.singletons[resource_name]
        self = super(HIVCodonPcnt, cls).__new__(cls)
        self.__init_resource(resource_name)
        cls.singletons[resource_name] = self
        return self

    def __init_resource(self, resource_name):
        with open(os.path.join(
                RES_DIRECTORY, 'codonpcnt', resource_name)) as fp:
            codonpcnts = self.__codonpcnts = json.load(fp)

        codonpcnts_dict = {}
        for codonpcnt in codonpcnts:
            genepos = (codonpcnt['gene'], codonpcnt['position'])
            (codonpcnts_dict.setdefault(genepos, {})
             [codonpcnt['codon']]) = codonpcnt

        self.__codonpcnts_dict = codonpcnts_dict

    def __get_codon(self, gene, position, codon):
        position_result = self.__codonpcnts_dict[(gene, position)]
        if codon in position_result:
            return position_result[codon]
        else:
            template_result = next(iter(position_result.values()))
            return {
                **template_result,
                'codon': codon,
                'aa': 'X',
                'percent': .0,
                'count': 0
            }

    def get(self, gene=None, position=None, codon=None, aa=None):
        if gene is None:
            # make a copy in case of any modification
            result = self.__codonpcnts
        elif position is None:
            result = [codonpcnt for codonpcnt in self.__codonpcnts
                      if codonpcnt['gene'] == gene]
        elif codon is None and aa is None:
            result = self.__codonpcnts_dict[(gene, position)]
        elif aa is not None:
            result = self.__codonpcnts_dict[(gene, position)]
            result = [r for r in result if r['aa'] == aa]
        else:
            result = self.__get_codon(gene, position, codon)
        return deepcopy(result)

    def get_highest_codon_percent_value(self, gene, position, mixture):
        """
        Returns the highest amino acid prevalence associated with each of
        the codon in a mixture.

        """
        pcntval = .0
        for codon in mixture:
            codon_pcntval = self.__get_codon(gene, position, codon)['percent']
            pcntval = max(pcntval, codon_pcntval)
        return pcntval

    # def contains_unusual_codon(self, gene, position, codons):
    #     """Returns True if the given mutation contains any unusual codon"""
    #     for codon in codons:
    #         codonpcnt = self.__get_codon(gene, position, codon)
    #         if codonpcnt['isUnusual']:
    #             return True
    #     return False
