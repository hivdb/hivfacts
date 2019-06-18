# -*- coding: utf-8 -*-

import os
import json
from copy import deepcopy

from . import data

RES_DIRECTORY = os.path.dirname(data.__file__)


class HIVAAPcnt:

    singletons = {}

    def __new__(cls, treatment, subtype):
        resource_name = '{}{}.json'.format(treatment, subtype)
        if resource_name in cls.singletons:
            return cls.singletons[resource_name]
        self = super(HIVAAPcnt, cls).__new__(cls)
        self.__init_resource(resource_name)
        cls.singletons[resource_name] = self
        return self

    def __init_resource(self, resource_name):
        with open(os.path.join(RES_DIRECTORY, resource_name)) as fp:
            aapcnts = self.__aapcnts = json.load(fp)

        aapcnts_dict = {}
        for aapcnt in aapcnts:
            genepos = (aapcnt['gene'], aapcnt['position'])
            aapcnts_dict.setdefault(genepos, {})[aapcnt['aa']] = aapcnt

        self.__aapcnts_dict = aapcnts_dict

    def get(self, gene=None, position=None, aa=None):
        if gene is None:
            # make a copy in case of any modification
            result = self.__aapcnts
        elif position is None:
            result = [aapcnt for aapcnt in self.__aapcnts
                      if aapcnt['gene'] == gene]
        elif aa is None:
            result = self.__aapcnts_dict[(gene, position)]
        else:
            result = self.__aapcnts_dict[(gene, position)][aa]
        return deepcopy(result)

    """
    Returns the highest amino acid prevalence associated with each of
    the AA in a mixture.
    """
    def get_highest_aa_percent_value(self, gene, position, mixture):
        pcntval = .0
        gpos = (gene, position)
        for aa in mixture:
            aa_pcntval = self.__aapcnts_dict[gpos][aa]['percent']
            pcntval = max(pcntval, aa_pcntval)
        return pcntval

    """Returns True if the given mutation contains any unusual AA"""
    def contains_unusual_aa(self, gene, position, aas):
        gpos = (gene, position)
        for aa in aas:
            aapcnt = self.__aapcnts_dict[gpos][aa]
            if aapcnt['isUnusual']:
                return True
        return False
