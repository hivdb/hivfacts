# -*- coding: utf-8 -*-

import os
import json
from copy import deepcopy

from . import data

RES_DIRECTORY = os.path.dirname(data.__file__)
ORDERED_AAS = 'ACDEFGHIKLMNPQRSTVWY_-*'
NUM_AAS = len(ORDERED_AAS)
GENE_INDEX_OFFSET = {
    'PR': 0,
    'RT': 99 * NUM_AAS,
    'IN': (99 + 560) * NUM_AAS
}
GENE_SIZE = {
    'PR': 99 * NUM_AAS,
    'RT': 560 * NUM_AAS,
    'IN': 288 * NUM_AAS
}
TOTAL_LEN = GENE_INDEX_OFFSET['IN'] + GENE_SIZE['IN']


class HIVAAPcnt:

    singletons = {}

    def __new__(cls, treatment, subtype):
        lrx = treatment.lower()
        resource_name = 'rx-{}_subtype-{}.json'.format(lrx, subtype)
        if resource_name in cls.singletons:
            return cls.singletons[resource_name]
        self = super(HIVAAPcnt, cls).__new__(cls)
        self.__init_resource(resource_name)
        cls.singletons[resource_name] = self
        return self

    @staticmethod
    def _get_slice(gene=None, position=None, aa=None):
        index = 0
        if gene is None:
            return slice(None)
        if isinstance(gene, int):
            return gene
        index += GENE_INDEX_OFFSET[gene]
        if position is None:
            return slice(index, index + GENE_SIZE[gene])
        index += (max(position, 1) - 1) * NUM_AAS
        if index >= GENE_INDEX_OFFSET[gene] + GENE_SIZE[gene]:
            return slice(TOTAL_LEN, TOTAL_LEN)
        if aa is None:
            return slice(index, index + 23)
        if len(aa) == 1:
            return index + ORDERED_AAS.index(aa)
        else:
            return tuple(index + ORDERED_AAS.index(a) for a in aa)

    def __init_resource(self, resource_name):
        with open(os.path.join(RES_DIRECTORY, 'aapcnt', resource_name)) as fp:
            aapcnts = self.__aapcnts = json.load(fp)

        aapcnts_dict = {}
        for aapcnt in aapcnts:
            genepos = (aapcnt['gene'], aapcnt['position'])
            aapcnts_dict.setdefault(genepos, {})[aapcnt['aa']] = aapcnt

        self.__aapcnts_dict = aapcnts_dict

    def __iter__(self):
        return iter(self.__aapcnts)

    def __getitem_by_index__(self, index):
        if isinstance(index, slice):
            result = self.__aapcnts[index]
        elif isinstance(index, tuple):
            result = []
            for i in index:
                partial = self.__getitem_by_index__(i)
                if isinstance(partial, list):
                    result.extend(partial)
                else:
                    result.append(partial)
        else:
            if isinstance(index, int) and (index < 0 or index > TOTAL_LEN - 1):
                raise IndexError('amino acid index out of range')
            result = self.__aapcnts[index]
        return deepcopy(result)

    def __getitem__(self, val):
        index = None
        if isinstance(val, slice):
            start = (val.start
                     if isinstance(val.start, int)
                     else self._get_slice(*val.start))
            if isinstance(start, slice):
                start = start.start
            elif isinstance(start, tuple):
                start = min(start)
            stop = (val.stop
                    if isinstance(val.stop, int)
                    else self._get_slice(*val.stop))
            if isinstance(stop, slice):
                stop = stop.stop
            elif isinstance(stop, tuple):
                stop = min(stop)
            index = slice(start, stop, val.step)
        elif isinstance(val, tuple):
            if val and isinstance(val[0], tuple):
                index = tuple(self._get_slice(*v) for v in val)
            else:
                index = self._get_slice(*val)
        else:
            index = self._get_slice(val)
        return self.__getitem_by_index__(index)

    def get(self, gene=None, position=None, aa=None):
        return self[gene, position, aa]

    def get_highest_aa_percent_value(self, gene, position, mixture):
        """
        Returns the highest amino acid prevalence associated with each of
        the AA in a mixture.
        """
        return max(self[gene, position, aa]['percent'] for aa in mixture)

    def contains_unusual_aa(self, gene, position, aas):
        """Returns True if the given mutation contains any unusual AA"""
        return any(self[gene, position, aa]['isUnusual'] for aa in aas)
