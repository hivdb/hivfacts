# -*- coding: utf-8 -*-

import os
import json

from . import data

RES_DIRECTORY = os.path.dirname(data.__file__)


class HIVSDRM:

    singleton = None

    @staticmethod
    def build_lookup(all_mutations):
        lookup = set()
        for _, mutations in all_mutations.items():
            for mut in mutations:
                for aa in mut['aa']:
                    lookup.add((mut['gene'], mut['position'], aa))
        return lookup

    def __new__(cls):
        if cls.singleton:
            return cls.singleton
        self = super(HIVSDRM, cls).__new__(cls)
        self.__init_resource()
        cls.singleton = self
        return self

    def __init_resource(self):
        with open(os.path.join(
                RES_DIRECTORY, 'sdrms_hiv1.json')) as fp:
            sdrms = self.__sdrms = json.load(fp)

        self.__sdrm_lookup = self.build_lookup(sdrms)

    def is_sdrm(self, gene, pos, aas):
        if '_' in aas:
            aas = '_'
        return any((gene, pos, aa) in self.__sdrm_lookup for aa in aas)

    def get_sdrms(self, gene, pos, aas):
        if '_' in aas:
            aas = '_'
        sdrms = set()
        for aa in aas:
            if (gene, pos, aa) in self.__sdrm_lookup:
                sdrms.add((gene, pos, aa))
        return sdrms

    def get_sdrm_lists(self):
        return self.__sdrms

    def get_sdrm_lookup_set(self):
        return self.__sdrm_lookup
