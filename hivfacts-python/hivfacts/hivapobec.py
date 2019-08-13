# -*- coding: utf-8 -*-

import os
import json

from . import data

RES_DIRECTORY = os.path.dirname(data.__file__)


class HIVAPOBEC:

    singleton = None

    @staticmethod
    def build_lookup(mutations):
        lookup = set()
        for mut in mutations:
            lookup.add((mut['gene'], mut['position'], mut['aa']))
        return lookup

    def __new__(cls):
        if cls.singleton:
            return cls.singleton
        self = super(HIVAPOBEC, cls).__new__(cls)
        self.__init_resource()
        cls.singleton = self
        return self

    def __init_resource(self):
        with open(os.path.join(
                RES_DIRECTORY, 'apobecs', 'apobecs.json')) as fp:
            apobecs = self.__apobecs = json.load(fp)

        with open(os.path.join(
                RES_DIRECTORY, 'apobecs', 'apobec_drms.json')) as fp:
            apobec_drms = self.__apobec_drms = json.load(fp)

        self.__apobec_lookup = self.build_lookup(apobecs)
        self.__apobec_drm_lookup = self.build_lookup(apobec_drms)

    def is_apobec_mutation(self, gene, pos, aas):
        aas = aas.split('_')[0]
        return any((gene, pos, aa) in self.__apobec_lookup for aa in aas)

    def is_apobec_drm(self, gene, pos, aas):
        aas = aas.split('_')[0]
        return any((gene, pos, aa) in self.__apobec_drm_lookup for aa in aas)

    def get_apobec_list(self):
        return self.__apobecs

    def get_apobec_drm_list(self):
        return self.__apobec_drms
