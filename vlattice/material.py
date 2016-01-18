import random

import numpy

from .lattice import Locator


class IncompleteMaterialInformationError(Exception):
    pass


class Anisotropy(object):
    '''
    Handles anisotropy parsing logic
    '''
    def __init__(self, data=None):
        self._data = data or {}
        self.vector = self._data.get('vector', [0.0, 0.0, 0.0])
        self.core = self._data.get('core', 0.0)
        self.surface = self._data.get('surface', self.core)

    @staticmethod
    def _axis_is_core(axis):
        return numpy.allclose(axis, 0.0)

    def from_axis(self, axis):
        _axis = numpy.array(axis)
        if Anisotropy._axis_is_core(_axis):
            return (self.core, self.vector)
        return (self.surface, _axis / numpy.linalg.norm(axis))


class Material(object):
    '''
    Represents a material descriptor, which is different to a lattice one,
    it contains information about the atoms in the sample, and the value of
    interactions and such.

    Furthermore, it can return an appropiate locator for the given material.
    '''

    def __init__(self, data):
        self.spins = data['spins']
        self.exchanges = data['exchanges']
        self.anisotropy = Anisotropy(data.get('anisotropy', {}))
        try:
            self.unitcell = numpy.array(data['unitcell']['vectors'])
        except KeyError:
            self.unitcell = numpy.eye(3)
        try:
            self.parameters = numpy.array(data['unitcell']['parameters'])
        except KeyError:
            self.parameters = numpy.array([1.0, ] * 3)
        try:
            self.crystal_coords = numpy.array(
                data['unitcell']['crystal_coords']
            )
        except KeyError:
            self.crystal_coords = False

    def locator(self):
        return Locator(self.unitcell * self.parameters,
                       crystal_coords=self.crystal_coords)

    def spin(self, kind):
        '''
        Yields the spin for the given kind of atom, raises an error if the
        information is not available.
        '''
        try:
            return random.choice(numpy.atleast_1d(self.spins[kind]))
        except KeyError:
            raise IncompleteMaterialInformationError(
                'You did not include the spin for the {} kind of atom'.format(
                    kind
                )
            )

    def exchange(self, kind):
        '''
        Yields the exchange constant for the given kind of link, raises an error
        if the information is not available.
        '''
        try:
            return random.choice(numpy.atleast_1d(self.exchanges[kind]))
        except KeyError:
            raise IncompleteMaterialInformationError(
                'You did not include the exchange for the {} exchange'.format(
                    kind
                )
            )

    def atom_kinds(self):
        return list(self.spins.keys())

    def exchange_kinds(self):
        return list(self.exchanges.keys())
