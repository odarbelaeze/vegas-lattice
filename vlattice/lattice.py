'''
A lattice generator scrap code, this is something that I often do, a whole
project with an example case of usage in a single file
'''
import collections
import itertools
import operator

# External dependencies
import numpy


Interaction = collections.namedtuple(
    'Interaction', ['source', 'target', 'exchange']
)


class Atom(object):
    '''
    Represents an atom
    '''
    def __init__(self, coords, kind, id):
        self.coords = coords
        self.kind = kind
        self.id = id

    def __repr__(self):
        return 'Atom({}, kind={}, id={})'.format(
            self.coords, self.kind, self.id)


class Site(object):
    '''
    Represents a site, marked with an atom and some coords
    '''
    def __init__(self, coords, atom):
        self.coords = coords
        self.atom = atom

    def __repr__(self):
        return 'Site({}, atom={})'.format(self.coords, self.atom)


class Vertex(object):
    '''
    Represents an interaction object
    '''
    def __init__(self, source, target, delta, kind):
        self.source = source
        self.target = target
        self.delta = delta
        self.kind = kind

    def target_for(self, site):
        '''
        Computes the target for a given site using our delta
        '''
        if site.atom.id != self.source:
            raise ValueError('This interaction is not meant for that atom')
        return (
            list(x + d for x, d in zip(site.coords, self.delta)),
            self.target)

    def suitable_for(self, site):
        if site.atom.id != self.source:
            return False
        return True


class Locator(object):
    '''
    A locator object that do coordinate system transforms
    '''
    def __init__(self, unitcell, crystal_coords=True):
        self.unitcell = unitcell
        self.crystal_coords = crystal_coords

    def locate(self, site):
        if self.crystal_coords:
            return self.unitcell.dot(
                numpy.array(site.coords) + site.atom.coords)
        return self.unitcell.dot(site.coords) + site.atom.coords


class Lattice(object):
    '''
    Represents a lattice with shape and periodicity
    '''
    def __init__(self, atoms, shape, pbc, vertices):
        self.atoms = atoms
        self.shape = shape
        self.pbc = pbc
        self.vertices = vertices

    @property
    def natoms(self):
        return len(self.atoms)

    def contains(self, site):
        '''
        Checks if a site is contained within the lattice
        '''
        if site.atom.id >= self.natoms:
            return False
        pbc, coo, shape = self.pbc, site.coords, self.shape
        if not all(p or (c >= 0 and c < l) for p, c, l in zip(pbc, coo, shape)):
            return False
        return True

    def bring_in(self, site, check_contains=True):
        '''
        Brings a site into the first replica of the lattice
        '''
        if check_contains and not self.contains(site):
            raise ValueError('Site not contained in the lattice')
        pbc, coo, shape = self.pbc, site.coords, self.shape
        return Site(
            list(c % l if p else c for p, c, l in zip(pbc, coo, shape)),
            atom=site.atom)

    def index(self, site, check_contains=True):
        '''
        Returns the index of a point in pbc
        '''
        if check_contains and not self.contains(site):
            raise ValueError('Site not contained in the lattice')
        coords = self.bring_in(site, check_contains=check_contains).coords
        weights = list(itertools.accumulate(self.shape, operator.mul))
        weights = [1] + weights[:-1]
        return site.atom.id + self.natoms * sum(
            map(operator.mul, weights, coords))

    def sites(self):
        '''
        Returns an iterator over all the sites of the lattice
        '''
        sequences = list(map(range, reversed(self.shape))) + [self.atoms]
        iterator = itertools.product(*sequences)
        for data in iterator:
            yield Site(tuple(reversed(data[:-1])), data[-1])

    def interactions_for(self, site):
        '''
        Yields the set of interactions and target sites for a given site
        '''
        for vertex in self.vertices:
            if not vertex.suitable_for(site):
                continue
            coord, iid = vertex.target_for(site)
            tgt = Site(coord, self.atoms[iid])
            if self.contains(tgt):
                yield Interaction(site, tgt, vertex)


class NanoParticle(Lattice):
    '''
    A spherical nanoparticle
    '''
    def __init__(self, locator, radius, *args):
        super().__init__(*args)
        self.locator = locator
        self.radius = radius
        self.centroid = numpy.mean(
            [locator.locate(site) for site in super().sites()],
            axis=0,
        )

    def contains(self, site):
        return super().contains(site) and numpy.linalg.norm(
            self.locator.locate(site) - self.centroid) <= self.radius

    def sites(self):
        for site in super().sites():
            if self.contains(site):
                yield site


class DepletedLattice(Lattice):
    '''
    A randomly depleted lattice
    '''
    def __init__(self, probability, *args):
        super().__init__(*args)
        self.probability = probability
        self._cache = {}

    def contains(self, site):
        if not super().contains(site):
            return False
        index = self.index(site, check_contains=False)
        if index not in self._cache:
            self._cache[index] = numpy.random.random() < self.probability
        return self._cache[index]

    def sites(self):
        for site in super().sites():
            if self.contains(site):
                yield site
