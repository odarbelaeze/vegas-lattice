'''
A lattice generator scrap code, this is something that I often do, a whole
project with an example case of usage in a single file
'''

import itertools
import json
import operator
import random

# External dependencies
import click
import numpy


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


class Interaction(object):
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

    def bring_in(self, site):
        '''
        Brings a site into the first replica of the lattice
        '''
        if not self.contains(site):
            raise ValueError('Site not contained in the lattice')
        pbc, coo, shape = self.pbc, site.coords, self.shape
        return Site(
            list(c % l if p else c for p, c, l in zip(pbc, coo, shape)),
            atom=site.atom)

    def index(self, site):
        '''
        Returns the index of a point in pbc
        '''
        if not self.contains(site):
            raise ValueError('Site not contained in the lattice')
        coords = self.bring_in(site).coords
        weights = list(itertools.accumulate(self.shape, operator.mul))
        weights = [1] + weights[:-1]
        return site.atom.id + self.natoms * sum(
            map(operator.mul, weights, coords))

    def sites(self):
        '''
        Returns an iterator over all the sites of the lattice
        '''
        iterator = itertools.product(
            *map(range, reversed(self.shape)), self.atoms)
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
                yield site, tgt, vertex


class Locator(object):
    '''
    A locator object that do coordinate system transforms
    '''
    def __init__(self, unitcell, crystal_coords=True):
        self.unitcell = unitcell
        self.crystal_coords = crystal_coords

    def locate(self, site):
        if self.crystal_coords:
            return self.unitcell @ (numpy.array(site.coords) + site.atom.coords)
        return (self.unitcell @ site.coords) + site.atom.coords


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


def get_spin(site):
    if site.atom.kind == 'A':
        return 2.5
    return random.choice([2.0, 2.5])


@click.group()
def cli():
    '''
    Some generators for magnetite based nanostructures
    '''
    pass


@cli.command()
@click.argument('descriptor', type=click.File('r'))
@click.option('--shape', default=(1, 1, 1),
              help='shape of the lattice')
@click.option('--pbc', default=(True, True, True),
              help='use periodic boundary conditions')
def bulk(descriptor, shape, pbc):

    data = json.load(descriptor)
    atoms = [Atom(**kw) for kw in data['atoms']]
    vertices = [Interaction(**kw) for kw in data['interactions']]
    exchanges = data['material']['exchanges']
    latt = Lattice(atoms, shape, pbc, vertices)

    lsites = []
    linteractions = []

    unitcell = 8.0 * numpy.eye(3)
    locator = Locator(unitcell, crystal_coords=False)

    for site in latt.sites():
        lsites.append(site)
        for interaction in latt.interactions_for(site):
            linteractions.append(interaction)

    click.echo("{}\t{}\t{}\t{}".format(len(lsites), len(linteractions), 2, 0))
    click.echo("A\nB")
    for site in lsites:
        px, py, pz = tuple(locator.locate(site))
        click.echo(
            "{id}\t{px}\t{py}\t{pz}\t{spin}\t0.0\t0.0\t0.0\t0.0\t{kind}".format(
                id=latt.index(site),
                px=px, py=py, pz=pz,
                spin=get_spin(site),
                kind=site.atom.kind
            ))

    for source, target, vertex in linteractions:
        click.echo('{source}\t{target}\t{exchange}'.format(
            source=latt.index(source),
            target=latt.index(target),
            exchange=1      # exchanges[vertex.kind]
            ))


@cli.command()
@click.argument('descriptor', type=click.File('r'))
@click.option('--diameter', default=5,
              help='Diameter of the nanoparticle')
def nanoparticle(descriptor, diameter):

    data = json.load(descriptor)
    atoms = [Atom(**kw) for kw in data['atoms']]
    vertices = [Interaction(**kw) for kw in data['interactions']]
    exchanges = data['material']['exchanges']
    unitcell = 8.0 * numpy.eye(3)
    locator = Locator(unitcell, crystal_coords=False)
    shape = (diameter, ) * 3
    pbc = (False, ) * 3

    latt = NanoParticle(locator, diameter*4, atoms, shape, pbc, vertices)

    lsites = []
    linteractions = []
    new_ids = {}

    for idx, site in enumerate(latt.sites()):
        lsites.append(site)
        new_ids[latt.index(site)] = idx
        for interaction in latt.interactions_for(site):
            linteractions.append(interaction)

    click.echo("{}\t{}\t{}\t{}".format(len(lsites), len(linteractions), 2, 0))
    click.echo("A\nB")

    for site in lsites:
        px, py, pz = tuple(locator.locate(site))
        click.echo(
            "{id}\t{px}\t{py}\t{pz}\t{spin}\t0.0\t0.0\t0.0\t0.0\t{kind}".format(
                id=new_ids[latt.index(site)],
                px=px, py=py, pz=pz,
                spin=get_spin(site),
                kind=site.atom.kind
            ))

    for source, target, vertex in linteractions:
        click.echo('{source}\t{target}\t{exchange}'.format(
            source=new_ids[latt.index(source)],
            target=new_ids[latt.index(target)],
            exchange=exchanges[vertex.kind]
            ))


if __name__ == '__main__':
    cli()
