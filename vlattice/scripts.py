import json
import random

import click
import numpy

from .lattice import Atom
from .lattice import Interaction
from .lattice import Lattice
from .lattice import Locator
from .lattice import NanoParticle


def get_spin(site, spins):
    return random.choice(numpy.atleast_1d(spins[site.atom.kind]))


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
    spins = data['material']['spins']
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
                spin=get_spin(site, spins),
                kind=site.atom.kind
            ))

    for source, target, vertex in linteractions:
        click.echo('{source}\t{target}\t{exchange}'.format(
            source=latt.index(source),
            target=latt.index(target),
            exchange=exchanges[vertex.kind]
            ))


@cli.command()
@click.argument('descriptor', type=click.File('r'))
@click.option('--diameter', default=5,
              help='Diameter of the nanoparticle')
def nanoparticle(descriptor, diameter):

    data = json.load(descriptor)
    atoms = [Atom(**kw) for kw in data['atoms']]
    vertices = [Interaction(**kw) for kw in data['interactions']]
    spins = data['material']['spins']
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
                spin=get_spin(site, spins),
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
