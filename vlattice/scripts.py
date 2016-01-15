import collections
import itertools
import json

import click
import numpy

from .lattice import Atom
from .lattice import Vertex
from .lattice import Lattice
from .lattice import NanoParticle
from .material import Material


Geometry = collections.namedtuple('Geometry', ['sites', 'links', 'axes', 'patch'])


def compute_geometry(lattice, material):
    locator = material.locator()

    lsites = []
    linteractions = []
    axes = []
    # Patch the ids according to the removed material
    new_ids = {}

    for idx, site in enumerate(lattice.sites()):
        lsites.append(site)
        new_ids[lattice.index(site)] = idx
        axis = numpy.array([0.0, 0.0, 0.0])
        for interaction in lattice.interactions_for(site):
            linteractions.append(interaction)
            axis += locator.locate(interaction.target) - locator.locate(site)
        axes.append(axis)

    return Geometry(lsites, linteractions, axes, new_ids)


def echo_header(material, sites, links):
    atom_kinds = material.atom_kinds()
    click.echo("{}\t{}\t{}\t{}".format(
        len(sites), len(links), len(atom_kinds), 0))
    click.echo('\n'.join(atom_kinds))


def echo_sites(material, lattice, sites, patch=None):
    patch = patch or {}
    locator = material.locator()
    for site in sites:
        px, py, pz = tuple(locator.locate(site))
        index = lattice.index(site)
        click.echo(
            "{id}\t{px}\t{py}\t{pz}\t{spin}\t0.0\t0.0\t0.0\t0.0\t{kind}".format(
                id=patch.get(index, index),
                px=px, py=py, pz=pz,
                spin=material.spin(site.atom.kind),
                kind=site.atom.kind
            ))


def echo_interactions(material, lattice, interactions, patch=None):
    patch = patch or {}
    for source, target, vertex in interactions:
        sindex = lattice.index(source)
        tindex = lattice.index(target)
        click.echo('{source}\t{target}\t{exchange}'.format(
            source=patch.get(sindex, sindex),
            target=patch.get(tindex, tindex),
            exchange=material.exchange(vertex.kind)
            ))


@click.group()
def cli():
    '''
    Some generators for grid graph based nanostructures, this program would
    generate samples from descriptor files using linear time algorithms which
    are much faster than those based on checking distances between points
    '''
    pass


@cli.command()
@click.argument('descriptor', type=click.File('r'))
@click.option('--shape', default=(1, 1, 1),
              help='shape of the lattice')
@click.option('--pbc', default=(True, True, True),
              help='use periodic boundary conditions')
@click.option('--lattice-params', default=(1.0, 1.0, 1.0),
              help='lattice parameters for the atoms in the descriptor')
def bulk(descriptor, shape, pbc, lattice_params):

    data = json.load(descriptor)
    atoms = [Atom(**kw) for kw in data['atoms']]
    vertices = [Vertex(**kw) for kw in data['interactions']]
    material = Material(data['material'])
    latt = Lattice(atoms, shape, pbc, vertices)

    lsites, linteractions, axes, _ = compute_geometry(latt, material)

    echo_header(material, lsites, linteractions)
    echo_sites(material, latt, lsites, axes)
    echo_interactions(material, latt, linteractions)


@cli.command()
@click.argument('descriptor', type=click.File('r'))
@click.option('--diameter', default=5,
              help='Diameter of the nanoparticle')
@click.option('--lattice-params', default=(1.0, 1.0, 1.0),
              help='lattice parameters for the atoms in the descriptor')
def nanoparticle(descriptor, diameter, lattice_params):

    data = json.load(descriptor)
    atoms = [Atom(**kw) for kw in data['atoms']]
    vertices = [Vertex(**kw) for kw in data['interactions']]
    material = Material(data['material'])
    shape = (diameter, ) * 3
    pbc = (False, ) * 3

    scale = numpy.amax(material.parameters)
    locator = material.locator()
    latt = NanoParticle(locator, diameter*scale/2, atoms, shape, pbc, vertices)

    lsites, linteractions, axes, patch = compute_geometry(latt, material)

    echo_header(material, lsites, linteractions)
    echo_sites(material, latt, lsites, patch=patch)
    echo_interactions(material, latt, linteractions, patch=patch)


@cli.command()
@click.argument('sites', type=click.File('r'))
@click.argument('descriptor', type=click.File('w'), default="-")
@click.option('--lattice-params', default=(1.0, 1.0, 1.0),
              help='Shape of the unitcell (default: a=b=c=1)')
@click.option('--cut', default=1.0,
              help='Neighbor cutoff radius (default: 1.0)')
def describe(sites, descriptor, lattice_params, cut):
    '''
    Generates a descriptor file from a site list file, the site list file format
    should be a plain text file with the following format:

        <x> <y> <z> <kind>

    the kind parameter is taked in as a string, and the identifiers of the sites
    are generated according to the order in which this file is presented,
    furthermore, you should include the lattice parameters in x y and z for
    the unit cell.

    If DESCRIPTOR is specified and writable, the json representation of the
    descriptor file will be written there, otherwise the standard output will
    be used.
    '''

    points = numpy.loadtxt(sites, usecols=(0, 1, 2))
    # Rewind the file handle
    sites.seek(0, 0)
    labels = numpy.loadtxt(sites, usecols=(3, ), dtype=str)
    labels = [l[2:-1] for l in labels]
    images = map(numpy.array, itertools.product((-1, 0, 1), repeat=3))

    spins = {label: None for label in labels}
    exchanges = {
        label + other: None
        for label, other in itertools.product(spins.keys(), repeat=2)}

    expanded = []
    for image in images:
        imaged = points + lattice_params * image
        for site, kind, uuid in zip(imaged, labels, itertools.count(0)):
            expanded.append({
                'site': site,
                'kind': kind,
                'image': image,
                'id': uuid,
            })

    norm = numpy.linalg.norm
    interactions = []
    for uuid, site in enumerate(points):
        for other in expanded:
            dis = norm(site - other['site'])
            is_real = (other['image'] == (0, 0, 0, )).all()
            if dis < cut and (not is_real or uuid != other['id']):
                interactions.append(collections.OrderedDict([
                    ('source', uuid),
                    ('target', other['id']),
                    ('delta', [str(i) for i in other['image']]),
                    ('type', labels[uuid] + other['kind']),
                ]))

    atoms = [{'coords': list(site), 'kind': kind, 'id': uuid}
             for site, kind, uuid in zip(points, labels, itertools.count(0))]

    data = collections.OrderedDict([
        ('material', {'spins': spins, 'exchanges': exchanges, }),
        ('atoms', atoms),
        ('interactions', interactions)
    ])
    json.dump(data, descriptor, indent=2)


if __name__ == '__main__':
    cli()
