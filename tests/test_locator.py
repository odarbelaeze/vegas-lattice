import numpy
import pytest

from vlattice.lattice import Atom
from vlattice.lattice import Locator
from vlattice.lattice import Site


@pytest.fixture(scope='module')
def cubic_locator(request):
    return Locator(numpy.eye(3), crystal_coords=False)


@pytest.fixture(scope='module')
def parallel_locator(request):
    return Locator(numpy.diag([1, 2, 3]), crystal_coords=True)


def test_basic_locator(cubic_locator):
    atom = Atom((0, 0, 0), 'A', 0)
    site = Site((0, 0, 0), atom)
    assert ([0, 0, 0] == cubic_locator.locate(site)).all()


def test_parallel_locator(parallel_locator):
    atom = Atom((0, 0, 0), 'A', 0)
    assert ([0, 0, 0] == parallel_locator.locate(Site((0, 0, 0), atom))).all()
    assert ([1, 0, 0] == parallel_locator.locate(Site((1, 0, 0), atom))).all()
    assert ([0, 2, 0] == parallel_locator.locate(Site((0, 1, 0), atom))).all()
    assert ([0, 0, 3] == parallel_locator.locate(Site((0, 0, 1), atom))).all()
    assert ([1, 2, 3] == parallel_locator.locate(Site((1, 1, 1), atom))).all()


def test_crystal_coords(parallel_locator):
    pl = parallel_locator
    atom = Atom((0.5, 0.5, 0.5), 'A', 0)
    assert ([0.5, 1.0, 1.5] == pl.locate(Site((0, 0, 0), atom))).all()
    assert ([1.5, 1.0, 1.5] == pl.locate(Site((1, 0, 0), atom))).all()
    assert ([0.5, 3.0, 1.5] == pl.locate(Site((0, 1, 0), atom))).all()
    assert ([0.5, 1.0, 4.5] == pl.locate(Site((0, 0, 1), atom))).all()
    assert ([1.5, 3.0, 4.5] == pl.locate(Site((1, 1, 1), atom))).all()
