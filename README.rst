.. image:: https://travis-ci.org/odarbelaeze/vegas-lattice.svg?branch=master
    :target: https://travis-ci.org/odarbelaeze/vegas-lattice

=============
vegas-lattice
=============

A lattice generator part of the **vegas** software initiative, it provides the
esential routines to generate regular graph lattices in linear time, as well as
some cuts of those lattices for nano particles and randomly depleted lattices.

Instalation
-----------

The **vegas-lattice** program has been tested in python 3.4 and 3.5 to work
properly, you can install the package from the python package index using the
following command:

.. code-block:: bash

    pip install vegas-lattice   # this might require sudo

After that the `vegas-lattice` command should be available to your command
prompt.

Feel free to use the descriptor files provided in the **docs** directory to do
some tests, and also to run

.. code-block:: bash

    vegas-lattice --help

to get some instructions in how to use the program.

Features
--------

The **vegas-lattice** package currently features three main routines:

- **bulk** creates a cubic lattice with the specified periodic boundary
  conditions.
- **nanoparticle** creates a spherical nanoparticle of the given diameter.
- **depleted** creates a bulk lattice with random site depletion.
- **describe** bootstraps a descriptor file using a simple site list.

