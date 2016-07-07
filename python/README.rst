lineagesim package
==================

This package contains the simulation model used for *Competition Between
Continuously Evolving Lineages in Asexual Populations* by Noah Ribeck,
Joseph S. Mulka, Luis Zaman, Brian D. Connelly, and Richard E. Lenski.

Installation
------------

``lineagesim`` is maintained in the `Python Package
Index <https://pypi.python.org/pypi>`__. The most recent version of the
model and its dependencies can be installed by running:

.. code:: python

    pip install lineagesim

Specific versions of the model can be given:

.. code:: python

    pip install lineagesim==0.9.0                                                     

Installing in a Virtual Environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To ensure consistency, each release of the model requires specific
versions of its dependencies. Because these versions may differ from
what you have on your machine, we recommend installing the model and its
dependencies into a virtual environment.

Using `virtualenv <https://virtualenv.pypa.io/en/latest/>`__, first
create the virtual environment into a folder (we'll use *lsenv*):

.. code:: sh

    virtualenv --python=python3.5 lsenv

We can then enter that environment by running the ``activate`` script:

.. code:: sh

    source lsenv/bin/activate

Finally, we'll install ``lineagesim`` into the virtual environment as
done above:

.. code:: python

    pip install lineagesim

Dependiencies
~~~~~~~~~~~~~

-  `Python <https://www.python.org>`__ 2.7 or 3.5 - 3.5.1 used for the
   paper
-  `six <https://pypi.python.org/pypi/six>`__ - 1.10.0 used for the
   paper
-  `NumPy <http://www.numpy.org>`__ - 1.11.0 used for the paper
-  `SciPy <http://www.scipy.org>`__ - 0.17.0 used for the paper
-  `iGraph <http://igraph.org/python/>`__ - 0.7.1 used for the paper

Running the Model
-----------------

Simulations are run using the ``lineagesim`` command. When no arguments
are supplied, the base conditions described in the paper are used. To
see what options can be provided, run ``lineagesim`` with the ``--help``
argument.

License
-------

This project is released under the `BSD 3-Clause
License <https://opensource.org/licenses/BSD-3-Clause>`__.
