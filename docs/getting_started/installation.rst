.. _installation_guide:

Installation
============

User install
------------

HydroMT is available from pypi and conda-forge, but we recommend installing with conda.
It is the same for the delwaq plugin.

To install HydroMT delwaq plugin using conda do:

.. code-block:: console

    $ conda install hydromt_delwaq -c conda-forge

To install a HydroMT environment with conda installed do:

.. code-block:: console

    $ conda create hydromt_delwaq -n hydromt-delwaq -c conda-forge

This will automatically install both HydroMT core library and dependencies as well as the delwaq plugin.


Alternatively, you can also install using pip.

.. code-block:: console

    $ pip install hydromt_delwaq


Developper install
------------------
If you want to download the DELWAQ plugin directly from git to easily have access to the latest developmemts or 
make changes to the code you can use the following steps.

First, clone hydromt's delwaq plugin ``git`` repo from
`github <https://github.com/Deltares/hydromt_delwaq>`_, then navigate into the 
the code folder (where the envs folder and pyproject.toml are located):

.. code-block:: console

    $ git clone https://github.com/Deltares/hydromt_delwaq.git
    $ cd hydromt_delwaq

Then, make and activate a new hydromt-delwaq conda environment based on the envs/hydromt-delwaq.yml
file contained in the repository:

.. code-block:: console

    $ conda env create -f envs/hydromt-delwaq.yml
    $ conda activate hydromt-delwaq

Finally, build and install hydromt_delwaq using pip.

.. code-block:: console

    $ pip install .

If you wish to make changes in hydromt_delwaq, you should make an editable install of hydromt. 
This is possible using the `flit <https://flit.readthedocs.io/en/latest/>`_ package and install command.

For Windows:

.. code-block:: console

    $ flit install --pth-file

For Linux:

.. code-block:: console

    $ flit install -s

For more information about how to contribute, see `HydroMT contributing guidelines <https://hydromt.readthedocs.io/en/latest/contributing.html>`_.
