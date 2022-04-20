.. _dev_env:

Developer's environment
-----------------------

Developing HydroMT requires Python >= 3.6. We prefer developing with the most recent 
version of Python. We strongly encourage you to develop in a separate conda environment.
All Python dependencies required to develop HydroMT can be found in `envs/hydromt-delwaq.yml <environment.yml>`__.

Developer installation guide
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First, clone hydromt-delwaq's ``git`` repo using `ssh <https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account>`_ from
`github <https://github.com/Deltares/hydromt_delwaq.git>`_.

.. code-block:: console

    $ git clone git@github.com:Deltares/hydromt_delwaq.git
    $ cd hydromt

Then, navigate into the the code folder (where the envs folder and pyproject.toml are located):
Make and activate a new hydromt-delwaq conda environment based on the envs/hydromt-delwaq.yml
file contained in the repository:

.. code-block:: console

    $ conda env create -f envs/hydromt-delwaq.yml
    $ conda activate hydromt-delwaq

Finally, create a developer installation of HydroMT:
This is possible using the `flit <https://flit.readthedocs.io/en/latest/>`_ package and install command.

For Windows:

.. code-block:: console

    $ flit install --pth-file

For Linux:

.. code-block:: console

    $ flit install -s