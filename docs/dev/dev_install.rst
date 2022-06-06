.. _dev_env:

Developer's environment
-----------------------

If you want to download the HydroMT-delwaq plugin directly from git to easily have access to the latest developments or make 
changes to the code you can use the following steps.

First, clone the HydroMT-delwaq plugin ``git`` repo from `github <https://github.com/Deltares/hydromt_delwaq.git>`_. 
If you need a github authentification key, you can use `ssh <https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account>`_

.. code-block:: console

    $ git clone git@github.com:Deltares/hydromt_delwaq.git

Then, navigate into the code folder (where the envs folder and pyproject.toml are located):
Make and activate a new Hydromt-delwaq conda environment based on the envs/hydromt-delwaq.yml
file contained in the repository:

.. code-block:: console

    $ cd hydromt_delwaq
    $ conda env create -f envs/hydromt-delwaq.yml
    $ conda activate hydromt-delwaq

Finally, create a developer installation of HydroMT-delwaq.
This is possible using the `flit <https://flit.readthedocs.io/en/latest/>`_ package and install command.

For Windows:

.. code-block:: console

    $ flit install --pth-file

For Linux:

.. code-block:: console

    $ flit install -s