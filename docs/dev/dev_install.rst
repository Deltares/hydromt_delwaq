.. _dev_env:

Developer's environment
-----------------------

If you want to download the HydroMT-DELWAQ plugin directly from git to easily have access to the latest developments or make
changes to the code you can use the following steps.

First, clone the HydroMT-DELWAQ plugin ``git`` repo from `github <https://github.com/Deltares/hydromt_delwaq.git>`_.
If you need a github authentification key, you can use `ssh <https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account>`_

.. code-block:: console

    $ git clone git@github.com:Deltares/hydromt_delwaq.git

Then, navigate into the code folder (where the envs folder and pyproject.toml are located).
Make and activate a new hydromt-delwaq conda environment based on a "full" installation
including all optional and developer dependencies as specified in the `pyproject.toml` and
available in PyPI. This will ensure you have all necessary packages to run examples and tests.

.. code-block:: console

    $ conda create -n hydromt-delwaq uv python=3.13
    $ conda activate hydromt-delwaq
    $ uv pip install "hydromt_delwaq[full]"

Finally, create a developer installation of hydromt_delwaq.
This can be done with:

.. code-block:: console

    $ pip install -e .
