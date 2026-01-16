.. image:: https://mybinder.org/badge_logo.svg
    :target: https://mybinder.org/v2/gh/Deltares/hydromt_delwaq/main?urlpath=lab/tree/examples

Several iPython notebook examples have been prepared for **HydroMT-DELWAQ** which you can
use as a HydroMT-DELWAQ tutorial.

These examples can be run online or on your local machine.
To run these examples online press the **binder** badge above.

To run these examples on your local machine you need a copy of the examples folder
of the repository and an installation of HydroMT-DELWAQ including some additional
packages required to run the notebooks.

1 - Install HydroMT-DELWAQ
--------------------------

The first step is to install HydroMT-DELWAQ and the other Python dependencies in a separate environment,
see the **Install HydroMT-DELWAQ in a new environment** section in the :ref:`installation guide <installation_guide>`.

To run the notebooks, you need to install additional dependencies that are jupyterlab to
run the notebooks and cartopy for some plots. To install these packages in your existing
hydromt-delwaq environment, do (these packages are also available in conda-forge):

.. code-block:: console

  $ conda activate hydromt-delwaq
  $(hydromt-delwaq) uv pip install "hydromt_delwaq[examples]"

2 - Download the content of the examples and notebooks
******************************************************
To run the examples locally, you will need to download the content of the hydromt_delwaq repository.
You have two options:

  1. Download and unzip the examples manually
  2. Clone the hydromt_delwaq GitHub repository

.. warning::

  Depending on your installed version of hydromt and hydromt_delwaq, you will need to download the correct versions of the examples.
  To check the version of hydromt_delwaq that you have installed, do:

  .. code-block:: console

    $(hydromt-delwaq) hydromt --models
        Model plugins:
          - model (hydromt 1.3.0)
          - delwaq (hydromt_delwaq 0.4.0)
          - demission (hydromt_delwaq 0.4.0)
          ...

In the examples above, we see version 0.4.0 of hydromt_delwaq is installed and version 1.3.0 of hydromt.

**Option 1: manual download and unzip**

To manually download the examples on Windows, do (!replace with your own hydromt_delwaq version!):

.. code-block:: console

  $ curl https://github.com/Deltares/hydromt_delwaq/archive/refs/tags/v0.4.0.zip -O -L
  $ tar -xf v0.4.0.zip
  $ ren hydromt_delwaq-0.4.0 hydromt_delwaq

You can also download, unzip and rename manually if you prefer, rather than using the windows command prompt.

**Option 2: cloning the hydromt_delwaq repository**
For git users, you can also get the examples by cloning the hydromt_delwaq github repository and checkout the version
you have installed:

.. code-block:: console

  $ git clone https://github.com/Deltares/hydromt_delwaq.git
  $ git checkout v0.4.0

3 - Running the examples
************************
Finally, start a jupyter lab server inside the **examples** folder
after activating the **hydromt-delwaq** environment, see below.

Alternatively, you can run the notebooks from `Visual Studio Code <https://code.visualstudio.com/download>`_.

.. code-block:: console

  $ conda activate hydromt-delwaq
  $ cd hydromt_delwaq/examples
  $ jupyter lab
