.. currentmodule:: hydromt_delwaq

.. _api_demissionmodel:

==============
DemissionModel
==============

The ``DemissionModel`` class represents the D-Emission model implementation.

.. autosummary::
   :toctree: _generated

   DemissionModel

Setup methods
-------------

.. autosummary::
   :toctree: _generated

   DemissionModel.setup_config
   DemissionModel.setup_basemaps
   DemissionModel.setup_monitoring
   DemissionModel.setup_emission_raster
   DemissionModel.setup_emission_vector
   DemissionModel.setup_emission_mapping
   DemissionModel.setup_roads
   DemissionModel.setup_hydrology_forcing
   DemissionModel.setup_climate_forcing

Components
----------

If you are using python, you can access and update the model data using the components such as ``config``, ``staticdata`` etc.
This is for example useful for computing statictics, plotting etc.
The components data are usually xarray, dictionary or geopandas objects and can be accessed via the data property: `model.staticdata.data`.
The components of the DemissionModel are;

+--------------------------+-----------------------------------------------------------------+
| **Model Component**      | **Component class**                                             |
+==========================+=================================================================+
| ``model.config``         | :class:`~hydromt_delwaq.components.DemissionConfigComponent`    |
+--------------------------+-----------------------------------------------------------------+
| ``model.staticdata``     | :class:`~hydromt_delwaq.components.DelwaqStaticdataComponent`   |
+--------------------------+-----------------------------------------------------------------+
| ``model.geometry``       | :class:`~hydromt_delwaq.components.DemissionGeometryComponent`  |
+--------------------------+-----------------------------------------------------------------+
| ``model.forcing``        | :class:`~hydromt_delwaq.components.DemissionForcingComponent`   |
+--------------------------+-----------------------------------------------------------------+
| ``model.hydromaps``      | :class:`~hydromt_delwaq.components.DelwaqHydromapsComponent`    |
+--------------------------+-----------------------------------------------------------------+
| ``model.geoms``          | :class:`~hydromt.model.components.GeomsComponent`               |
+--------------------------+-----------------------------------------------------------------+

Component-level API
-------------------

The table below summarizes the important methods and attributes for each component.
These allow for fine-grained reading, writing, modification, and inspection of component data.
They are particularly useful when working interactively in Python, for example when updating
specific configuration parameters, clipping static maps, or inspecting the forcing data.

Each component exposes a ``data`` attribute, which holds the underlying model data
(e.g. :class:`dict`, :class:`xarray.Dataset`, or :class:`geopandas.GeoDataFrame`),
and supports a common set of I/O and manipulation methods such as
:meth:`read`, :meth:`write`, and :meth:`set`.

For general I/O at the model level, refer to:
:class:`~hydromt_wflow.model.WflowSbmModel` and its
:meth:`~hydromt_wflow.model.WflowSbmModel.read` and
:meth:`~hydromt_wflow.model.WflowSbmModel.write` methods.

The following table provides a detailed overview of the component-level APIs.

+-----------------------------------------------------------------+-----------------------------------------------------------------------------+--------------------------------------------------------------------+
| **Component class**                                             | **Methods**                                                                 | **Attributes**                                                     |
+=================================================================+=============================================================================+====================================================================+
| :class:`~hydromt_delwaq.components.DemissionConfigComponent`    | :meth:`~hydromt_delwaq.components.DemissionConfigComponent.read`,           | :attr:`~hydromt_delwaq.components.DemissionConfigComponent.data`   |
|                                                                 | :meth:`~hydromt_delwaq.components.DemissionConfigComponent.write`,          |                                                                    |
|                                                                 | :meth:`~hydromt_delwaq.components.DemissionConfigComponent.get_value`,      |                                                                    |
|                                                                 | :meth:`~hydromt_delwaq.components.DemissionConfigComponent.set`,            |                                                                    |
|                                                                 | :meth:`~hydromt_delwaq.components.DemissionConfigComponent.update`,         |                                                                    |
+-----------------------------------------------------------------+-----------------------------------------------------------------------------+--------------------------------------------------------------------+
| :class:`~hydromt_delwaq.components.DelwaqStaticdataComponent`   | :meth:`~hydromt_delwaq.components.DelwaqStaticdataComponent.read`,          | :attr:`~hydromt_delwaq.components.DelwaqStaticdataComponent.data`  |
|                                                                 | :meth:`~hydromt_delwaq.components.DelwaqStaticdataComponent.write`,         |                                                                    |
|                                                                 | :meth:`~hydromt_delwaq.components.DelwaqStaticdataComponent.set`,           |                                                                    |
+-----------------------------------------------------------------+-----------------------------------------------------------------------------+--------------------------------------------------------------------+
| :class:`~hydromt_delwaq.components.DemissionGeometryComponent`  | :meth:`~hydromt_delwaq.components.DemissionGeometryComponent.read`,         | :attr:`~hydromt_delwaq.components.DemissionGeometryComponent.data` |
|                                                                 | :meth:`~hydromt_delwaq.components.DemissionGeometryComponent.write`,        |                                                                    |
|                                                                 | :meth:`~hydromt_delwaq.components.DemissionGeometryComponent.set`,          |                                                                    |
+-----------------------------------------------------------------+-----------------------------------------------------------------------------+--------------------------------------------------------------------+
| :class:`~hydromt_delwaq.components.DemissionForcingComponent`   | :meth:`~hydromt_delwaq.components.DemissionForcingComponent.read`,          | :attr:`~hydromt_delwaq.components.DemissionForcingComponent.data`  |
|                                                                 | :meth:`~hydromt_delwaq.components.DemissionForcingComponent.write`,         |                                                                    |
|                                                                 | :meth:`~hydromt_delwaq.components.DemissionForcingComponent.set`,           |                                                                    |
+-----------------------------------------------------------------+-----------------------------------------------------------------------------+--------------------------------------------------------------------+
| :class:`~hydromt_delwaq.components.DelwaqHydromapsComponent`    | :meth:`~hydromt_delwaq.components.DelwaqHydromapsComponent.read`,           | :attr:`~hydromt_delwaq.components.DelwaqHydromapsComponent.data`   |
|                                                                 | :meth:`~hydromt_delwaq.components.DelwaqHydromapsComponent.write`,          |                                                                    |
|                                                                 | :meth:`~hydromt_delwaq.components.DelwaqHydromapsComponent.set`,            |                                                                    |
+-----------------------------------------------------------------+-----------------------------------------------------------------------------+--------------------------------------------------------------------+

The ``GeomsComponent`` is inherited from HydroMT core and documented in `hydromt.model.components.GeomsComponent <https://deltares.github.io/hydromt/latest/api/model_components.html#geomscomponent>`_


I/O methods
-------------

If you are using python, you can read and write the different model components using the methods below.

.. autosummary::
   :toctree: _generated

   DelwaqModel.build
   DelwaqModel.update

   DemissionModel.read
   DemissionModel.write
   DemissionModel.write_waqgeom

Attributes
----------

.. autosummary::
   :toctree: _generated

   DemissionModel.crs
   DemissionModel.root
   DemissionModel.basins
   DemissionModel.nrofseg
   DemissionModel.nrofexch
   DemissionModel.fluxes
