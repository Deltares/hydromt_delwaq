.. currentmodule:: hydromt_delwaq

.. _api_delwaqmodel:

===========
DelwaqModel
===========

The ``DelwaqModel`` class represents the D-Water Quality fate and transport model implementation.

.. autosummary::
   :toctree: _generated

   DelwaqModel

Setup methods
-------------

.. autosummary::
   :toctree: _generated

   DelwaqModel.setup_config
   DelwaqModel.setup_basemaps
   DelwaqModel.setup_monitoring
   DelwaqModel.setup_staticdata_from_rasterdataset
   DelwaqModel.setup_hydrology_forcing
   DelwaqModel.setup_sediment_forcing
   DelwaqModel.setup_climate_forcing

Components
----------

If you are using python, you can access and update the model data using the components such as ``config``, ``staticdata`` etc.
This is for example useful for computing statictics, plotting etc.
The components data are usually xarray, dictionary or geopandas objects and can be accessed via the data property: `model.staticdata.data`.
The components of the DelwaqModel are;

+--------------------------+-----------------------------------------------------------------+
| **Model Component**      | **Component class**                                             |
+==========================+=================================================================+
| ``model.config``         | :class:`~hydromt_delwaq.components.DelwaqConfigComponent`       |
+--------------------------+-----------------------------------------------------------------+
| ``model.staticdata``     | :class:`~hydromt_delwaq.components.DelwaqStaticdataComponent`   |
+--------------------------+-----------------------------------------------------------------+
| ``model.pointer``        | :class:`~hydromt_delwaq.components.DelwaqPointerComponent`      |
+--------------------------+-----------------------------------------------------------------+
| ``model.forcing``        | :class:`~hydromt_delwaq.components.DelwaqForcingComponent`      |
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

+-----------------------------------------------------------------+-----------------------------------------------------------------------------+-------------------------------------------------------------------+
| **Component class**                                             | **Methods**                                                                 | **Attributes**                                                    |
+=================================================================+=============================================================================+===================================================================+
| :class:`~hydromt_delwaq.components.DelwaqConfigComponent`       | :meth:`~hydromt_delwaq.components.DelwaqConfigComponent.read`,              | :attr:`~hydromt_delwaq.components.DelwaqConfigComponent.data`     |
|                                                                 | :meth:`~hydromt_delwaq.components.DelwaqConfigComponent.write`,             |                                                                   |
|                                                                 | :meth:`~hydromt_delwaq.components.DelwaqConfigComponent.get_value`,         |                                                                   |
|                                                                 | :meth:`~hydromt_delwaq.components.DelwaqConfigComponent.set`,               |                                                                   |
|                                                                 | :meth:`~hydromt_delwaq.components.DelwaqConfigComponent.update`,            |                                                                   |
+-----------------------------------------------------------------+-----------------------------------------------------------------------------+-------------------------------------------------------------------+
| :class:`~hydromt_delwaq.components.DelwaqStaticdataComponent`   | :meth:`~hydromt_delwaq.components.DelwaqStaticdataComponent.read`,          | :attr:`~hydromt_delwaq.components.DelwaqStaticdataComponent.data` |
|                                                                 | :meth:`~hydromt_delwaq.components.DelwaqStaticdataComponent.write`,         |                                                                   |
|                                                                 | :meth:`~hydromt_delwaq.components.DelwaqStaticdataComponent.set`,           |                                                                   |
+-----------------------------------------------------------------+-----------------------------------------------------------------------------+-------------------------------------------------------------------+
| :class:`~hydromt_delwaq.components.DelwaqPointerComponent`      | :meth:`~hydromt_delwaq.components.DelwaqPointerComponent.read`,             | :attr:`~hydromt_delwaq.components.DelwaqPointerComponent.data`    |
|                                                                 | :meth:`~hydromt_delwaq.components.DelwaqPointerComponent.write`,            |                                                                   |
|                                                                 | :meth:`~hydromt_delwaq.components.DelwaqPointerComponent.set`,              |                                                                   |
+-----------------------------------------------------------------+-----------------------------------------------------------------------------+-------------------------------------------------------------------+
| :class:`~hydromt_delwaq.components.DelwaqForcingComponent`      | :meth:`~hydromt_delwaq.components.DelwaqForcingComponent.read`,             | :attr:`~hydromt_delwaq.components.DelwaqForcingComponent.data`    |
|                                                                 | :meth:`~hydromt_delwaq.components.DelwaqForcingComponent.write`,            |                                                                   |
|                                                                 | :meth:`~hydromt_delwaq.components.DelwaqForcingComponent.set`,              |                                                                   |
+-----------------------------------------------------------------+-----------------------------------------------------------------------------+-------------------------------------------------------------------+
| :class:`~hydromt_delwaq.components.DelwaqHydromapsComponent`    | :meth:`~hydromt_delwaq.components.DelwaqHydromapsComponent.read`,           | :attr:`~hydromt_delwaq.components.DelwaqHydromapsComponent.data`  |
|                                                                 | :meth:`~hydromt_delwaq.components.DelwaqHydromapsComponent.write`,          |                                                                   |
|                                                                 | :meth:`~hydromt_delwaq.components.DelwaqHydromapsComponent.set`,            |                                                                   |
+-----------------------------------------------------------------+-----------------------------------------------------------------------------+-------------------------------------------------------------------+

The ``GeomsComponent`` is inherited from HydroMT core and documented in `hydromt.model.components.GeomsComponent <https://deltares.github.io/hydromt/latest/api/model_components.html#geomscomponent>`_

I/O methods
-------------

If you are using python, you can read and write the different model components using the methods below.

.. autosummary::
   :toctree: _generated

   DelwaqModel.build
   DelwaqModel.update

   DelwaqModel.read
   DelwaqModel.write
   DelwaqModel.write_waqgeom

Attributes
----------

.. autosummary::
   :toctree: _generated

   DelwaqModel.crs
   DelwaqModel.root
   DelwaqModel.basins
   DelwaqModel.nrofseg
   DelwaqModel.nrofexch
   DelwaqModel.surface_water
   DelwaqModel.fluxes
