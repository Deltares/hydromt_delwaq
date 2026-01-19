.. currentmodule:: hydromt_delwaq.components

.. _api_components:

==========
Components
==========

The ``hydromt_delwaq.components`` module defines reusable data container classes
that represent configuration, static maps, states, outputs, and other model data.

Each component exposes a ``data`` attribute, which holds the underlying model data
(e.g. :class:`dict`, :class:`xarray.Dataset`, or :class:`geopandas.GeoDataFrame`),
and supports a common set of I/O and manipulation methods such as
:meth:`read`, :meth:`write`, and :meth:`set`.

Component-level API
-------------------

The table below summarizes the important methods and attributes for each component.
These allow for fine-grained reading, writing, modification, and inspection of component data.
They are particularly useful when working interactively in Python, for example when updating
specific configuration parameters, or inspecting data.


.. autosummary::
   :toctree: _generated
   :nosignatures:
   :template: component-class.rst

   DelwaqConfigComponent
   DemissionConfigComponent
   DelwaqStaticdataComponent
   DelwaqPointerComponent
   DemissionGeometryComponent
   DelwaqForcingComponent
   DemissionForcingComponent
   DelwaqHydromapsComponent
