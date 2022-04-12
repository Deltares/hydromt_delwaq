.. _intro_getting_started:

Getting Started
===============

.. grid:: 2
    :gutter: 1 

    .. grid-item-card:: 
        :text-align: center
        :link: installation_guide
        :link-type: ref
        
        :octicon:`settings;10em`
        +++
        Installation guide

    .. grid-item-card:: 
        :text-align: center
        :link: intro_user_guide
        :link-type: ref
        
        :octicon:`book;10em`
        +++
        User guide

    .. grid-item-card:: 
        :text-align: center
        :link: api_reference
        :link-type: ref
        
        :octicon:`list-unordered;10em`
        +++
        API reference

    .. grid-item-card:: 
        :text-align: center
        :link: examples
        :link-type: ref
                
        :octicon:`graph;10em`
        +++
        Examples

    .. grid-item-card:: 
        :text-align: center
        :link: user_stories
        :link-type: ref
        
        :octicon:`file-media;10em`
        +++
        User stories

With the **hydromt_delwaq plugin**, users can easily benefit from the rich set of tools of the 
`HydroMT package <https://github.com/Deltares/hydromt>`_ to build and update 
DELWAQ models from available global and local data.

This plugin assists the modeller in:

- Quickly setting up a model and default parameter values
- Making maximum use of the best available global or local data
- Adjusting and updating components of a delwaq model and their associated parameters in a consistent way

Water quality and emission, fate and transport of pollutants through the landscape, rivers and oceans is very much linked to the water movements. 
For this reason, hydroMT strongly links DELWAQ modelbuilding to an underlying hydrology / hydraulic or hydrodynamic model. As of now, only the link 
between ** :ref:`Wflow and DELWAQ <coupling_wflow>` ** is supported.

This plugin supports both:

- emission models, ** `D-Emission <https://www.deltares.nl/en/software/module/D-Emissions/>`_ ** (also referred to as EM in hydroMT), that predicts the quantity (fluxes) of substances being released from sources of pollutants to the surface waters (mass/time).
- fate and transport models, ** `D-Water Quality <https://www.deltares.nl/en/software/module/d-water-quality/>`_ ** (also referred to as WQ in hydroMT), that predicts the fate and transport of substances in the surface waters (concentrations in mass / volume).

.. image:: ../img/D-emissions.png

License
-------

Copyright (c) 2021, Deltares

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public 
License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You can find the full terms of the GNU General Public License at <https://www.gnu.org/licenses/>.

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Table of Contents

   installation.rst
   ../examples/index.rst
   user_stories.rst