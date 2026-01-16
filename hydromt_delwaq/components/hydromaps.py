"""Custom Delwaq hydromaps component module."""

import glob
import logging
from os.path import join

import numpy as np
from hydromt import hydromt_step, readers
from hydromt.model import Model
from hydromt.model.components import GridComponent, ModelComponent
from hydromt_wflow.components.utils import test_equal_grid_data

__all__ = ["DelwaqHydromapsComponent"]

logger = logging.getLogger(f"hydromt.{__name__}")


class DelwaqHydromapsComponent(GridComponent):
    """Delwaq hydromaps component.

    Inherits from the HydroMT-core GridComponent model-component.
    It is used to store metadata from the gridded hydrological or hydraulic model Delwaq
    is linked to. Contains for example the segment ID, mask of active cells...
    """

    def __init__(
        self,
        model: Model,
        *,
        filename: str = "hydromodel/{name}.tif",
        region_component: str | None = None,
    ):
        """Initialize a DelwaqHydromapsComponent.

        Parameters
        ----------
        model : Model
            HydroMT model instance
        filename : str
            The path to use for reading and writing of component data by default.
            By default "hydromodel/{name}.tif".
        region_component : str, optional
            The name of the region component to use as reference for this component's
            region. If None, the region will be set to the grid extent.
        """
        super().__init__(
            model,
            filename=filename,
            region_component=region_component,
            region_filename=None,
        )

    ## I/O methods
    @hydromt_step
    def read(self, **kwargs):
        """
        Read hydromaps at <root/hydromodel> and parse to xarray.

        Keywords arguments are passed to the underlying hydromt.readers.open_mfraster
        function.
        """
        self.root._assert_read_mode()
        self._initialize_grid(skip_read=True)

        if "chunks" not in kwargs:
            kwargs.update(chunks={"y": -1, "x": -1})
        fns = glob.glob(join(self.root.path, "hydromodel", "*.tif"))
        if len(fns) > 0:
            ds_hydromaps = readers.open_mfraster(fns, **kwargs)

        # Load grid data in r+ mode to allow overwriting netcdf files
        if self.root.is_reading_mode() and self.root.is_writing_mode():
            ds_hydromaps.load()
        if ds_hydromaps.raster.crs is None and self.model.crs is not None:
            ds_hydromaps.raster.set_crs(self.model.crs)
        ds_hydromaps.coords["mask"] = ds_hydromaps["modelmap"].astype(bool)

        # When reading tif files, default lat/lon names are y/x
        # Rename to name used in grid for consistency
        if self._region_component is not None:
            region_grid = self._get_grid_data()
            if ds_hydromaps.raster.x_dim != region_grid.raster.x_dim:
                ds_hydromaps = ds_hydromaps.rename(
                    {ds_hydromaps.raster.x_dim: region_grid.raster.x_dim}
                )
            if ds_hydromaps.raster.y_dim != region_grid.raster.y_dim:
                ds_hydromaps = ds_hydromaps.rename(
                    {ds_hydromaps.raster.y_dim: region_grid.raster.y_dim}
                )

        # Set hydromaps data
        self.set(ds_hydromaps)

    @hydromt_step
    def write(self):
        """Write hydromaps at <root/hydromodel> in tif format."""
        self.root._assert_write_mode()
        if len(self.data) == 0:
            logger.info("No hydromaps data found, skip writing.")
            return

        ds_out = self.data
        logger.info("Writing hydromap files.")
        # Convert bool dtype before writing
        for var in ds_out.raster.vars:
            if ds_out[var].dtype == "bool":
                ds_out[var] = ds_out[var].astype(np.int32)
        ds_out.raster.to_mapstack(
            root=join(self.root.path, "hydromodel"),
            mask=True,
        )

    def test_equal(self, other: ModelComponent) -> tuple[bool, dict[str, str]]:
        """Test if two hydromaps components are equal.

        Checks the model component type as well as the data variables and their values.

        Parameters
        ----------
        other : ModelComponent
            The component to compare against.

        Returns
        -------
        tuple[bool, Dict[str, str]]
            True if the components are equal, and a dict with the associated errors per
            property checked.
        """
        errors: dict[str, str] = {}
        if not isinstance(other, self.__class__):
            errors["__class__"] = f"other does not inherit from {self.__class__}."

        # Check data
        _, data_errors = test_equal_grid_data(self.data, other.data)
        if len(data_errors) > 0:
            errors.update({"hydromaps": data_errors})

        return len(errors) == 0, errors
