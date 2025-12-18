"""Component for D-Emission geometry file."""

import logging
import os
from os.path import join

import numpy as np
import pandas as pd
from hydromt import hydromt_step
from hydromt.model import Model
from hydromt.model.components import ModelComponent

__all__ = ["DemissionGeometryComponent"]

logger = logging.getLogger(f"hydromt.{__name__}")


class DemissionGeometryComponent(ModelComponent):
    """Demission geometry component.

    Inherits from the HydroMT-core ModelComponent.
    Handles the Demission geometry file as a Pandas DataFrame.

    Columns of the DataFrame are:

    * TotArea: area of the segment
    * fPaved: fraction of the segment that is paved
    * fUnpaved: fraction of the segment that is unpaved
    * fOpenWater: fraction of the segment that is open water

    """

    def __init__(
        self,
        model: Model,
        *,
        filename: str = "config/B7_geometry.inc",
    ):
        """Initialize a DemissionGeometryComponent.

        Parameters
        ----------
        model : Model
            HydroMT model instance
        filename : str
            The path to use for reading and writing of component data by default.
            By default "config/B7_geometry.inc".
        """
        self._data: pd.DataFrame | None = None
        self._filename: str = filename

        super().__init__(model=model)

    @property
    def data(self) -> pd.DataFrame:
        """Return the geometry DataFrame."""
        if self._data is None:
            self._initialize()
        return self._data

    def _initialize(self, skip_read=False) -> None:
        """Initialize the model geometry."""
        if self._data is None:
            self._data = pd.DataFrame()
            if self.root.is_reading_mode() and not skip_read:
                self.read()

    def set(self, data: pd.DataFrame) -> None:
        """
        Add geometry data.

        Parameters
        ----------
        data : pd.DataFrame
            DataFrame with geometry data.

        """
        self._initialize()

        # Check that required columns are present
        required_columns = {"TotArea", "fPaved", "fUnpaved", "fOpenWater"}
        if not required_columns.issubset(data.columns):
            logger.warning(
                f"Geometry data must contain the following columns: {required_columns}"
            )
            return

        # Check if _data is empty
        if not self._data.empty:
            logger.info("Overwriting existing geometry data in the component.")
        self._data = data

    @hydromt_step
    def read(self) -> None:
        """Read geometry file in ASCII format at <root>/filename."""
        self.root._assert_read_mode()
        self._initialize(skip_read=True)

        filepath = join(self.root.path, self._filename)

        if not os.path.isfile(filepath):
            logger.warning(f"Geometry file {filepath} does not exist.")
            return

        # Read the geometry file into a DataFrame
        logger.info(f"Reading geometry data from {filepath}.")

        # Find where to start reading the file (skip header lines)
        with open(filepath, "r") as f:
            lines = f.readlines()

        for i, line in enumerate(lines):
            if line == "DATA\n":
                start_line = i + 1
                break

        data = pd.read_csv(
            filepath,
            skiprows=start_line,
            sep="\s+",
            names=["TotArea", "fPaved", "fUnpaved", "fOpenWater"],
            dtype={
                "TotArea": float,
                "fPaved": float,
                "fUnpaved": float,
                "fOpenWater": float,
            },
        )

        self.set(data)

    @hydromt_step
    def write(self):
        """Write geometry at <root/config> in ASCII and binary format."""
        self.root._assert_write_mode()
        if not self.data.empty:
            logger.info("No geometry data found, skip writing.")
            return

        logger.info("Writing geometry file in root/config")
        fname = os.path.splitext(join(self.root.path, self._filename))[0]

        # Write ASCII file
        exfile = open(fname + ".inc", "w")
        print(";Geometry of the EM compartment", file=exfile)
        print("PARAMETERS TotArea fPaved fUnpaved fOpenWater ALL", file=exfile)
        print("DATA", file=exfile)
        np.savetxt(exfile, self._geometry.values, fmt="%10.4f")
        exfile.close()

        # Write binary file
        # Flatten the geometry data and repeat them for each compartment
        geometry_data = np.tile(self._geometry.values.flatten(), 1)
        artow = np.array(geometry_data, dtype=np.float32).copy()
        # Define dummy time
        timear = np.array(0, dtype=np.int32)
        # Open and write the data
        fp = open(fname + ".bin", "wb")
        tstr = timear.tobytes() + artow.tobytes()
        fp.write(tstr)
        fp.close()

        # Write corresponding def file of the geometry
        fpa = open(fname + "-parameters.inc", "w")
        print("PARAMETERS TotArea fPaved fUnpaved fOpenWater", file=fpa)
        fpa.close()
