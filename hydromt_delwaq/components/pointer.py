"""Component for Delwaq pointer data."""

import logging
import os
import struct
from os.path import join
from typing import Any

import numpy as np
from hydromt import hydromt_step
from hydromt.model import Model
from hydromt.model.components import ModelComponent

__all__ = ["DelwaqPointerComponent"]

logger = logging.getLogger(f"hydromt.{__name__}")


class DelwaqPointerComponent(ModelComponent):
    """Delwaq pointer component.

    Inherits from the HydroMT-core ModelComponent.
    Used to create the Delwaq pointer file and access some of its properties.

    Keys of the dictionary and their content are:

    - ``pointer``: pointer data as np.ndarray with 4 columns (ID, x, y, z)
    - ``nrofseg``: number of segments as int
    - ``nrofexch``: number of exchanges as int
    - ``surface_water``: name of the surface water compartment as str
    - ``boundaries``: list of the boundaries names as list of str
    - ``fluxes``: list of the flux names as list of str

    """

    def __init__(
        self,
        model: Model,
        *,
        filename: str = "config/B4_pointer.inc",
    ):
        """Initialize a DelwaqPointerComponent.

        Parameters
        ----------
        model : Model
            HydroMT model instance
        filename : str
            The path to use for reading and writing of component data by default.
            By default "config/B4_pointer.inc".
        """
        self._data: dict[str, Any] | None = None
        self._filename: str = filename

        super().__init__(model=model)

    @property
    def data(self) -> dict[str, Any]:
        """Return the pointer data dictionary."""
        if self._data is None:
            self._initialize()
        return self._data

    def _initialize(self, skip_read=False) -> None:
        """Initialize the model pointer."""
        if self._data is None:
            self._data = {}
            if self.root.is_reading_mode() and not skip_read:
                self.read()

    def set(self, name: str, value: Any):
        """
        Add model attribute property to pointer.

        Parameters
        ----------
        name : str
            Name of the pointer attribute to set.
        value : Any
            Value of the pointer attribute to set.

        """
        self._initialize()

        # Check that pointer attr is a four column np.array
        if name == "pointer":
            if not isinstance(value, np.ndarray) and value.shape[1] == 4:
                logger.warning(
                    "pointer values in self.pointer should be a"
                    "np.ndarray with four columns."
                )
                return
        elif np.isin(name, ["boundaries", "fluxes"]):
            if not isinstance(value, list):
                logger.warning(
                    f"{name} object in self.pointer should be a list of names."
                )
                return
        elif np.isin(name, ["surface_water"]):
            if not isinstance(value, str):
                logger.warning(f"{name} object in self.pointer should be a string.")
                return
        elif np.isin(name, ["nrofseg", "nrofexch"]):
            if not isinstance(value, int):
                logger.warning(f"{name} object in self.pointer should be an integer.")
                return

        # Add to pointer
        self._data[name] = value

    @hydromt_step
    def read(self):
        """Read pointer file in ASCII format at <root/config/B4_pointer.inc>."""
        self.root._assert_read_mode()
        self._initialize(skip_read=True)

        fname = self._filename.format()
        with open(join(self.root.path, fname), "r") as f:
            lines = f.readlines()

        pointer = []
        for line in lines:
            # skip commented lines
            if line.startswith(";"):
                continue
            parts = line.split()
            if len(parts) == 4:
                pointer.append(
                    [int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3])]
                )
            else:
                logger.warning(f"Skipping invalid line in pointer file: {line.strip()}")

        self.set("pointer", value=np.array(pointer))

    @hydromt_step
    def write(self):
        """Write pointer at <root/dynamicdata> in ASCII and binary format."""
        self.root._assert_write_mode()
        if not self.data or "pointer" not in self.data:
            logger.info("No pointer data found, skip writing.")
            return

        pointer = self.data["pointer"]

        logger.info("Writing pointer file in root/config")
        fname = join(self.root.path, "config", "B4_pointer")

        if not os.path.isdir(join(self.root.path, "config")):
            os.makedirs(join(self.root.path, "config"))

        # Write ASCII file
        exfile = open((fname + ".inc"), "w")
        print(";Pointer for WAQ simulation in Surface Water", file=exfile)
        print(";nr of pointers is: ", str(pointer.shape[0]), file=exfile)
        np.savetxt(exfile, pointer, fmt="%10.0f")
        exfile.close()

        # Write binary file
        f = open((fname + ".poi"), "wb")
        for i in range(pointer.shape[0]):
            f.write(struct.pack("4i", *np.int_(pointer[i, :])))
        f.close()

    def test_equal(self, other: ModelComponent) -> tuple[bool, dict[str, str]]:
        """Test if two pointer components are equal.

        Parameters
        ----------
        other : ModelComponent
            The component to compare against.

        Returns
        -------
        tuple[bool, dict[str, str]]
            True if the components are equal, and a dict with the associated errors per
            property checked.
        """
        errors: dict[str, str] = {}
        if not isinstance(other, self.__class__):
            errors["__class__"] = f"other does not inherit from {self.__class__}."

        # for once python does the recursion for us
        if not self.data == other.data:
            errors["config"] = "Configs are not equal"

        return len(errors) == 0, errors
