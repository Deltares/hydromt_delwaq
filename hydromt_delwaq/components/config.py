"""Custom Delwaq config component module."""

import glob
import logging
import os
from os.path import basename, join, splitext
from pathlib import Path
from typing import List

from hydromt import hydromt_step
from hydromt.model import Model
from hydromt.model.components import ConfigComponent, ModelComponent

__all__ = ["DelwaqConfigComponent", "DemissionConfigComponent"]

logger = logging.getLogger(f"hydromt.{__name__}")


class DelwaqConfigComponent(ConfigComponent):
    """Manage the delwaq configuration input file for model simulations/settings.

    ``DelwaqConfigComponent`` data is stored in a dictionary. The component
    is used to prepare and update model simulations/settings of the delwaq model.
    """

    def __init__(
        self,
        model: Model,
        *,
        filename: str = "delwaq.inp",
    ):
        """Manage configuration files for delwaq simulations/settings.

        Parameters
        ----------
        model : Model
            HydroMT model instance
        filename : str
            A path relative to the root where the configuration files will
            be read and written if user does not provide a path themselves.
            By default 'delwaq.inp'
        """
        super().__init__(
            model,
            filename=filename,
        )

    ## I/O methods
    @hydromt_step
    def read(
        self,
        filename: str = None,
        skip: List[str] = [
            "B4_pointer",
            "B2_stations",
            "B2_stations-balance",
            "B2_monareas",
        ],
    ):
        """Read config files in ASCII format at <root/config>."""
        # Because of template config can be read also in write mode.
        # Use function from hydromt core to read the main config file
        super().read(path=filename)

        # Add the other files in the config folder
        config_fns = glob.glob(join(self.root.path, "config", "*.inc"))
        for fn in config_fns:
            name = splitext(basename(fn))[0]
            if name in skip:
                # Skip pointer file (should be read with read_pointer())
                # Skip monitoring files (should be read with read_monitoring())
                continue
            self.data[name] = dict()
            with open(fn) as f:
                for i, line in enumerate(f):
                    # Remove line breaks
                    self.data[name][f"l{i+1}"] = line.replace("\n", "")

    @hydromt_step
    def write(self):
        """Write config files in ASCII format at <root/config>."""
        self.root._assert_write_mode()

        if not self.data:
            logger.info(
                f"{self.model.name}.{self.name_in_model}: "
                "No config data found, skip writing."
            )
            return

        logger.info("Writing model config to file.")

        if not os.path.isdir(join(self.root.path, "config")):
            os.makedirs(join(self.root.path, "config"))

        for name, lines in self.data.items():
            fn_out = join(self.root.path, "config", f"{name}.inc")
            Path(fn_out).parent.mkdir(parents=True, exist_ok=True)

            with open(fn_out, "w") as exfile:
                for _, line in lines.items():
                    print(line, file=exfile)

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
        for name in other.data.keys():
            if name not in self.data:
                errors[name] = "Config file missing in self."
        for name in self.data.keys():
            if name not in other.data:
                errors[name] = "Config file not found in other."
            elif self.data[name] != other.data[name]:
                errors[name] = "Config file data not equal."

        if len(errors) == 0:
            errors = {}
        else:
            errors = {"config": errors}

        return len(errors) == 0, errors


class DemissionConfigComponent(DelwaqConfigComponent):
    """Manage the D-EMmission configuration input file for model simulations/settings.

    ``DemissionConfigComponent`` data is stored in a dictionary. The component
    is used to prepare and update model simulations/settings of the D-Emission model.
    """

    ## I/O methods
    @hydromt_step
    def read(
        self,
        filename: str = None,
        skip: List[str] = [
            "B7_geometry",
            "B2_stations",
            "B2_stations-balance",
            "B2_monareas",
        ],
    ):
        """Read config files in ASCII format at <root/config>."""
        # Because of template config can be read also in write mode.
        # Use function from hydromt core to read the main config file
        super().read(filename=filename, skip=skip)
