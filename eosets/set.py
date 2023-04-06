import os
import shutil
import tempfile
from abc import abstractmethod
from typing import Tuple, Union

from cloudpathlib import AnyPath, CloudPath
from sertit import files

from eosets.env_vars import CI_EOSETS_BAND_FOLDER
from eosets.utils import AnyPathType


class Set:
    """Abstract class of set. Basically implementing output management"""

    def __init__(
        self,
        paths: Union[list, str, AnyPathType],
        output_path: Union[str, AnyPathType] = None,
        id: str = None,
        remove_tmp: bool = True,
        **kwargs,
    ):
        # Manage output
        # TODO : create a temp folder for the pairs ?
        """Output path of the pairs."""

        # Remove temporary files
        self._tmp_output = None
        self._output = None
        self._remove_tmp = remove_tmp
        """ Remove temporary files, propagated to EOReader's Products. """

        # -- Other parameters --
        # Full name
        self.full_name: str = ""
        """ Mosaic full name. """

        # Condensed name
        self.condensed_name = ""
        """ Mosaic condensed name, a mix based on the dates and constellations of the components of the mosaic. """

        # Manage output path
        if output_path:
            self._tmp_output = None
            self._output = AnyPath(output_path)
        else:
            self._tmp_output = tempfile.TemporaryDirectory()
            self._output = AnyPath(self._tmp_output.name)

        self._tmp_process = self.output.joinpath("tmp_mosaic")
        os.makedirs(self._tmp_process, exist_ok=True)

        self.id: str = id
        """ ID of the reference product, given by the creator of the mosaic. If not, a mix based on the dates and constellations of its components. """

        # Nodata (by default use EOReader's)
        self.nodata = None
        """ Nodata of the mosaic. If not provided in kwargs, using the first product's nodata. """

        # Resolution (by default use EOReader's)
        self.resolution = None
        """ Resolution of the mosaic. If not provided in kwargs, using the first product's resolution. """

        self.crs = None
        """ CRS of the mosaic. If not provided in kwargs, using the first product's crs. """

        self.same_constellation = None
        """ Is the mosaic constituted of the same constellation? """

        self.same_crs = None
        """ Is the mosaic constituted of the same sensor type? """

        self.constellations = None
        """ List of unique constellations constituting the set """

    @abstractmethod
    def clean_tmp(self):
        """
        Clean the temporary directory of the current product
        """
        raise NotImplementedError

    @abstractmethod
    def clear(self):
        """
        Clear this product's cache
        """
        raise NotImplementedError

    @abstractmethod
    def _manage_output(self):
        """
        Manage the output specifically for this child class
        """
        raise NotImplementedError

    def __del__(self):
        """Cleaning up _tmp directory"""
        self.clear()

        # -- Remove temp folders
        if self._tmp_output:
            self._tmp_output.cleanup()

        elif self._remove_tmp:
            files.remove(self._tmp_process)
            self.clean_tmp()

    @property
    def output(self) -> AnyPathType:
        """
        Output directory of the mosaic

        Returns:
            AnyPathType: Output path ofthe mosaic
        """
        return self._output

    @output.setter
    def output(self, value: Union[str, AnyPathType]) -> None:
        """
        Output directory of the mosaic

        Args:
            value (Union[str, AnyPathType]): Output path ofthe mosaic
        """
        # Set the new output
        self._output = AnyPath(value)
        if not isinstance(self._output, CloudPath):
            self._output = self._output.resolve()

        # Create temporary process folder
        old_tmp_process = self._tmp_process
        self._tmp_process = self._output.joinpath(f"tmp_{self.condensed_name}")
        os.makedirs(self._tmp_process, exist_ok=True)

        # Update for prods
        self._manage_output()

        # Move all files from old process folder into the new one
        for file in files.listdir_abspath(old_tmp_process):
            try:
                shutil.move(str(file), self._tmp_process)
            except shutil.Error:
                # Don't overwrite file
                pass

        # Remove old output if existing into the new output
        if self._tmp_output:
            self._tmp_output.cleanup()
            self._tmp_output = None

    def _get_tmp_folder(self, writable: bool = False) -> AnyPathType:
        """
        Manage the case of CI bands

        Returns:
            AnyPathType : Band folder
        """
        tmp_folder = self._tmp_process

        # Manage CI bands (when we do not write anything, read only)
        if not writable:
            ci_tmp_folder = os.environ.get(CI_EOSETS_BAND_FOLDER)
            if ci_tmp_folder:
                ci_tmp_folder = AnyPath(ci_tmp_folder)
                if ci_tmp_folder.is_dir():
                    # If we need a writable directory, check it
                    tmp_folder = ci_tmp_folder

        return tmp_folder

    def _get_out_path(self, filename: str) -> Tuple[AnyPathType, bool]:
        """
        Returns the output path of a file to be written, depending on if it already exists or not (manages CI folders)

        Args:
            filename (str): Filename

        Returns:
            Tuple[AnyPathType , bool]: Output path and if the file already exists or not
        """
        out = self._get_tmp_folder() / filename
        exists = True
        if not out.exists():
            exists = False
            out = self._get_tmp_folder(writable=True) / filename

        return out, exists
