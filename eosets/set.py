# Copyright 2025, SERTIT-ICube - France, https://sertit.unistra.fr/
# This file is part of eosets project
#     https://github.com/sertit/eosets
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Super class of all sets implemented in EOSets"""

import contextlib
import logging
import os
import shutil
import tempfile
from abc import abstractmethod
from enum import unique
from typing import Any, Union

import geopandas as gpd
import xarray as xr
from eoreader.bands import BandType, to_str
from eoreader.env_vars import CI_EOREADER_BAND_FOLDER
from eoreader.products import Product, SensorType
from sertit import AnyPath, files, path
from sertit.misc import ListEnum
from sertit.types import AnyPathStrType, AnyPathType, AnyXrDataStructure

from eosets import EOSETS_NAME
from eosets.env_vars import CI_EOSETS_BAND_FOLDER
from eosets.utils import BandsType

LOGGER = logging.getLogger(EOSETS_NAME)


@unique
class GeometryCheck(ListEnum):
    """Available geometry checks."""

    FOOTPRINT = "footprint"
    """ Ensure the checks are done regarding the footprints."""

    EXTENT = "extent"
    """ Ensure the checks are done regarding the extents."""

    NONE = "none"
    """ No geometric check will be applied."""


GeometryCheckType = Union[GeometryCheck, str]


class Set:
    """Abstract class of set. Basically implementing output management"""

    def __init__(
        self,
        output_path: AnyPathStrType = None,
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
        # Don't use output properties here, subsets are not yet initialized!
        if output_path:
            self._tmp_output = None
            self._output = AnyPath(output_path)
        else:
            self._tmp_output = tempfile.TemporaryDirectory()
            self._output = AnyPath(self._tmp_output.name)

        self.id: str = id
        """ ID of the reference product, given by the creator of the mosaic. If not, a mix based on the dates and constellations of its components. """

        # Nodata (by default use EOReader's)
        self.nodata = None
        """ Nodata of the mosaic. If not provided in kwargs, using the first product's nodata. """

        # Pixel size
        self.pixel_size = None
        """ Pixel size of the set. If not provided in kwargs, using the first product's pixel size. """

        self.crs = None
        """ CRS of the mosaic. If not provided in kwargs, using the first product's crs. """

        self.same_constellation = None
        """ Is the mosaic constituted of the same constellation? """

        self.same_crs = None
        """ Is the mosaic constituted of the same sensor type? """

        self.constellations = None
        """ List of unique constellations constituting the set """

        self.nof_prods: int = 0
        """ Number of products. """

        self.is_sar: bool = False
        """ All products of this set are SAR data. """

        self.is_optical: bool = False
        """ All products of this set are Optical data. """

        self.nof_prods: int = 0
        """ Number of products. """

        # Set tmp process
        self._set_tmp_process()

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
        Output directory of the set

        Returns:
            AnyPathType: Output path of the set
        """
        return self._output

    def _set_tmp_process(self):
        """Set temporary folder avoiding recursive tmps."""
        if self._output.name == "tmp":
            # Avoid nested "tmp" folders
            self._tmp_process = self._output
        else:
            self._tmp_process = self._output.joinpath("tmp")

        os.makedirs(self._tmp_process, exist_ok=True)

    @output.setter
    def output(self, value: AnyPathStrType) -> None:
        """
        Output directory of the set

        Args:
            value (AnyPathStrType): Output path of the set
        """
        # Set the new output
        self._output = AnyPath(value)
        if not path.is_cloud_path(self._output):
            self._output = self._output.resolve()

        # Create temporary process folder
        old_tmp_process = self._tmp_process
        self._set_tmp_process()

        # Update for every sets
        self._manage_output()

        # Move all files from old process folder into the new one
        for file in path.listdir_abspath(old_tmp_process):
            # Don't overwrite file
            with contextlib.suppress(shutil.Error):
                shutil.move(str(file), self._tmp_process)

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
            ci_tmp_folder = os.environ.get(
                CI_EOSETS_BAND_FOLDER, os.environ.get(CI_EOREADER_BAND_FOLDER)
            )
            if ci_tmp_folder:
                ci_tmp_folder = AnyPath(ci_tmp_folder)
                if ci_tmp_folder.is_dir():
                    # If we need a writable directory, check it
                    tmp_folder = ci_tmp_folder

        return tmp_folder

    def _get_out_path(self, filename: str) -> tuple[AnyPathType, bool]:
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

    @abstractmethod
    def get_prods(self) -> list:
        """
        Get all the products as a list.

        Returns:
            list: Products list
        """
        raise NotImplementedError

    def get_attr(self, attr: str, **kwargs) -> Any:
        """
        Get attribute, either from kwargs or from the first product (default)

        Args:
            attr (str): Wanted attribute
            **kwargs: Other args

        Returns:
            Any: Attribute result
        """
        attr = kwargs.pop(attr, getattr(self.get_first_prod(), attr))
        if callable(attr):
            attr = attr()

        return attr

    def is_homogeneous(self, attr: str) -> bool:
        """
        Check if the given attribute is the same for all products constituting the mosaic.

        Args:
            attr (str): Attribute to be checked. Must be available in EOReader's Product

        Returns:
            bool: True if this attribute is the same for all products constituting the mosaic.
        """
        ref_attr = getattr(self.get_first_prod(), attr)

        if self.nof_prods > 1:
            if callable(ref_attr):
                is_homogeneous = all(
                    ref_attr() == getattr(sec, attr)() for sec in self.get_prods()[1:]
                )
            else:
                is_homogeneous = all(
                    ref_attr == getattr(sec, attr) for sec in self.get_prods()[1:]
                )
        else:
            is_homogeneous = True

        return is_homogeneous

    def _has_only_sar(self, **kwargs):
        """Check if the set has only SAR products."""
        return (
            self.same_constellation
            and self.get_attr("sensor_type", **kwargs) == SensorType.SAR
        )

    def _has_only_optical(self, **kwargs):
        """Check if the set has only Optical products."""
        return (
            self.same_constellation
            and self.get_attr("sensor_type", **kwargs) == SensorType.OPTICAL
        )

    def post_init(self, **kwargs):
        """Post initialization as the set level."""
        # Constellations
        self.same_constellation = self.is_homogeneous("constellation")
        self.constellations = list(set(prod.constellation for prod in self.get_prods()))

        # CRS
        self.crs = self.get_attr("crs", **kwargs)
        self.same_crs = self.is_homogeneous("crs")

        # Other
        self.nodata = self.get_attr("nodata", **kwargs)
        self.pixel_size = self.get_attr("pixel_size", **kwargs)

        # Product types
        self.is_sar = self._has_only_sar(**kwargs)
        self.is_optical = self._has_only_optical(**kwargs)

    def get_first_prod(self) -> Product:
        """
        Get first product, which should be coherent with all others

        Returns:
            Product: First reference product
        """
        return self.get_prods()[0]

    @abstractmethod
    def read_mtd(self):
        """"""
        # TODO: how ? Just return the fields that are shared between set's components ? Or create a XML from scratch ?
        raise NotImplementedError

    @abstractmethod
    def footprint(self) -> gpd.GeoDataFrame:
        """
        Get the footprint of the set.

        Returns:
            gpd.GeoDataFrame: Footprint of the set
        """
        raise NotImplementedError

    @abstractmethod
    def extent(self) -> gpd.GeoDataFrame:
        """
        Get the extent of the set.

        Returns:
            gpd.GeoDataFrame: Extent of the set

        """
        raise NotImplementedError

    def has_band(self, band: BandType) -> bool:
        """
        Does this moasic have products with the specified band ?

        By band, we mean:

        - satellite band
        - index
        - DEM band
        - cloud band

        Args:
            band (BandType): EOReader band (optical, SAR, clouds, DEM)

        Returns:
            bool: True if the products has the specified band
        """
        return all(prod.has_band(band) for prod in self.get_prods())

    def has_bands(self, bands: BandsType) -> bool:
        """
        Does this moasic have products with the specified bands ?

        By band, we mean:

        - satellite band
        - index
        - DEM band
        - cloud band

        See :code:`has_band` for a code example.

        Args:
            bands (BandsType): EOReader bands (optical, SAR, clouds, DEM)

        Returns:
            bool: True if the products has the specified band
        """

        return all(prod.has_bands(bands) for prod in self.get_prods())

    def _update_attrs(
        self, xarr: AnyXrDataStructure, bands: list, **kwargs
    ) -> AnyXrDataStructure:
        """
        Update attributes of the given array
        Args:
            xarr (AnyXrDataStructure): Array whose attributes need an update
            bands (list): Bands
        Returns:
            AnyXrDataStructure: Updated DataArray/Dataset
        """
        # Clean attributes, we don't want to pollute our attributes by default ones (not deterministic)
        # Are we sure of that ?
        xarr.attrs = {}

        if not isinstance(bands, list):
            bands = [bands]

        long_name = to_str(bands)
        xr_name = "_".join(long_name)
        attr_name = " ".join(long_name)

        with contextlib.suppress(ValueError):
            xarr = xarr.rename(xr_name)

        xarr.attrs["long_name"] = attr_name
        xarr.attrs["condensed_name"] = self.condensed_name

        # TODO: complete that
        return self._update_attrs_constellation_specific(xarr, bands, **kwargs)

    def _update_attrs_constellation_specific(
        self, xarr: xr.DataArray, bands: list, **kwargs
    ) -> xr.DataArray:
        """
        Update attributes of the given array (constellation specific)

        Args:
            xarr (xr.DataArray): Array whose attributes need an update
            bands (list): Array name (as a str or a list)

        Returns:
            xr.DataArray: Updated array/dataset
        """
        raise NotImplementedError

    def _update_xds_attrs(
        self,
        xds: xr.Dataset,
        bands: BandsType,
        **kwargs,
    ) -> xr.Dataset:
        """ """
        # Rename all bands and add attributes
        for key, val in xds.items():
            xds[key] = self._update_attrs(val, key, **kwargs)

        # Update stack's attributes
        if len(xds) > 0:
            xds = self._update_attrs(xds, bands, **kwargs)

        return xds

    def get_bands_to_load(self, bands, out_suffix="tif") -> (list, dict):
        # Get the bands to be loaded
        bands_path = {}
        bands_to_load = []
        for band in bands:
            band_path, exists = self._get_out_path(
                f"{self.id}_{to_str(band)[0]}.{out_suffix}"
            )
            bands_path[band] = band_path
            if not exists:
                bands_to_load.append(band)

        return bands_to_load, bands_path
