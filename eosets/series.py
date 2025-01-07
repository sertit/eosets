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
"""Class implementing the series object"""

import logging
import os
from collections import defaultdict
from enum import unique
from typing import Union

import geopandas as gpd
import numpy as np
import xarray as xr
from eoreader import cache
from eoreader.bands import to_band, to_str
from eoreader.reader import Reader
from eoreader.utils import UINT16_NODATA
from rasterio.enums import Resampling
from sertit import AnyPath, types
from sertit.misc import ListEnum
from sertit.types import AnyPathStrType

from eosets import EOSETS_NAME
from eosets.exceptions import IncompatibleProducts
from eosets.mosaic import Mosaic
from eosets.set import GeometryCheck, GeometryCheckType, Set
from eosets.utils import BandsType, read, stack, write

LOGGER = logging.getLogger(EOSETS_NAME)
READER = Reader()


@unique
class Alignment(ListEnum):
    """Available alignment methods for series."""

    FIRST = "first"
    """
    Align the series on the first mosaic.
    """

    LAST = "last"
    """
    Align the series on the last mosaic.
    """

    MEAN = "mean"
    """
    Align the series on an hypothetical mean mosaic. Not yet implemented.
    """

    CUSTOM = "custom"
    """
    Align the series on a custom mosaic. Not yet implemented.
    """

    EXTERNAL = "external"
    """
    Align the series on an external mosaic. Not yet implemented.
    """


class Series(Set):
    """Class of series"""

    def __init__(
        self,
        paths: list,
        id: str = None,
        output_path: AnyPathStrType = None,
        remove_tmp: bool = True,
        overlap_check: GeometryCheckType = GeometryCheck.EXTENT,
        contiguity_check: GeometryCheckType = GeometryCheck.EXTENT,
        alignement: Union[Alignment, str] = Alignment.FIRST,
        coregister: bool = False,
        reference_mosaic: Union[Mosaic, int, str] = None,
        **kwargs,
    ):
        # Manage mosaics
        self.mosaics = None
        """ Mosaics composing the series. """

        self.id = None
        """ ID of the series """

        self.alignment = Alignment.convert_from(alignement)[0]
        """ Series alignment (on first mosaic, on last mosaic, on a hypothetical mean product). """

        self.coregister = coregister
        """ Do we need to coregister the time series? """

        self.reference_mosaic = reference_mosaic
        """ Reference mosaic """

        self._unique_mosaic = None

        contiguity_check = GeometryCheck.convert_from(contiguity_check)[0]
        overlap_check = GeometryCheck.convert_from(overlap_check)[0]

        # Init the base class
        super().__init__(
            output_path,
            id,
            remove_tmp,
            **kwargs,
        )

        # Update mosaics of the pair
        self._manage_mosaics(paths, contiguity_check, overlap_check)

        # Post init at the set level
        self.post_init(**kwargs)

    def clean_tmp(self):
        """
        Clean the temporary directory of the current series
        """
        for mos in self.mosaics:
            mos.clean_tmp()

    def clear(self):
        """
        Clear this series' cache
        """
        # Delete all cached properties and functions
        for mos in self.mosaics:
            mos.clear()

    def _manage_output(self):
        """
        Manage the output specifically for this child class
        """
        for mos in self.mosaics:
            mos.output = self.output

    def get_prods(self) -> list:
        """
        Get all the products as a list.

        Returns:
            list: Products list
        """
        prods = []
        for mos in self.mosaics:
            prods += mos.get_prods()

        return prods

    def _manage_mosaics(
        self,
        mosaic_paths: list,
        contiguity_check: GeometryCheck = GeometryCheck.EXTENT,
        overlap_check: GeometryCheck = GeometryCheck.EXTENT,
    ) -> None:
        """

        Check if all the mosaics are overlapping and have the same CRS. Each mosaic of a series should have a different datetime.
        If not, throws a IncompatibleProducts error.
        Args:
            mosaic_paths (list): List of paths of the mosaics composing the series
            contiguity_check (GeometryCheck): Check regarding the contiguity of the products of the mosaics
            overlap_check (GeometryCheck): Check regarding the overlapping of the mosaics of the series

        Raises:
            IncompatibleProducts: Incompatible products if not contiguous or not the same date

        """
        # Information regarding the pair composition
        assert types.is_iterable(
            mosaic_paths
        ), "You should give a list of paths or Mosaics to create your Series!"

        self.mosaics: list[Mosaic] = []
        for paths in mosaic_paths:
            # Open the mosaic
            if isinstance(paths, Mosaic):
                mos = paths
            else:
                mos = Mosaic(
                    paths,
                    output_path=self.output,
                    remove_tmp=self._remove_tmp,
                    contiguity_check=contiguity_check,
                )

            # Check datetime compatibility
            if all(other_mos.datetime != mos.datetime for other_mos in self.mosaics):
                self.mosaics.append(mos)
            else:
                raise IncompatibleProducts(
                    "All mosaics of a series should have a different datetime! If not, please regroup them in a unique mosaic."
                )

        # Check geometry
        if overlap_check != GeometryCheck.NONE:
            reference_geom: gpd.GeoDataFrame = getattr(
                self.reference_mosaic, str(overlap_check.value)
            )()
            for mos in self.mosaics:
                if mos.id != self.reference_mosaic.id:
                    mos_geom: gpd.GeoDataFrame = getattr(
                        mos, str(overlap_check.value)
                    )()
                    if not reference_geom.intersects(
                        mos_geom.to_crs(self.reference_mosaic.crs)
                    ).all():
                        raise IncompatibleProducts("All mosaics should overlap!")

        # Fill other attributes
        self.nof_prods = len(self.get_prods())
        self.id = "_".join(mos.id for mos in self.mosaics)
        self.full_name = "_".join(mos.full_name for mos in self.mosaics)
        self.condensed_name = self.full_name

        self._unique_mosaic = len(self.mosaics) == 1
        # TODO (how to name series ???)

    @property
    def reference_mosaic(self) -> Mosaic:
        """
        Get the reference mosaic of the series

        Returns:
            Mosaic: Output path ofthe mosaic
        """
        if self.alignment in [Alignment.FIRST, Alignment.LAST]:
            return self.mosaics[self._reference_mosaic]
        elif self.alignment == Alignment.MEAN:
            raise NotImplementedError
        elif self.alignment == Alignment.EXTERNAL:
            return self._reference_mosaic
        elif self.alignment == Alignment.CUSTOM and isinstance(
            self._reference_mosaic, int
        ):
            return self.mosaics[self._reference_mosaic]

    @reference_mosaic.setter
    def reference_mosaic(self, mosaic: Union[Mosaic, int] = None):
        if self.alignment == Alignment.FIRST:
            self._reference_mosaic = 0
        elif self.alignment == Alignment.LAST:
            self._reference_mosaic = -1
        elif self.alignment == Alignment.MEAN:
            raise NotImplementedError
        elif self.alignment == Alignment.EXTERNAL:
            assert mosaic is not None
            self._reference_mosaic = mosaic
        elif self.alignment == Alignment.CUSTOM:
            assert isinstance(mosaic, (Mosaic, int))
            if isinstance(mosaic, int):
                assert abs(mosaic) <= len(self.mosaics)
            self._reference_mosaic = mosaic

    def read_mtd(self):
        """Read the pair's metadata, but not implemented for now."""
        # TODO: how ? Just return the fields that are shared between series' components ? Or create a XML from scratch ?
        raise NotImplementedError

    @cache
    def footprint(self) -> gpd.GeoDataFrame:
        """
        Get the footprint of the series, i.e. the intersection between all mosaics in the reference mosaic's CRS.

        Returns:
            gpd.GeoDataFrame: Footprint of the mosaic
        """
        footprint = None
        for mos in self.mosaics:
            geom: gpd.GeoDataFrame = mos.footprint().to_crs(self.reference_mosaic.crs)
            if footprint is None:
                footprint = geom
            else:
                footprint = footprint.overlay(geom, "intersection")

        return footprint

    @cache
    def extent(self) -> gpd.GeoDataFrame:
        """
        Get the extent of the series, i.e. the intersection between all mosaics in the reference mosaic's CRS.

        Returns:
            gpd.GeoDataFrame: Extent of the mosaic

        """
        extent = None
        for mos in self.mosaics:
            geom: gpd.GeoDataFrame = mos.footprint().to_crs(self.reference_mosaic.crs)
            extent = geom if extent is None else extent.overlay(geom, "intersection")

        return extent

    def load(
        self,
        bands: BandsType = None,
        pixel_size: float = None,
        resampling: Resampling = Resampling.bilinear,
        **kwargs,
    ) -> xr.Dataset:
        """

        Load the bands and compute the wanted spectral indices.

        Args:
            bands (BandsType): Wanted bands
            pixel_size (float): Pixel size of the returned Dataset. If not specified, use the mosaic's pixel size.
            resampling (Resampling): Resampling method
            **kwargs: Other arguments used to load bands

        Returns:
            xr.Dataset: Wanted bands as xr.Datasets
        """
        bands = to_band(bands)

        # Load bands
        window = kwargs.pop("window", self.footprint())

        # Load reference mosaic bands
        reference_ds = self.reference_mosaic.load(
            bands, pixel_size=pixel_size, window=window, **kwargs
        ).expand_dims({"time": [self.reference_mosaic.datetime]}, axis=-1)

        # Load mosaic bands
        arr_dict = defaultdict(list)
        for mos in self.mosaics:
            # Get mosaic datetime
            dt = mos.datetime

            # Manage if band to be loaded or already on memory
            if mos.id != self.reference_mosaic.id:
                # Load bands for new mosaic
                mos_ds = mos.load(bands, pixel_size=pixel_size, window=window, **kwargs)

                # Add the bands to the dataset
                for band in bands:
                    mos_arr = mos_ds[band]
                    reference_arr = reference_ds[band]

                    # To be sure, always collocate arrays, even if the size is the same
                    # Indeed, a small difference in the coordinates will lead to empy arrays
                    # So the bands MUST BE exactly aligned
                    mos_arr = (
                        mos_arr.rio.reproject_match(
                            reference_arr, resampling=resampling, **kwargs
                        )
                        .expand_dims({"time": [dt]}, axis=-1)
                        .assign_coords({"time": [dt]})
                    )

                    # Save mosaics
                    arr_dict[band].append(mos_arr)
            else:
                # We already have the bands in memory, just assign time to every arrays and save it
                for band in bands:
                    mos_arr = reference_ds[band].assign_coords({"time": [dt]})
                    arr_dict[band].append(mos_arr)

        # Create empty dataset with correct coordinates
        series_ds = xr.Dataset(
            coords={
                "x": reference_ds.x,
                "y": reference_ds.y,
                "band": np.arange(1, len(self.mosaics) + 1),  # One array per band
                "time": [mos.datetime for mos in self.mosaics],
            }
        )

        # Make sure the dataset has the bands in the right order -> re-order the input dict
        series_ds = xr.merge([xr.merge(arr_dict[band]) for band in bands])

        # Update attributes
        series_ds = self._update_xds_attrs(series_ds, bands)

        return series_ds

    def stack(
        self,
        bands: list,
        pixel_size: float = None,
        stack_path: AnyPathStrType = None,
        save_as_int: bool = False,
        **kwargs,
    ) -> xr.DataArray:
        """
        Stack bands and index of a series.

        Args:
            bands (list): Bands and index combination
            pixel_size (float): Stack pixel size. If not specified, use the product pixel size.
            stack_path (AnyPathStrType): Stack path
            save_as_int (bool): Convert stack to uint16 to save disk space (and therefore multiply the values by 10.000)
            **kwargs: Other arguments passed to :code:`load` or :code:`rioxarray.to_raster()` (such as :code:`compress`)

        Returns:
            xr.DataArray: Stack as a DataArray
        """
        # Convert just in case
        if bands is None:
            bands = []
        bands = to_band(bands)

        if stack_path:
            stack_path = AnyPath(stack_path)
            if stack_path.is_file():
                return read(stack_path, pixel_size=pixel_size)
            else:
                os.makedirs(str(stack_path.parent), exist_ok=True)

        # Load all bands
        band_ds = self.load(bands, pixel_size=pixel_size, **kwargs)

        # Rename bands and remove time variable
        new_bands = []
        for band in bands:
            for dt in band_ds.time.values:
                new_bands.append(
                    f"{np.datetime_as_string(dt, unit='s')}_{to_str(band)[0]}"
                )

        # Stack bands
        if save_as_int:
            nodata = kwargs.get("nodata", UINT16_NODATA)
        else:
            nodata = kwargs.get("nodata", self.nodata)
        stk, dtype = stack(band_ds, save_as_int, nodata, **kwargs)

        # Update stack's attributes
        stk = self._update_attrs(stk, new_bands, **kwargs)

        # Write on disk
        if stack_path:
            LOGGER.debug("Saving stack")
            write(stk, stack_path, dtype=dtype, **kwargs)

        return stk

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
        # TODO: complete that

        return xarr
