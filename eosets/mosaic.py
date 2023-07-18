# -*- coding: utf-8 -*-
# Copyright 2023, SERTIT-ICube - France, https://sertit.unistra.fr/
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
""" Class implementing the mosaic object """
import logging
import os
import shutil
from collections import defaultdict
from enum import unique
from glob import glob
from pathlib import Path
from typing import Union

import geopandas as gpd
import xarray as xr
from cloudpathlib import AnyPath, CloudPath
from eoreader import cache, utils
from eoreader.bands import BandNames, is_spectral_band, to_band, to_str
from eoreader.products import Product
from eoreader.reader import Reader
from eoreader.utils import UINT16_NODATA
from sertit import rasters
from sertit.misc import ListEnum

from eosets import EOSETS_NAME
from eosets.exceptions import IncompatibleProducts
from eosets.set import GeometryCheck, Set
from eosets.utils import AnyPathType

READER = Reader()

LOGGER = logging.getLogger(EOSETS_NAME)


@unique
class MosaicMethod(ListEnum):
    """Available mosaicing methods."""

    GTIFF = "merge_gtiff"
    VRT = "merge_vrt"


class Mosaic(Set):
    """Class of mosaic objetcs, composed of several contiguous EOReader's products acquired the same day."""

    def __init__(
        self,
        paths: Union[list, str, AnyPathType],
        output_path: Union[str, AnyPathType] = None,
        id: str = None,
        remove_tmp: bool = True,
        contiguity_check: Union[GeometryCheck, str] = GeometryCheck.EXTENT,
        mosaic_method: Union[MosaicMethod, str] = MosaicMethod.VRT,
        **kwargs,
    ):
        # Manage reference product
        self.prods: dict = {}
        """ Products (contiguous and acquired the same day). """

        # We need the date in _manage_prods
        self.date = None
        """ Date of the mosaic. If not provided in kwargs, using the first product's date. """

        self.datetime = None
        """ Datetime of the mosaic. If not provided in kwargs, using the first product's datetime. """

        self.mosaic_method = MosaicMethod.convert_from(mosaic_method)[0]
        """ Mosaicing method. If GTIFF is specified, the temporary files from every products will be removed, if VRT is spoecified, they will not."""

        contiguity_check = GeometryCheck.convert_from(contiguity_check)[0]

        # Init the base class
        super().__init__(
            output_path,
            id,
            remove_tmp,
            **kwargs,
        )
        # Update products of the mosaic
        self._manage_prods(paths, contiguity_check, **kwargs)

        # Fill attributes
        self.nodata = self.get_attr("nodata", **kwargs)
        self.pixel_size = self.get_attr("pixel_size", **kwargs)
        self.crs = self.get_attr("crs", **kwargs)
        self.same_constellation: bool = self.is_homogeneous("constellation")
        self.same_crs: bool = self.is_homogeneous("crs")
        self.constellations = list(set(prod.constellation for prod in self.get_prods()))

    def clean_tmp(self):
        """
        Clean the temporary directory of the current mosaic
        """
        for prod in self.get_prods():
            prod.clean_tmp()

    def clear(self):
        """
        Clear this mosaic's cache
        """
        # Delete all cached properties and functions
        for prod in self.get_prods():
            prod.clear()

    def _manage_output(self):
        """
        Manage the output specifically for this child class
        """
        for prod in self.get_prods():
            try:
                prod.output = self._get_tmp_folder(writable=True)
            except FileNotFoundError:
                # Never mind for non-existing files: they have already been copied :)
                pass

    def _manage_prods(
        self,
        paths: Union[list, str, Path, CloudPath],
        contiguity_check: GeometryCheck,
        **kwargs,
    ):
        """
        Manage products attributes and check the compatibility of the mosaic's components

        Args:
            paths (Union[list, str, Path, CloudPath]): Paths of the mosaic
            contiguity_check (GeometryCheck): Method to check the contiguity of the mosaic
            **kwargs: Other arguments

        Raises:
            IncompatibleProducts: Incompatible products if not contiguous or not the same date
        """
        # Manage reference product
        assert len(paths) > 0
        if not isinstance(paths, list):
            paths = [paths]

        # Remove EOReader tmp, the mosaic will save files in its tmp when needed
        remove_tmp = True

        # Nof prods (before checks)
        self.nof_prods = len(paths)

        # Open first product as a reference
        first_prod: Product = READER.open(
            paths[0],
            remove_tmp=remove_tmp,
            output_path=self._get_tmp_folder(writable=True),
            **kwargs,
        )
        if first_prod is None:
            raise ValueError(
                f"There is no existing products in EOReader corresponding to {paths[0]}"
            )

        self.prods[first_prod.condensed_name] = first_prod

        # Open others
        for path in paths[1:]:
            prod: Product = READER.open(
                path,
                remove_tmp=remove_tmp,
                output_path=self._get_tmp_folder(writable=True),
                **kwargs,
            )
            if prod is None:
                raise ValueError(
                    f"There is no existing products in EOReader corresponding to {path}"
                )

            # Ensure compatibility of the mosaic component, i.e. unique date and contiguous product
            self.prods[prod.condensed_name] = prod
            self.check_compatibility(first_prod, prod)
        self.check_contiguity(contiguity_check)

        # Create full_name
        self.date = kwargs.pop("date", first_prod.date)
        self.datetime = kwargs.pop("datetime", first_prod.datetime)
        self.full_name = (
            f"{'-'.join([prod.condensed_name for prod in self.get_prods()])}"
        )

        # Create condensed_name: [{date}-{sat_id}_]{???}
        # TODO: is it OK ?
        # TODO: if fixed date, change that
        # TODO: if all same constellation, set it only once
        # TODO: add sth ?
        self.condensed_name = f"{self.date.strftime('%Y%m%d')}_{'-'.join(list(set([prod.constellation_id for prod in self.get_prods()])))}"

        if self.id is None:
            self.id = self.condensed_name

    def get_prods(self) -> list:
        """
        Get all the products as a list.

        Returns:
            list: Products list
        """
        return list(self.prods.values())

    def check_compatibility(self, first_prod: Product, prod: Product) -> None:
        """
        Check if the mosaic products are coherent between each other.
        - Same sensor type
        - Same date

        TODO: same constellation ?

        If not, throws a IncompatibleProducts error.

        Args:
            first_prod(Product): First product, to be checked against
            prod (Product): Product to check

        Raises:
            IncompatibleProducts: Incompatible products if not contiguous or not the same date
        """
        # Check same sensor_type
        if first_prod.sensor_type != prod.sensor_type:
            raise IncompatibleProducts(
                f"Components of a mosaic should have the same sensor type! {first_prod.sensor_type.name=} != {prod.sensor_type.name=}"
            )

        # Check same date
        if first_prod.date != prod.date:
            raise IncompatibleProducts(
                f"Components of a mosaic should have the same date! {first_prod.date=} != {prod.date=}"
            )

    def check_contiguity(self, check_contiguity: GeometryCheck):
        """
        Check the contiguity of the mosaic

        Args:
            check_contiguity (GeometryCheck): Contiguity checking method

        Raises:
            IncompatibleProducts: Incompatible products if not contiguous according to the given method
        """
        if check_contiguity == GeometryCheck.EXTENT:
            union_extent = self.extent()
            if len(union_extent) > 1:
                raise IncompatibleProducts(
                    "The mosaic should have a contiguous extent!"
                )
        elif check_contiguity == GeometryCheck.FOOTPRINT:
            union_footprint = self.footprint()
            if len(union_footprint) > 1:
                raise IncompatibleProducts(
                    "The mosaic should have a contiguous footprint!"
                )
        else:
            LOGGER.warning("The contiguity of your mosaic won't be checked!")
            pass

    def read_mtd(self):
        """Read the pair's metadata, but not implemented for now."""
        # TODO: how ? Just return the fields that are shared between mosaic's components ? Or create a XML from scratch ?
        raise NotImplementedError

    @cache
    def footprint(self) -> gpd.GeoDataFrame:
        """
        Get the footprint of the mosaic.

        Returns:
            gpd.GeoDataFrame: Footprint of the mosaic
        """
        ref_prod = self.get_first_prod()
        footprint: gpd.GeoDataFrame = self.get_first_prod().footprint()

        if self.nof_prods > 1:
            for prod in self.get_prods()[1:]:
                footprint = footprint.overlay(
                    prod.footprint().to_crs(ref_prod.crs()), how="union"
                )

            # Dissolve and explode the footprint
            footprint = footprint.dissolve().explode(index_parts=True)

        return footprint

    @cache
    def extent(self) -> gpd.GeoDataFrame:
        """
        Get the extent of the mosaic.

        Returns:
            gpd.GeoDataFrame: Extent of the mosaic

        """
        ref_prod = self.get_first_prod()
        extent: gpd.GeoDataFrame = self.get_first_prod().extent()

        if self.nof_prods > 1:
            for prod in self.get_prods()[1:]:
                extent = extent.overlay(
                    prod.extent().to_crs(ref_prod.crs()), how="union"
                )

            # Dissolve and explode the extent
            extent = extent.dissolve().explode(index_parts=True)

        return extent

    def _get_band_suffix(self):
        """Get the band suffix"""
        return f"{self.mosaic_method.name.lower()}"

    def load(
        self,
        bands: Union[list, BandNames, str],
        pixel_size: float = None,
        **kwargs,
    ) -> xr.Dataset:
        """
        Load the bands and compute the wanted spectral indices.

        Args:
            bands (Union[list, BandNames, str]): Wanted bands
            pixel_size (float): Pixel size of the returned Dataset. If not specified, use the mosaic's pixel size.
            **kwargs: Other arguments used to load bands

        Returns:
            xr.Dataset: Wanted bands as xr.Datasets
        """
        # Get merge function and extension
        merge_fct = getattr(rasters, self.mosaic_method.value)

        # Convert just in case
        bands = to_band(bands)

        # Get the bands to be loaded
        bands_to_load, bands_path = self.get_bands_to_load(
            bands, self._get_band_suffix()
        )

        # Check validity of the bands
        for prod in self.get_prods():
            for band in bands_to_load:
                assert prod.has_band(
                    band
                ), f"{prod.condensed_name} has not a {to_str(band)[0]} band."

        # Load and reorganize bands
        prod_band_paths = defaultdict(list)
        for prod in self.get_prods():
            prod: Product
            LOGGER.debug(
                f"*** Loading {to_str(bands_to_load)} for {prod.condensed_name} ***"
            )

            # Load bands
            prod.load(bands_to_load, pixel_size, **kwargs).keys()

            # Store paths
            for band in bands_to_load:
                if is_spectral_band(band):
                    band_path = prod.get_band_paths([band], pixel_size, **kwargs)[band]
                else:
                    # Use glob fct as _get_band_folder is a tmpDirectory
                    band_path = glob(
                        os.path.join(prod._get_band_folder(), f"*{to_str(band)[0]}*")
                    )[0]

                prod_band_paths[band].append(str(band_path))

        # Merge
        merged_dict = {}
        for band in bands_path:
            output_path = bands_path[band]
            if not output_path.is_file():
                LOGGER.debug(f"Merging bands {to_str(band)[0]}")
                if self.mosaic_method == MosaicMethod.VRT:
                    prod_path = []
                    for path in prod_band_paths[band]:
                        out_path, exists = self._get_out_path(os.path.basename(path))
                        if not exists:
                            shutil.move(path, out_path)
                        prod_path.append(out_path)
                else:
                    prod_path = prod_band_paths[band]

                # Don't pass kwargs here because of unwanted errors
                merge_fct(prod_path, output_path)

            # Load in memory and update attribute
            merged_dict[band] = self._update_attrs(
                rasters.read(output_path), bands, **kwargs
            )

        # Collocate VRTs
        LOGGER.debug("Collocating bands")
        merged_dict = self._collocate_bands(merged_dict)

        # Create a dataset (only after collocation)
        coords = None
        if merged_dict:
            coords = merged_dict[bands[0]].coords

        # Make sure the dataset has the bands in the right order -> re-order the input dict
        mos_ds = xr.Dataset({key: merged_dict[key] for key in bands}, coords=coords)

        # Update attributes
        mos_ds = self._update_xds_attrs(mos_ds, bands)

        return mos_ds

    def stack(
        self,
        bands: list,
        pixel_size: float = None,
        stack_path: Union[str, AnyPathType] = None,
        save_as_int: bool = False,
        **kwargs,
    ) -> xr.DataArray:
        """
        Stack bands and index of a mosaic.

        Args:
            bands (list): Bands and index combination
            pixel_size (float): Stack pixel size. . If not specified, use the product pixel size.
            stack_path (Union[str, AnyPathType]): Stack path
            save_as_int (bool): Convert stack to uint16 to save disk space (and therefore multiply the values by 10.000)
            **kwargs: Other arguments passed to :code:`load` or :code:`rioxarray.to_raster()` (such as :code:`compress`)

        Returns:
            xr.DataArray: Stack as a DataArray
        """
        if stack_path:
            stack_path = AnyPath(stack_path)
            if stack_path.is_file():
                return utils.read(stack_path, pixel_size=pixel_size)
            else:
                os.makedirs(str(stack_path.parent), exist_ok=True)

        # Create the analysis stack
        band_ds = self.load(bands, pixel_size=pixel_size, **kwargs)

        # Stack bands
        if save_as_int:
            nodata = kwargs.get("nodata", UINT16_NODATA)
        else:
            nodata = kwargs.get("nodata", self.nodata)
        stack, dtype = utils.stack_dict(bands, band_ds, save_as_int, nodata, **kwargs)

        # Update stack's attributes
        stack = self._update_attrs(stack, bands, **kwargs)

        # Write on disk
        if stack_path:
            LOGGER.debug("Saving stack")
            utils.write(stack, stack_path, dtype=dtype, **kwargs)

        return stack

    def _collocate_bands(self, bands: dict, reference: xr.DataArray = None) -> dict:
        """
        Collocate all bands from a dict

        Args:
            bands (dict): Dict of bands to collocate if needed
            reference (xr.DataArray): Reference array

        Returns:
            dict: Collocated bands
        """
        return self.get_first_prod()._collocate_bands(bands, reference)

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
        xarr.attrs["acquisition_date"] = self.date

        return xarr
