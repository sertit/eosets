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
"""Class implementing a two-products pair"""

import logging
import os
from enum import unique

import geopandas as gpd
import xarray as xr
from eoreader import cache
from eoreader.bands import to_band, to_str
from eoreader.utils import UINT16_NODATA
from rasterio.enums import Resampling
from sertit import AnyPath, rasters
from sertit.misc import ListEnum
from sertit.types import AnyPathStrType

from eosets import EOSETS_NAME
from eosets.exceptions import IncompatibleProducts
from eosets.mosaic import AnyMosaicType, Mosaic
from eosets.set import GeometryCheck, GeometryCheckType, Set
from eosets.utils import BandsType, read, stack, write

LOGGER = logging.getLogger(EOSETS_NAME)


@unique
class DiffMethod(ListEnum):
    """Available difference methods."""

    REFERENCE_SECONDARY = "reference-secondary"
    SECONDARY_REFERENCE = "secondary-reference"


class Pair(Set):
    """Class of two-products pair"""

    def __init__(
        self,
        reference_paths: AnyMosaicType,
        secondary_paths: AnyMosaicType = None,
        id: str = None,
        output_path: AnyPathStrType = None,
        remove_tmp: bool = True,
        overlap_check: GeometryCheckType = GeometryCheck.EXTENT,
        contiguity_check: GeometryCheckType = GeometryCheck.EXTENT,
        **kwargs,
    ):
        # Manage reference mosaic
        self.reference_mosaic = None
        """ Reference mosaic (unique date and contiguous). The one on which the secondary will be aligned. """

        self.reference_id = None
        """ ID of the reference product """

        # Manage secondary mosaic
        self.secondary_mosaic = None
        """ Secondary mosaic (unique date and contiguous). The one which will be aligned on the reference. """

        self.secondary_id = None
        """ ID of the secondary product """

        # Information regarding the pair composition
        self.has_secondary = None
        """ Does the pair have a secondary mosaic? (Pair with only one reference is allowed) """

        # Convert the checks to the corresponding enums
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
        if secondary_paths is None:
            secondary_paths = []
        self._manage_mosaics(
            reference_paths, secondary_paths, contiguity_check, overlap_check
        )

        # Fill attributes
        self.full_name = f"{self.reference_id}"
        if self.has_secondary:
            self.full_name += f"_{self.secondary_id}"

        self.condensed_name = self.full_name
        # TODO (how to name pair ???)

        # Post init at the set level
        self.post_init(**kwargs)

    def clean_tmp(self):
        """
        Clean the temporary directory of the current pair
        """
        self.reference_mosaic.clean_tmp()

        if self.has_secondary:
            self.secondary_mosaic.clean_tmp()

    def clear(self):
        """
        Clear this pair's cache
        """
        # Delete all cached properties and functions
        self.reference_mosaic.clear()

        if self.has_secondary:
            self.secondary_mosaic.clear()

    def _manage_output(self):
        """
        Manage the output specifically for this child class
        """
        self.reference_mosaic.output = self.output
        try:
            if self.has_secondary:
                self.secondary_mosaic.output = self.output
        except FileNotFoundError:
            # Never mind for non-existing files: they have already been copied :)
            pass

    def get_prods(self) -> list:
        """
        Get all the products as a list.

        Returns:
            list: Products list
        """
        prods = self.reference_mosaic.get_prods()

        if self.has_secondary:
            prods += self.secondary_mosaic.get_prods()
        return prods

    def _manage_mosaics(
        self,
        reference_paths: AnyMosaicType,
        secondary_paths: AnyMosaicType = None,
        contiguity_check: GeometryCheck = GeometryCheck.EXTENT,
        overlap_check: GeometryCheck = GeometryCheck.EXTENT,
    ) -> None:
        """
        Check if the reference and secondary mosaics are overlapping.

        TODO: check if same constellation ?

        If not, throws a IncompatibleProducts error.

        Args:
            reference_paths (AnyMosaicType): Paths corresponding to the reference mosaic
            secondary_paths (AnyMosaicType): Paths corresponding to the secondary mosaic
            contiguity_check (GeometryCheck): Check regarding the contiguity of the products of the mosaics
            overlap_check (GeometryCheck): Check regarding the overlapping of the two mosaics

        Raises:
            IncompatibleProducts: Incompatible products if not contiguous or not the same date
        """
        # Manage reference product
        if isinstance(reference_paths, Mosaic):
            self.reference_mosaic = reference_paths
        else:
            self.reference_mosaic: Mosaic = Mosaic(
                reference_paths,
                output_path=self._get_tmp_folder(writable=True),
                remove_tmp=self._remove_tmp,
                contiguity_check=contiguity_check,
            )
        self.reference_id: str = self.reference_mosaic.id

        # Information regarding the pair composition
        self.has_secondary: bool = len(secondary_paths) > 0

        if self.has_secondary:
            if isinstance(reference_paths, Mosaic):
                self.reference_mosaic = reference_paths
            else:
                self.secondary_mosaic: Mosaic = Mosaic(
                    secondary_paths,
                    output_path=self._get_tmp_folder(writable=True),
                    remove_tmp=self._remove_tmp,
                    contiguity_check=contiguity_check,
                )
            self.secondary_id: str = self.secondary_mosaic.id

            # Check Geometry
            if overlap_check != GeometryCheck.NONE:
                ref_geom: gpd.GeoDataFrame = getattr(
                    self.reference_mosaic, str(overlap_check.value)
                )()
                sec_geom: gpd.GeoDataFrame = getattr(
                    self.secondary_mosaic, str(overlap_check.value)
                )()
                if not ref_geom.intersects(
                    sec_geom.to_crs(self.reference_mosaic.crs)
                ).all():
                    raise IncompatibleProducts(
                        "Reference and secondary mosaics should overlap!"
                    )

        self.nof_prods = len(self.get_prods())

    def read_mtd(self):
        """Read the pair's metadata, but not implemented for now."""
        # TODO: how ? Just return the fields that are shared between pair's components ? Or create a XML from scratch ?
        raise NotImplementedError

    @cache
    def footprint(self) -> gpd.GeoDataFrame:
        """
        Get the footprint of the pair, i.e. the intersection between reference and secondary footprints.

        Returns:
            gpd.GeoDataFrame: Footprint of the pair
        """
        ref_geom: gpd.GeoDataFrame = self.reference_mosaic.footprint()

        if self.has_secondary:
            second_geom: gpd.GeoDataFrame = self.secondary_mosaic.footprint().to_crs(
                self.reference_mosaic.crs
            )
            footprint = ref_geom.overlay(second_geom, "intersection")
        else:
            footprint = ref_geom

        return footprint

    @cache
    def extent(self) -> gpd.GeoDataFrame:
        """
        Get the extent of the pair, i.e. the intersection between reference and secondary extents.

        Returns:
            gpd.GeoDataFrame: Extent of the pair

        """
        ref_geom: gpd.GeoDataFrame = self.reference_mosaic.extent()

        if self.has_secondary:
            second_geom: gpd.GeoDataFrame = self.secondary_mosaic.extent().to_crs(
                self.reference_mosaic.crs
            )
            extent = ref_geom.overlay(second_geom, "intersection")
        else:
            extent = ref_geom
        return extent

    def load(
        self,
        reference_bands: BandsType = None,
        secondary_bands: BandsType = None,
        diff_bands: BandsType = None,
        pixel_size: float = None,
        diff_method: DiffMethod = DiffMethod.REFERENCE_SECONDARY,
        resampling: Resampling = Resampling.bilinear,
        **kwargs,
    ) -> (xr.Dataset, xr.Dataset, xr.Dataset):
        """
        Load the bands and compute the wanted spectral indices for reference, secondary and diff.

        Args:
            reference_bands (BandsType): Wanted reference bands
            secondary_bands (BandsType): Wanted secondary bands
            diff_bands (BandsType): Wanted diff bands
            pixel_size (float): Pixel size of the returned Datasets. If not specified, use the pair's pixel size.
            diff_method (DiffMethod): Difference method for the computation of diff_bands
            resampling (Resampling): Resampling method
            kwargs: Other arguments used to load bands

        Returns:
            (xr.Dataset, xr.Dataset, xr.Dataset): Reference, secondary and diff wanted bands as xr.Datasets
        """
        assert any(
            [
                reference_bands is not None,
                secondary_bands is not None,
                diff_bands is not None,
            ]
        )

        # Convert just in case
        if reference_bands is None:
            reference_bands = []
        if secondary_bands is None:
            secondary_bands = []
        if diff_bands is None:
            diff_bands = []

        reference_bands = to_band(reference_bands)
        secondary_bands = to_band(secondary_bands)
        diff_bands = to_band(diff_bands)

        # Check existing diff paths
        diff_bands_to_load, diff_bands_path = self.get_bands_to_load(diff_bands)

        # Overload reference and secondary bands with diff bands
        ref_bands_to_load = reference_bands.copy()
        sec_bands_to_load = secondary_bands.copy()
        for band in diff_bands_to_load:
            if band not in reference_bands:
                ref_bands_to_load.append(band)
            if band not in secondary_bands:
                sec_bands_to_load.append(band)

        # -- Load bands
        window = kwargs.pop("window", self.footprint())

        # Load reference bands
        ref_ds: xr.Dataset = self.reference_mosaic.load(
            ref_bands_to_load, pixel_size=pixel_size, window=window, **kwargs
        )

        # Load secondary bands
        if self.has_secondary:
            sec_ds: xr.Dataset = self.secondary_mosaic.load(
                sec_bands_to_load, pixel_size=pixel_size, window=window, **kwargs
            )

            # Load diff bands
            diff_dict = {}
            for band in diff_bands:
                diff_path, exists = self._get_out_path(
                    f"{self.condensed_name}_d{to_str(band)[0]}.tif"
                )
                if exists:
                    diff_arr = read(
                        diff_path,
                        pixel_size=pixel_size,
                        resampling=resampling,
                        **kwargs,
                    )
                else:
                    f"*** Loading d{to_str(band)} for {self.condensed_name} ***"
                    ref_arr = ref_ds[band]
                    sec_arr = sec_ds[band]

                    # To be sure, always collocate arrays, even if the size is the same
                    # Indeed, a small difference in the coordinates will lead to empy arrays
                    # So the bands MUST BE exactly aligned
                    sec_arr = sec_arr.rio.reproject_match(
                        ref_arr, resampling=resampling, **kwargs
                    )

                    # Nans are conserved with +/-
                    # So only the overlapping extent WITH nodata of both reference and secondary is loaded
                    if diff_method == DiffMethod.REFERENCE_SECONDARY:
                        diff_arr = ref_arr - sec_arr
                    else:
                        diff_arr = sec_arr - ref_arr

                    # Save diff band
                    diff_name = f"d{to_str(band)[0]}"
                    diff_arr = ref_arr.copy(data=diff_arr).rename(diff_name)
                    diff_arr.attrs["long_name"] = diff_name

                    # Write on disk
                    write(diff_arr, diff_path)

                diff_dict[band] = diff_arr

            # Collocate diff bands
            diff_dict = self._collocate_bands(diff_dict)

            # Drop not wanted bands from secondary dataset
            sec_ds = sec_ds.drop_vars(
                [band for band in sec_ds if band not in secondary_bands]
            )

            # Create diff dataset
            coords = None
            if diff_dict:
                coords = diff_dict[diff_bands[0]].coords

            # Make sure the dataset has the bands in the right order -> re-order the input dict
            diff_ds = xr.Dataset(
                {key: diff_dict[key] for key in diff_bands}, coords=coords
            )

            # Update attributes
            sec_ds = self._update_xds_attrs(sec_ds, secondary_bands)
            diff_ds = self._update_xds_attrs(diff_ds, diff_bands)
        else:
            sec_ds = xr.Dataset()
            diff_ds = xr.Dataset()

        # Update reference dataset
        # Drop not wanted bands from reference dataset
        ref_ds = ref_ds.drop_vars(
            [band for band in ref_ds if band not in reference_bands]
        )

        ref_ds = self._update_xds_attrs(ref_ds, reference_bands)

        return ref_ds, sec_ds, diff_ds

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

    def stack(
        self,
        reference_bands: BandsType = None,
        secondary_bands: BandsType = None,
        diff_bands: BandsType = None,
        pixel_size: float = None,
        diff_method: DiffMethod = DiffMethod.REFERENCE_SECONDARY,
        stack_path: AnyPathStrType = None,
        save_as_int: bool = False,
        **kwargs,
    ) -> xr.DataArray:
        """
        Stack bands and index of a pair.

        Args:
            reference_bands (BandsType): Bands and index combination for the reference mosaic
            secondary_bands (BandsType): Bands and index combination for the secondary mosaic
            diff_bands (BandsType): Bands and index combination for the difference between reference and secondary mosaic
            pixel_size (float): Stack pixel size. If not specified, use the pair's pixel size.
            stack_path (AnyPathStrType): Stack path
            save_as_int (bool): Convert stack to uint16 to save disk space (and therefore multiply the values by 10.000)
            **kwargs: Other arguments passed to :code:`load` or :code:`rioxarray.to_raster()` (such as :code:`compress`)

        Returns:
            xr.DataArray: Stack as a DataArray
        """
        assert any(
            [
                reference_bands is not None,
                secondary_bands is not None,
                diff_bands is not None,
            ]
        )

        # Convert just in case
        if reference_bands is None:
            reference_bands = []
        if secondary_bands is None:
            secondary_bands = []
        if diff_bands is None:
            diff_bands = []

        reference_bands = to_band(reference_bands)
        secondary_bands = to_band(secondary_bands)
        diff_bands = to_band(diff_bands)

        if stack_path:
            stack_path = AnyPath(stack_path)
            if stack_path.is_file():
                return read(stack_path, pixel_size=pixel_size)
            else:
                os.makedirs(str(stack_path.parent), exist_ok=True)

        # Load all bands
        ref_ds, sec_ds, diff_ds = self.load(
            reference_bands,
            secondary_bands,
            diff_bands,
            pixel_size=pixel_size,
            **kwargs,
        )

        # Rename bands
        ref_band_mapping = {
            band: f"Reference_{to_str(band)[0]}" for band in reference_bands
        }
        sec_band_mapping = {
            band: f"Secondary_{to_str(band)[0]}" for band in secondary_bands
        }
        diff_band_mapping = {band: f"d{to_str(band)[0]}" for band in diff_bands}
        all_bands = (
            list(ref_band_mapping.values())
            + list(sec_band_mapping.values())
            + list(diff_band_mapping.values())
        )
        ref_ds = ref_ds.rename_vars(ref_band_mapping)

        if self.has_secondary:
            sec_ds = rasters.collocate(ref_ds, sec_ds.rename_vars(sec_band_mapping))
            diff_ds = rasters.collocate(ref_ds, diff_ds.rename_vars(diff_band_mapping))

        # Merge datasets
        band_ds = xr.merge([ref_ds, sec_ds, diff_ds])

        # Stack bands
        if save_as_int:
            nodata = kwargs.get("nodata", UINT16_NODATA)
        else:
            nodata = kwargs.get("nodata", self.nodata)
        stk, dtype = stack(band_ds, save_as_int, nodata, **kwargs)

        # Update stack's attributes
        stk = self._update_attrs(stk, all_bands, **kwargs)

        # Write on disk
        if stack_path:
            LOGGER.debug("Saving stack")
            write(stk, stack_path, dtype=dtype, **kwargs)

        return stk

    def _collocate_bands(
        self, bands: dict, reference: xr.DataArray = None
    ) -> xr.Dataset:
        """
        Collocate all bands from a dict

        Args:
            bands (dict): Dict of bands to collocate if needed
            reference (xr.DataArray): Reference array

        Returns:
            xr.Dataset: Collocated bands
        """
        return self.reference_mosaic._collocate_bands(bands, reference)
