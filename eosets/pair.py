""" Class implementing the pairs """
import logging
import os
from enum import unique
from pathlib import Path
from typing import Union

import geopandas as gpd
import xarray as xr
from cloudpathlib import AnyPath, CloudPath
from eoreader import cache
from eoreader.bands import BandNames, to_band, to_str
from eoreader.utils import UINT16_NODATA
from rasterio.enums import Resampling
from sertit import rasters
from sertit.misc import ListEnum

from eosets import utils
from eosets.exceptions import IncompatibleProducts
from eosets.mosaic import Mosaic
from eosets.set import GeometryCheck, Set
from eosets.utils import EOSETS_NAME, AnyPathType

LOGGER = logging.getLogger(EOSETS_NAME)


@unique
class DiffMethod(ListEnum):
    """Available difference methods."""

    PIVOT_CHILD = "pivot-child"
    CHILD_PIVOT = "child-pivot"


class Pair(Set):
    """Class of two-products pair"""

    def __init__(
        self,
        pivot_paths: Union[list, str, Path, CloudPath],
        child_paths: Union[list, str, Path, CloudPath] = None,
        id: str = None,
        output_path: Union[str, Path, CloudPath] = None,
        remove_tmp: bool = True,
        overlap_check: Union[GeometryCheck, str] = GeometryCheck.EXTENT,
        contiguity_check: Union[GeometryCheck, str] = GeometryCheck.EXTENT,
        **kwargs,
    ):
        # Manage reference product
        self.pivot_mosaic = None
        """ Pivot mosaic (unique date and contiguous). The one on which the child will be aligned. """

        self.pivot_id = None
        """ ID of the pivot product """

        self.child_mosaic = None
        """ Child mosaic (unique date and contiguous). The one which will be aligned on the pivot. """

        self.child_id = None
        """ ID of the child product """

        # Information regarding the pair composition
        self.has_child = None
        """ Does the pair have a child? (Pair with only one reference is allowed) """

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
        self._manage_mosaics(pivot_paths, child_paths, contiguity_check, overlap_check)

        # Fill attributes
        self.full_name = f"{self.pivot_id}_{self.child_id}"
        self.condensed_name = self.full_name
        # TODO (how to name pair ???)

        self.nodata = self.get_attr("nodata", **kwargs)
        self.pixel_size = self.get_attr("pixel_size", **kwargs)
        self.crs = self.get_attr("crs", **kwargs)
        self.same_constellation = self.is_homogeneous("constellation")
        self.constellations = list(set(prod.constellation for prod in self.get_prods()))

    def clean_tmp(self):
        """
        Clean the temporary directory of the current product
        """
        self.pivot_mosaic.clean_tmp()
        self.child_mosaic.clean_tmp()

    def clear(self):
        """
        Clear this product's cache
        """
        # Delete all cached properties and functions
        self.pivot_mosaic.clear()
        self.child_mosaic.clear()

    def _manage_output(self):
        """
        Manage the output specifically for this child class
        """
        self.pivot_mosaic._manage_output()
        self.child_mosaic._manage_output()

    def get_prods(self) -> list:
        """
        Get all the products as a list.

        Returns:
            list: Products list
        """
        return self.pivot_mosaic.get_prods() + self.child_mosaic.get_prods()

    def _manage_mosaics(
        self,
        pivot_paths: Union[list, str, Path, CloudPath],
        child_paths: Union[list, str, Path, CloudPath] = None,
        contiguity_check: GeometryCheck = GeometryCheck.EXTENT,
        overlap_check: GeometryCheck = GeometryCheck.EXTENT,
    ) -> None:
        """
        Check if the pivot and child mosaics are overlapping

        TODO: same constellation ? same CRS ?...

        If not, throws a IncompatibleProducts error.

        Raises:
            IncompatibleProducts: Incompatible products if not contiguous or not the same date
        """
        # Manage reference product
        self.pivot_mosaic: Mosaic = Mosaic(
            pivot_paths,
            output_path=self._get_tmp_folder(writable=True),
            remove_tmp=self._remove_tmp,
            contiguity_check=contiguity_check,
        )
        self.pivot_id: str = self.pivot_mosaic.id

        # Information regarding the pair composition
        self.has_child: bool = len(child_paths) > 0

        if self.has_child:
            self.child_mosaic: Mosaic = Mosaic(
                child_paths,
                output_path=self._get_tmp_folder(writable=True),
                remove_tmp=self._remove_tmp,
                contiguity_check=contiguity_check,
            )
            self.child_id: str = self.child_mosaic.id

            # Make the checks
            # CRS
            if self.pivot_mosaic.crs != self.pivot_mosaic.crs:
                raise IncompatibleProducts(
                    f"Pivot and child mosaics should have the same CRS! {self.pivot_mosaic.crs=} != {self.pivot_mosaic.crs=}"
                )

            # Geometry
            if overlap_check != GeometryCheck.NONE:
                pivot_geom: gpd.GeoDataFrame = getattr(
                    self.pivot_mosaic, str(overlap_check.value)
                )()
                child_geom: gpd.GeoDataFrame = getattr(
                    self.child_mosaic, str(overlap_check.value)
                )()
                if not pivot_geom.intersects(
                    child_geom.to_crs(self.pivot_mosaic.crs)
                ).all():
                    raise IncompatibleProducts(
                        "Pivot and child mosaics should overlap!"
                    )

        self.nof_prods = len(self.get_prods())

    def read_mtd(self):
        """"""
        # TODO: how ? Just return the fields that are shared between mosaic's components ? Or create a XML from scratch ?
        raise NotImplementedError

    @cache
    def footprint(self) -> gpd.GeoDataFrame:
        """
        Get the footprint of the pair, i.e. the intersection between pivot and child footprints.

        Returns:
            gpd.GeoDataFrame: Footprint of the mosaic
        """
        pivot_geom: gpd.GeoDataFrame = self.pivot_mosaic.footprint()
        child_geom: gpd.GeoDataFrame = self.child_mosaic.footprint().to_crs(
            self.pivot_mosaic.crs
        )
        footprint = pivot_geom.overlay(child_geom, "intersection")
        return footprint

    @cache
    def extent(self) -> gpd.GeoDataFrame:
        """
        Get the extent of the pair, i.e. the intersection between pivot and child extents.

        Returns:
            gpd.GeoDataFrame: Extent of the mosaic

        """
        pivot_geom: gpd.GeoDataFrame = self.pivot_mosaic.extent()
        child_geom: gpd.GeoDataFrame = self.child_mosaic.extent().to_crs(
            self.pivot_mosaic.crs
        )
        extent = pivot_geom.overlay(child_geom, "intersection")
        return extent

    def load(
        self,
        pivot_bands: Union[list, BandNames, str] = None,
        child_bands: Union[list, BandNames, str] = None,
        diff_bands: Union[list, BandNames, str] = None,
        pixel_size: float = None,
        diff_method: DiffMethod = DiffMethod.PIVOT_CHILD,
        resampling: Resampling = Resampling.bilinear,
        **kwargs,
    ) -> (xr.Dataset, xr.Dataset, xr.Dataset):
        """

        Args:
            pivot_bands:
            child_bands:
            diff_bands:
            pixel_size:
            diff_method:
            resampling:
            **kwargs:

        Returns:
            (xr.Dataset, xr.Dataset, xr.Dataset):
        """
        assert any(
            [pivot_bands is not None, child_bands is not None, diff_bands is not None]
        )

        # Convert just in case
        if pivot_bands is None:
            pivot_bands = []
        if child_bands is None:
            child_bands = []
        if diff_bands is None:
            diff_bands = []

        pivot_bands = to_band(pivot_bands)
        child_bands = to_band(child_bands)
        diff_bands = to_band(diff_bands)

        # Check existing diff paths
        diff_bands_to_load, diff_bands_path = self.get_bands_to_load(diff_bands)

        # Overload pivot and child bands with diff bands
        pivot_bands_to_load = pivot_bands.copy()
        child_bands_to_load = child_bands.copy()
        for band in diff_bands_to_load:
            if band not in pivot_bands:
                pivot_bands_to_load.append(band)
            if band not in child_bands:
                child_bands_to_load.append(band)

        # Load bands
        window = kwargs.pop("window", self.footprint())

        # Load pivot bands
        pivot_ds: xr.Dataset = self.pivot_mosaic.load(
            pivot_bands_to_load, pixel_size=pixel_size, window=window, **kwargs
        )

        # Load pivot bands
        child_ds: xr.Dataset = self.child_mosaic.load(
            child_bands_to_load, pixel_size=pixel_size, window=window, **kwargs
        )

        # Load pivot bands
        diff_dict = {}
        for band in diff_bands:
            diff_path, exists = self._get_out_path(
                f"{self.condensed_name}_{to_str(band)[0]}.tif"
            )
            if exists:
                diff_arr = utils.read(
                    path=diff_path,
                    pixel_size=pixel_size,
                    resampling=resampling,
                    **kwargs,
                )
            else:
                f"*** Loading d{to_str(band)} for {self.condensed_name} ***"
                pivot_arr = pivot_ds[band]
                child_arr = child_ds[band]

                # To be sure, always collocate arrays, even if the size is the same
                # Indeed, a small difference in the coordinates will lead to empy arrays
                # So the bands MUST BE exactly aligned
                child_arr = child_arr.rio.reproject_match(
                    pivot_arr, resampling=resampling, **kwargs
                )

                # Nans are conserved with +/-
                # So only the overlapping extent WITH nodata of both pivot and child is loaded
                if diff_method == DiffMethod.PIVOT_CHILD:
                    diff_arr = pivot_arr - child_arr
                else:
                    diff_arr = child_arr - pivot_arr

                # Save diff band
                diff_name = f"d{to_str(band)[0]}"
                diff_arr = pivot_arr.copy(data=diff_arr).rename(diff_name)
                diff_arr.attrs["long_name"] = diff_name

                # Write on disk
                utils.write(diff_arr, diff_path)

            diff_dict[band] = diff_arr

        # Collocate diff bands
        diff_dict = self._collocate_bands(diff_dict)

        # Drop not wanted bands from pivot and child datasets
        pivot_ds = pivot_ds.drop_vars(
            [band for band in pivot_ds.keys() if band not in pivot_bands]
        )
        child_ds = child_ds.drop_vars(
            [band for band in child_ds.keys() if band not in child_bands]
        )

        # Create diff dataset
        coords = None
        if diff_dict:
            coords = diff_dict[diff_bands[0]].coords

        # Make sure the dataset has the bands in the right order -> re-order the input dict
        diff_ds = xr.Dataset({key: diff_dict[key] for key in diff_bands}, coords=coords)

        return pivot_ds, child_ds, diff_ds

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
        pivot_bands: Union[list, BandNames, str] = None,
        child_bands: Union[list, BandNames, str] = None,
        diff_bands: Union[list, BandNames, str] = None,
        pixel_size: float = None,
        diff_method: DiffMethod = DiffMethod.PIVOT_CHILD,
        stack_path: Union[str, AnyPathType] = None,
        save_as_int: bool = False,
        **kwargs,
    ) -> xr.DataArray:
        """
        Stack bands and index of a pair.

        Args:
            pivot_bands (list): Bands and index combination for the pivot mosaic
            child_bands (list): Bands and index combination for the child mosaic
            diff_bands (list): Bands and index combination for the difference between pivot and child mosaic
            pixel_size (float): Stack pixel size. . If not specified, use the product pixel size.
            stack_path (Union[str, AnyPathType]): Stack path
            save_as_int (bool): Convert stack to uint16 to save disk space (and therefore multiply the values by 10.000)
            **kwargs: Other arguments passed to :code:`load` or :code:`rioxarray.to_raster()` (such as :code:`compress`)

        Returns:
            xr.DataArray: Stack as a DataArray
        """
        assert any(
            [pivot_bands is not None, child_bands is not None, diff_bands is not None]
        )

        # Convert just in case
        if pivot_bands is None:
            pivot_bands = []
        if child_bands is None:
            child_bands = []
        if diff_bands is None:
            diff_bands = []

        pivot_bands = to_band(pivot_bands)
        child_bands = to_band(child_bands)
        diff_bands = to_band(diff_bands)

        if stack_path:
            stack_path = AnyPath(stack_path)
            if stack_path.is_file():
                return utils.read(stack_path, pixel_size=pixel_size)
            else:
                os.makedirs(str(stack_path.parent), exist_ok=True)

        # Load all bands
        pivot_ds, child_ds, diff_ds = self.load(
            pivot_bands, child_bands, diff_bands, pixel_size=pixel_size, **kwargs
        )

        # Rename bands
        pivot_band_mapping = {band: f"Pivot_{band}" for band in pivot_bands}
        child_band_mapping = {band: f"Child_{band}" for band in child_bands}
        diff_band_mapping = {band: f"d{band}" for band in diff_bands}
        all_bands = (
            list(pivot_band_mapping.values())
            + list(child_band_mapping.values())
            + list(diff_band_mapping.values())
        )
        pivot_ds = pivot_ds.rename_vars(pivot_band_mapping)
        child_ds = rasters.collocate(pivot_ds, child_ds.rename_vars(child_band_mapping))
        diff_ds = rasters.collocate(pivot_ds, diff_ds.rename_vars(diff_band_mapping))

        # Merge datasets
        band_ds = xr.merge([pivot_ds, child_ds, diff_ds])

        # Stack bands
        if save_as_int:
            nodata = kwargs.get("nodata", UINT16_NODATA)
        else:
            nodata = kwargs.get("nodata", self.nodata)
        stack, dtype = utils.stack_dict(
            all_bands, band_ds, save_as_int, nodata, **kwargs
        )

        # Update stack's attributes
        stack = self._update_attrs(stack, all_bands, **kwargs)

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
        return self.pivot_mosaic._collocate_bands(bands, reference)
