""" Class implementing the pairs """
from enum import unique
from pathlib import Path
from typing import Union

import geopandas as gpd
import xarray as xr
from cloudpathlib import CloudPath
from eoreader import cache
from eoreader.bands import BandNames, to_band, to_str
from eoreader.reader import Reader
from rasterio.enums import Resampling
from sertit.misc import ListEnum

from eosets.exceptions import IncompatibleProducts
from eosets.mosaic import Mosaic
from eosets.set import GeometryCheck, Set

READER = Reader


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
            output_path=self.output,
            remove_tmp=self._remove_tmp,
            contiguity_check=contiguity_check,
        )
        self.pivot_id: str = self.pivot_mosaic.id

        # Information regarding the pair composition
        self.has_child: bool = len(child_paths) > 0

        if self.has_child:
            self.child_mosaic: Mosaic = Mosaic(
                child_paths,
                output_path=self.output,
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
    ) -> (dict, dict, dict):
        """"""
        assert any(
            [pivot_bands is not None, child_bands is not None, diff_bands is not None]
        )

        pivot_dict = {}
        child_dict = {}
        diff_dict = {}

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

        # Overload pivot and child bands with diff bands
        pivot_bands_to_load = pivot_bands.copy()
        child_bands_to_load = child_bands.copy()
        for band in diff_bands:
            if band not in pivot_bands:
                pivot_bands_to_load.append(band)
            if band not in child_bands:
                child_bands_to_load.append(band)

        # Load bands
        window = kwargs.pop("window", self.footprint())

        # Load pivot bands
        pivot_dict = self.pivot_mosaic.load(
            pivot_bands_to_load, pixel_size=pixel_size, window=window, **kwargs
        )
        # TODO: crop to footprint ?

        # Load pivot bands
        child_dict = self.child_mosaic.load(
            child_bands_to_load, pixel_size=pixel_size, window=window, **kwargs
        )
        # TODO: crop to footprint ?

        # Load pivot bands
        diff_dict = {}
        for band in diff_bands:
            pivot_arr = pivot_dict[band]
            child_arr = child_dict[band]
            if pivot_arr.shape != child_arr.shape:
                child_arr = child_arr.rio.reproject_match(
                    pivot_arr, resampling=resampling, **kwargs
                )

            # Nans are conserved with +/-
            # So only the overlapping footprint is loaded
            """
            a = xr.DataArray(np.ones([2, 2]))
            a[0, 0] = np.nan
            b = xr.DataArray(2*np.ones([2, 2]))
            b[1, 0] = np.nan
            a + b
                <xarray.DataArray (dim_0: 2, dim_1: 2)>
                array([[nan,  3.],
                       [nan,  3.]])
            a - b
                <xarray.DataArray (dim_0: 2, dim_1: 2)>
                array([[nan, -1.],
                       [nan, -1.]])
            """
            if diff_method == DiffMethod.PIVOT_CHILD:
                diff_arr = pivot_arr - child_arr
            else:
                diff_arr = child_arr - pivot_arr

            # Save diff band
            diff_name = f"d{to_str(band)[0]}"
            diff_arr = pivot_arr.copy(data=diff_arr).rename(diff_name)
            diff_arr.attrs["long_name"] = diff_name
            diff_dict[band] = diff_arr

        return (
            {band: pivot_dict[band] for band in pivot_bands},
            {band: child_dict[band] for band in child_bands},
            diff_dict,
        )

    def _update_attrs(self, xarr: xr.DataArray, bands: list, **kwargs) -> xr.DataArray:
        """
        Update attributes of the given array
        Args:
            xarr (xr.DataArray): Array whose attributes need an update
            bands (list): Bands
        Returns:
            xr.DataArray: Updated array
        """
        # Clean attributes, we don't want to pollute our attributes by default ones (not deterministic)
        # Are we sure of that ?
        xarr.attrs = {}

        if not isinstance(bands, list):
            bands = [bands]
        long_name = to_str(bands)
        xr_name = "_".join(long_name)
        attr_name = " ".join(long_name)

        xarr = xarr.rename(xr_name)
        xarr.attrs["long_name"] = attr_name
        xarr.attrs["condensed_name"] = self.condensed_name

        # TODO: complete that

        return xarr
