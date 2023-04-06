""" Class implementing the pairs """
from enum import unique
from pathlib import Path
from typing import Union

import geopandas as gpd
import xarray as xr
from cloudpathlib import CloudPath
from eoreader.bands import BandNames, to_band, to_str
from eoreader.reader import Reader
from sertit.misc import ListEnum

from eosets.exceptions import IncompatibleProducts
from eosets.mosaic import Mosaic
from eosets.set import GeometryCheck, Set

READER = Reader


@unique
class DiffMethod(ListEnum):
    """Available difference methods."""

    PRE_POST = "pre-post"
    POST_PRE = "post-pre"


class Pair(Set):
    """Class of multiple pairs"""

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
        super().__init__(output_path, id, remove_tmp, **kwargs)

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
        self._manage_mosaics(
            pivot_paths,
            child_paths,
            output_path,
            remove_tmp,
            contiguity_check,
            overlap_check,
        )

        # Full name
        self.full_name = f"{self.pivot_id}_{self.child_id}"
        """ Pair full name. """

        # Condensed name
        self.condensed_name = self.full_name
        # TODO (how to name pairs ???)

        # Nodata (by default use EOReader's)
        self.nodata = kwargs.get("nodata")
        """ Nodata of the pairs. """

        # Pixel size (by default use EOReader's)
        self.pixel_size = kwargs.get("pixel_size")
        """ Pixel size of the pair. """

        self.same_constellation = self.is_homogeneous("constellation")
        """ Are the pairs constituted of the same constellation? """

        self.constellations = list(set(prod.constellation for prod in self.get_prods()))
        """ List of unique constellations constituting the pair """

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
        # -- Delete all cached properties and functions
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
        output_path: Union[str, Path, CloudPath] = None,
        remove_tmp: bool = True,
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
            output_path=output_path,
            remove_tmp=remove_tmp,
            contiguity_check=contiguity_check,
        )
        self.pivot_id: str = self.pivot_mosaic.id

        # Information regarding the pair composition
        self.has_child: bool = len(child_paths) > 0

        if self.has_child:
            self.child_mosaic: Mosaic = Mosaic(
                child_paths,
                output_path=output_path,
                remove_tmp=remove_tmp,
                contiguity_check=contiguity_check,
            )
            self.child_id: str = self.child_mosaic.id

            # Make the checks
            # CRS
            if self.pivot_mosaic.crs() != self.pivot_mosaic.crs():
                raise IncompatibleProducts(
                    f"Pivot and child mosaics should have the same CRS! {self.pivot_mosaic.crs()=} != {self.pivot_mosaic.crs()=}"
                )

            # Geometry
            if overlap_check != GeometryCheck.NONE:
                pivot_geom: gpd.GeoDataFrame = getattr(
                    self.pivot_mosaic, str(overlap_check.value)
                )()
                child_geom: gpd.GeoDataFrame = getattr(
                    self.child_mosaic, str(overlap_check.value)
                )()
                if not pivot_geom.overlaps(child_geom).all():
                    raise IncompatibleProducts(
                        "Pivot and child mosaics should overlap!"
                    )

    def read_mtd(self):
        """"""
        # TODO: how ? Just return the fields that are shared between mosaic's components ? Or create a XML from scratch ?
        raise NotImplementedError

    def footprint(self) -> gpd.GeoDataFrame:
        """
        Get the footprint of the pair, i.e. the intersection between pivot and child footprints.

        Returns:
            gpd.GeoDataFrame: Footprint of the mosaic
        """
        pivot_geom: gpd.GeoDataFrame = self.pivot_mosaic.footprint()
        child_geom: gpd.GeoDataFrame = self.child_mosaic.footprint()
        footprint = pivot_geom.overlay(child_geom, "intersection")
        return footprint

    def extent(self) -> gpd.GeoDataFrame:
        """
        Get the extent of the pair, i.e. the intersection between pivot and child extents.

        Returns:
            gpd.GeoDataFrame: Extent of the mosaic

        """
        pivot_geom: gpd.GeoDataFrame = self.pivot_mosaic.extent()
        child_geom: gpd.GeoDataFrame = self.child_mosaic.extent()
        extent = pivot_geom.overlay(child_geom, "intersection")
        return extent

    def load(
        self,
        pivot_bands: Union[list, BandNames, str],
        child_bands: Union[list, BandNames, str],
        diff_bands: Union[list, BandNames, str],
        pixel_size: float = None,
        diff_method: DiffMethod = DiffMethod.PRE_POST,  # TODO: check that (dnbr = pre - post, idem ndvi)
        **kwargs,
    ) -> (dict, dict, dict):
        """"""
        pivot_dict = {}
        child_dict = {}
        diff_dict = {}

        # Convert just in case
        pivot_bands = to_band(pivot_bands)
        child_bands = to_band(child_bands)
        diff_bands = to_band(diff_bands)

        # Load bands

        # out = pair.load(pre_bands=[NBR, ...], post_bands=[NDVI, ...], diff_bands=[NBR, ...])
        # out: (out_pre, out_post, out_diff)
        # >>> out_pre:{NBR: nbr_arr_from_pre_over_pair_footprint}
        # >>> out_post: {NDVI: nbr_arr_from_pos_over_pair_footprint}
        # >>> out_diff: {NBR: dnbr_arr}

        # TODO

        # TODO: mange difference automatically

        return diff_dict, pivot_dict, child_dict

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
