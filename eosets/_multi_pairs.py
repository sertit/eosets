""" Class implementing the multi-pairs. Load and stack raise NotImplementedError for now. Please write an issue if you see any usecase. """
import logging
from pathlib import Path
from typing import Union

import geopandas as gpd
import xarray as xr
from cloudpathlib import CloudPath
from eoreader import cache
from eoreader.bands import BandNames
from rasterio.enums import Resampling

from eosets import EOSETS_NAME
from eosets.mosaic import Mosaic
from eosets.pair import DiffMethod, Pair
from eosets.set import GeometryCheck, Set
from eosets.utils import AnyPathType

LOGGER = logging.getLogger(EOSETS_NAME)


class MultiPairs(Set):
    """Class of multi-products pairs, but with the same pivot product. Load and stack raise NotImplementedError for now. Please write an issue if you see any usecase."""

    def __init__(
        self,
        pivot_paths: Union[list, str, Path, CloudPath],
        children_paths: list = None,
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

        self.children_mosaics = None
        """ Child mosaics (each mosaic should have a unique date and be contiguous). These mosaics will be aligned on the pivot. """

        self.children_id = None
        """ ID of the children products """

        # Information regarding the pair composition
        self.has_children = None
        """ Does the pair have a children? (Pair with only one pivot is allowed) """

        self.pairs = None
        """ All pairs contained in this multi-pairs object. """

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
        self._manage_mosaics(
            pivot_paths, children_paths, contiguity_check, overlap_check
        )

        # Fill attributes
        self.full_name = f"{self.pivot_id}_{self.children_id}"
        self.condensed_name = self.full_name
        # TODO (how to name multi_pairs ???)

        self.nodata = self.get_attr("nodata", **kwargs)
        self.pixel_size = self.get_attr("pixel_size", **kwargs)
        self.crs = self.get_attr("crs", **kwargs)
        self.same_constellation = self.is_homogeneous("constellation")
        self.constellations = list(set(prod.constellation for prod in self.get_prods()))

    def clean_tmp(self):
        """
        Clean the temporary directory of the current multi-pairs
        """
        self.pivot_mosaic.clean_tmp()
        for mos in self.children_mosaics:
            mos.clean_tmp()

    def clear(self):
        """
        Clear this multi-pairs' cache
        """
        # Delete all cached properties and functions
        self.pivot_mosaic.clear()
        for mos in self.children_mosaics:
            mos.clear()

    def _manage_output(self):
        """
        Manage the output specifically for this child class
        """
        self.pivot_mosaic._manage_output()
        for mos in self.children_mosaics:
            mos._manage_output()

    def get_prods(self) -> list:
        """
        Get all the products as a list.

        Returns:
            list: Products list
        """
        prods = self.pivot_mosaic.get_prods()
        for mos in self.children_mosaics:
            prods += mos.get_prods()

        return prods

    def _manage_mosaics(
        self,
        pivot_paths: Union[list, str, Path, CloudPath],
        children_paths: Union[list, str, Path, CloudPath] = None,
        contiguity_check: GeometryCheck = GeometryCheck.EXTENT,
        overlap_check: GeometryCheck = GeometryCheck.EXTENT,
    ) -> None:
        """
        Check if the pivot and children mosaics are overlapping and if their CRS are the same.

        TODO: check if same constellation ?

        If not, throws a IncompatibleProducts error.

        Args:
            pivot_paths (Union[list, str, Path, CloudPath, Mosaic]): Paths corresponding to the pivot mosaic
            children_paths (Union[list, str, Path, CloudPath, Mosaic]): Paths corresponding to the children mosaics
            contiguity_check (GeometryCheck): Check regarding the contiguity of the products of the mosaics
            overlap_check (GeometryCheck): Check regarding the overlapping of the two mosaics

        Raises:
            IncompatibleProducts: Incompatible products if not contiguous or not the same date
        """
        # Manage reference product
        if isinstance(pivot_paths, Mosaic):
            self.pivot_mosaic = pivot_paths
        else:
            self.pivot_mosaic: Mosaic = Mosaic(
                pivot_paths,
                output_path=self._get_tmp_folder(writable=True),
                remove_tmp=self._remove_tmp,
                contiguity_check=contiguity_check,
            )
        self.pivot_id: str = self.pivot_mosaic.id

        # Information regarding the pair composition
        self.has_children: bool = len(children_paths) > 0

        if self.has_children:
            self.children_mosaics = []
            self.pairs = []
            for child_paths in children_paths:
                if isinstance(child_paths, Mosaic):
                    mos = child_paths
                else:
                    mos: Mosaic = Mosaic(
                        child_paths,
                        output_path=self._get_tmp_folder(writable=True),
                        remove_tmp=self._remove_tmp,
                        contiguity_check=contiguity_check,
                    )
                self.pairs.append(
                    Pair(
                        self.pivot_mosaic,
                        mos,
                        contiguity_check=contiguity_check,
                        overlap_check=overlap_check,
                    )
                )

                # Everything is OK, store the mosaic as a child
                self.children_mosaics.append(mos)

            # Create children ID
            self.children_id: str = "_".join(mos.id for mos in self.children_mosaics)

        # Get number of products
        self.nof_prods = len(self.get_prods())

    def read_mtd(self):
        """Read the pair's metadata, but not implemented for now."""
        # TODO: how ? Just return the fields that are shared between multi_pairs' components ? Or create a XML from scratch ?
        raise NotImplementedError

    @cache
    def footprint(self) -> gpd.GeoDataFrame:
        """
        Get the footprint of the multi-pairs, i.e. the intersection between pivot and children footprints.

        Returns:
            gpd.GeoDataFrame: Footprint of the multi-pairs
        """
        footprint: gpd.GeoDataFrame = self.pivot_mosaic.footprint()

        for mos in self.children_mosaics:
            footprint = footprint.overlay(
                mos.footprint().to_crs(footprint.crs), "intersection"
            )

        return footprint

    @cache
    def extent(self) -> gpd.GeoDataFrame:
        """
        Get the extent of the multi-pairs, i.e. the intersection between pivot and children extents.

        Returns:
            gpd.GeoDataFrame: Extent of the multi-pairs

        """
        extent: gpd.GeoDataFrame = self.pivot_mosaic.footprint()

        for mos in self.children_mosaics:
            extent = extent.overlay(mos.extent().to_crs(extent.crs), "intersection")

        return extent

    def load(
        self,
        pivot_bands: Union[list, BandNames, str] = None,
        children_bands: Union[list, BandNames, str] = None,
        diff_bands: Union[list, BandNames, str] = None,
        pixel_size: float = None,
        diff_method: DiffMethod = DiffMethod.PIVOT_CHILD,
        resampling: Resampling = Resampling.bilinear,
        **kwargs,
    ) -> (xr.Dataset, xr.Dataset, xr.Dataset):
        """
        **NotImplemented**

        Load the bands and compute the wanted spectral indices for pivot, child and diff.

        Args:
            pivot_bands (Union[list, BandNames, str]): Wanted pivot bands
            children_bands (Union[list, BandNames, str]): Wanted child bands
            diff_bands (Union[list, BandNames, str]): Wanted diff bands
            pixel_size (float): Pixel size of the returned Datasets. If not specified, use the pair's pixel size.
            diff_method (DiffMethod): Difference method for the computation of diff_bands
            resampling (Resampling): Resampling method
            kwargs: Other arguments used to load bands

        Returns:
            (xr.Dataset, xr.Dataset, xr.Dataset): Pivot, child and diff wanted bands as xr.Datasets
        """
        raise NotImplementedError

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
        children_bands: Union[list, BandNames, str] = None,
        diff_bands: Union[list, BandNames, str] = None,
        pixel_size: float = None,
        diff_method: DiffMethod = DiffMethod.PIVOT_CHILD,
        stack_path: Union[str, AnyPathType] = None,
        save_as_int: bool = False,
        **kwargs,
    ) -> xr.DataArray:
        """
        **NotImplemented**

        Stack bands and index of a pair.

        Args:
            pivot_bands (list): Bands and index combination for the pivot mosaic
            children_bands (list): Bands and index combination for the child mosaic
            diff_bands (list): Bands and index combination for the difference between pivot and child mosaic
            pixel_size (float): Stack pixel size. . If not specified, use the product pixel size.
            stack_path (Union[str, AnyPathType]): Stack path
            save_as_int (bool): Convert stack to uint16 to save disk space (and therefore multiply the values by 10.000)
            **kwargs: Other arguments passed to :code:`load` or :code:`rioxarray.to_raster()` (such as :code:`compress`)

        Returns:
            xr.DataArray: Stack as a DataArray
        """
        raise NotImplementedError

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
