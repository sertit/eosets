""" Class implementing the series object """
from enum import unique
from pathlib import Path
from typing import Union

import geopandas as gpd
import xarray as xr
from cloudpathlib import CloudPath
from eoreader import cache
from eoreader.bands import BandNames, to_band
from eoreader.reader import Reader
from rasterio.enums import Resampling
from sertit.misc import ListEnum

from eosets.exceptions import IncompatibleProducts
from eosets.mosaic import Mosaic
from eosets.set import GeometryCheck, Set
from eosets.utils import AnyPathType

READER = Reader()


@unique
class Alignment(ListEnum):
    """Available alignment methods for time-series."""

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
        paths: Union[list, str, Path, CloudPath],
        id: str = None,
        output_path: Union[str, Path, CloudPath] = None,
        remove_tmp: bool = True,
        overlap_check: Union[GeometryCheck, str] = GeometryCheck.EXTENT,
        contiguity_check: Union[GeometryCheck, str] = GeometryCheck.EXTENT,
        alignement: Union[Alignment, str] = Alignment.FIRST,
        coregister: bool = False,
        ruling_mosaic: Union[Mosaic, int, str] = None,
        **kwargs,
    ):
        # Manage reference product
        self.mosaics = None
        """ Pivot mosaic (unique date and contiguous). The one on which the child will be aligned. """

        self.id = None
        """ ID of the pivot product """

        self.alignment = Alignment.convert_from(alignement)[0]
        """ Time Series alignment (on first mosaic, on last mosaic, on a hypothetical mean product). """

        self.coregister = coregister
        """ Do we need to coregister the time series? """

        self._ruling_mosaic = ruling_mosaic
        """ Ruling mosaic """

        self._unique_mosaic = len(paths) == 1

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

        # Fill attributes
        self.nodata = self.get_attr("nodata", **kwargs)
        self.pixel_size = self.get_attr("pixel_size", **kwargs)
        self.crs = self.get_attr("crs", **kwargs)
        self.same_constellation = self.is_homogeneous("constellation")
        self.constellations = list(set(prod.constellation for prod in self.get_prods()))
        self._unique_mosaic = len(self.get_prods()) == 1

    def clean_tmp(self):
        """
        Clean the temporary directory of the current product
        """
        for mos in self.mosaics:
            mos.clean_tmp()

    def clear(self):
        """
        Clear this product's cache
        """
        # Delete all cached properties and functions
        for mos in self.mosaics:
            mos.clear()

    def _manage_output(self):
        """
        Manage the output specifically for this child class
        """
        for mos in self.mosaics:
            mos._manage_output()

    def get_prods(self) -> list:
        """
        Get all the products as a list.

        Returns:
            list: Products list
        """
        return [mos.get_prods() for mos in self.mosaics]

    def _manage_mosaics(
        self,
        mosaic_paths: list,
        contiguity_check: GeometryCheck = GeometryCheck.EXTENT,
        overlap_check: GeometryCheck = GeometryCheck.EXTENT,
    ) -> None:
        """
        Check if all the mosaics are overlapping and have the same CRS. Each mosaic of a series should have a different datetime.
        If not, throws a IncompatibleProducts error.

        Raises:
            IncompatibleProducts: Incompatible products if not contiguous or not the same date
        """
        # Information regarding the pair composition
        assert isinstance(mosaic_paths, list)

        self.mosaics: list[Mosaic] = []
        for paths in mosaic_paths:
            # Open the mosaic
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

        # Make the checks
        # CRS
        if any(self.ruling_mosaic.crs != mos.crs for mos in self.mosaics):
            raise IncompatibleProducts(
                f"All mosaics should have the same CRS! (All mosaics should be aligned on {self.ruling_mosaic.crs=})"
            )

        # Geometry
        if overlap_check != GeometryCheck.NONE:
            ruling_geom: gpd.GeoDataFrame = getattr(
                self.ruling_mosaic, str(overlap_check.value)
            )()
            for mos in self.mosaics:
                if mos.id != self.ruling_mosaic.id:
                    mos_geom: gpd.GeoDataFrame = getattr(
                        mos, str(overlap_check.value)
                    )()
                    if not ruling_geom.intersects(
                        mos_geom.to_crs(self.ruling_mosaic.crs)
                    ).all():
                        raise IncompatibleProducts("All mosaics should overlap!")

        # Fill other attributes
        self.nof_prods = len(self.get_prods())
        self.id = "_".join(mos.id for mos in self.mosaics)
        self.full_name = "_".join(mos.full_name for mos in self.mosaics)
        self.condensed_name = self.full_name
        # TODO (how to name series ???)

    @property
    def ruling_mosaic(self) -> Mosaic:
        """
        Get the ruling mosaic of the series

        Returns:
            Mosaic: Output path ofthe mosaic
        """
        if self.alignment in [Alignment.FIRST, Alignment.LAST]:
            return self.mosaics[self._ruling_mosaic]
        elif self.alignment == Alignment.MEAN:
            raise NotImplementedError
        elif self.alignment == Alignment.EXTERNAL:
            return self._ruling_mosaic
        elif self.alignment == Alignment.CUSTOM:
            if isinstance(self._ruling_mosaic, int):
                return self.mosaics[self._ruling_mosaic]

    @ruling_mosaic.setter
    def ruling_mosaic(self, mosaic: Union[Mosaic, int] = None):
        if self.alignment == Alignment.FIRST:
            self._ruling_mosaic = 0
        elif self.alignment == Alignment.LAST:
            self._ruling_mosaic = -1
        elif self.alignment == Alignment.MEAN:
            raise NotImplementedError
        elif self.alignment == Alignment.EXTERNAL:
            assert mosaic is not None
            self._ruling_mosaic = mosaic
        elif self.alignment == Alignment.CUSTOM:
            assert isinstance(mosaic, (Mosaic, int))
            if isinstance(mosaic, int):
                assert abs(mosaic) <= len(self.mosaics)
            self._ruling_mosaic = mosaic

    def read_mtd(self):
        """"""
        # TODO: how ? Just return the fields that are shared between mosaic's components ? Or create a XML from scratch ?
        raise NotImplementedError

    @cache
    def footprint(self) -> gpd.GeoDataFrame:
        """
        Get the footprint of the series, i.e. the intersection between all mosaics in the ruling mosaic's CRS.

        Returns:
            gpd.GeoDataFrame: Footprint of the mosaic
        """
        footprint = None
        for mos in self.mosaics:
            geom: gpd.GeoDataFrame = mos.footprint().to_crs(self.ruling_mosaic.crs)
            if not footprint:
                footprint = geom
            else:
                footprint = footprint.overlay(geom, "intersection")

        return footprint

    @cache
    def extent(self) -> gpd.GeoDataFrame:
        """
        Get the extent of the series, i.e. the intersection between all mosaics in the ruling mosaic's CRS.

        Returns:
            gpd.GeoDataFrame: Extent of the mosaic

        """
        extent = None
        for mos in self.mosaics:
            geom: gpd.GeoDataFrame = mos.footprint().to_crs(self.ruling_mosaic.crs)
            if not extent:
                extent = geom
            else:
                extent = extent.overlay(geom, "intersection")

        return extent

    def load(
        self,
        bands: Union[list, BandNames, str] = None,
        pixel_size: float = None,
        resampling: Resampling = Resampling.bilinear,
        **kwargs,
    ) -> (dict, dict, dict):
        """"""
        bands_dict = {}
        bands = to_band(bands)

        # Load bands
        window = kwargs.pop("window", self.footprint())

        # Load pivot bands
        ruling_ds = self.ruling_mosaic.load(
            bands, pixel_size=pixel_size, window=window, **kwargs
        )

        # Load mosaic bands
        for mos in self.mosaics:
            if mos.id != self.ruling_mosaic.id:
                mos_ds = mos.load(bands, pixel_size=pixel_size, window=window, **kwargs)

                # Add the bands to the dataset
                for band in bands:
                    mos_arr = mos_ds[band]
                    ruling_arr = ruling_ds[band]

                    # To be sure, always collocate arrays, even if the size is the same
                    # Indeed, a small difference in the coordinates will lead to empy arrays
                    # So the bands MUST BE exactly aligned
                    mos_ds[band] = mos_arr.rio.reproject_match(
                        ruling_arr, resampling=resampling, **kwargs
                    )
            else:
                mos_ds = ruling_ds

            # Set a time coordinate to distinguish the mosaics
            dt = mos.datetime()
            bands_dict[dt] = mos_ds.assign_coords({"time": dt})

        # Create a dataset (only after collocation)
        coords = None
        if bands_dict:
            coords = ruling_ds.coords
            coords["time"] = bands_dict.keys()

        # Make sure the dataset has the bands in the right order -> re-order the input dict
        return xr.Dataset({key: bands_dict[key] for key in bands}, coords=coords)

    def stack(
        self,
        bands: list,
        pixel_size: float = None,
        stack_path: Union[str, AnyPathType] = None,
        save_as_int: bool = False,
        **kwargs,
    ) -> xr.DataArray:
        """
        Stack bands and index of a series.

        Args:
            bands (list): Bands and index combination
            pixel_size (float): Stack pixel size. . If not specified, use the product pixel size.
            stack_path (Union[str, AnyPathType]): Stack path
            save_as_int (bool): Convert stack to uint16 to save disk space (and therefore multiply the values by 10.000)
            **kwargs: Other arguments passed to :code:`load` or :code:`rioxarray.to_raster()` (such as :code:`compress`)

        Returns:
            xr.DataArray: Stack as a DataArray
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
