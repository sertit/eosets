""" Class implementing the pairs """
import logging
import os
import shutil
import tempfile
from collections import defaultdict
from enum import unique
from glob import glob
from pathlib import Path
from typing import Any, Tuple, Union

import geopandas as gpd
import xarray as xr
from cloudpathlib import AnyPath, CloudPath
from eoreader import utils
from eoreader.bands import BandNames, is_spectral_band, to_band, to_str
from eoreader.products import Product
from eoreader.reader import Reader
from eoreader.utils import UINT16_NODATA
from sertit import files, rasters
from sertit.misc import ListEnum

from eosets.env_vars import CI_EOSETS_BAND_FOLDER
from eosets.exceptions import IncompatibleProducts
from eosets.utils import EOPAIRS_NAME, AnyPathType

READER = Reader()

LOGGER = logging.getLogger(EOPAIRS_NAME)


@unique
class ContiguityCheck(ListEnum):
    """Available contiguity checks."""

    FOOTPRINT = "footprint"
    EXTENT = "extent"
    NONE = "none"


@unique
class MosaicMethod(ListEnum):
    """Available mosaicing methods."""

    GTIFF = "merge_gtiff"
    VRT = "merge_vrt"


class Mosaic:
    """Class of mosaic, composed by several contiguous EOReader's products acquired the same day"""

    def __init__(
        self,
        paths: Union[list, str, AnyPathType],
        output_path: Union[str, AnyPathType] = None,
        id: str = None,
        remove_tmp: bool = True,
        contiguity_check: Union[ContiguityCheck, str] = ContiguityCheck.EXTENT,
        mosaic_method: Union[MosaicMethod, str] = MosaicMethod.VRT,
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

        # Manage output path
        if output_path:
            self._tmp_output = None
            self._output = AnyPath(output_path)
        else:
            self._tmp_output = tempfile.TemporaryDirectory()
            self._output = AnyPath(self._tmp_output.name)

        self._tmp_process = self.output.joinpath("tmp_mosaic")
        os.makedirs(self._tmp_process, exist_ok=True)

        # Manage reference product
        self.prods: dict = {}
        """ Products (contiguous and acquired the same day). """

        self.id: str = id
        """ ID of the reference product, given by the creator of the mosaic. If not, a mix based on the dates and constellations of its components. """

        self.nof_prods: int = 0
        """ Number of products. """

        # -- Other parameters --
        # Full name
        self.full_name: str = ""
        """ Mosaic full name. """

        # Condensed name
        self.condensed_name = ""
        """ Mosaic condensed name, a mix based on the dates and constellations of the components of the mosaic. """

        # We need the date in _manage_prods
        self.date = None
        """ Date of the mosaic. If not provided in kwargs, using the first product's date. """

        contiguity_check = ContiguityCheck.convert_from(contiguity_check)[0]
        self._manage_prods(paths, contiguity_check, **kwargs)

        # Nodata (by default use EOReader's)

        # Nodata (by default use EOReader's)
        self.nodata = self.get_attr("nodata", **kwargs)
        """ Nodata of the mosaic. If not provided in kwargs, using the first product's nodata. """

        # Resolution (by default use EOReader's)
        self.resolution = self.get_attr("resolution", **kwargs)
        """ Resolution of the mosaic. If not provided in kwargs, using the first product's resolution. """

        self.crs = self.get_attr("crs", **kwargs)
        """ CRS of the mosaic. If not provided in kwargs, using the first product's crs. """

        self.same_constellation: bool = self.is_homogeneous("constellation")
        """ Is the mosaic constituted of the same constellation? """

        self.same_crs: bool = self.is_homogeneous("crs")
        """ Is the mosaic constituted of the same sensor type? """

        self.constellations = list(set(prod.constellation for prod in self.get_prods()))
        """ List of unique constellations constituting the pairs """

        self.mosaic_method = MosaicMethod.convert_from(mosaic_method)[0]
        """ Mosaicing method. If GTIFF is specified, the temporary files from every products will be removed, if VRT is spoecified, they will not."""

    def clean_tmp(self):
        """
        Clean the temporary directory of the current product
        """
        for prod in self.get_prods():
            prod.clean_tmp()

    def clear(self):
        """
        Clear this product's cache
        """
        # -- Delete all cached properties and functions
        for prod in self.get_prods():
            prod.clear()

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
        for prod in self.get_prods():
            prod.output = self._get_tmp_folder(writable=True)

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

    def _manage_prods(
        self,
        paths: Union[list, str, Path, CloudPath],
        contiguity_check: ContiguityCheck,
        **kwargs,
    ):
        """
        Manage products attributes and check the compatibility of the mosaic's components

        Args:
            paths (Union[list, str, Path, CloudPath]): Paths of the mosaic
            contiguity_check (ContiguityCheck): Method to check the contiguity of the mosaic
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

        # Open first product as a reference
        first_prod: Product = READER.open(paths[0], remove_tmp=remove_tmp, **kwargs)
        if first_prod is None:
            raise ValueError(
                f"There is no existing products in EOReader corresponding to {paths[0]}"
            )

        self.prods[first_prod.condensed_name] = first_prod

        # Open others
        for path in paths[1:]:
            prod: Product = READER.open(path, remove_tmp=remove_tmp, **kwargs)
            if prod is None:
                raise ValueError(
                    f"There is no existing products in EOReader corresponding to {path}"
                )

            # Ensure compatibility of the mosaic component, i.e. unique date and contiguous product
            self.check_compatibility(first_prod, prod)
            self.prods[prod.condensed_name] = prod

        self.check_contiguity(contiguity_check)

        # Create full_name
        self.nof_prods = len(self.prods)
        self.date = first_prod.date.date()
        self.full_name = (
            f"{'-'.join([prod.condensed_name for prod in self.get_prods()])}"
        )

        # Create condensed_name: [{date}-{sat_id}_]{???}
        # TODO: is it OK ?
        # TODO: if fixed date, change that
        # TODO: if all same constellation, set it only once
        # TODO: add sth ?
        self.condensed_name = f"{self.date.strftime('%Y%m%d')}_{'-'.join(list(set([prod.constellation_id for prod in self.get_prods()])))}"

        # Rename tmp_process and set product outputs
        self._tmp_process = AnyPath(
            shutil.move(
                str(self._tmp_process),
                str(self.output.joinpath(f"tmp_{self.condensed_name}")),
            )
        )
        for prod in self.get_prods():
            prod.output = self._get_tmp_folder(writable=True)

        if self.id is None:
            self.id = self.condensed_name

    def get_prods(self) -> list:
        """
        Get all the products as a list.

        Returns:
            list: Products list
        """
        return list(self.prods.values())

    def get_first_prod(self) -> Product:
        """
        Get first product, which should be coherent with all others

        Returns:
            Product: First reference product
        """
        return self.get_prods()[0]

    def get_attr(self, attr: str, **kwargs) -> Any:
        """
        Get attribute, either from kwargs or from the first product (default)

        Args:
            attr (str): Wanted attribute
            **kwargs: Other args

        Returns:
            Any: Attribute result
        """
        attr = kwargs.pop("crs", getattr(self.get_first_prod(), attr))
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
                    ref_attr() == getattr(child, attr)()
                    for child in self.get_prods()[1:]
                )
            else:
                is_homogeneous = all(
                    ref_attr == getattr(child, attr) for child in self.get_prods()[1:]
                )
        else:
            is_homogeneous = True

        return is_homogeneous

    def check_compatibility(self, first_prod, prod) -> None:
        """
        Check if the mosaic products are coherent between each other.
        - Same sensor type
        - Same date

        TODO: same constellation ? same CRS ?...

        If not, throws a IncompatibleProducts error.

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

    def check_contiguity(self, check_contiguity: ContiguityCheck):
        """
        Check the contiguity of the mosaic

        Args:
            check_contiguity (ContiguityCheck): Contiguity checking method

        Raises:
            IncompatibleProducts: Incompatible products if not contiguous according to the given method
        """
        if check_contiguity == ContiguityCheck.EXTENT:
            union_extent = self.extent()
            if len(union_extent) > 1:
                raise IncompatibleProducts(
                    "The mosaic should have a contiguous extent!"
                )
        elif check_contiguity == ContiguityCheck.FOOTPRINT:
            union_footprint = self.footprint()
            if len(union_footprint) > 1:
                raise IncompatibleProducts(
                    "The mosaic should have a contiguous footprint!"
                )
        else:
            LOGGER.warning("The contiguity of your mosaic won't be checked!")
            pass

    def read_mtd(self):
        """"""
        # TODO: how ? Just return the fields that are shared between mosaic's components ? Or create a XML from scratch ?
        raise NotImplementedError

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
                footprint = footprint.overlay(prod.footprint().to_crs(ref_prod.crs()))

            # Dissolve and explode the footprint
            footprint = footprint.dissolve().explode()

        return footprint

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
                extent = extent.overlay(prod.extent().to_crs(ref_prod.crs()))

            # Dissolve and explode the extent
            extent = extent.dissolve().explode()

        return extent

    def has_band(self, band: Union[BandNames, str]) -> bool:
        """
        Does this moasic have products with the specified band ?

        By band, we mean:

        - satellite band
        - index
        - DEM band
        - cloud band

        Args:
            band (Union[BandNames, str]): EOReader band (optical, SAR, clouds, DEM)

        Returns:
            bool: True if the products has the specified band
        """
        return all(prod.has_band(band) for prod in self.get_prods())

    def has_bands(self, bands: Union[list, BandNames, str]) -> bool:
        """
        Does this moasic have products with the specified bands ?

        By band, we mean:

        - satellite band
        - index
        - DEM band
        - cloud band

        See :code:`has_band` for a code example.

        Args:
            bands (Union[list, BandNames, str]): EOReader bands (optical, SAR, clouds, DEM)

        Returns:
            bool: True if the products has the specified band
        """

        return all(prod.has_bands(bands) for prod in self.get_prods())

    def load(
        self,
        bands: Union[list, BandNames, str],
        resolution: float = None,
        **kwargs,
    ) -> dict:
        """"""
        # Get merge function and extension
        merge_fct = getattr(rasters, self.mosaic_method.value)
        out_suffix = f"{self.mosaic_method.name.lower()}"

        # Convert just in case
        bands = to_band(bands)

        # Get the bands to be loaded
        band_paths = {}
        bands_to_load = []
        for band in bands:
            band_path, exists = self._get_out_path(
                f"{self.id}_{to_str(band)[0]}.{out_suffix}"
            )
            band_paths[band] = band_path
            if not exists:
                bands_to_load.append(band)

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
            prod.load(bands_to_load, resolution, **kwargs).keys()

            # Store paths
            for band in bands_to_load:
                if is_spectral_band(band):
                    band_path = prod.get_band_paths([band], resolution, **kwargs)[band]
                else:
                    # Use glob fct as _get_band_folder is a tmpDirectory
                    band_path = glob(
                        os.path.join(prod._get_band_folder(), f"*{to_str(band)[0]}*")
                    )[0]

                prod_band_paths[band].append(str(band_path))

        # Merge
        merged_dict = {}
        for band in band_paths:
            output_path = band_paths[band]
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
                merge_fct(prod_path, output_path, **kwargs)

            # Load in memory and update attribute
            merged_dict[band] = self._update_attrs(
                rasters.read(output_path), bands, **kwargs
            )

        # Collocate VRTs
        LOGGER.debug("Collocating bands")
        merged_dict = self._collocate_bands(merged_dict)

        return merged_dict

    def _collocate_bands(self, bands: dict, master_xds: xr.DataArray = None) -> dict:
        """
        Collocate all bands from a dict if needed (if a raster shape is different)

        Args:
            bands (dict): Dict of bands to collocate if needed
            master_xds (xr.DataArray): Master array

        Returns:
            dict: Collocated bands
        """
        return self.get_first_prod()._collocate_bands(bands, master_xds)

    def stack(
        self,
        bands: list,
        resolution: float = None,
        stack_path: Union[str, AnyPathType] = None,
        save_as_int: bool = False,
        **kwargs,
    ) -> xr.DataArray:
        """
        Stack bands and index of a mosaic.

        Args:
            bands (list): Bands and index combination
            resolution (float): Stack resolution. . If not specified, use the product resolution.
            stack_path (Union[str, AnyPathType]): Stack path
            save_as_int (bool): Convert stack to uint16 to save disk space (and therefore multiply the values by 10.000)
            **kwargs: Other arguments passed to :code:`load` or :code:`rioxarray.to_raster()` (such as :code:`compress`)

        Returns:
            xr.DataArray: Stack as a DataArray
        """
        # Create the analysis stack
        band_dict = self.load(bands, resolution=resolution, **kwargs)

        # Stack bands
        if save_as_int:
            nodata = kwargs.get("nodata", UINT16_NODATA)
        else:
            nodata = kwargs.get("nodata", self.nodata)
        stack, dtype = utils.stack_dict(bands, band_dict, save_as_int, nodata, **kwargs)

        # Update stack's attributes
        stack = self._update_attrs(stack, bands, **kwargs)

        # Write on disk
        LOGGER.debug("Saving stack")
        if stack_path:
            stack_path = AnyPath(stack_path)
            if not stack_path.parent.exists():
                os.makedirs(str(stack_path.parent), exist_ok=True)

            utils.write(stack, stack_path, dtype=dtype, **kwargs)

        return stack

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
        xarr.attrs["acquisition_date"] = self.date
        xarr.attrs["condensed_name"] = self.condensed_name

        # TODO: complete that

        return xarr
