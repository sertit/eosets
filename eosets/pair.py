""" Class implementing the pairs """
from pathlib import Path
from typing import Union

import geopandas as gpd
from cloudpathlib import CloudPath
from eoreader.products import Product
from eoreader.reader import Reader

# from eosets.exceptions import IncompatibleProducts

READER = Reader()


class Pair:
    """Class of multiple pairs"""

    def __init__(
        self,
        pivot_paths: Union[list, str, Path, CloudPath],
        child_paths: Union[list, str, Path, CloudPath] = None,
        pivot_id: str = None,
        child_id: str = None,
        output_path: Union[str, Path, CloudPath] = None,
        remove_tmp: bool = True,
        **kwargs,
    ):
        """
        Initialization of a Pair

        Args:
            pivot_paths:
            child_paths:
            pivot_id:
            child_id:
            output_path:
            remove_tmp:
            **kwargs:
        """
        # Manage output
        # TODO : create a temp folder for the pairs ?
        self.output_path = output_path
        """ Output path of the pairs. """

        # Remove temporary files
        self._remove_tmp = remove_tmp
        """ Remove temporary files, propagated to EOReader's Products. """

        # Manage reference product
        self.reference_prods: dict = {}
        """ Pivot products (unique date and contiguous). The ones on which everything will be aligned. """

        self.pivot_id: str = pivot_id
        """ ID of the pivot product """

        self.is_pivot_single: bool = True
        """ If the pivot products are composed of one product only. """

        self._manage_pivot(pivot_paths, pivot_id, **kwargs)

        # Open child products
        self.child_prods: dict = {}
        """ Child products (unique date and contiguous), to be aligned on the reference one. """

        self.child_id: str = child_id
        """ ID of the child product """

        self.is_child_single: bool = True
        """ If the child products are composed of one product only. """

        # Information regarding the pair composition
        self.has_child: bool = True
        """ Does the pair have a child? (Pair with only one reference is allowed) """

        self._manage_child(child_paths, child_id, **kwargs)

        # -- Other parameters --
        # Full name
        self.full_name = f"{self.pivot_id}_{self.child_id}"
        """ Pairs full name. """

        # Condensed name
        self.condensed_name = self.full_name
        # TODO (how to name pairs ???)

        # Nodata (by default use EOReader's)
        self.nodata = kwargs.get("nodata")
        """ Nodata of the pairs. """

        # Resolution (by default use EOReader's)
        self.resolution = kwargs.get("resolution")
        """ Resolution of the pairs. """

        self.same_constellation = self.homogeneous_attribute("constellation")
        """ Are the pairs constituted of the same constellation? """

        self.same_sensor_type = self.homogeneous_attribute("sensor_type")
        """ Are the pairs constituted of the same sensor type? """

        self.same_crs = self.homogeneous_method("crs")
        """ Are the pairs constituted of the same sensor type? """

        self.constellations = list(
            set(prod.constellation for prod in self.get_products_list())
        )
        """ List of unique constellations constitutig the pairs """

        # if self.same_constellation:
        #     self.constellation = self.ref_prod.constellation
        # else:
        #     self.constellation = None

    def _manage_pivot(self, reference_paths, reference_id, **kwargs):
        # Manage reference product

        if reference_paths is None:
            reference_paths = []
        elif not isinstance(reference_paths, list):
            reference_paths = [reference_paths]
        for path in reference_paths:
            ref_prod: Product = READER.open(
                path,
                remove_tmp=self._remove_tmp,
                output_path=self.output_path,
                **kwargs,
            )

            # TODO: ensure unique date and contiguous product
            self.check_reference_compatibility()
            self.reference_prods[ref_prod.condensed_name] = ref_prod

        if self.pivot_id is None:
            # TODO
            pass

    def _manage_child(self, child_paths, child_id, **kwargs):

        if child_paths is None:
            child_paths = []
        elif not isinstance(child_paths, list):
            child_paths = [child_paths]
        for path in child_paths:
            child_prod: Product = READER.open(
                path,
                remove_tmp=self._remove_tmp,
                output_path=self.output_path,
                **kwargs,
            )

            # Check the pair compatibility (if incompatible, the function throws a IncompatibleProducts error)
            self.check_child_compatibility()
            self.child_prods[child_prod.condensed_name] = child_prod

        if self.child_id is None:
            # TODO
            pass

    def get_products_list(self) -> list:
        """
        Get all the products as a list. Reference is the first one.

        Returns:
            list: Products list
        """
        return list(self.reference_prods.values()) + list(self.child_prods.values())

    def get_first_reference_prod(self) -> Product:
        """
        Get first reference product, that should be coherent with all others

        Returns:
            Product: First reference product
        """
        return list(self.reference_prods.values())[0]

    def homogeneous_attribute(self, attr: str) -> bool:
        """
        Check if the given attribute is the same for all products constituting the pairs.

        Args:
            attr (str): Attribute to be cecked. Must be available in EOReader's Product

        Returns:
            bool: True if this attribute is the same for all products constituting the pairs.
        """
        ref_attr = getattr(self.get_first_reference_prod(), attr)

        return all(
            ref_attr == getattr(child, attr) for child in self.child_prods.values()
        )

    def homogeneous_method(self, attr: str) -> bool:
        """
        Check if the given method (with empty arguments) is the same for all products constituting the pairs.

        Args:
            attr (str): Method to be cecked. Must be available in EOReader's Product

        Returns:
            bool: True if this method is the same for all products constituting the pairs.
        """
        ref_method = getattr(self.get_first_reference_prod(), attr)()

        return all(
            ref_method == getattr(child, attr)() for child in self.child_prods.values()
        )

    def check_reference_compatibility(self) -> None:
        """
        Check if the reference products are coherent between each other, in order to create a coherent reference product

        TODO: same constellation ? same CRS ?...

        If not, throws a IncompatibleProducts error.

        Raises:
            IncompatibleProducts: Products are incompatible and therefore pairs cannot be created.
        """
        # Check contiguous

        # Check date ? necessary ?
        pass

    def check_child_compatibility(self) -> None:
        """
        Check if the products are coherent between each other, in order to create pairs
        If not, throws a IncompatibleProducts error.

        Raises:
            IncompatibleProducts: Products are incompatible and therefore pairs cannot be created.
        """
        # Check overlap

        # Check...
        pass

    def get_reference_mtd(self):
        """"""
        raise NotImplementedError

    def overlapping_footprint(self) -> gpd.GeoDataFrame:
        """
        Get the footprint of the overlapping area between every product of the pairs.

        Returns:
            gpd.GeoDataFrame: Footprint of the overlapping area
        """
        overlapping_footprint: gpd.GeoDataFrame = self.ref_prod.footprint()
        for prod in self.child_prods.values():
            overlapping_footprint = overlapping_footprint.overlay(
                prod.footprint().to_crs(self.ref_prod.crs())
            )

        return overlapping_footprint

    def overlapping_extent(self) -> gpd.GeoDataFrame:
        """
        Get the extent of the overlapping area between every product of the pairs.

        Returns:
            gpd.GeoDataFrame: Extent of the overlapping area

        """
        overlapping_extent: gpd.GeoDataFrame = self.ref_prod.extent()
        for prod in self.child_prods.values():
            overlapping_extent = overlapping_extent.overlay(
                prod.extent().to_crs(self.ref_prod.crs())
            )

        return overlapping_extent

    def load(self) -> bool:
        """"""
        # TODO
        pass

    def stack(self) -> bool:
        """"""
        # TODO
        pass
