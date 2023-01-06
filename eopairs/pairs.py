""" Class implementing the pairs """
from pathlib import Path
from typing import Union

import geopandas as gpd
from cloudpathlib import CloudPath
from eoreader.products import Product
from eoreader.reader import Reader

# from eopairs.exceptions import IncompatibleProducts

READER = Reader()


class Pairs:
    """Class of multiple pairs"""

    def __init__(
        self,
        reference_path: Union[str, Path, CloudPath],
        children_paths: Union[list, str, Path, CloudPath] = None,
        output_path: Union[str, Path, CloudPath] = None,
        remove_tmp: bool = None,
        **kwargs,
    ):

        # Open reference product
        self.ref_prod: Product = READER.open(
            reference_path, remove_tmp, output_path, **kwargs
        )
        """ Reference product (unique). The one on which everything will be aligned. """

        # Open children products
        self.children_prods = {}
        """ Children products, to be aligned on the reference one. """
        if children_paths is None:
            children_paths = []
        if not isinstance(children_paths, list):
            children_paths = [children_paths]
        for path in children_paths:
            child_prod: Product = READER.open(path, remove_tmp, output_path, **kwargs)

            # Check the pair compatibility (if incompatible, the function throws a IncompatibleProducts error)
            self.check_compatibility()
            self.children_prods[child_prod.condensed_name] = child_prod

        # -- Other parameters --
        # TODO : create a temp folder for the pairs ?
        self.output_path = output_path
        """ Output path of the pairs. """

        # Full name
        # TODO (how to name pairs ???)
        self.full_name = f"{self.ref_prod.condensed_name}__{'_'.join(prod.condensed_name for prod in self.children_prods.values())}"
        """ Pairs full name. """

        # Condensed name
        # TODO (how to name pairs ???)

        # Nodata (by default use EOReader's)
        self.nodata = kwargs.get("nodata")
        """ Nodata of the pairs. """

        # Resolution (by default use EOReader's)
        self.resolution = kwargs.get("resolution")
        """ Resolution of the pairs. """

        # Information regarding the pair composition
        self.has_child = len(self.children_prods) > 0
        """ Does the pairs have at least one child? """

        self.has_children = len(self.children_prods) > 1
        """ Does the pairs have children? """

        self.has_unique_child = len(self.children_prods) == 1
        """ Does the pairs have a unique child? """

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

    def get_products_list(self) -> list:
        """
        Get all the products as a list. Reference is the first one.

        Returns:
            list: Products list
        """
        return [self.ref_prod] + list(self.children_prods.values())

    def homogeneous_attribute(self, attr: str) -> bool:
        """
        Check if the given attribute is the same for all products constituting the pairs.

        Args:
            attr (str): Attribute to be cecked. Must be available in EOReader's Product

        Returns:
            bool: True if this attribute is the same for all products constituting the pairs.
        """
        return all(
            getattr(self.ref_prod, attr) == getattr(child, attr)
            for child in self.children_prods.values()
        )

    def homogeneous_method(self, attr: str) -> bool:
        """
        Check if the given attribute is the same for all products constituting the pairs.

        Args:
            attr (str): Attribute to be cecked. Must be available in EOReader's Product

        Returns:
            bool: True if this attribute is the same for all products constituting the pairs.
        """
        return all(
            getattr(self.ref_prod, attr)() == getattr(child, attr)()
            for child in self.children_prods.values()
        )

    def check_compatibility(self) -> None:
        """
        Check if the products are coherent between each other, in order to create pairs
        If not, throws a IncompatibleProducts error.

        Raises:
            IncompatibleProducts: Products are incompatible and therefore pairs cannot be created.
        """
        # Check overlap

        # Check...
        pass

    def overlapping_footprint(self) -> gpd.GeoDataFrame:
        """
        Get the footprint of the overlapping area between every product of the pairs.

        Returns:
            gpd.GeoDataFrame: Footprint of the overlapping area
        """
        overlapping_footprint: gpd.GeoDataFrame = self.ref_prod.footprint()
        for prod in self.children_prods.values():
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
        for prod in self.children_prods.values():
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
