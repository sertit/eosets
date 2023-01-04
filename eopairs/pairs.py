""" Class implementing the pairs """
from pathlib import Path
from typing import Union

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
        children_paths: Union[list, str, Path, CloudPath],
        output_path: Union[str, Path, CloudPath] = None,
        remove_tmp: bool = None,
        **kwargs,
    ):

        # Open reference product
        self.ref_prod: Product = READER.open(
            reference_path, remove_tmp, output_path, **kwargs
        )

        # Open children products
        self.children_prods = {}
        for path in children_paths:
            child_prod: Product = READER.open(path, remove_tmp, output_path, **kwargs)

            # Check the pair compatibility (if incompatible, the function throws a IncompatibleProducts error)
            self.check_compatibility()
            self.children_prods[child_prod.condensed_name] = child_prod

        # -- Other parameters --
        # TODO : create a temp folder for the pairs ?
        self.output_path = output_path

        # Full name
        # TODO (how to name pairs ???)
        self.full_name = f"{self.ref_prod.condensed_name}__{'_'.join(prod.condensed_name for prod in self.children_prods.values())}"

        # Condensed name
        # TODO (how to name pairs ???)

        # Nodata (by default use EOReader's)
        self.nodata = kwargs.get("nodata")

        # Resolution (by default use EOReader's)
        self.resolution = kwargs.get("resolution")

        # Information regarding the pair composition
        self.has_child = len(self.children_prods) > 0
        self.has_children = len(self.children_prods) > 1
        self.has_unique_child = len(self.children_prods) == 1
        self.same_constellation = self.homogeneous_attribute("constellation")
        self.same_sensor_type = self.homogeneous_attribute("sensor_type")
        self.constellations = list(
            set(
                [self.ref_prod.constellation]
                + [prod.constellation for prod in self.children_prods.values()]
            )
        )

        # if self.same_constellation:
        #     self.constellation = self.ref_prod.constellation
        # else:
        #     self.constellation = None

    def homogeneous_attribute(self, attr) -> bool:
        """"""
        return all(
            getattr(self.ref_prod, attr) == getattr(child, attr)
            for child in self.children_prods.values()
        )

    def check_compatibility(self) -> bool:
        """"""
        # TODO
        pass

    def overlapping_footprint(self) -> bool:
        """"""
        # TODO
        pass

    def overlapping_extent(self) -> bool:
        """"""
        # TODO
        pass

    def load(self) -> bool:
        """"""
        # TODO
        pass

    def stack(self) -> bool:
        """"""
        # TODO
        pass
