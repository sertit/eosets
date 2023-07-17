"""
**EOSets** library
"""
# flake8: noqa
from .__meta__ import (
    __author__,
    __author_email__,
    __copyright__,
    __description__,
    __documentation__,
    __license__,
    __title__,
    __url__,
    __version__,
)

EOSETS_NAME = "eosets"

__all__ = ["Mosaic", "Pair", "Series"]

from .mosaic import Mosaic
from .pair import Pair
from .series import Series
