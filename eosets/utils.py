# Copyright 2026, SERTIT-ICube - France, https://sertit.unistra.fr/
# This file is part of eosets project
#     https://github.com/sertit/eosets
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Utils file"""

import contextlib
import logging

from eoreader import utils
from eoreader.products import Product
from sertit import AnyPath
from sertit.types import AnyPathStrType, AnyPathType

from eosets import EOSETS_NAME

read = utils.read
write = utils.write
stack = utils.stack
convert_to_uint16 = utils.convert_to_uint16
write_path_in_attrs = utils.write_path_in_attrs

AnyProductType = AnyPathStrType | Product
""" Any Product Type, either a path or an eoreader.Product"""

LOGGER = logging.getLogger(EOSETS_NAME)

# Band(s) Type
try:
    from eoreader.bands import (
        BandNames,
        BandsType,
        CloudsBandNames,
        DemBandNames,
        SarBandNames,
        SpectralBandNames,
        is_spectral_band,
        to_str,
    )

    BandType = (
        str
        | SpectralBandNames
        | SarBandNames
        | CloudsBandNames
        | DemBandNames
        | BandNames
    )
    BandsType = BandsType
except ImportError:
    from eoreader.bands import BandType

    BandsType = list | BandType


def look_for_prod_band_file(prod: Product, band: BandType, pixel_size: float, **kwargs):
    """
    Look for a product's band file

    Args:
        prod (Product): Product to look in
        band (BandType): Band to look for
        pixel_size (float): Pixel size in meters (if needed)
        **kwargs: Other args

    Returns:
        AnyPathType: Band file path
    """
    band_path = _look_for_prod_band_file(
        prod, band, pixel_size, writable=False, **kwargs
    )

    if band_path is None:
        band_path = _look_for_prod_band_file(
            prod, band, pixel_size, writable=True, **kwargs
        )

    if band_path is None:
        raise FileNotFoundError(
            f"Non-existing processed band {to_str(band)[0]} in {prod.condensed_name}!"
        )

    return band_path


def _look_for_prod_band_file(
    prod: Product, band: BandType, pixel_size: float, writable: bool, **kwargs
) -> AnyPathType:
    """
    Look for a product's band file

    Args:
        prod (Product): Product to look in
        band (BandType): Band to look for
        pixel_size (float): Pixel size in meters (if needed)
        writable (bool): Whether to force look in writable folder or not
        **kwargs: Other args

    Returns:
        AnyPathType: Band file path
    """
    band_path = None

    # Get the band name
    band_name = to_str(band)[0]

    # Spectral band case: use a dedicated function
    if is_spectral_band(band):
        band_path = prod.get_band_paths(
            [band], pixel_size, writable=writable, **kwargs
        )[band]

        # Check if the band exists in a writable directory if not existing in the default one
        if not AnyPath(band_path).is_file():
            band_path = None

    else:
        # No window in non-spectral paths (for now ?)
        with contextlib.suppress(StopIteration):
            # Check if the band exists in a non-writable directory
            band_regex = f"*{prod.condensed_name}*_{band_name}_*"
            band_path = next(prod._get_band_folder(writable=writable).glob(band_regex))

    return band_path
