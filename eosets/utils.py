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

import logging
import os

from eoreader import utils
from eoreader.bands import DEM, HILLSHADE, SLOPE, BandType, is_dem, to_str
from eoreader.env_vars import DEM_PATH
from eoreader.keywords import DEM_KW, HILLSHADE_KW, SLOPE_KW
from eoreader.products import Product
from sertit import path
from sertit.types import AnyPathStrType

from eosets import EOSETS_NAME

read = utils.read
write = utils.write
stack = utils.stack
convert_to_uint16 = utils.convert_to_uint16
write_path_in_attrs = utils.write_path_in_attrs

AnyProductType = AnyPathStrType | Product
""" Any Product Type, either a path or an eoreader.Product"""

LOGGER = logging.getLogger(EOSETS_NAME)


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
    dem_name = ""
    if is_dem(band):
        # We already checked if it exists when loading the DEM band
        dem_path = os.environ.get(DEM_PATH)
        if band == DEM:
            dem_path = kwargs.get(DEM_KW, dem_path)
        elif band == SLOPE:
            dem_path = kwargs.get(SLOPE_KW, dem_path)
        elif band == HILLSHADE:
            dem_path = kwargs.get(HILLSHADE_KW, dem_path)

        dem_name = path.get_filename(dem_path)
    band_path = prod.get_band_path(
        band, pixel_size, dem_name=dem_name, writable=False, **kwargs
    )

    if not band_path.exists():
        band_path = prod.get_band_path(
            band, pixel_size, dem_name=dem_name, writable=True, **kwargs
        )

    if not band_path.exists():
        raise FileNotFoundError(
            f"Non-existing processed band {to_str(band)[0]} in {prod.condensed_name}!"
        )

    return band_path
