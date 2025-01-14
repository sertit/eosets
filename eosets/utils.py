# Copyright 2025, SERTIT-ICube - France, https://sertit.unistra.fr/
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

from eoreader import utils
from eoreader.products import Product
from sertit.types import AnyPathStrType

read = utils.read
write = utils.write
stack = utils.stack

# Bands Type
try:
    from eoreader.bands import BandsType

    BandsType = BandsType
except ImportError:
    from typing import Union

    from eoreader.bands import BandType

    BandsType = Union[list, BandType]

AnyProductType = Union[AnyPathStrType, Product]
