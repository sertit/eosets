# -*- coding: utf-8 -*-
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
"""
**EOSets** library
"""

# flake8: noqa
from .__meta__ import __version__

EOSETS_NAME = "eosets"

__all__ = ["Mosaic", "Pair", "Series"]

try:
    from .mosaic import Mosaic
    from .pair import Pair
    from .series import Series
except (ModuleNotFoundError, ImportError):
    # Don't fail for check with empty environment
    pass
