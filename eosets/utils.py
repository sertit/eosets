""" Utils file """
from pathlib import Path
from typing import Union

from cloudpathlib import CloudPath
from eoreader import utils

EOPAIRS_NAME = "eosets"

AnyPathType = Union[CloudPath, Path]

read = utils.read
write = utils.write
