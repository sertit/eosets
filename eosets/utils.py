""" Utils file """
from pathlib import Path
from typing import Union

from cloudpathlib import CloudPath
from eoreader import utils

AnyPathType = Union[CloudPath, Path]

read = utils.read
write = utils.write
stack_dict = utils.stack_dict
