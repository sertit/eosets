""" Utils module for scripts """
import logging
import os
from pathlib import Path
from typing import Union

from cloudpathlib import AnyPath, CloudPath
from eoreader.reader import Reader
from sertit import ci

from eopairs.utils import EOPAIRS_NAME

LOGGER = logging.getLogger(EOPAIRS_NAME)
READER = Reader()

CI_PAIRS_S3 = "CI_PAIRS_USE_S3"


def get_ci_dir() -> Union[CloudPath, Path]:
    """
    Get CI DATA directory
    Returns:
        str: CI DATA directory
    """
    return AnyPath(__file__).parent.parent


def get_ci_db_dir() -> Union[CloudPath, Path]:
    """
    Get CI database directory (S3 bucket)
    Returns:
        str: CI database directory
    """
    if int(os.getenv(CI_PAIRS_S3, 0)):
        # ON S3
        ci.define_s3_client()
        return AnyPath("s3://sertit-eopairs-ci")
    else:
        # ON DISK
        try:
            # CI
            return AnyPath(ci.get_db3_path(), "CI", "eopairs")
        except NotADirectoryError:
            # Windows
            path = AnyPath(r"//ds2/database03/CI/eopairs")
            if not path.is_dir():
                raise NotADirectoryError("Impossible to find get_ci_db_dir")

            return path


def get_ci_data_dir() -> Union[CloudPath, Path]:
    """
    Get CI DATA directory (S3 bucket)
    Returns:
        str: CI DATA directory
    """
    if len(os.getenv(ci.AWS_ACCESS_KEY_ID, "")) > 0:
        return get_ci_db_dir().joinpath("DATA")
    else:
        return get_ci_dir().joinpath("DATA")


def data_path():
    return get_ci_db_dir().joinpath("DATA")
