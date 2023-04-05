""" Utils module for scripts """
import logging
import os

from cloudpathlib import AnyPath
from eoreader.reader import Reader
from sertit import ci

from eosets.utils import EOPAIRS_NAME, AnyPathType

LOGGER = logging.getLogger(EOPAIRS_NAME)
READER = Reader()

CI_EOSETS_S3 = "CI_EOSETS_USE_S3"


def get_ci_dir() -> AnyPathType:
    """
    Get CI DATA directory
    Returns:
        AnyPathType: CI DATA directory
    """
    return AnyPath(__file__).parent.parent


def get_ci_db_dir() -> AnyPathType:
    """
    Get CI database directory (S3 bucket)
    Returns:
        AnyPathType: CI database directory
    """
    if int(os.getenv(CI_EOSETS_S3, 0)):
        # ON S3
        ci.define_s3_client()
        return AnyPath("s3://sertit-eosets-ci")
    else:
        # ON DISK
        try:
            # CI
            return AnyPath(ci.get_db3_path(), "CI", "eosets")
        except NotADirectoryError:
            # Windows
            path = AnyPath(r"//ds2/database03/CI/eosets")
            if not path.is_dir():
                raise NotADirectoryError("Impossible to find get_ci_db_dir")

            return path


def get_ci_data_dir() -> AnyPathType:
    """
    Get CI DATA directory (S3 bucket)
    Returns:
        AnyPathType: CI DATA directory
    """
    if len(os.getenv(ci.AWS_ACCESS_KEY_ID, "")) > 0:
        return get_ci_db_dir().joinpath("DATA")
    else:
        return get_ci_dir().joinpath("DATA")


def data_path():
    return get_ci_db_dir().joinpath("DATA")


def get_db_dir_on_disk() -> AnyPathType:
    """
    Get database directory in the DS2

    Returns:
        AnyPathType: Database directory
    """
    # ON DISK
    db_dir = AnyPath(r"//ds2/database02/BASES_DE_DONNEES")

    if not db_dir.is_dir():
        try:
            db_dir = AnyPath(ci.get_db2_path(), "BASES_DE_DONNEES")
        except NotADirectoryError:
            db_dir = AnyPath("/home", "ds2_db2", "BASES_DE_DONNEES")

    if not db_dir.is_dir():
        raise NotADirectoryError("Impossible to open database directory !")

    return db_dir


def get_db_dir() -> AnyPathType:
    """
    Get database directory in the DS2

    Returns:
        AnyPathType: Database directory
    """

    if int(os.getenv(CI_EOSETS_S3, 0)):
        # ON S3
        ci.define_s3_client()
        return AnyPath("s3://sertit-geodatastore")
    else:
        # ON DISK
        db_dir = get_db_dir_on_disk()

    return db_dir
