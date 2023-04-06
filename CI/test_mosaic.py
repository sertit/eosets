""" Testing pairs """
import os
import tempfile

import pytest
from eoreader.bands import CLOUDS, MNDWI, RED, SLOPE
from eoreader.env_vars import DEM_PATH
from sertit import ci

from CI.scripts_utils import data_path, get_db_dir
from eosets.exceptions import IncompatibleProducts
from eosets.mosaic import Mosaic

ci.reduce_verbosity()

ON_DISK = False


def test_s2_mosaic():
    dem_sub_dir_path = ["GLOBAL", "COPDEM_30m", "COPDEM_30m.vrt"]
    os.environ[DEM_PATH] = str(get_db_dir().joinpath(*dem_sub_dir_path))

    if ON_DISK:
        from cloudpathlib import AnyPath

        mosaic_path = AnyPath("D:\_EXTRACTEO\DS3\CI\eosets\MOSAIC")
    else:
        mosaic_path = data_path() / "MOSAIC"
    # Get some Sentinel-2 paths
    s2_32umu = (
        mosaic_path
        / "S2B_MSIL2A_20220330T102619_N0400_R108_T32UMU_20220330T141833.SAFE"
    )
    s2_32ulu = (
        mosaic_path
        / "S2B_MSIL2A_20220228T102849_N0400_R108_T32ULU_20220228T134712.SAFE"
    )
    s2_32ulv = (
        mosaic_path
        / "S2B_MSIL2A_20220228T102849_N0400_R108_T32ULV_20220228T134712.SAFE"
    )

    with tempfile.TemporaryDirectory() as output:
        if ON_DISK:
            output = r"D:\_EXTRACTEO\OUTPUT\eosets\Mosaic"

        # First try with incompatible products
        with pytest.raises(IncompatibleProducts):
            mosaic = Mosaic(
                [s2_32umu, s2_32ulu], output_path=output, mosaic_method="VRT"
            )

        # Then with compatible
        mosaic = Mosaic([s2_32ulv, s2_32ulu], mosaic_method="VRT")
        mosaic.output = os.path.join(output, mosaic.condensed_name)

        # Stack with a pixel_size of 60m
        mosaic.stack(
            [MNDWI, RED, CLOUDS, SLOPE],
            stack_path=mosaic._output / "water_stack.tif",
            pixel_size=60,
        )
