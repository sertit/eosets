""" Testing mosaics """
import os
import tempfile

import pytest
from eoreader.bands import RED
from eoreader.env_vars import DEM_PATH
from sertit import ci

from CI.scripts_utils import compare_geom, data_folder, get_db_dir, mosaic_folder
from eosets import Mosaic
from eosets.exceptions import IncompatibleProducts

ci.reduce_verbosity()

ON_DISK = False


def test_s2_mosaic():
    """Test mosaic object with Sentinel-2 products"""
    dem_sub_dir_path = ["GLOBAL", "COPDEM_30m", "COPDEM_30m.vrt"]
    os.environ[DEM_PATH] = str(get_db_dir().joinpath(*dem_sub_dir_path))

    # Get some Sentinel-2 paths
    s2_32umu = (
        data_folder()
        / "S2B_MSIL2A_20220330T102619_N0400_R108_T32UMU_20220330T141833.SAFE"
    )
    s2_32ulu = (
        data_folder()
        / "S2B_MSIL2A_20220228T102849_N0400_R108_T32ULU_20220228T134712.SAFE"
    )
    s2_32ulv = (
        data_folder()
        / "S2B_MSIL2A_20220228T102849_N0400_R108_T32ULV_20220228T134712.SAFE"
    )

    with tempfile.TemporaryDirectory() as output:
        if ON_DISK:
            output = r"/mnt/ds2_db3/CI/eosets/MOSAIC"

        # First try with incompatible products
        with pytest.raises(IncompatibleProducts):
            mosaic = Mosaic(
                [s2_32umu, s2_32ulu], output_path=output, mosaic_method="VRT"
            )

        # Create object
        mosaic = Mosaic([s2_32ulv, s2_32ulu], mosaic_method="VRT")
        mosaic.output = os.path.join(output, mosaic.condensed_name)

        # Check extent
        compare_geom("extent", mosaic, mosaic_folder(), ON_DISK)

        # Check footprint
        compare_geom("footprint", mosaic, mosaic_folder(), ON_DISK)

        # Stack with a pixel_size of 600m
        mosaic_out = mosaic.output / "red_stack.tif"
        assert mosaic.has_bands(RED)
        mosaic.stack(
            [RED],
            stack_path=mosaic_out,
            pixel_size=600,
        )

        # Test it
        if ON_DISK:
            ci_path = mosaic_out
        else:
            ci_path = mosaic_folder() / mosaic.condensed_name / mosaic_out.name

        ci.assert_raster_equal(mosaic_out, ci_path)

        # Not implemented
        with pytest.raises(NotImplementedError):
            mosaic.read_mtd()

        # Clean everything
        mosaic.clear()
        mosaic.clean_tmp()
