# """ Testing pair """
import os
import tempfile

import pytest
from eoreader.bands import RED
from sertit import ci

from CI.scripts_utils import compare_geom, data_folder, pair_folder
from eosets.exceptions import IncompatibleProducts
from eosets.pair import Pair

ci.reduce_verbosity()

ON_DISK = False


def _test_pair_core(paths: dict) -> None:
    """
    Test pair core

    Args:
        paths (dict): Pair paths
    """

    aoi_path = data_folder() / "Fire_Spain.geojson"

    with tempfile.TemporaryDirectory() as output:
        if ON_DISK:
            output = r"/mnt/ds2_db3/CI/eosets/PAIR"

        # Create object
        pair = Pair(**paths)
        pair.output = os.path.join(output, pair.condensed_name)

        # Check extent
        compare_geom("extent", pair, pair_folder(), ON_DISK)

        # Check footprint
        compare_geom("footprint", pair, pair_folder(), ON_DISK)

        # TODO: check with input mosaic, check child-pivot

        # Stack with a pixel_size of 60m
        pair_out = pair.output / "red_stack.tif"
        assert pair.has_bands(RED)
        pair.stack(
            pivot_bands=RED,
            child_bands=RED,
            diff_bands=RED,
            window=aoi_path,
            pixel_size=60,
            stack_path=pair_out,
        )

        # Test it
        if ON_DISK:
            ci_path = pair_out
        else:
            ci_path = pair_folder() / pair.condensed_name / "red_stack.tif"

        ci.assert_raster_equal(pair_out, ci_path)

        # Not implemented
        with pytest.raises(NotImplementedError):
            pair.read_mtd()

        # Clean everything
        pair.clear()
        pair.clean_tmp()


def test_s2_pair():
    """Test pair object with Sentinel-2 products"""
    s2_paths = {
        "pivot_paths": [
            data_folder()
            / "S2A_MSIL1C_20200824T110631_N0209_R137_T29TQE_20200824T150432.zip"
        ],
        "child_paths": [
            data_folder()
            / "S2B_MSIL1C_20200908T110619_N0209_R137_T29TQE_20200908T132324.zip"
        ],
    }
    _test_pair_core(s2_paths)


def test_s3_pair():
    """Test pair object with Sentinel-3 products"""
    s3_paths = {
        "pivot_paths": [
            data_folder()
            / "S3B_SL_1_RBT____20200824T105515_20200824T105815_20200825T151744_0179_042_322_2340_LN2_O_NT_004.SEN3"
        ],
        "child_paths": [
            data_folder()
            / "S3B_SL_1_RBT____20200909T104016_20200909T104316_20200910T161910_0179_043_165_2340_LN2_O_NT_004.SEN3"
        ],
    }
    _test_pair_core(s3_paths)


def test_l8_pair():
    """Test pair object with Landsat-8 products"""
    l8_paths = {
        "pivot_paths": [data_folder() / "LC08_L1TP_202032_20200828_20200906_02_T1.tar"],
        "child_paths": [data_folder() / "LC08_L1TP_202032_20200929_20201006_02_T1.tar"],
    }
    _test_pair_core(l8_paths)


def test_pair_fail():
    """Test failure for pair objects"""
    paths = {
        "pivot_paths": [
            data_folder()
            / "S2A_MSIL1C_20200824T110631_N0209_R137_T29TQE_20200824T150432.zip"
        ],
        "child_paths": [
            data_folder()
            / "S2B_MSIL2A_20220228T102849_N0400_R108_T32ULU_20220228T134712.SAFE"
        ],
    }

    # Fails with not overlapping products
    with pytest.raises(IncompatibleProducts):
        Pair(**paths)
