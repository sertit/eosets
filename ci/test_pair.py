# """ Testing pair """
import os

import pytest
from eoreader.bands import BLUE, GREEN, NIR, PAN, RED
from eoreader.env_vars import CI_EOREADER_BAND_FOLDER, DEM_PATH
from eoreader.reader import Reader
from sertit import ci
from tempenv import tempenv

from ci.scripts_utils import (
    compare_geom,
    data_folder,
    get_ci_data_dir,
    get_copdem_30,
    get_output,
    pair_folder,
    s3_env,
)
from eosets.exceptions import IncompatibleProducts
from eosets.pair import Pair

ci.reduce_verbosity()

ON_DISK = False


def get_ci_pair_data_dir():
    return str(get_ci_data_dir() / "PAIR")


@s3_env
def test_pair_non_ortho_with_window(tmp_path):
    """Test pair non-ortho native with a window"""
    # DO NOT REPROJECT BANDS --> WAY TOO SLOW
    with tempenv.TemporaryEnvironment(
        {DEM_PATH: get_copdem_30(), CI_EOREADER_BAND_FOLDER: get_ci_pair_data_dir()}
    ):
        aoi = data_folder() / "psh_pld_aoi.shp"
        pld_paths = {
            "reference_paths": [data_folder() / "IMG_PHR1A_P_001"],
            "secondary_paths": [data_folder() / "IMG_PHR1A_MS_004"],
        }

        output = get_output(tmp_path, "PAIR", ON_DISK)

        pld_pair = Pair(**pld_paths, remove_tmp=not ON_DISK)
        pld_pair.output = os.path.join(output, pld_pair.condensed_name)

        ms_path = pld_pair.output / "rgbn_stack.tif"
        pan_path = pld_pair.output / "pan_stack.tif"

        # RGBN
        rgbn_stck = pld_pair.secondary_mosaic.stack(
            [RED, GREEN, BLUE, NIR], stack_path=ms_path, window=aoi
        )
        ci.assert_val(rgbn_stck.rio.resolution()[0], 2.0, "RGBN resolution")
        ci.assert_val(rgbn_stck.rio.count, 4, "RGBN number of bands")

        # PAN
        pan_stck = pld_pair.reference_mosaic.stack(
            [PAN], stack_path=pan_path, window=aoi
        )
        ci.assert_val(pan_stck.rio.resolution()[0], 0.5, "PAN resolution")
        ci.assert_val(pan_stck.rio.count, 1, "PAN number of bands")


def _test_pair_core(paths: dict, tmp_path) -> None:
    """
    Test pair core

    Args:
        paths (dict): Pair paths
    """

    with tempenv.TemporaryEnvironment(
        {CI_EOREADER_BAND_FOLDER: get_ci_pair_data_dir()}
    ):
        aoi_path = data_folder() / "Fire_Spain.geojson"

        output = get_output(tmp_path, "PAIR", ON_DISK)

        # Create object
        pair = Pair(**paths, remove_tmp=not ON_DISK)
        pair.output = os.path.join(output, pair.condensed_name)

        # Check extent
        compare_geom("extent", pair, pair_folder(), ON_DISK)

        # Check footprint
        compare_geom("footprint", pair, pair_folder(), ON_DISK)

        # Check some properties
        assert pair.is_homogeneous

        # TODO: check with input mosaic, check secondary-reference

        # Test to see if there is an error
        pair.load(
            diff_bands=RED,
            window=aoi_path,
            pixel_size=60,
        )

        # Stack with a pixel_size of 60m
        pair_out = pair.output / "red_stack.tif"
        assert pair.has_bands(RED)
        pair.stack(
            reference_bands=RED,
            secondary_bands=RED,
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

        ci.assert_raster_almost_equal(pair_out, ci_path)

        # Not implemented
        with pytest.raises(NotImplementedError):
            pair.read_mtd()

        # Clean everything
        pair.clear()
        pair.clean_tmp()


@s3_env
def test_s2_pair(tmp_path):
    """Test pair object with Sentinel-2 products"""
    s2_paths = {
        "reference_paths": [
            data_folder()
            / "S2A_MSIL1C_20200824T110631_N0209_R137_T29TQE_20200824T150432.zip"
        ],
        "secondary_paths": [
            data_folder()
            / "S2B_MSIL1C_20200908T110619_N0209_R137_T29TQE_20200908T132324.zip"
        ],
    }
    _test_pair_core(s2_paths, tmp_path)


@s3_env
def test_s3_pair(tmp_path):
    """Test pair object with Sentinel-3 products"""
    s3_paths = {
        "reference_paths": [
            data_folder()
            / "S3B_SL_1_RBT____20200824T105515_20200824T105815_20200825T151744_0179_042_322_2340_LN2_O_NT_004.SEN3"
        ],
        "secondary_paths": [
            data_folder()
            / "S3B_SL_1_RBT____20200909T104016_20200909T104316_20200910T161910_0179_043_165_2340_LN2_O_NT_004.SEN3"
        ],
    }
    _test_pair_core(s3_paths, tmp_path)


@s3_env
def test_l8_pair(tmp_path):
    """Test pair object with Landsat-8 products"""
    l8_paths = {
        "reference_paths": [
            data_folder() / "LC08_L1TP_202032_20200828_20200906_02_T1.tar"
        ],
        "secondary_paths": [
            data_folder() / "LC08_L1TP_202032_20200929_20201006_02_T1.tar"
        ],
    }
    _test_pair_core(l8_paths, tmp_path)


@s3_env
def test_pair_no_secondary(tmp_path):
    """Test pair object with Landsat-8 products"""
    l8_paths = {
        "reference_paths": [
            data_folder() / "LC08_L1TP_202032_20200828_20200906_02_T1.tar"
        ],
    }
    _test_pair_core(l8_paths, tmp_path)


def test_pair_from_custom_prod(tmp_path):
    with tempenv.TemporaryEnvironment(
        {CI_EOREADER_BAND_FOLDER: get_ci_pair_data_dir()}
    ):
        output = get_output(tmp_path, "PAIR", ON_DISK)

        # Get a custom stack path
        pld_psh_path = data_folder() / "pld_psh.tif"

        # Create object
        pld_psh = Reader().open(
            pld_psh_path,
            custom=True,
            sensor_type="OPTICAL",
            band_map={"BLUE": 1, "GREEN": 2, "RED": 3, "NIR": 4},
        )
        pair = Pair(**{"reference_paths": pld_psh}, remove_tmp=not ON_DISK)
        pair.output = os.path.join(output, pair.condensed_name)
        pair.stack(
            ["NDVI"],
            pixel_size=60,
        )


@s3_env
def test_pair_fail():
    """Test failure for pair objects"""
    paths = {
        "reference_paths": [
            data_folder()
            / "S2A_MSIL1C_20200824T110631_N0209_R137_T29TQE_20200824T150432.zip"
        ],
        "secondary_paths": [
            data_folder()
            / "S2B_MSIL2A_20220228T102849_N0400_R108_T32ULU_20220228T134712.SAFE"
        ],
    }

    # Fails with not overlapping products
    with pytest.raises(IncompatibleProducts):
        Pair(**paths)


# @s3_env
# def test_pair_multi_res():
#     """ Test mosaic multi res with a window """
#     # DO NOT REPROJECT BANDS --> WAY TOO SLOW
#     with tempenv.TemporaryEnvironment(
#             {
#                 CI_EOREADER_BAND_FOLDER: str(get_ci_data_dir())
#             }
#     ):
#         l9_paths = {
#             "reference_paths": [
#                 data_folder() / "LC09_L1TP_200030_20220201_20220201_02_T1.tar",
#             ],
#         }
#
#         with tempfile.TemporaryDirectory() as output:
#             if ON_DISK:
#                 output = r"/mnt/ds2_db3/CI/eosets/PAIR"
#
#             l9_pair = Pair(**l9_paths, remove_tmp=not ON_DISK)
#             l9_pair.output = os.path.join(output, l9_pair.condensed_name)
