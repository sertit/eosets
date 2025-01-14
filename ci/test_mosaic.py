"""Testing mosaics"""

import os

import pytest
from eoreader.bands import NBR, NDVI, RED
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
    mosaic_folder,
    s3_env,
)
from eosets import Mosaic
from eosets.exceptions import IncompatibleProducts

ci.reduce_verbosity()

ON_DISK = False


def get_ci_mosaic_data_dir():
    return str(get_ci_data_dir() / "MOSAIC")


@s3_env
def test_s2_mosaic(tmp_path):
    """Test mosaic object with Sentinel-2 products"""
    with tempenv.TemporaryEnvironment(
        {DEM_PATH: get_copdem_30(), CI_EOREADER_BAND_FOLDER: get_ci_mosaic_data_dir()}
    ):
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

        output = get_output(tmp_path, "MOSAIC", ON_DISK)

        # First try with incompatible products
        with pytest.raises(IncompatibleProducts):
            mosaic = Mosaic(
                [s2_32umu, s2_32ulu], output_path=output, mosaic_method="VRT"
            )

        # Create object
        mosaic = Mosaic(
            [s2_32ulv, s2_32ulu], mosaic_method="VRT", remove_tmp=not ON_DISK
        )
        mosaic.output = os.path.join(output, mosaic.condensed_name)

        # Some checks
        # TODO: add more
        ci.assert_val(mosaic.is_optical, True, "Is Optical?")
        ci.assert_val(mosaic.is_sar, False, "Is SAR?")

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


@s3_env
def test_mono_mosaic(tmp_path):
    """Test mosaic object with Sentinel-2 products (only one product)"""

    with tempenv.TemporaryEnvironment(
        {CI_EOREADER_BAND_FOLDER: get_ci_mosaic_data_dir()}
    ):
        output = get_output(tmp_path, "MOSAIC", ON_DISK)

        # Get some Sentinel-2 paths
        s2_32umu = (
            data_folder()
            / "S2B_MSIL2A_20220330T102619_N0400_R108_T32UMU_20220330T141833.SAFE"
        )

        # Create object
        mosaic = Mosaic([s2_32umu], mosaic_method="VRT", remove_tmp=not ON_DISK)
        mosaic.output = os.path.join(output, mosaic.condensed_name)
        mosaic.stack(
            [NDVI],
            pixel_size=600,
        )

        # Some checks
        # TODO: add more
        ci.assert_val(mosaic.is_optical, True, "Is Optical?")
        ci.assert_val(mosaic.is_sar, False, "Is SAR?")

        # Just see if this doesn't fail


def test_mosaic_from_custom_prod(tmp_path):
    with tempenv.TemporaryEnvironment(
        {CI_EOREADER_BAND_FOLDER: get_ci_mosaic_data_dir()}
    ):
        output = get_output(tmp_path, "MOSAIC", ON_DISK)

        # Get a custom stack path
        pld_psh_path = data_folder() / "pld_psh.tif"

        # Create object
        pld_psh = Reader().open(
            pld_psh_path,
            custom=True,
            sensor_type="OPTICAL",
            band_map={"BLUE": 1, "GREEN": 2, "RED": 3, "NIR": 4},
        )
        mosaic = Mosaic(pld_psh, mosaic_method="VRT", remove_tmp=not ON_DISK)
        mosaic.output = os.path.join(output, mosaic.condensed_name)
        mosaic.stack(
            [NDVI],
            pixel_size=60,
        )

        # Some checks
        # TODO: add more
        ci.assert_val(mosaic.is_optical, True, "Is Optical?")
        ci.assert_val(mosaic.is_sar, False, "Is SAR?")


def test_ci_eoreader_band_folder(tmp_path):
    """Test mosaic with CI_EOREADER_BAND_FOLDER set to an arbitrary diretcory."""
    with tempenv.TemporaryEnvironment(
        {CI_EOREADER_BAND_FOLDER: get_ci_mosaic_data_dir()}
    ):
        output = get_output(tmp_path, "MOSAIC", ON_DISK)

        # Get some Sentinel-2 paths
        s2_32umu = (
            data_folder()
            / "S2B_MSIL2A_20220330T102619_N0400_R108_T32UMU_20220330T141833.SAFE"
        )

        # Create object
        mosaic = Mosaic([s2_32umu], mosaic_method="VRT", remove_tmp=not ON_DISK)
        mosaic.output = os.path.join(output, mosaic.condensed_name)
        mosaic.stack(
            [NBR],
            pixel_size=600,
        )

        # Some checks
        # TODO: add more
        ci.assert_val(mosaic.is_optical, True, "Is Optical?")
        ci.assert_val(mosaic.is_sar, False, "Is SAR?")


# TODO: Add tests for SAR mosaics
