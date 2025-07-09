# """ Testing series """
import os

import pytest
from eoreader.bands import RED
from eoreader.env_vars import CI_EOREADER_BAND_FOLDER
from eoreader.reader import Reader
from sertit import ci
from tempenv import tempenv

from ci.scripts_utils import (
    compare_geom,
    data_folder,
    get_ci_data_dir,
    get_output,
    s3_env,
    series_folder,
)
from eosets import Series
from eosets.exceptions import IncompatibleProducts

ci.reduce_verbosity()

ON_DISK = False


def get_ci_series_data_dir():
    return str(get_ci_data_dir() / "SERIES")


@s3_env
def test_s2_series(tmp_path):
    """Test series object with Sentinel-2 products"""
    s2_paths = [
        [
            data_folder()
            / "S2A_MSIL1C_20200824T110631_N0209_R137_T29TQE_20200824T150432.SAFE"
        ],
        [
            data_folder()
            / "S2B_MSIL1C_20200908T110619_N0209_R137_T29TQE_20200908T132324.SAFE"
        ],
    ]
    aoi_path = data_folder() / "Fire_Spain.geojson"

    with tempenv.TemporaryEnvironment(
        {CI_EOREADER_BAND_FOLDER: get_ci_series_data_dir()}
    ):
        output = get_output(tmp_path, "SERIES", ON_DISK)

        # Try with same datetime
        with pytest.raises(IncompatibleProducts):
            Series([s2_paths[0], s2_paths[0]])

        # Try with non-overlapping products
        with pytest.raises(IncompatibleProducts):
            Series(
                [
                    s2_paths[0],
                    [
                        data_folder()
                        / "S2B_MSIL2A_20220228T102849_N0400_R108_T32ULU_20220228T134712.SAFE"
                    ],
                ]
            )

        # Create object
        series = Series(paths=s2_paths, remove_tmp=not ON_DISK)
        series.output = os.path.join(output, series.condensed_name)

        # Check extent
        compare_geom("extent", series, series_folder(), ON_DISK)

        # Check footprint
        compare_geom("footprint", series, series_folder(), ON_DISK)

        # TODO: check with input mosaic, different reference mosaic

        # Stack with a pixel_size of 60m
        series_path = series.output / "red_stack.tif"
        assert series.has_bands(RED)
        series.stack(RED, window=aoi_path, pixel_size=60, stack_path=series_path)

        # Test it
        if ON_DISK:
            ci_path = series_path
        else:
            ci_path = series_folder() / series.condensed_name / "red_stack.tif"

        ci.assert_raster_equal(series_path, ci_path)

        # Not implemented
        with pytest.raises(NotImplementedError):
            series.read_mtd()

        # Clean everything
        if not ON_DISK:
            series.clear()
            series.clean_tmp()


@s3_env
def test_mono_series(tmp_path):
    """Test series object with Sentinel-2 products (only one mosaic)"""

    with tempenv.TemporaryEnvironment(
        {CI_EOREADER_BAND_FOLDER: get_ci_series_data_dir()}
    ):
        output = get_output(tmp_path, "SERIES", ON_DISK)
        s2_paths = [
            [
                data_folder()
                / "S2A_MSIL1C_20200824T110631_N0209_R137_T29TQE_20200824T150432.SAFE"
            ],
        ]
        aoi_path = data_folder() / "Fire_Spain.geojson"

        # Create object
        series = Series(paths=s2_paths, remove_tmp=not ON_DISK)
        series.output = os.path.join(output, series.condensed_name)

        # Stack with a pixel_size of 60m
        series.stack(RED, window=aoi_path, pixel_size=60)

        # Just see if this doesn't fail


@s3_env
def test_series_from_custom_prod(tmp_path):
    with tempenv.TemporaryEnvironment(
        {CI_EOREADER_BAND_FOLDER: get_ci_series_data_dir()}
    ):
        output = get_output(tmp_path, "SERIES", ON_DISK)

        # Get a custom stack path
        pld_psh_path = data_folder() / "pld_psh.tif"

        # Create object
        pld_psh = Reader().open(
            pld_psh_path,
            custom=True,
            datetime="20250117T030415",
            sensor_type="OPTICAL",
            band_map={"BLUE": 1, "GREEN": 2, "RED": 3, "NIR": 4},
        )
        pld_psh_2 = Reader().open(
            pld_psh_path,
            custom=True,
            datetime="20250118T030415",
            sensor_type="OPTICAL",
            band_map={"BLUE": 1, "GREEN": 2, "RED": 3, "NIR": 4},
        )
        series = Series([pld_psh, pld_psh_2], remove_tmp=not ON_DISK)
        series.output = os.path.join(output, series.condensed_name)
        series.stack(
            ["NDVI"],
            pixel_size=60,
        )

        ci.assert_val(len(series), 2, "Number of products")


@s3_env
@pytest.mark.parametrize(
    ("default_pixel_size", "pixel_size"),
    [
        pytest.param(None, 30),
        pytest.param("coarsest", 30),
        pytest.param("most_resolute", 10),
        pytest.param(50, 50),
    ],
)
def test_default_pixel_size(default_pixel_size, pixel_size, tmp_path):
    output = get_output(tmp_path, "SERIES", ON_DISK)
    paths = [
        [
            data_folder()
            / "S2A_MSIL1C_20200824T110631_N0209_R137_T29TQE_20200824T150432.SAFE",
        ],
        [
            data_folder() / "LC08_L1TP_202032_20200929_20201006_02_T1",
        ],
    ]
    aoi_path = data_folder() / "Fire_Spain.geojson"

    # Create series by default
    series = Series(
        paths=paths,
        remove_tmp=not ON_DISK,
        default_pixel_size=default_pixel_size,
        overlap_check="none",
        contiguity_check="none",
    )
    series.output = os.path.join(output, series.condensed_name)
    red = series.load(RED, window=aoi_path)[RED]
    ci.assert_val(series.default_pixel_size, pixel_size, "Default pixel size")
    ci.assert_val(round(red.rio.resolution()[0]), pixel_size, "Band pixel size")
