# """ Testing series """
import os
import tempfile

import pytest
from eoreader.bands import RED
from sertit import ci

from CI.scripts_utils import compare_geom, data_folder, series_folder
from eosets import Series
from eosets.exceptions import IncompatibleProducts

ci.reduce_verbosity()

ON_DISK = False


def test_s2_series():
    """Test series object with Sentinel-2 products"""
    s2_paths = [
        [
            data_folder()
            / "S2A_MSIL1C_20200824T110631_N0209_R137_T29TQE_20200824T150432.zip"
        ],
        [
            data_folder()
            / "S2B_MSIL1C_20200908T110619_N0209_R137_T29TQE_20200908T132324.zip"
        ],
    ]
    aoi_path = data_folder() / "Fire_Spain.geojson"

    with tempfile.TemporaryDirectory() as output:
        if ON_DISK:
            output = r"/mnt/ds2_db3/CI/eosets/SERIES"

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
        series = Series(paths=s2_paths)
        series.output = os.path.join(output, series.condensed_name)

        # Check extent
        compare_geom("extent", series, series_folder(), ON_DISK)

        # Check footprint
        compare_geom("footprint", series, series_folder(), ON_DISK)

        # TODO: check with input mosaic, different ruling mosaic

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
        series.clear()
        series.clean_tmp()
