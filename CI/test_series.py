# """ Testing series """
import os
import tempfile

from eoreader.bands import RED
from sertit import ci

from CI.scripts_utils import data_folder, series_folder
from eosets import Series

ci.reduce_verbosity()

ON_DISK = False


def test_s2_series():
    import rasterio

    print(rasterio.__version__)

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

        series = Series(paths=s2_paths)
        series.output = os.path.join(output, series.condensed_name)

        # %%
        series_path = series.output / "red_stack.tif"
        series.stack(RED, window=aoi_path, pixel_size=60, stack_path=series_path)

        # Test it
        if ON_DISK:
            ci_path = series_path
        else:
            ci_path = series_folder() / series.condensed_name / "red_stack.tif"

        ci.assert_raster_equal(series_path, ci_path)
