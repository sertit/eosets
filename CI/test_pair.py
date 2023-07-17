# """ Testing pair """
import os
import tempfile

from eoreader.bands import RED
from sertit import ci

from CI.scripts_utils import data_folder, pair_folder
from eosets.pair import Pair

ci.reduce_verbosity()

ON_DISK = False


def _test_core(paths):

    aoi_path = data_folder() / "Fire_Spain.geojson"

    with tempfile.TemporaryDirectory() as output:
        if ON_DISK:
            output = r"/mnt/ds2_db3/CI/eosets/PAIR"

        pair = Pair(**paths)
        pair.output = os.path.join(output, pair.condensed_name)

        # Stack
        pair_out = pair.output / "red_stack.tif"
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


def test_s2_pair():
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
    _test_core(s2_paths)


def test_s3_pair():
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
    _test_core(s3_paths)


def test_l8_pair():
    l8_paths = {
        "pivot_paths": [data_folder() / "LC08_L1TP_202032_20200828_20200906_02_T1.tar"],
        "child_paths": [data_folder() / "LC08_L1TP_202032_20200929_20201006_02_T1.tar"],
    }
    _test_core(l8_paths)
