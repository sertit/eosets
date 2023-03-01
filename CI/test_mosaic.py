""" Testing pairs """

from eoreader.bands import MNDWI, NDWI
from sertit import ci

from CI.scripts_utils import data_path
from eosets.mosaic import Mosaic

ci.reduce_verbosity()


def test_s2_mosaic():

    s2_32umu = (
        data_path()
        / "S2B_MSIL2A_20210623T102559_N0300_R108_T32UMU_20210623T133440.SAFE.zip"
    )
    s2_32umv = (
        data_path()
        / "S2B_MSIL2A_20210623T102559_N0300_R108_T32UMV_20210623T133440.SAFE.zip"
    )

    mosaic = Mosaic(
        [s2_32umu, s2_32umv], output_path=r"D:\_EXTRACTEO\OUTPUT\EOPairs\Mosaic"
    )

    mosaic.stack([NDWI, MNDWI], stack_path=mosaic._output / "water_stack.tif")
