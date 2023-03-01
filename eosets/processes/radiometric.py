""" Radiometric processes for Pairs """
import logging
from enum import unique
from pathlib import Path
from typing import Union

from cloudpathlib import AnyPath, CloudPath
from eoreader.bands import BLUE, GREEN, NIR, PAN, RED
from eoreader.products import LandsatInstrument, LandsatProduct, Product
from eoreader.reader import Constellation
from sertit import misc, strings
from sertit.misc import ListEnum
from sertit.rasters import MAX_CORES

from eosets.multi_pairs import MultiPairs
from eosets.utils import EOPAIRS_NAME

LOGGER = logging.getLogger(EOPAIRS_NAME)


@unique
class PshMethod(ListEnum):
    """
    Pansharpening methods
    """

    GDAL = "GDAL"
    """
    Uses GDAL's weighted Brovey algorithm. The only method usable if :code:`arcpy` is not available.
    """

    IHS = "Intensity, Hue, and Saturation"
    """
    Uses Intensity, Hue, and Saturation color space for data fusion. ⚠ Needs :code:`arcpy`!
    """

    BROVEY = "Brovey"
    """
    Uses the Brovey algorithm based on spectral modeling for data fusion. ⚠ Needs :code:`arcpy`!
    """

    SM = "SimpleMean"
    """
    Uses the averaged value between the red, green, and blue values and the panchromatic pixel value. ⚠ Needs :code:`arcpy`!
    """

    GS = "Gram-Schmidt"
    """
    Uses the Gram-Schmidt spectral-sharpening algorithm to sharpen multispectral data. ⚠ Needs :code:`arcpy`!
    """


def pansharpen(
    pairs: MultiPairs,
    method: PshMethod = PshMethod.GDAL,
    output_path: Union[str, Path, CloudPath] = None,
) -> Union[Path, CloudPath]:
    """
    Pansharpening a pair (pairs with only one child), with the same constellation.
    This process uses weights derived from ArcGis Pro.

    ⚠ Some products have the panchromatic band inside themselves and therefore the pair can be composed of only one reference product (Landsat, Superview...)

    ⚠ Only a RGBNIR stack will be pansharpen!

    ⚠ Other methods than ~PshMethod.GDAL need :code:`arcpy` to be installed!

    Args:
        pairs (MultiPairs) : Pair which will be pansharpened
        method (PshMethod): Pansharpening method (GDAL by default)
        output_path (Union[str, Path, CloudPath]): Output path, where to write the pansharpened stack

    Returns:

    """
    # Some checks
    assert len(pairs.children_prods) <= 1
    assert pairs.same_constellation

    # Manage the products to have a PAN and a MS one
    pan_prod = pairs.ref_prod

    if pairs.has_child:
        ms_prod = pairs.children_prods[0]
        if not pairs.has_unique_child:
            LOGGER.warning(
                "Multiple children found in the pair. Only considering the first one."
            )
    else:
        ms_prod = pairs.ref_prod

    # Convert method if needed
    method = PshMethod.from_value(method)

    # Create output path
    if not output_path:
        output_path = (
            pairs.output_path / f"{pairs.children_prods[0].condensed_name}_PSH.tif"
        )

    # Ensure output_path is a string here
    output_path = str(output_path)

    # Get sensor and sensor pansharpening weights (from ArcGis Pro)
    sensor, sensor_weights = _get_sensor_psh_weights(pan_prod)

    # Use arcpy (for now)
    if method != PshMethod.GDAL:
        try:
            from arcpy.sa import Pansharpen

        # If arcpy is not implemented
        except ImportError:
            raise ImportError(
                "You must install `arcpy` in your environment to use other methods than `GDAL`."
            )

        # Manage Landsat data that hasn't a stack to pansharpen and the PAN band in the same product
        if isinstance(pan_prod, LandsatProduct):
            ms_path = ms_prod._output / f"{ms_prod.condensed_name}_MS_stack.tif"
            ms_prod.stack([BLUE, GREEN, RED, NIR], stack_path=ms_path)
            pan_path = pan_prod.get_band_paths([PAN])[PAN]
        else:
            pan_path = pan_prod.get_default_band_path()
            ms_path = ms_prod.get_default_band_path()

        # Pansharpen
        pansharpen_raster = Pansharpen(
            ms_raster=ms_path,
            pan_raster=pan_path,
            fourth_band_of_ms_is_ir=True,
            weights=sensor_weights,
            type=method.value,
            sensor=sensor,
        )

        # Save output
        pansharpen_raster.save(output_path)
    else:
        # Manage Landsat data that haven't a stack to pansharpen and the PAN band in the same product
        if isinstance(pan_prod, LandsatProduct):
            pan_path = pan_prod.get_band_paths([PAN])[PAN]
            band_paths = ms_prod.get_band_paths([RED, GREEN, BLUE, NIR])
            ms_cli = [
                strings.to_cmd_string(band_paths[RED]),
                strings.to_cmd_string(band_paths[GREEN]),
                strings.to_cmd_string(band_paths[BLUE]),
                strings.to_cmd_string(band_paths[NIR]),
            ]
        else:
            pan_path = pan_prod.get_default_band_path()
            ms_path = strings.to_cmd_string(ms_prod.get_default_band_path())
            ms_cli = [
                f"{ms_path}, band={ms_prod.bands[RED].id}",
                f"{ms_path}, band={ms_prod.bands[GREEN].id}",
                f"{ms_path}, band={ms_prod.bands[BLUE].id}",
                f"{ms_path}, band={ms_prod.bands[NIR].id}",
            ]

        # BGRNIR weights
        weight_cli = ["-w".join(str(weight) for weight in sensor_weights)]

        # Run CLI
        misc.run_cli(
            [
                "gdal_pansharpen.py",
                *weight_cli,
                strings.to_cmd_string(pan_path),
                *ms_cli,
                "-threads",
                str(MAX_CORES),
                "-nodata",
                str(pairs.nodata),
            ]
        )

    return AnyPath(output_path)


def _get_sensor_psh_weights(prod: Product) -> (str, list):
    """
    Get sensor and sensor pansharpening weights, from ArcGis Pro.

    Args:
        prod (Product): Product from which to derive the weights

    Returns:
        (str, list): ArcGis Pro sensor's name and the weight list (RGBNIR)
    """
    weights = {
        # "Jilin-1": [0.166, 0.167, 0.167, 0.5],
        "Landsat 1-5 MSS": [0.166, 0.167, 0.167, 0],
        "Landsat 7 ETM+": [0.11, 0.14, 0.14, 0.61],
        "Landsat 8": [0.42, 0.51, 0.07, 0],
        "Pleiades-1": [0.9, 0.75, 0.5, 0.5],
        "Pléiades Neo": [0.9, 0.75, 0.5, 0.5],
        "SPOT 5": [0.166, 0.167, 0.167, 0],
        "SPOT 6": [0.9, 0.75, 0.5, 0.5],
        "SPOT 7": [0.9, 0.75, 0.5, 0.5],
        "Quickbird": [0.85, 0.7, 0.35, 1],
        "GeoEye-1": [0.41, 0.16, 0.13, 0.3],
        "WorldView-2": [0.39, 0.23, 0.21, 0.17],
        "WorldView-3": [0.38, 0.25, 0.2, 0.16],
        "WorldView-4": [0.38, 0.25, 0.2, 0.16],
        "SkySat-C": [0.378, 0.211, 0, 0.411],
        "SuperView-1": [0.85, 0.7, 0.35, 1],
        "Unknown": [0.166, 0.167, 0.167, 0],
    }

    sensor = None
    if prod.constellation.name in weights:
        # Works for some
        sensor = prod.constellation.name
    elif prod.constellation.name.replace("-", " ") in weights:
        # Works for some
        sensor = prod.constellation.name.replace("-", " ")
    else:
        # Won't work for Landsats
        if isinstance(prod, LandsatProduct):
            if prod.instrument == LandsatInstrument.MSS:
                sensor = "Landsat 1-5 MSS"
            elif prod.instrument == LandsatInstrument.ETM:
                sensor = "Landsat 7 ETM+"
            elif prod.instrument in [LandsatInstrument.OLI, LandsatInstrument.OLI_TIRS]:
                sensor = "Landsat 8"
        else:
            if prod.constellation == Constellation.PLD:
                sensor = "Pleiades-1"
            elif prod.constellation == Constellation.PNEO:
                sensor = "Pléiades Neo"
            elif prod.constellation == Constellation.QB:
                sensor = "Quickbird"
            elif prod.constellation == Constellation.SKY:
                sensor = "SkySat-C"

    if not sensor:
        LOGGER.warning(f"Unknown sensor: {prod.constellation.name}")
        sensor = "Unknown"

    return sensor, weights[sensor]
