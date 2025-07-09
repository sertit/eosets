# EOSets

This library aims to simplify any process working with sets of EO data handled by [EOReader](https://github.com/sertit/eoreader).

## Content

EOSets implement these objects:

![eosets_objects](docs/_static/eosets_objects.png)

NB: multi-pairs are not really implemented for now. We are waiting for true usecases. For now please use a list of pairs to process your data.

## Features

Implemented features are similar to those implemented in EOReader:
- `load`
- `stack`
- `extent`
- `footprint`
- ...

### Default pixel size

In case of heterogeneous Set, it is not straightforward to determine its default pixel size.
A keyword `default_pixel_size` has been created, set to `coarsest` by default which means that the default pixel size of a heterogeneous set is the coarsest pixel size.
This keyword can be set to either None, `coarsest`, `most_resolute` or any floating point number as a custom default resolution.


For example, in the case of a Series of Sentinel-2 and Landsat-8 data:
- Nothing specified: the default pixel size is 30 meters
- `default_pixel_size="coarsest"`: the default pixel size is 30 meters
- `default_pixel_size="most_resolute"`: the default pixel size is 10 meters
- `default_pixel_size="20"`: the default pixel size is 20 meters

## Examples

Available notebooks provided as examples:

- [Mosaic](https://eosets.readthedocs.io/latest/notebooks/mosaic.html#)
- [Pair](https://eosets.readthedocs.io/latest/notebooks/pair.html)
- [Series](https://eosets.readthedocs.io/latest/notebooks/series.html)