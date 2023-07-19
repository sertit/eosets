import os

import setuptools

from eosets.__meta__ import (
    __author__,
    __author_email__,
    __description__,
    __documentation__,
    __title__,
    __url__,
    __version__,
)

BASEDIR = os.path.dirname(os.path.abspath(os.path.realpath(__file__)))
with open(os.path.join(BASEDIR, "README.md"), "r", encoding="utf8") as f:
    readme = f.read()

setuptools.setup(
    name=__title__,
    version=__version__,
    author=__author__,
    author_email=__author_email__,
    description=__description__,
    long_description=readme,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    install_requires=[
        "lxml",
        "h5netcdf",
        "scipy",
        "rasterio>=1.3.0",
        "xarray>=0.18.0",
        "rioxarray>=0.10.0",
        "geopandas>=0.11.0",
        "spyndex>=0.3.0",
        "pyresample",
        "zarr",
        "rtree",
        "validators",
        "methodtools",
        "dicttoxml",
        "eoreader>0.20.0",
        "sertit[full]>=1.27.3",
    ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: GIS",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    python_requires=">=3.9",
    project_urls={
        "Bug Tracker": f"{__url__}/issues/",
        "Documentation": __documentation__,
        "Source Code": __url__,
    },
)
