# Release History

## 0.4.1 (2025-01-07)

- CI: Change CI bucket
- CI: Update test mosaic with custom product: use several bands and write stack on disk
- CI: Add weekly tests on Python 3.13 and 3.14 ([#30](https://github.com/sertit/eosets/issues/30))
- DOC: Update copyright to 2026 ([#33](https://github.com/sertit/eosets/issues/33))
- DOC: Improve documentation
- DEPS: Drop Python 3.9 support and support Python 3.13 and 3.14 ([#30](https://github.com/sertit/eosets/issues/30))
- DEPS: Update deps to align on `sertit` and `eoreader` ([#34](https://github.com/sertit/eosets/issues/34))

## 0.4.0 (2025-07-09)

- **ENH: Better manage and advertize default resolution in case of heterogeneous Set** ([#20](https://github.com/sertit/eosets/issues/20))
- **ENH: Add the `len()` property to any Set, giving the number of EOReader Products contained into the Set** 
- FIX: Fix pair creation with paths and Mosaic instead of list of paths
- FIX: Better management of `condensed_name` / `full_name` and `id`
- CI: Fix scheduled pipeline

## 0.3.3 (2025-04-08)

- FIX: Adapt the code to `eoreader>=0.22.0`
- DEPS: Update `sertit` and `eoreader`

## 0.3.2 (2025-01-22)

- FIX: Correctly delete old temporary process folder if a new output is given.
- FIX: Correctly look for band files in product's temporary directory
- FIX: Better manage pairs with no secondary product
- CI: Enhance tests, loading both an index and a spectral band

## 0.3.1 (2025-01-14)

- FIX: Fixing import in utils

## 0.3.0 (2025-01-14)

- **ENH: Allow to give `eoreader.Product` instead of paths to create any Set** ([#8](https://github.com/sertit/eosets/issues/8))
- **ENH: Drop `isort`, `black` and `flake8` and use `ruff`**
- **ENH: Use `pyproject.toml` instead of `setup.py`**
- FIX: Don't fail in case of a string is given as Mosaic path. 
- FIX: Update and refactor types 

## 0.2.5 (2024-10-21)

- FIX: Fix retrieval of `is_optical` and `is_sar` Mosaic members
- FIX: Fix changes looked for to run CI

## 0.2.4 (2024-10-18)

- FIX: Fix band retrieving when the env variable `CI_EOREADER_BAND_FOLDER` is set, in case of multiple files of the same band from different satellite data are present in the directory
- FIX: Don't force `remove_tmp` to `True` for `eoreader.Product` in `Mosaic`.
- FIX: Don't try to mosaic bands of mono-mosaics
- CI: Add more tests, speed up and refactor

## 0.2.3 (2024-10-16)

- FIX: Fix band retrieving when the env variable `CI_EOREADER_BAND_FOLDER` is set

## 0.2.2 (2024-10-08)

- ENH: Add `is_sar` and `is_optical` member to any `Set`
- FIX: Fix Mosaic stacking with multi-resolution constellations ([#6](https://github.com/sertit/eosets/issues/6))
- CI: Test with S3-stored products
- PUBLISH: Use PyPI's Trusted Publisher Management mechanism

## 0.2.1 (2024-04-25)

- FIX: Fix loading difference of bands with a Pair without a reference band
- FIX: Remove deprecation warnings from other libs (`sertit`, `eoreader`)
- CI: rename `CI` in `ci`
- DOC: Fix caching notebooks in readthedocs

## 0.2.0 (2024-04-24)

- **BREAKING CHANGE**: Switching from `pivot`/`child` to `reference`/`secondary`
- **BREAKING CHANGE**: Switching from `ruling` to `reference`
- FIX: Add a `d` before the band in the difference between reference and secondary saved on disk (i.e. `20200824_S2_20200908_S2_dBAI` for difference of BAI, instead of `20200824_S2_20200908_S2_BAI`)
- FIX: With mosaics, fix paths retrieval for similar bands (i.e. BAI and BAIS2)
- FIX: With VRT mosaics, move EOReader's bands and copy raw data bands to mosaic output folder (instead of always moving)
- FIX: Fix Pair without secondary
- FIX: Mosaic - Don't try to load sth if bands_to_load is empty 
- CI: Test Pair without secondary and Mosaic and Series with only one product
- CI: Update pre-commit hooks
- CI: Enabling pre-commit.ci and dependabot bots
- DOC: Update doc and notebooks

## 0.1.0 (2023-07-19)

- First release. Going opensource 🚀