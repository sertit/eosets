# Release History

## 0.3.0 (2025-01-14)

- **ENH: Allow to give `eoreader.Product` instead of paths to create any Set** ([#8](https://github.com/sertit/eoreader/issues/8))
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
- FIX: Fix Mosaic stacking with multi-resolution constellations ([#6](https://github.com/sertit/eoreader/issues/6))
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