# Release History


## 0.2.0 (2023-MM-DD)

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

## 0.1.0 (2023-07-19)

- First release. Going opensource 🚀