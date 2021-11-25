# Changelog

### [0.3.1](https://www.github.com/RIVM-bioinformatics/AmpliGone/compare/v0.3.0...v0.3.1) (2021-11-25)


### Bug Fixes

* Distinguish (mis)match and cut soft clips ([d4132c1](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/d4132c1d648b9a2cd5c4331dfc6998c2a99da376))


### Performance Improvements

* Add minimap scoring to presets ([11158e5](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/11158e5184ea3761f15bf7bcb3594237a5f4049d))

## [0.3.0](https://www.github.com/RIVM-bioinformatics/AmpliGone/compare/v0.2.3...v0.3.0) (2021-11-17)


### Features

* Use CIGAR strings to cut end-to-end ([ebd2684](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/ebd2684eedcd4db6649f19b3c0c92abbe020e9b9))


### Bug Fixes

* Fix typo in 75e792fa8d69746 ([eca004b](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/eca004b4945305b83169826ea21c97d7fd5d03a4))
* only return data that got an alignment-hit ([f9e64e8](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/f9e64e8f4c1f939dcae7fcccb552a5d4686b70e5))


### Performance Improvements

* Improved (fastq)data loading speed ([9d1ae23](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/9d1ae2306d11f92a326130ad3b1df87243dd033f))


### Documentation

* Fix comment in cut_read function ([b920738](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/b9207388e20675a90f8c13868f5bb63414af8de6))
* update docs to match new primer removal process ([035be92](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/035be9299384dd3e0d225ef06f0618744132ed17))
* update readme ([95ff984](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/95ff9840550e25480b7f7d35cabcd070e5428d32))

### [0.2.3](https://www.github.com/RIVM-bioinformatics/AmpliGone/compare/v0.2.2...v0.2.3) (2021-10-28)


### Bug Fixes

* set correct exit code when `-to` flag is given ([5a93ff0](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/5a93ff00ac3c4d19f635a8564efd21a5e9166e01))

### [0.2.2](https://www.github.com/RIVM-bioinformatics/AmpliGone/compare/v0.2.1...v0.2.2) (2021-10-27)


### Bug Fixes

* add missing primer-orientation keywords ([41a84e6](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/41a84e6688432cf195c6c777986d57473c89dd05))
* exit cleanly when primers aren't found ([#11](https://www.github.com/RIVM-bioinformatics/AmpliGone/issues/11)) ([e927591](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/e9275912ffb8a181afb54bb890eb5e4e39b2c28b))
* Exit cleanly with empty input ([#12](https://www.github.com/RIVM-bioinformatics/AmpliGone/issues/12)) ([35aae65](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/35aae655d8833afececd857a9e2557b4c8b9a98d))


### Documentation

* update user-guide to include newly added flag ([4d28ac1](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/4d28ac180e97ae0088a891d5704da0fe8ad35ebe))

### [0.2.1](https://www.github.com/RIVM-bioinformatics/AmpliGone/compare/v0.2.0...v0.2.1) (2021-10-25)


### Documentation

* add changelog link ([4b62a92](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/4b62a92abd245e286a02de9e278f88f6cbd66589))
