# Changelog

## [2.0.0](https://github.com/RIVM-bioinformatics/AmpliGone/compare/v1.3.1...v2.0.0) (2025-03-31)


### ⚠ BREAKING CHANGES

* change fasta2bed method from fuzzyregex to semi-global alignment based search with parasail

### Features

* .github.workflows.tests - created automatic testing workflow ([41e6222](https://github.com/RIVM-bioinformatics/AmpliGone/commit/41e6222ac1f686839a7d354805702c75f66e7b3c))
* add quiet mode, only prints WARNING and ERROR log messages ([16967e7](https://github.com/RIVM-bioinformatics/AmpliGone/commit/16967e76294f446d512137ab256944610d8ced52))
* add support for gzipped fastq output ([5b16301](https://github.com/RIVM-bioinformatics/AmpliGone/commit/5b163013cc2c1632cfab9a26e4d1d3906128c22f))
* add the --verbose flag and debugging log messages ([5b16301](https://github.com/RIVM-bioinformatics/AmpliGone/commit/5b163013cc2c1632cfab9a26e4d1d3906128c22f))
* add the 'virtual primer support' to link closely positioned primers to each other in the index for accurate removal ([3620f2f](https://github.com/RIVM-bioinformatics/AmpliGone/commit/3620f2f1b111c2c1e2da763ed657e51407f258de))
* change fasta2bed method from fuzzyregex to semi-global alignment based search with parasail ([36c4479](https://github.com/RIVM-bioinformatics/AmpliGone/commit/36c44798d5981cf33d4358159271918188f01845))
* sonar-project.properties - added branch analysis config ([41e6222](https://github.com/RIVM-bioinformatics/AmpliGone/commit/41e6222ac1f686839a7d354805702c75f66e7b3c))
* test-requirements- added seperate test reqs ([41e6222](https://github.com/RIVM-bioinformatics/AmpliGone/commit/41e6222ac1f686839a7d354805702c75f66e7b3c))


### Bug Fixes

* .gitignore - added *:Zone.Identifier ([2f24ebe](https://github.com/RIVM-bioinformatics/AmpliGone/commit/2f24ebe82e98757be251fd5f089a61a2a265d06b))
* add an explicit import for IUPACData from Bio.Data (instead of the implicit import that didn't always work depending on the interpreter) ([5abae9a](https://github.com/RIVM-bioinformatics/AmpliGone/commit/5abae9a8e0c22849158a5816b9f46bddf2511f2b))
* Add extra checks to cleanly exit when dealing with corrupted or incomplete inputfiles ([2155999](https://github.com/RIVM-bioinformatics/AmpliGone/commit/2155999c305fd41dedc0bdbf086e85cad5039229))
* add UTF-8 encoding when opening export primers file when `-to` is given ([c0bfe5f](https://github.com/RIVM-bioinformatics/AmpliGone/commit/c0bfe5fae875222e58ab020a28b4f27c85cd89ea))
* change the table of IUPAC data to actually reference the IUPACData table in fasta2bed ([daed988](https://github.com/RIVM-bioinformatics/AmpliGone/commit/daed988a5c1bef3e1785612dcdef71af6de95cdc))
* disable the progress bar when in debug mode (verbose) ([5b16301](https://github.com/RIVM-bioinformatics/AmpliGone/commit/5b163013cc2c1632cfab9a26e4d1d3906128c22f))
* do not update per-thread progression when processing very low read counts (will throw zerodivisionerror) ([d3bb8b5](https://github.com/RIVM-bioinformatics/AmpliGone/commit/d3bb8b5a7bca2d15dd70d8a1d8669447010b7bc7))
* fixed a typo in fasta2bed ([7968429](https://github.com/RIVM-bioinformatics/AmpliGone/commit/7968429f2a2dac8320998ee3f4c66e689f91d0cc))
* improve error message when too little threads are given ([57e74b5](https://github.com/RIVM-bioinformatics/AmpliGone/commit/57e74b55acff38684b3094ffd6aca7efaa20c40f))
* improved static type checking, linting, and formatting for all files ([2f24ebe](https://github.com/RIVM-bioinformatics/AmpliGone/commit/2f24ebe82e98757be251fd5f089a61a2a265d06b))
* proper grouping of the regex to clean cigar strings ([abe1f69](https://github.com/RIVM-bioinformatics/AmpliGone/commit/abe1f69323d7fac1129d34d75a9f1ce81a70b91b))
* **tests:** add conditional to handle boolean arguments in test ([fc786fe](https://github.com/RIVM-bioinformatics/AmpliGone/commit/fc786fe41e9a80002b8a443a6e54465724cdf5a9))
* **tests:** too short e2e test is expected behaviour. Updated test configuration to include comparison file and correct failure status ([c7ab8a7](https://github.com/RIVM-bioinformatics/AmpliGone/commit/c7ab8a70655285756076bd3a88f220b128f29869))
* **tests:** use relative path for unittest reference file path ([03a2ddf](https://github.com/RIVM-bioinformatics/AmpliGone/commit/03a2ddf663e8e0f593c60bdf932ce589c5f8107a))
* **tests:** use the correct referenceID in the synthetic bedfile ([a5c4d8c](https://github.com/RIVM-bioinformatics/AmpliGone/commit/a5c4d8c32f375a9424cf6cecf39256c687919f80))
* typo in github workflow ([fed3c2e](https://github.com/RIVM-bioinformatics/AmpliGone/commit/fed3c2eeea024107445ffca42a372be87c218ecd))


### Performance Improvements

* perform most alignment-preset related calculations in parallel (multiprocessing with concurrent.futures) ([edf23ac](https://github.com/RIVM-bioinformatics/AmpliGone/commit/edf23aca19d21e0cceb169e1d239c42963f0ac39))
* use pgzip for multithreaded writing of gzip file ([5b16301](https://github.com/RIVM-bioinformatics/AmpliGone/commit/5b163013cc2c1632cfab9a26e4d1d3906128c22f))


### Dependencies

* add flake8-pyproject ([06cb0cb](https://github.com/RIVM-bioinformatics/AmpliGone/commit/06cb0cb1c8e2fe9e4a7c295380252ccb4181f111))
* add pgzip to the dependency list ([cfd25e3](https://github.com/RIVM-bioinformatics/AmpliGone/commit/cfd25e3288b7ee7dfb6816b889b1cb53fa0e7127))
* remove intel channel from conda recipe file ([ed916b8](https://github.com/RIVM-bioinformatics/AmpliGone/commit/ed916b89041b5bd34e299fac2a2f9912d1efd353))
* update dependencies in env.yml ([ed916b8](https://github.com/RIVM-bioinformatics/AmpliGone/commit/ed916b89041b5bd34e299fac2a2f9912d1efd353))
* update dependencies in setup.py ([ed916b8](https://github.com/RIVM-bioinformatics/AmpliGone/commit/ed916b89041b5bd34e299fac2a2f9912d1efd353))
* update dependencies in setup.py and env.yml ([e09cb10](https://github.com/RIVM-bioinformatics/AmpliGone/commit/e09cb109da8cce9a4722a38893c1fbe6f8a2db85))
* update minimum required python version in setup.py ([0fe1beb](https://github.com/RIVM-bioinformatics/AmpliGone/commit/0fe1bebe165f258573f2de9ffb11abeb752f46ab))


### Documentation

* add docstring for percentage_representation function in fasta2bed ([a625b8b](https://github.com/RIVM-bioinformatics/AmpliGone/commit/a625b8b63283281fd61212a96582df50277baad0))
* add virtual primer mode section and related image to user guide ([5da3690](https://github.com/RIVM-bioinformatics/AmpliGone/commit/5da3690b5bbbcf7141ff8ea1407e3f526abd2288))
* expand function docstring to better explain the edge cases ([910f75a](https://github.com/RIVM-bioinformatics/AmpliGone/commit/910f75a7058187e298bf5c43067d04b8c0db042a))
* improves the docstring for this find_or_read_primers() method ([e670d2b](https://github.com/RIVM-bioinformatics/AmpliGone/commit/e670d2bf7ce6783bde9e91434c5ffd23c715602c))
* update citation information in CITATION.cff and README.md ([7e8e0c6](https://github.com/RIVM-bioinformatics/AmpliGone/commit/7e8e0c6d35cd1f6033a8b52ccab5a4378cb43f0c))

## [1.3.1](https://github.com/RIVM-bioinformatics/AmpliGone/compare/v1.3.0...v1.3.1) (2024-04-03)


### Bug Fixes

* don't replace NA name values when reading BED file ([6d0c192](https://github.com/RIVM-bioinformatics/AmpliGone/commit/6d0c192a0b997d98778c2a2f0281e19208aedcac))

## [1.3.0](https://github.com/RIVM-bioinformatics/AmpliGone/compare/v1.2.1...v1.3.0) (2023-08-08)


### Features

* allow multi-fasta reference input ([cd4f413](https://github.com/RIVM-bioinformatics/AmpliGone/commit/cd4f4138118482dea683fde4d201a6e4159b37bb))
* primer removal process now works with reference-specific primer coordinates ([7a1a899](https://github.com/RIVM-bioinformatics/AmpliGone/commit/7a1a8991d455d8eb8767f09b3ea80356ddf42ca3))


### Bug Fixes

* add a default value of 0 to min() to avoid a ValueError when an empty list is provided ([c23301f](https://github.com/RIVM-bioinformatics/AmpliGone/commit/c23301f74d355857f73b79e3ab4629bead26192e))

## [1.2.1](https://github.com/RIVM-bioinformatics/AmpliGone/compare/v1.2.0...v1.2.1) (2023-02-28)


### Bug Fixes

* Fix extension checking for cli arguments ([015e10d](https://github.com/RIVM-bioinformatics/AmpliGone/commit/015e10d9828bc7f5974d914401454c278efde897))
* fix pandas version string in setup.py install_requires ([d0a40f0](https://github.com/RIVM-bioinformatics/AmpliGone/commit/d0a40f04d3e7302fc6065df5f0751812cef32955))
* split file extension checking for input and output files ([c5cc43a](https://github.com/RIVM-bioinformatics/AmpliGone/commit/c5cc43acd501b707580889acb8d399bcffda10fa))
* update permissions for GH-actions workflows ([59fead2](https://github.com/RIVM-bioinformatics/AmpliGone/commit/59fead2015808db8f80026317db3b32d0c5349c0))
* use log.warning() instead of the deprecated log.warn() ([87685f6](https://github.com/RIVM-bioinformatics/AmpliGone/commit/87685f6d0109d108eab6da507ddd8f9159ed06c0))


### Documentation

* add citation file ([59fead2](https://github.com/RIVM-bioinformatics/AmpliGone/commit/59fead2015808db8f80026317db3b32d0c5349c0))
* update conda installation badge ([d88f4f4](https://github.com/RIVM-bioinformatics/AmpliGone/commit/d88f4f4fe0544bea58f4a85ac1d73305532dbe26))

## [1.2.0](https://github.com/RIVM-bioinformatics/AmpliGone/compare/v1.1.0...v1.2.0) (2022-12-20)


### Features

* add support for fragmented amplicon reads with the `--amplicon-type fragmented` mode ([9d5559b](https://github.com/RIVM-bioinformatics/AmpliGone/commit/9d5559bd261af37350ad760bb9212f82685b4caf))


### Bug Fixes

* only check the last file extension to see for a fastq or bam file (fixes check for file with multiple dots in the filename) ([386fb1f](https://github.com/RIVM-bioinformatics/AmpliGone/commit/386fb1ff9a6cd405166fc75fd5b0cd2f2c1fd89a))
* print the correct coordinates when a primer is found multiple times on the reverse strand ([420438c](https://github.com/RIVM-bioinformatics/AmpliGone/commit/420438cb18d870f5807998996fde8c7ddc2e362c))


### Documentation

* update docstrings ([c8af27a](https://github.com/RIVM-bioinformatics/AmpliGone/commit/c8af27a995bc04296e3553e437642af4297ce4f9))
* update readme ([0287408](https://github.com/RIVM-bioinformatics/AmpliGone/commit/02874081733dc5dcd1aa06a9ae266755a93596db))
* update user documentation to include the new "fragmented" amplicon-type ([159366e](https://github.com/RIVM-bioinformatics/AmpliGone/commit/159366e0a6832dc4dacfcd484fbe2ca2debd2806))

## [1.1.0](https://www.github.com/RIVM-bioinformatics/AmpliGone/compare/v1.0.3...v1.1.0) (2022-08-11)


### Features

* output some basic stats in the log ([c5ad032](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/c5ad032e38bf86c341be8f685ad70f327c1ecdfb))
* use `rich` as a logging handler to improve logging with timestamps and various information levels. ([c5ad032](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/c5ad032e38bf86c341be8f685ad70f327c1ecdfb))


### Dependencies

* add `rich` as a dependency ([c5ad032](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/c5ad032e38bf86c341be8f685ad70f327c1ecdfb))
* change parmap version to lenient version 1.5.x ([0c4d2a5](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/0c4d2a54974e98542b2b64565eb645898fa36710))
* remove `tqdm` as a dependency ([c5ad032](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/c5ad032e38bf86c341be8f685ad70f327c1ecdfb))
* update mappy to latest version 2.24 ([efed0a1](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/efed0a162746afb306351d42139a474e3854b873))
* update pandas to >=1.3 ([efed0a1](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/efed0a162746afb306351d42139a474e3854b873))
* update parmap to latest version 1.5.3 ([efed0a1](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/efed0a162746afb306351d42139a474e3854b873))
* update pysam to lenient version 0.19.* ([efed0a1](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/efed0a162746afb306351d42139a474e3854b873))
* updated environment recipe and setup.py to newest dependencies ([efed0a1](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/efed0a162746afb306351d42139a474e3854b873))

### [1.0.3](https://www.github.com/RIVM-bioinformatics/AmpliGone/compare/v1.0.2...v1.0.3) (2022-07-12)


### Bug Fixes

* **ci:** Use other access key in order to trigger publishing workflow ([384625e](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/384625e5171da3e9fc674f0dd3174dfd7ba54ecf))

### [1.0.2](https://www.github.com/RIVM-bioinformatics/AmpliGone/compare/v1.0.1...v1.0.2) (2022-07-12)


### Bug Fixes

* Fix cutting reverse primers too short ([a01c69d](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/a01c69ddd21d69bd9c098f6215be518c5455dbed))
* FWSet and RVSet are now different objects ([a28d629](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/a28d6292aa8c08e0bdd2c030455112d846626803))


### Documentation

* **ci:** update mkdocs config to support tabbed content ([bdfb52d](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/bdfb52d00318378e9dacc918b7e8ce855ca9c319))
* restructured user guide page, clarified some sections and added many more example commands for reference ([d723b97](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/d723b975a6378b1174dac2b011c6b4b7382cf78f))
* rewritten installation page with clarified instructions and added conda instructions ([271a0a8](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/271a0a8d238de0a9ad492a829f8597a2e87a8d66))
* update documentation site landing page ([45a3fca](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/45a3fca2c5a9ed7e719957e9f5f2c29a25baf300))
* update readme with quick installation instructions ([bdfb52d](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/bdfb52d00318378e9dacc918b7e8ce855ca9c319))

### [1.0.1](https://www.github.com/RIVM-bioinformatics/AmpliGone/compare/v1.0.0...v1.0.1) (2022-04-14)


### Continuous Integration

* update github actions workflow, split package upload to separate workflow ([77c13a1](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/77c13a1a84323070541ef65bad57922cf29a0b77))

## [1.0.0](https://www.github.com/RIVM-bioinformatics/AmpliGone/compare/v0.4.3...v1.0.0) (2022-04-14)


### ⚠ BREAKING CHANGES

* Export BED file of cut primers

### Features

* Add support for BED file input ([eaa60a0](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/eaa60a048e0361e8b1ed8ba0f55ec995ab8314e8))
* Export BED file of cut primers ([f51d222](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/f51d222f9bc6f0e6ca0e55267b128cf7a45722eb))


### Bug Fixes

* Add seq and revseq to coordinatelist ([eb27331](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/eb27331ec3110b8924eee0c6dc1bc4f8201516bf))
* Fix parsing error rate from cli ([28976f2](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/28976f2260385481e67d4232dbfbf7f643e054b6))
* Fix parsing of "." in score column ([cc6913b](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/cc6913b515d549394b2c900356930ef4be79d7b8))
* Fix running fasta2bed script ([f8545a4](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/f8545a432097ea27c39d0941177b34988ca52845))
* Output empty bed if no primers are found ([ca67b0b](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/ca67b0bfa79d3fcd6d759f91abf4809505f72fbf))


### Dependencies

* update all python versions mentions to 3.8 ([89ced9c](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/89ced9cb1d0070b051300e987c7edc5cfe62ea23))


### Documentation

* add docstrings to `AmpliGone.mappreset` ([6d9c32c](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/6d9c32c3d25f84299cf0432c8d3aeb2c64dd5c13))
* add docstrings to functions in `AmpliGone.AmpliGone` ([54a0a01](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/54a0a01531953cc1e38b2e6a7eecb5d8fd21357c))
* add docstrings to functions in AmpliGone.fasta2bed ([ea2d5c6](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/ea2d5c6b09bf976b052e3c4d716a0ca216c785e8))
* add docstrings to functions in AmpliGone.io_ops ([d9c60ae](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/d9c60aeaf10f51245b89593c4f4f1a08c3980b92))
* Update installation instructions ([fbb7e0e](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/fbb7e0ec47008301f2a303044abe20790f0864d4))
* Update user guide ([071e71f](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/071e71fecf972731ccb8db83715a032c7fdc4406))

### [0.4.3](https://www.github.com/RIVM-bioinformatics/AmpliGone/compare/v0.4.2...v0.4.3) (2022-02-03)


### Bug Fixes

* remove unnecessary print statement ([c7bd712](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/c7bd71241f036791517d1d884f2458a45698e1a6))


### Documentation

* add `--error-rate` flag to user guide documentation ([9db3848](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/9db3848c6fcb817bd730273ec82eec5eeab7e835))

### [0.4.2](https://www.github.com/RIVM-bioinformatics/AmpliGone/compare/v0.4.1...v0.4.2) (2022-02-03)


### Features

* allow the primer-search error rate to be set by the user ([fddd94a](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/fddd94a555b57a35d3bb064c1d53b03777c72f24))

### [0.4.1](https://www.github.com/RIVM-bioinformatics/AmpliGone/compare/v0.4.0...v0.4.1) (2022-02-02)


### Bug Fixes

* don't report duplicate coordinates of same primer when ambiguity-primers are given ([f5449e2](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/f5449e2a5d5c9e502192ace4f836c90a45060afa))
* use Nonetypes instead of string "None" ([f5449e2](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/f5449e2a5d5c9e502192ace4f836c90a45060afa))

## [0.4.0](https://www.github.com/RIVM-bioinformatics/AmpliGone/compare/v0.3.3...v0.4.0) (2022-01-27)


### Features

* find coordinates of primers with mismatches to reference ([d6bf925](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/d6bf925087a904daa2568bedda01e9b3a1504011))


### Bug Fixes

* consistently remove reverse primers in end-to-end mode ([0c4f621](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/0c4f621022f3d58820aa125da97a18b345d5e649))
* Fix cutting too many bases ([78bd9ea](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/78bd9eac8bb4b37d03c043fde35164f0d3507fc5))
* Make number of iterations variable ([62ac2d7](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/62ac2d792b36a7b1206ba069b05f7d045f0c35ec))


### Documentation

* update installation instructions and include primer-search changes ([5fdc6fa](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/5fdc6fa4a6371a9ea86c2771129182d7a98af9e2))

### [0.3.3](https://www.github.com/RIVM-bioinformatics/AmpliGone/compare/v0.3.2...v0.3.3) (2021-12-08)


### Bug Fixes

* solve preset-sampling crash when less than 3 input reads are provided ([29d5624](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/29d5624e0b2871663a7c8e6cee9f3553c4d1debb))
* solve reference before assignment error in end-to-mid mode ([14c9d24](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/14c9d241fc9fd806fd018049ff67118c7b59517b))

### [0.3.2](https://www.github.com/RIVM-bioinformatics/AmpliGone/compare/v0.3.1...v0.3.2) (2021-11-26)


### Bug Fixes

* report removed coordinates from all iterations ([f726e9a](https://www.github.com/RIVM-bioinformatics/AmpliGone/commit/f726e9a9e740f9591496c8722aba2ed10aba8b80))

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
