# Change Log
All notable changes to this project will be documented in this file, following the suggestions of [Keep a CHANGELOG](http://keepachangelog.com/). This project adheres to [Semantic Versioning](http://semver.org/) for its most widely used - and defacto - public interfaces.

Note that since we don't clearly distinguish between a public and private interfaces there will be changes in non-major versions that are potentially breaking. If we make breaking changes to less used interfaces we will highlight it in here.


## [Unreleased]
### Changed
- [Breaking] The `zip` function is now asynchronous and expects a `RuntimeContext`. Also added `Zip()` returning a `Task`.
- Add support for ``ColorTheme.palette`` designed for providing gradient-like coloring.


## [v2.0.2] - 2021-03-29
### Added
- `Canvas3D.getRenderObjects`.
- [WIP] Animate state interpolating, including model trajectories

### Changed
- Recognise MSE, SEP, TPO, PTR and PCA as non-standard amino-acids.

### Fixed
- VolumeFromDensityServerCif transform label


## [v2.0.1] - 2021-03-23
### Fixed
- Exclude tsconfig.commonjs.tsbuildinfo from npm bundle


## [v2.0.0] - 2021-03-23
Too many changes to list as this is the start of the changelog... Notably, default exports are now forbidden.
