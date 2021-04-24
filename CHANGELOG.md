# Change Log
All notable changes to this project will be documented in this file, following the suggestions of [Keep a CHANGELOG](http://keepachangelog.com/). This project adheres to [Semantic Versioning](http://semver.org/) for its most widely used - and defacto - public interfaces.

Note that since we don't clearly distinguish between a public and private interfaces there will be changes in non-major versions that are potentially breaking. If we make breaking changes to less used interfaces we will highlight it in here.

## [Unreleased]

- Ability to pass ``Canvas3DContext`` to ``PluginContext.fromCanvas``.
- Relative frame support for ``Canvas3D`` viewport.
- Fix bug in screenshot copy UI.
- Add ability to select residues from a list of identifiers to the Selection UI.
- Fix SSAO bugs when used with ``Canvas3D`` viewport.
- Support for  full pausing (no draw) rendering: ``Canvas3D.pause(true)``.
- Add `MeshBuilder.addMesh`.
- Add `Torus` primitive.

## [v2.0.4] - 2021-04-20

- [WIP] Mesh export extension
- ``Structure.eachAtomicHierarchyElement`` (#161)
- Fixed reading multi-line values in SDF format
- Fixed Measurements UI labels (#166)

## [v2.0.3] - 2021-04-09
### Added
- Support for ``ColorTheme.palette`` designed for providing gradient-like coloring.

### Changed
- [Breaking] The `zip` function is now asynchronous and expects a `RuntimeContext`. Also added `Zip()` returning a `Task`.
- [Breaking] Add ``CubeGridFormat`` in ``alpha-orbitals`` extension.

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
