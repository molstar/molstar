# Change Log
All notable changes to this project will be documented in this file, following the suggestions of [Keep a CHANGELOG](http://keepachangelog.com/). This project adheres to [Semantic Versioning](http://semver.org/) for its most widely used - and defacto - public interfaces.

Note that since we don't clearly distinguish between a public and private interfaces there will be changes in non-major versions that are potentially breaking. If we make breaking changes to less used interfaces we will highlight it in here.


## [Unreleased]

- Add `tubularHelices` parameter to Cartoon representation
- Add USDZ support to ``geo-export`` extension.

## [v2.1.0] - 2021-07-05

- Add parameter for to display aromatic bonds as dashes next to solid cylinder/line.
- Add backbone representation
- Fix outline in orthographic mode and set default scale to 2.

## [v2.0.7] - 2021-06-23

- Add ability to specify ``volumeIndex`` in ``Viewer.loadVolumeFromUrl`` to better support Volume Server inputs.
- Support in-place reordering for trajectory ``Frame.x/y/z`` arrays for better memory efficiency.
- Fixed text CIF encoder edge cases (most notably single whitespace not being escaped).

## [v2.0.6] - 2021-06-01

- Add glTF (GLB) and STL support to ``geo-export`` extension.
- Protein crosslink improvements
    - Change O-S bond distance to allow for NOS bridges (doi:10.1038/s41586-021-03513-3)
    - Added NOS-bridges query & improved disulfide-bridges query
- Fix #178: ``IndexPairBonds`` for non-single residue structures (bug due to atom reordering).
- Add volumetric color smoothing for MolecularSurface and GaussianSurface representations (#173)
- Fix nested 3d grid lookup that caused results being overwritten in non-covalent interactions computation.
- Basic implementation of ``BestDatabaseSequenceMapping`` (parse from CIF, color theme, superposition).
- Add atom id ranges support to Selection UI.

## [v2.0.5] - 2021-04-26

- Ability to pass ``Canvas3DContext`` to ``PluginContext.fromCanvas``.
- Relative frame support for ``Canvas3D`` viewport.
- Fix bug in screenshot copy UI.
- Add ability to select residues from a list of identifiers to the Selection UI.
- Fix SSAO bugs when used with ``Canvas3D`` viewport.
- Support for  full pausing (no draw) rendering: ``Canvas3D.pause(true)``.
- Add `MeshBuilder.addMesh`.
- Add `Torus` primitive.
- Lazy volume loading support.
- [Breaking] ``Viewer.loadVolumeFromUrl`` signature change.
    - ``loadVolumeFromUrl(url, format, isBinary, isovalues, entryId)`` => ``loadVolumeFromUrl({ url, format, isBinary }, isovalues, { entryId, isLazy })``
- Add ``TextureMesh`` support to ``geo-export`` extension.

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
