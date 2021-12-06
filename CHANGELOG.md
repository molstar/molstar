# Change Log
All notable changes to this project will be documented in this file, following the suggestions of [Keep a CHANGELOG](http://keepachangelog.com/). This project adheres to [Semantic Versioning](http://semver.org/) for its most widely used - and defacto - public interfaces.

Note that since we don't clearly distinguish between a public and private interfaces there will be changes in non-major versions that are potentially breaking. If we make breaking changes to less used interfaces we will highlight it in here.


## [Unreleased]

- Add ``bumpiness`` (per-object and per-group), ``bumpFrequency`` & ``bumpAmplitude`` (per-object) render parameters (#299)

## [v3.0.0-dev.3] - 2021-12-4

- Fix OBJ and USDZ export

## [v3.0.0-dev.2] - 2021-12-1

- Do not include tests and source maps in NPM package

## [v3.0.0-dev.0] - 2021-11-28

- Add multiple lights support (with color, intensity, and direction parameters)
- [Breaking] Add per-object material rendering properties
  - ``SimpleSettingsParams.lighting.renderStyle`` and ``RendererParams.style`` were removed
- Add substance theme with per-group material rendering properties
- ``StructureComponentManager.Options`` state saving support
- ``ParamDefinition.Group.presets`` support

## [v2.4.1] - 2021-11-28

- Fix: allow atoms in aromatic rings to do hydrogen bonds

## [v2.4.0] - 2021-11-25

- Fix secondary-structure property handling
    - StructureElement.Property was incorrectly resolving type & key
    - StructureSelectionQuery helpers 'helix' & 'beta' were not ensuring property availability
- Re-enable VAO with better workaround (bind null elements buffer before deleting)
- Add ``Representation.geometryVersion`` (increments whenever the geometry of any of its visuals changes)
- Add support for grid-based smoothing of Overpaint and Transparency visual state for surfaces

## [v2.3.9] - 2021-11-20

- Workaround: switch off VAO support for now

## [v2.3.8] - 2021-11-20

- Fix double canvas context creation (in plugin context)
- Fix unused vertex attribute handling (track which are used, disable the rest)
- Workaround for VAO issue in Chrome 96 (can cause WebGL to crash on geometry updates)

## [v2.3.7] - 2021-11-15

- Added ``ViewerOptions.collapseRightPanel``
- Added ``Viewer.loadTrajectory`` to support loading "composed" trajectories (e.g. from gro + xtc)
- Fix: handle parent in Structure.remapModel
- Add ``rounded`` and ``square`` helix profile options to Cartoon representation (in addition to the default ``elliptical``)

## [v2.3.6] - 2021-11-8

- Add additional measurement controls: orientation (box, axes, ellipsoid) & plane (best fit)
- Improve aromatic bond visuals (add ``aromaticScale``, ``aromaticSpacing``, ``aromaticDashCount`` params)
- [Breaking] Change ``adjustCylinderLength`` default to ``false`` (set to true for focus representation)
- Fix marker highlight color overriding select color
- CellPack extension update
    - add binary model support
    - add compartment (including membrane) geometry support
    - add latest mycoplasma model example
- Prefer WebGL1 in Safari 15.1.

## [v2.3.5] - 2021-10-19

- Fix sequence viewer for PDB files with COMPND record and multichain entities.
- Fix index pair bonds order assignment

## [v2.3.4] - 2021-10-12

- Fix pickScale not taken into account in line/point shader
- Add pixel-scale, pick-scale & pick-padding GET params to Viewer app
- Fix selecting bonds not adding their atoms in selection manager
- Add ``preferAtoms`` option to SelectLoci/HighlightLoci behaviors
- Make the implicit atoms of bond visuals pickable
    - Add ``preferAtomPixelPadding`` to Canvas3dInteractionHelper
- Add points & crosses visuals to Line representation
- Add ``pickPadding`` config option (look around in case target pixel is empty)
- Add ``multipleBonds`` param to bond visuals with options: off, symmetric, offset
- Fix ``argparse`` config in servers.

## [v2.3.3] - 2021-10-01

- Fix direct volume shader

## [v2.3.2] - 2021-10-01

- Prefer WebGL1 on iOS devices until WebGL2 support has stabilized.

## [v2.3.1] - 2021-09-28

- Add Charmm saccharide names
- Treat missing occupancy column as occupancy of 1
- Fix line shader not accounting for aspect ratio
- [Breaking] Fix point repr & shader
    - Was unusable with ``wboit``
    - Replaced ``pointFilledCircle`` & ``pointEdgeBleach`` params by ``pointStyle`` (square, circle, fuzzy)
    - Set ``pointSizeAttenuation`` to false by default
    - Set ``sizeTheme`` to ``uniform`` by default
- Add ``markerPriority`` option to Renderer (useful in combination with edges of marking pass)
- Add support support for ``chem_comp_bond`` and ``struct_conn`` categories (fixes ModelServer behavior where these categories should have been present)
- Model and VolumeServer: fix argparse config

## [v2.3.0] - 2021-09-06

- Take include/exclude flags into account when displaying aromatic bonds
- Improve marking performance
    - Avoid unnecessary draw calls/ui updates when marking
    - Check if loci is superset of visual
    - Check if loci overlaps with unit visual
    - Ensure ``Interval`` is used for ranges instead of ``SortedArray``
    - Add uniform marker type
    - Special case for reversing previous mark
- Add optional marking pass
    - Outlines visible and hidden parts of highlighted/selected groups
    - Add highlightStrength/selectStrength renderer params

## [v2.2.3] - 2021-08-25

- Add ``invertCantorPairing`` helper function
- Add ``Mesh`` processing helper ``.smoothEdges``
- Smooth border of molecular-surface with ``includeParent`` enabled
- Hide ``includeParent`` option from gaussian-surface visuals (not particularly useful)
- Improved ``StructureElement.Loci.size`` performance (for marking large cellpack models)
- Fix new ``TransformData`` issues (camera/bounding helper not showing up)
- Improve marking performance (avoid superfluous calls to ``StructureElement.Loci.isWholeStructure``)

## [v2.2.2] - 2021-08-11

- Fix ``TransformData`` issues [#133](https://github.com/molstar/molstar/issues/133)
- Fix ``mol-script`` query compiler const expression recognition.

## [v2.2.1] - 2021-08-02

- Add surrounding atoms (5 Angstrom) structure selection query
- [Breaking] Add maxDistance prop to ``IndexPairBonds``
- Fix coordinateSystem not handled in ``Structure.asParent``
- Add ``dynamicBonds`` to ``Structure`` props (force re-calc on model change)
    - Expose as optional param in root structure transform helper
- Add overpaint support to geometry exporters
- ``InputObserver`` improvements
  - normalize wheel speed across browsers/platforms
  - support Safari gestures (used by ``TrackballControls``)
  - ``PinchInput.fractionDelta`` and use it in ``TrackballControls``

## [v2.2.0] - 2021-07-31

- Add ``tubularHelices`` parameter to Cartoon representation
- Add ``SdfFormat`` and update SDF parser to be able to parse data headers according to spec (hopefully :)) #230
- Fix mononucleotides detected as polymer components (#229)
- Set default outline scale back to 1
- Improved DCD reader cell angle handling (interpret near 0 angles as 90 deg)
- Handle more residue/atom names commonly used in force-fields
- Add USDZ support to ``geo-export`` extension.
- Fix ``includeParent`` support for multi-instance bond visuals.
- Add ``operator`` Loci granularity, selecting everything with the same operator name.
- Prefer ``_label_seq_id`` fields in secondary structure assignment.
- Support new EMDB API (https://www.ebi.ac.uk/emdb/api/entry/map/[EMBD-ID]) for EM volume contour levels.
- ``Canvas3D`` tweaks:
    - Update ``forceDraw`` logic.
    - Ensure the scene is re-rendered when viewport size changes.
    - Support ``noDraw`` mode in ``PluginAnimationLoop``.

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
- Add ``MeshBuilder.addMesh``.
- Add ``Torus`` primitive.
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

- Add support for ``ColorTheme.palette`` designed for providing gradient-like coloring.
- [Breaking] The ``zip`` function is now asynchronous and expects a ``RuntimeContext``. Also added ``Zip()`` returning a ``Task``.
- [Breaking] Add ``CubeGridFormat`` in ``alpha-orbitals`` extension.

## [v2.0.2] - 2021-03-29

- Add ``Canvas3D.getRenderObjects``.
- [WIP] Animate state interpolating, including model trajectories
- Recognise MSE, SEP, TPO, PTR and PCA as non-standard amino-acids.
- Fix VolumeFromDensityServerCif transform label

## [v2.0.1] - 2021-03-23

- Exclude tsconfig.commonjs.tsbuildinfo from npm bundle

## [v2.0.0] - 2021-03-23

Too many changes to list as this is the start of the changelog... Notably, default exports are now forbidden.
