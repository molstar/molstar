# Change Log
All notable changes to this project will be documented in this file, following the suggestions of [Keep a CHANGELOG](http://keepachangelog.com/). This project adheres to [Semantic Versioning](http://semver.org/) for its most widely used - and defacto - public interfaces.

Note that since we don't clearly distinguish between a public and private interfaces there will be changes in non-major versions that are potentially breaking. If we make breaking changes to less used interfaces we will highlight it in here.

## [Unreleased]
- [Breaking] Renamed some color schemes ('inferno' -> 'inferno-no-black', 'magma' -> 'magma-no-black', 'turbo' -> 'turbo-no-black', 'rainbow' -> 'simple-rainbow')
- [Breaking] `Box3D.nearestIntersectionWithRay` -> `Ray3D.intersectBox3D`
- [Breaking] `Plane3D.distanceToSpher3D` -> `distanceToSphere3D` (fix spelling)
- [Breaking] fix typo `MarchinCubes` -> `MarchingCubes`
- [Breaking] `PluginContext.initViewer/initContainer/mount` are now async and have been renamed to include `Async` postfix
- [Breaking] Add `Volume.instances` support and a `VolumeInstances` transform to dynamically assign it
  - This change is breaking because all volume objects require the `instances` field now.
- [Breaking] `Canvas3D.identify` now expects `Vec2` or `Ray3D`
- [Breaking] `TrackballControlsParams.animate.spin.speed` now means "Number of rotations per second" instead of "radians per second"
- Update production build to use `esbuild`
- Emit explicit paths in `import`s in `lib/`
- Fix outlines on opaque elements using illumination mode
- Change `Representation.Empty` to a lazy property to avoid issue with some bundlers
- MolViewSpec extension:
  - Generic color schemes (`palette` parameter for color_from_* nodes)
  - Annotation field remapping (`field_remapping` parameter for color_from_* nodes)
  - `representation` node: support custom property `molstar_reprepresentation_params`
  - Add `backbone` and `line` representation types
  - `primitives` node: support custom property `molstar_mesh/label/line_params`
  - `canvas` node: support custom property `molstar_postprocessing` with the ability to customize outline, depth of field, bloom, shadow, occlusion (SSAO), and fog
  - `clip` node support for structure and volume representations
  - `grid_slice` representation support for volumes
  - Support tethers and background for primitive labels
  - Support `snapshot_key` parameter on primitives that enables transition between states via clicking on 3D objects
  - Inline selectors and MVS annotations support `instance_id`
  - Support `matrix` on transform params
  - Support `surface_type` (`molecular` / `gaussian`) on for `surface` representation nodes
  - Add `instance` node type
  - Add `transform.rotation_center` property that enables rotating an object around its centroid or a specific point
  - Support transforming and instancing of structures, components, and volumes
  - Use params hash for node version for more performant tree diffs
  - Add `Snapshot.animation` support that enables animating almost every property in a given tree
  - Add `createMVSX` helper function
  - Support Mol* trackball animation via `animation.custom.molstar_trackball`
  - MVSX - use Murmur hash instead of FNV in archive URI
  - Support additional file formats (pdbqt, gro, xyz, mol, sdf, mol2, xtc, lammpstrj)
  - Support loading trajectory coordinates from separate nodes
  - Trigger markdown commands from primitives using `molstar_markdown_commands` custom extensions
  - Support `molstar_on_load_markdown_commands` custom state on the `root` node
- Added new color schemes, synchronized with D3.js ('inferno', 'magma', 'turbo', 'rainbow', 'sinebow', 'warm', 'cool', 'cubehelix-default', 'category-10', 'observable-10', 'tableau-10')
- Snapshot Markdown improvements
  - Add `MarkdownExtensionManager` (`PluginContext.managers.markdownExtensions`)
  - Support custom markdown commands to control the plugin via the `[link](!command)` pattern
  - Support rendering custom elements via the `![alt](!parameters)` pattern
  - Support tables
  - Support loading images and audio from MVSX files
  - Indicate external links with ⤴
  - Audio support
  - Add `PluginState.Snapshot.onLoadMarkdownCommands`
- Avoid calculating rings for coarse-grained structures
- Fix isosurface compute shader normals when transformation matrix is applied to volume
- Symmetry operator naming for spacegroup symmetry - parenthesize multi-character indices (1_111-1 -> 1_(11)1(-1))
- Add `SymmetryOperator.instanceId` that corresponds to a canonical operator name (e.g. ASM-1, ASM-X0-1 for assemblies, 1_555, 1_(11)1(-1) for crystals)
- Mol2 Reader
    - Fix column count parsing
    - Add support for substructure
- Fix shader error when clipping flags are set without clip objects present
- Fix wrong group count calculation on geometry update (#1562)
- Fix wrong instance index in `calcMeshColorSmoothing`
- Add `Ray3D` object and helpers
- Volume slice representation: add `relativeX/Y/Z` options for dimension
- Add `StructureInstances` transform
- `mvs-stories` app
  - Add `story-id` URL arg support
  - Add `story-session-url` URL arg support
  - Add "Download MVS State" link
  - Add "Open in Mol*" link
  - Add "Edit in MolViewStories" link for story states
- Add ray-based picking
    - Render narrow view of scene scene from ray origin & direction to a few pixel sized viewport
    - Cast ray on every input as opposed to the standard "whole screen" picking
    - Can be enabled with new `Canvas3dInteractionHelperParams.convertCoordsToRay` param
    - Allows to have input methods that are 3D pointers in the scene
    - Add `ray: Ray3D` property to `DragInput`, `ClickInput`, and `MoveInput`
- Add async, non-blocking picking (only WebGL2)
    - Refactor `Canvas3dInteractionHelper` internals to use async picking for move events
- Add `enable` param for post-processing effects. If false, no effects are applied.
- Dot volume representation improvements
    - Add positional perturbation to avoid camera artifacts
    - Fix handling of negative isoValues by considering only volume cells with values lower than isoValue (#1559)
    - Fix volume-value size theme
- Change the parsing of residue names in PDB files from 3-letter to 4-letter.
- Support versioning transform using a hash function in `mol-state`
- Support for "state snapshot transitions"
    - Add `PluginState.Snapshot.transition` that enables associating a state snapshot with a list states that can be animated
    - Add `AnimateStateSnapshotTransition` animation
    - Update the snapshots UI to support this feature
- Use "proper time" in the animation loop to prevent animation skips during blocking operations (e.g., shader complication)
- Add `Hsl` and (normalized) `Rgb` color spaces
- Add `Color.interpolateHsl`
- Add `rotationCenter` property to `TransformParam`
- Add Monolayer transparency (exploiting dpoit).
- Add plugin config item ShowReset (shows/hides "Reset Zoom" button)
- Fix transform params not being normalized when used together with param hash version
- Replace `immer` with `mutative`

## [v4.18.0] - 2025-06-08
- MolViewSpec extension:
  - Support for label_comp_id and auth_comp_id in annotations
  - Geometric primitives - do not render if position refers to empty substructure
  - Primitive arrow - nicer default cap size (relative to tube_radius)
  - Primitive angle_measurement - added vector_radius param
  - Fix MVSX file assets being disposed in multi-snapshot states
- Add `mol-utils/camera.ts` with `fovAdjustedPosition` and `fovNormalizedCameraPosition`
- Show FOV normalized position in `CameraInfo` UI and use it in "Copy MVS State"
- Support static resources in `AssetManager`
- General:
  - Use `isolatedModules` tsconfig flag
  - Fix TurboPack build when using ES6 modules
- Support `pickingAlphaThreshold` when `xrayShaded` is enabled
- Support sampling from arbitrary planes for structure plane and volume slice representations
- Refactor SCSS to not use `@import` (fixes deprecation warnings)

## [v4.17.0] - 2025-05-22
- Remove `xhr2` dependency for NodeJS, use `fetch`
- Add `mvs-stories` app included in the `molstar` NPM package
  - Use the app in the corresponding example
- Interactions extension: remove `salt-bridge` interaction kind (since `ionic` is supported too)

## [v4.16.0] - 2025-05-20
- Load potentially big text files as `StringLike` to bypass string size limit
- MolViewSpec extension:
  - Load single-state MVS as if it were multi-state with one state
  - Merged `loadMVS` options `keepCamera` and `keepSnapshotCamera` -> `keepCamera`
  - Removed `loadMVS` option `replaceExisting` (is now default)
  - Added `loadMVS` option `appendSnapshots`
- Fix camera not being interpolated in MP4 export due to updates in WebGL ContextLost handling

## [v4.15.0] - 2025-05-19
- IHM improvements:
  - Disable volume streaming
  - Disable validation report visualization
  - Enable assembly symmetry for integrative models
- Fix transparency rendering with occlusion in NodeJS
- mmCIF Support
  - Add custom `molstar_bond_site` category that enables serializing explicit bonds by referencing `atom_site.id`
  - Add `includeCategoryNames`, `keepAtomSiteId`, `exportExplicitBonds`, `encoder` properties to `to_mmCIF` exporter
- Add support for attachment points property (`M APO`) to the MOL V2000 parser
- Add `json-cif` extension that should pave way towards structure editing capabilities in Mol\*
  - JSON-based encoding of the CIF data format
  - `JSONCifLigandGraph` that enables editing of small molecules via modifying `atom_site` and `molstar_bond_site` categories
- Add `ligand-editor` example that showcases possible use-cases of the `json-cif` extension
- Breaking (minor): Changed `atom_site.id` indexing to 1-based in `mol-model-formats/structure/mol.ts::getMolModels`.
- WebGL ContextLost handling
    - Fix missing framebuffer & drawbuffer re-attachments
    - Fix missing cube texture re-initialization
    - Fix missing extensions reset
    - Fix timer clearing edge case
    - Add reset support for geometry generated on the GPU

## [v4.14.1] - 2025-05-09
- Do not raise error when creating duplicate state transformers and print console warning instead

## [v4.14.0] - 2025-05-07
- Fix `Viewer.loadTrajectory` when loading a topology file
- Fix `StructConn.residueCantorPairs` to not include identity pairs
- Add format selection option to image export UI (PNG, WebP, JPEG)
- Add `StateBuilder.To.updateState`
- MVS:
  - Support updating transform states
  - Add support for `is_hidden` custom state as an extension
  - Add `queryMVSRef` and `createMVSRefMap` utility functions
- Adjust max resolution of surfaces for auto quality (#1501)
- Fix switching representation type in Volume UI
- VolumeServer: Avoid grid expansion when requiring unit cell (avoids including an extra layer of cells outside the unit cell query box)

## [v4.13.0] - 2025-04-14
- Support `--host` option for build-dev.mjs script
- Add `Viewer.loadFiles` to open supported files
- Support installing the viewer as a Progressive Web App (PWA)
- `ihm-restraints` example: show entity labels
- Fix `element-point` visual not using child unit
- Ignore `renderables` with empty draw count
- Add experimental support for `esbuild` for development
  - Use `npm run dev` for faster development builds
- Use `StructureElement.Bundle` instead of expressions to serialize measurement elements
  - Fixes measurements not being supported for coarse models
- Implementation of `ColorScale.createDiscrete` (#1458)
- Add `ColorScale.createDiscrete` to the `uncertainty` color theme
- Fix color palette shown in the UI (for non-gradient palettes)
- Fix colors description in the UI (when using custom thresholds)
- Fix an edge case in the UI when the user deletes all colors from the color list
- Add `interactions` extension and a corresponding example that utilizes it
- Add element source index to default atomic granularity hover labels
- Add `StructureElement.Schema` based on corresponding MolViewSpec implementation that allows data-driven selection of structural elements
- Add `StructureElement.Loci/Bundle.fromExpression/Query/Schema` helper functions
- Add `addLinkCylinderMesh` (from `createLinkCylinderMesh`)
- Add `Unit.transientCache` and `Unit.getCopy`
- Fix `ElementBondIterator` indices mapping logic for inter-unit bonds
- Fix `pickPadding` and `pickScale` not updating `PickHelper`
- MolViewSpec extension: support loading extensions when loading multistate files
- Do not add bonds for pairs of residues that have a `struct_conn` entry
- Improved `ma_qa_metric` support
  - Parse all local metrics
  - Ability to select alternate metrics in the pLDDT/qmean themes
  - Do not assume PAE plot is symmetric
- Added `PluginConfig.Viewport.ShowScreenshotControls` to control visibility of screenshot controls
- Fix MolViewSpec builder for volumes.
- Generalize `mvs-kinase-story` example to `mvs-stories`
  - Add TATA-binding protein story
  - Improve the Kinase story
- Fix alpha orbitals example

## [v4.12.0] - 2025-02-28

- Fix PDBj structure data URL
- Improve logic when to cull in renderer
- Add `atom.ihm.has-seq-id` and `atom.ihm.overlaps-seq-id-range` symbol to the query language
- MolViewSpec extension:
  - Add box, arrow, ellipse, ellipsoid, angle primitives
  - Add basic support for volumetric data (map, Volume Server)
  - Add support for `molstar_color_theme_name` custom extension
  - Better IH/M support:
    - Support `coarse` components
    - Support `spacefill` representation
    - Support `carbohydrate` representation
    - Support for `custom.molstar_use_default_coloring` property on Color node.
    - Use `atom.ihm.has-seq-id` and `atom.ihm.overlaps-seq-id-range` for matching `label_seq_id` locations to support querying coarse elements.
    - Add ihm-restraints example
- Add `mvs-kinase-story` example
- Remove static uses of `ColorTheme` and `SizeTheme` fields. Should resolvent "undefined" errors in certain builds
- Add `transform` property to clip objects
- Add support for trimming `image` geometry to a box
- Improve/fix iso-level support of `slice` representation
- Add support for rotating `slice` representation around an axis
- Add default color support for palette based themes
- Add `plane` structure representation
    - Can be colored with any structure theme
    - Can be colored with the `external-volume` theme
    - Can show atoms as a cutout
    - Supports principal axes and bounding box as a reference frame
- Add `Camera` section to "Screenshot / State" controls
- Add `CoarseIndex` for fast lookup of coarse elements

## [v4.11.0] - 2025-01-26

- Fix for tubular helices issue (Fixes #1422)
- Volume UI improvements
    - Render all volume entries instead of selecting them one-by-one
    - Toggle visibility of all volumes
    - More accessible iso value control
- Support wheel event on sliders
- MolViewSpec extension:
    - Add validation for discriminated union params
    - Primitives: remove triangle_colors, line_colors, have implicit grouping instead; rename many parameters
- UI configuration options
    - Support removal of independent selection controls in the viewport
    - Support custom selection controls
    - Support for custom granularity dropdown options
    - Support for custom Sequence Viewer mode options
- Add `external-structure` theme that colors any geometry by structure properties
- Support float and half-float data type for direct-volume rendering and GPU isosurface extraction
- Minor documentation updates
- Add support for position-location to `volume-value` color theme
- Add support for color themes to `slice` representation
- Improve/fix palette support in volume color themes
- Fix `Plane3D.projectPoint`
- Fix marking related `image` rendering issues
    - Handle pixels without a group
    - Take fog into account
- MolViewSpec extension: Initial support for customizable representation parameters
- Quick Styles section reorganized
- UI color improvements (scrollbar contrast, toggle button hover color)
- Add `overrideWater` param for entity-id color theme
- Renames PDB-Dev to PDB-IHM and adjusts data source
- Fix vertex based themes for spheres shader
- Add volume dot representation
- Add volume-value size theme
- Sequence panel: Mark focused loci (bold+underline)
- Change modifier key behavior in Normal Mode (default = select only, Ctrl/Cmd = add to selection, Shift = extend last selected range)
- Handle Firefox's limit on vertex ids per draw (#1116)
- Fix behavior of `Vec3.makeRotation(out, a, b)` when `a ≈ -b`

## [v4.10.0] - 2024-12-15

- Add `ModelWithCoordinates` decorator transform.
- Fix outlines on transparent background using illumination mode (#1364)
- Fix transparent depth texture artifacts using illumination mode
- Fix marking of consecutive gap elements (#876)
- Allow React 19 in dependencies
- Fix missing deflate header if `CompressionStream` is available
- Fix is_iOS check for NodeJS
- Added PluginCommands.Camera.FocusObject
- Plugin state snapshot can have instructions to focus objects (PluginState.Snapshot.camera.focus)
- MolViewSpec extension: Support for multi-state files (animations)
- Fix units transform data not fully updated when structure child changes
- Fix `addIndexPairBonds` quadratic runtime case
- Use adjoint matrix to transform normals in shaders
- Fix resize handling in `tests/browser`

## [v4.9.1] - 2024-12-05

- Fix iOS check when running on Node

## [v4.9.0] - 2024-12-01

- Fix artifacts when using xray shading with high xrayEdgeFalloff values
- Enable double rounded capping on tubular helices
- Fix single residue tubular helices not showing up
- Fix outlines on volume and surface reps that do not disappear (#1326)
- Add example `glb-export`
- Membrane orientation: Improve `isApplicable` check and error handling (#1316)
- Fix set fenceSync to null after deleteSync.
- Fix operator key-based `IndexPairBonds` assignment
    - Don't add bonds twice
    - Add `IndexPairs.bySameOperator` to avoid looping over all bonds for each unit
- Add `Structure.intraUnitBondMapping`
- Add more structure-based visuals to avoid too many (small) render-objects
    - `structure-intra-bond`, `structure-ellipsoid-mesh`, `structure-element-point`, `structure-element-cross`
- Upgrade to express v5 (#1311)
- Fix occupancy check using wrong index for inter-unit bond computation (@rxht, #1321)
- Fix transparent SSAO for image rendering, e.g., volumne slices (#1332)
- Fix bonds not shown with `ignoreHydrogens` on (#1315)
    - Better handle mmCIF files with no entities defined by using `label_asym_id`
    - Show bonds in water chains when `ignoreHydorgensVariant` is `non-polar`
- Add MembraneServer API, generating data to be consumed in the context of MolViewSpec
- Fix `StructConn.isExhaustive` for partial models (e.g., returned by the model server)
- Refactor value swapping in molstar-math to fix SWC (Next.js) build (#1345)
- Fix transform data not updated when structure child changes
- Fix `PluginStateSnapshotManager.syncCurrent` to work as expected on re-loaded states.
- Fix do not compute implicit hydrogens when unit is explicitly protonated (#1257)
- ModelServer and VolumeServer: support for input files from Google Cloud Storage (gs://)
- Fix color of missing partial charges for SB partial charges extension

## [v4.8.0] - 2024-10-27

- Add SSAO support for transparent geometry
- Fix SSAO color not updating
- Improve blending of overlapping outlines from transparent & opaque geometries
- Default to `blended` transparency on iOS due to `wboit` not being supported.
- Fix direct-volume with fog off (and on with `dpoit`) and transparent background on (#1286)
- Fix missing pre-multiplied alpha for `blended` & `wboit` with no fog (#1284)
- Fix backfaces visible using blended transparency on impostors (#1285)
- Fix StructureElement.Loci.isSubset() only considers common units (#1292)
- Fix `Scene.opacityAverage` calculation never 1
- Fix bloom in illumination mode
- Fix `findPredecessorIndex` bug when repeating values
- MolViewSpec: Support for transparency and custom properties
- MolViewSpec: MVP Support for geometrical primitives (mesh, lines, line, label, distance measurement)
- Mesoscale Explorer: Add support for 4-character PDB IDs (e.g., 8ZZC) in PDB-IHM/PDB-Dev loader
- Fix Sequence View in Safari 18
- Improve performance of `IndexPairBonds` assignment when operator keys are available
- ModelArchive QualityAssessment extension:
    - Add support for ma_qa_metric_local_pairwise mmCIF category
    - Add PAE plot component
- Add new AlphaFoldDB-PAE example app
- Add support for LAMMPS data and dump formats
- Remove extra anti-aliasing from text shader (fixes #1208 & #1306)

## [v4.7.1] - 2024-09-30

- Improve `resolutionMode` (#1279)
    - Add `auto` that picks `scaled` for mobile devices and `native` elsewhere
    - Add `resolution-mode` Viewer GET param
    - Add `PluginConfig.General.ResolutionMode` config item

## [v4.7.0] - 2024-09-29

- Add illumination mode
    - Path-traced SSGI
    - Automatic thickness (estimate)
        - Base thickness as max(backface depth) - min(frontface depth)
        - Per object density factor to adjust thickness
    - Progressively trace samples to keep viewport interactive
    - Toggle on/off by pressing "G"
    - `illumination` Viewer GET param
- Enables dXrayShaded define when rendering depth
- Fix handling of PDB files that have chains with same id separated by TER record (#1245)
- Sequence Panel: Improve visuals of unmodeled sequence positions (#1248)
- Fix no-compression xtc parser (#1258)
- Mol2 Reader: Fix mol2 status_bit read error (#1251)
- Fix shadows with multiple lights
- Fix impostor sphere interior normal when using orthographic projection
- Add `resolutionMode` parameter to `Canvas3DContext`
    - `scaled`, divides by `devicePixelRatio`
    - `native`, no changes
- Add `CustomProperty.Context.errorContext` to support reporting errors during loading of custom properties (#1254)
    - Use in MolViewSpec extension
- Mesoscale Explorer: fix color & style issues
- Remove use of deprecated SASS explicit color functions
- Allow "Components" section to display nested components created by "Apply Action > Selection".

## [v4.6.0] - 2024-08-28

- Add round-caps option on tubular alpha helices
- Fix missing Sequence UI update on state object removal (#1219)
- Improved prmtop format support (CTITLE, %COMMENT)
- Avoid calculating bonds for water units when `ignoreHydrogens` is on
- Add `Water` trait to `Unit`
- Improve entity-id coloring for structures with multiple models from the same source (#1221)
- Wrap screenshot & image generation in a `Task`
- AlphaFold DB: Add BinaryCIF support when fetching data
- PDB-IHM/PDB-Dev: Add support for 4-character PDB IDs (e.g., 8ZZC)
- Fix polymer-gap visual coloring with cartoon theme
- Add formal-charge color theme (#328)
- Add more coloring options to cartoon theme
- Use `CompressionStream` Browser API when available
- Add `pdbx_structure_determination_methodology` mmcif field and `Model` helpers
- Fix cartoon representation not updated when secondary structure changes
- Add Zhang-Skolnick secondary-structure assignment method which handles coarse-grained models (#49)
- Calculate bonds for coarse-grained models
- VolumeServer: Add `health-check` endpoint + `healthCheckPath` config prop to report service health
- ModelServer: Add `health-check` endpoint + `healthCheckPath` config prop to report service health

## [v4.5.0] - 2024-07-28

- Separated postprocessing passes
- Take into account explicit hydrogens when computing hydrogen bonds
- Fix DoF with pixel ratios =! 1
- Fix DoF missing transparent depth
- Fix trackball pinch zoom and add pan
- Fix aromatic link rendering when `adjustCylinderLength` is true
- Change trackball animate spin speed unit to radians per second
- Fix `mol-plugin-ui/skin/base/components/misc.scss` syntax to be in line with latest Sass syntax
- Handle missing theme updates
    - Fix trajectory-index color-theme not always updated (#896)
    - Fix bond cylinders not updated on size-theme change with `adjustCylinderLength` enabled (#1215)
- Use `OES_texture_float_linear` for SSAO when available

## [v4.4.1] - 2024-06-30

- Clean `solidInterior` transparent cylinders
- Create a transformer to deflate compressed data
- Adjust Quick Styles panel button labels
- Improve camera interpolation code (interpolate camera rotation instead of just position)
- Mesoscale Explorer
    - Add `illustrative` coloring option
    - Press 'C' to toggle between center and zoom & center on click
    - Add entities selection description
    - Clicking a leaf node in the right panel tree will center each instance in turn
    - Add measurement controls to right panel
    - Mouse left click on label with snapshot key will load snapshot
    - Mouse hover over label with protein name highlight entities with the same name
    - Custom ViewportSnapshotDescription with custom MarkdowAnchor
        - \# other snapshots with a given key \[...](#key)
        - i highlight a protein with a given NAME \[...](iNAME)
        - g highlight a group with a given group type and group name \[...](ggrouptype.groupname)
        - h URLs with a given link \[...](http...)
    - Snapshot description panel window size and text can be resized and hidden with new icons
    - Add styles controls to right panel
    - Add viewport settings to left panel
    - Add app info component to left panel with interactive tour and doc link
- Fixes SSAO edge artifacts (#1122)
    - Add `reuseOcclusion` parameter to multi-sample pass
    - Add `blurDepthBias` parameter to occlusion pass
    - Handle near clip in SSAO blur
- Support reading score from B-factor in pLDDT color theme
- Add Cel-shading support
    - `celShaded` geometry parameter
    - `celSteps` renderer parameter
- Add the ability to customize the Snapshot Description component via `PluginUISpec.components.viewport.snapshotDescription`
- Add `doNotDisposeCanvas3DContext` option to `PluginContext.dispose`
- Remove support for density data from edmaps.rcsb.org

## [v4.3.0] - 2024-05-26

- Fix State Snapshots export animation (#1140)
- Add depth of field (dof) postprocessing effect
- Add `SbNcbrTunnels` extension for for visualizing tunnels in molecular structures from ChannelsDB (more info in [tunnels.md](./docs/docs/extensions/tunnels.md))
- Fix edge case in minimizing RMSD transform computation

## [v4.2.0] - 2024-05-04

- Add emissive material support
- Add bloom post-processing
- MolViewSpec extension: `loadMVS` supports `keepCamera` parameter
- Return StateTransform selectors from measurements API (addDistance, addAngle, etc.)
- Refactor transparency rendering
    - More uniform behavior for blended, wboit, dpoit
    - Fix issues with text & image geometry
- Fix render-spheres example (#1100)
    - Wrong step size in sphere geometry boundingSphere & groupmapping
    - Handle empty `instanceGrid` in renderer & renderable
- Fix bond assignment from `IndexPairBonds`
    - Can not always be cached in `ElementSetIntraBondCache`
    - Wrong operator checks in `findPairBonds`
- Fix SSAO artifacts (@corredD, #1082)
- Fix bumpiness artifacts (#1107, #1084)

## [v4.1.0] - 2024-03-31

- Add `VolumeTransform` to translate/rotate a volume like in a structure superposition
- Fix BinaryCIF encoder edge cases caused by re-encoding an existing BinaryCIF file
- Fix edge-case where width/height in InputObserver are not correct
- Fix transparency rendering fallback (#1058)
- Fix SSAO broken when `OES_texture_float_linear` is unavailable
- Add `normalOffset` to `external-volume` color theme
    - This can give results similar to pymol's surface_ramp_above_mode=1
- Add `rotation` parameter to skybox background

## [v4.0.1] - 2024-02-19

- Fix BinaryCIF decoder edge cases. Fixes mmCIF model export from data provided by ModelServer.
- MolViewSpec extension: support for MVSX file format
- Revert "require WEBGL_depth_texture extension" & "remove renderbuffer use"

## [v4.0.0] - 2024-02-04

- Add Mesoscale Explorer app for investigating large systems
- [Breaking] Remove `cellpack` extension (superseded by Mesoscale Explorer app)
- [Breaking] Set minimal node.js version to 18
- [Breaking] Generalize rcsb/assembly-symmetry/ extension
    - Move to assembly-symmetry/
    - Remove RCSB specific dependencies and prefixes
- [Breaking] Require `WEBGL_depth_texture` webgl extension
    - Remove `renderbuffer` use
- [Breaking] Change build target to ES2018
    - Custom builds only require ES6 for dependencies like immer.js
- [Breaking] Changed `createPluginUI`
    - The function now takes a single `options` argument
    - The caller must specify a `render` method that mounts the Mol* react component to DOM
        - A default `renderReact18` method is provided, but needs to be imported separately
        - To support React 16 and 17, `ReactDOM.render` can be passed
- Improve `SetUtils` performance using ES6 features
- [Breaking] Reduce memory usage of `SymmetryOperator.ArrayMapping`
    - Requires calling methods from instance
- [Breaking] Fix `mol-model/structure/model/properties/seconday-structure.ts` file name (#938)
- [Breaking] Add `Canvas3DContext` runtime props
    - Props: pixelScale, pickScale, transparency (blended, wboit, dpoit)
    - Replaces instantiation-time attribs
- [Breaking] Change default compile target to ES2018
- [Breaking] Add culling & LOD support
    - Cull per-object and per-instance
    - Cull based on frustum and camera distance
    - LOD visibility based on camera distance
    - Special LOD mode for spheres with automatic levels
    - Occlusion culling (only WebGL2)
        - Hi-Z pass
        - Cull based on previous frame's Hi-Z buffer
- Add stochastic/dithered transparency to fade overlapping LODs in and out
- Add "Automatic Detail" preset that shows surface/cartoon/ball & stick based on camera distance

## [v3.45.0] - 2024-02-03

- Add color interpolation to impostor cylinders
- MolViewSpec components are applicable only when the model has been loaded from MolViewSpec
- Add `snapshotKey` and `tooltip` params to loci `LabelRepresentation`
- Update `FocusLoci` behavior to support `snapshotKey` param
  - Clicking a visual with `snapshotKey` will trigger that snapshot
- Render multiline loci label tooltips as Markdown
- `ParamDefinition.Text` updates:
  - Support `multiline` inputs
  - Support `placeholder` parameter
  - Support `disableInteractiveUpdates` to only trigger updates once the control loses focus
- Move dependencies related to the headless context from optional deps to optional peer deps

## [v3.44.0] - 2024-01-06

- Add new `cartoon` visuals to support atomic nucleotide base with sugar
- Add `thicknessFactor` to `cartoon` representation for scaling nucleotide block/ring/atomic-fill visuals
- Use bonds from `_struct_conn` in mmCIF files that use `label_seq_id`
- Fix measurement label `offsetZ` default: not needed when `scaleByRadius` is enbaled
- Support for label rendering in HeadlessPluginContext
- MolViewSpec extension
  - Support all X11 colors
  - Support relative URIs
  - CLI tools: mvs-validate, mvs-render, mvs-print-schema
  - Labels applied in one node
- ModelServer SDF/MOL2 ligand export: fix atom indices when additional atoms are present
- Avoid showing (and calculating) inter-unit bonds for huge structures
- Fixed `DragOverlay` on WebKit/Safari browsers

## [v3.43.1] - 2023-12-04

- Fix `react-markdown` dependency

## [v3.43.0] - 2023-12-02

- Fix `State.tryGetCellData` (return type & data check)
- Don't change camera.target unless flyMode or pointerLock are enabled
- Handle empty CIF files
- Snapshot improvements:
    - Add `key` property
    - Ability to existing snapshot name, key, and description
    - Support markdown in descriptions (ignores all HTML tags)
    - Ability to link to snapshots by key from descriptions
    - Separate UI control showing description of the current snapshot
- Do not activate drag overlay for non-file content
- Add `structure-element-sphere` visual to `spacefill` representation
- Fix missing `await` in `HeadlessPluginContext.saveStateSnapshot`
- Added support for providing custom sequence viewers to the plugin spec
- MolViewSpec extension (MVS)
- Add URL parameters `mvs-url`, `mvs-data`, `mvs-format`
- Add drag&drop for `.mvsj` files
- Fix `bumpiness` scaling with `ignoreLight` enabled
- Add `transforms` & `label` params to `ShapeFromPly`
- Optimize `LociSelectManager.selectOnly` to avoid superfluous loci set operations
- Dispose of viewer on `unload` event to aid GC

## [v3.42.0] - 2023-11-05

- Fix handling of PDB files with insertion codes (#945)
- Fix de-/saturate of colors with no hue
- Improve `distinctColors` function
    - Add `sort` and `sampleCountFactor` parameters
    - Fix clustering issues
- Add `clipPrimitive` option to spheres geometry, clipping whole spheres instead of cutting them
- Add `DragAndDropManager`
- Add `options` support for default bond labels

## [v3.41.0] - 2023-10-15

- Add `PluginContext.initialized` promise & support for it in the `Plugin` UI component.
- Fix undesired interaction between settings panel and the panel on the right.
- Add ability to customize server parameters for `RCSBAssemblySymmetry`.

## [v3.40.1] - 2023-09-30

- Do not call `updateFocusRepr` if default `StructureFocusRepresentation` isn't present.
- Treat "tap" as a click in `InputObserver`
- ModelServer ligand queries: fix atom count reported by SDF/MOL/MOL2 export
- CCD extension: Make visuals for aromatic bonds configurable
- Add optional `file?: CifFile` to `MmcifFormat.data`
- Add support for webgl extensions
    - `WEBGL_clip_cull_distance`
    - `EXT_conservative_depth`
    - `WEBGL_stencil_texturing`
    - `EXT_clip_control`
- Add `MultiSampleParams.reduceFlicker` (to be able to switch it off)
- Add `alphaThickness` parameter to adjust alpha of spheres for radius
- Ability to hide "right" panel from simplified viewport controls
- Add `blockIndex` parameter to TrajectoryFromMmCif
- Fix bounding sphere calculation for "element-like" visuals
- Fix RCSB PDB validation report URL
- Add sharpening postprocessing option
- Take pixel-ratio into account for outline scale
- Gracefully handle missing HTMLImageElement
- Fix pixel-ratio changes not applied to all render passes

## [v3.39.0] - 2023-09-02

- Add some elements support for `guessElementSymbolString` function
- Faster bounding rectangle calculation for imposter spheres
- Allow toggling of hydrogens as part of `LabelTextVisual`

## [v3.38.3] - 2023-07-29

- Fix imposter spheres not updating, e.g. in trajectories (broke in v3.38.0)

## [v3.38.2] - 2023-07-24

- Don't rely solely on `chem_comp_atom` when detecting CCD files (#877)
- Actually support non-physical keys in `Bindings.Trigger.code`

## [v3.38.1] - 2023-07-22

- Fix pixel-scale not updated in SSAO pass

## [v3.38.0] - 2023-07-18

- Fix display issue with SIFTS mapping
- Support non-physical keys in `Bindings.Trigger.code`
- Update `getStateSnapshot` to only overwrite current snapshot if it was created automatically
- Fix distinct palette's `getSamples` infinite loop
- Add 'NH2', 'FOR', 'FMT' to `CommonProteinCaps`
- Add `opened` event to `PluginStateSnapshotManager`
- Properly switch-off fog
- Add `approximate` option for spheres rendering
- Reduce `Spheres` memory usage
    - Derive mapping from VertexID
    - Pull position and group from texture
- Add `Euler` math primitive
- Add stride option to element sphere & point visuals
- Add `disabledExtensions` field to default viewer's options
- Add `LRUCache.remove`
- Add 'Chain Instance' and 'Uniform' options for 'Carbon Color' param (in Color Theme: Element Symbol)

## [v3.37.1] - 2023-06-20

- Fix issues with wboit/dpoit in large scenes
- Fix lines, text, points rendering (broken in v3.37.0)

## [v3.37.0] - 2023-06-17

- Add `inverted` option to `xrayShaded` parameter
- Model-export extension: Add ability to set a file name for structures
- Add `contextHash` to `SizeTheme`
- Add mipmap-based blur for image backgrounds

## [v3.36.1] - 2023-06-11

- Allow parsing of CCD ligand files
- Add dedicated wwPDB CCD extension to align and visualize ideal & model CCD coordinates
- Make operators in `IndexPairBonds` a directed property
- Remove erroneous bounding-box overlap test in `Structure.eachUnitPair`
- Fix `EdgeBuilder.addNextEdge` for loop edges
- Optimize inter unit bond compute
- Ensure consistent state for volume representation (#210)
- Improve SSAO for thin geometry (e.g. lines)
- Add snapshot support for structure selections
- Add `nucleicProfile` parameter to cartoon representation
- Add `cartoon` theme with separate colorings for for mainchain and sidechain visuals

## [v3.35.0] - 2023-05-14

- Enable odd dash count (1,3,5)
- Add principal axes spec and fix edge cases
- Add a uniform color theme for NtC tube that still paints residue and segment dividers in a different color
- Mesh exporter improvements
    - Support points & lines in glTF export
    - Set alphaMode and doubleSided in glTF export
    - Fix flipped cylinder caps
- Fix bond assignments `struct_conn` records referencing waters
- Add StructConn extension providing functions for inspecting struct_conns
- Fix `PluginState.setSnapshot` triggering unnecessary state updates
- Fix an edge case in the `mol-state`'s `State` when trying to apply a transform to an existing Null object
- Add `SbNcbrPartialCharges` extension for coloring and labeling atoms and residues by partial atomic charges
  - uses custom mmcif categories `_sb_ncbr_partial_atomic_charges_meta` and `_sb_ncbr_partial_atomic_charges` (more info in [README.md](./src/extensions/sb-ncbr/README.md))
- Parse HEADER record when reading PDB file
- Support `ignoreHydrogens` in interactions representation
- Add hydroxyproline (HYP) commonly present in collagen molecules to the list of amino acids
- Fix assemblies for Archive PDB files (do not generate unique `label_asym_id` if `REMARK 350` is present)
- Add additional functions to `core.math` in `mol-script`
    - `cantorPairing`, `sortedCantorPairing`, `invertCantorPairing`,
    - `trunc`, `sign`

## [v3.34.0] - 2023-04-16

- Avoid `renderMarkingDepth` for fully transparent renderables
- Remove `camera.far` doubling workaround
- Add `ModifiersKeys.areNone` helper function
- Do not render NtC tube segments unless all required atoms are present in the structure
- Fix rendering issues caused by VAO reuse
- Add "Zoom All", "Orient Axes", "Reset Axes" buttons to the "Reset Camera" button
- Improve trackball move-state handling when key bindings use modifiers
- Fix rendering with very small viewport and SSAO enabled
- Fix `.getAllLoci` for structure representations with `structure.child`
- Fix `readAllLinesAsync` refering to dom length property
- Make mol-util/file-info node compatible
- Add `eachLocation` to representation/visual interface

## [v3.33.0] - 2023-04-02

- Handle resizes of viewer element even when window remains the same size
- Throttle canvas resize events
- Selection toggle buttons hidden if selection mode is off
- Camera focus loci bindings allow reset on click-away to be overridden
- Input/controls improvements
    - Move or fly around the scene using keys
    - Pointer lock to look around scene
    - Toggle spin/rock animation using keys
- Apply bumpiness as lightness variation with `ignoreLight`
- Remove `JSX` reference from `loci-labels.ts`
- Fix overpaint/transparency/substance smoothing not updated when geometry changes
- Fix camera project/unproject when using offset viewport
- Add support for loading all blocks from a mmcif file as a trajectory
- Add `Frustum3D` and `Plane3D` math primitives
- Include `occupancy` and `B_iso_or_equiv` when creating `Conformation` from `Model`
- Remove LazyImports (introduced in v3.31.1)

## [v3.32.0] - 2023-03-20

- Avoid rendering of fully transparent renderables
- Add occlusion color parameter
- Fix issue with outlines and orthographic camera
- Reduce over-blurring occlusion at larger view distances
- Fix occlusion artefact with non-canvas viewport and pixel-ratio > 1
- Update nodejs-shims conditionals to handle polyfilled document object in NodeJS environment.
- Ensure marking edges are at least one pixel wide
- Add exposure parameter to renderer
- Only trigger marking when mouse is directly over canvas
- Fix blurry occlusion in screenshots
- [Breaking] Add `setFSModule` to `mol-util/data-source` instead of trying to trick WebPack

## [v3.31.4] - 2023-02-24

- Allow link cylinder/line `dashCount` set to '0'
- Stop animation loop when disposing `PluginContext` (thanks @gfrn for identifying the issue)

## [v3.31.3] - 2023-02-22

- Fix impostor bond visuals not correctly updating on `sizeFactor` changes
- Fix degenerate case in PCA
- Fix near clipping avoidance in impostor shaders
- Update `fs` import in `data-source.ts`

## [v3.31.2] - 2023-02-12

- Fix exit code of volume pack executable (pack.ts). Now exits with non-0 status when an error happens
- Remove pca transform from components ui focus (too distracting)
- Fix artefacts with opaque outlines behind transparent objects
- Fix polymer trace visual not updating
- Fix use of `WEBGL_provoking_vertex`

## [v3.31.1] - 2023-02-05

- Improve Component camera focus based on the PCA of the structure and the following rules:
    - The first residue should be in first quadrant if there is only one chain
    - The average position of the residues of the first chain should be in the first quadrant if there is more than one chain
- Add `HeadlessPluginContext` and `HeadlessScreenshotHelper` to be used in Node.js
- Add example `image-renderer`
- Fix wrong offset when rendering text with orthographic projection
- Update camera/handle helper when `devicePixelRatio` changes
- Add various options to customize the axes camera-helper
- Fix issue with texture-mesh color smoothing when changing themes
- Add fast boundary helper and corresponding unit trait
- Add Observable for Canvas3D commits

## [v3.30.0] - 2023-01-29

- Improve `Dnatco` extension
    - Factor out common code in `Dnatco` extension
    - Add `NtC tube` visual. Applicable for structures with NtC annotation
    - [Breaking] Rename `DnatcoConfalPyramids` to `DnatcoNtCs`
- Improve boundary calculation performance
- Add option to create & include images in state snapshots
- Fix SSAO artefacts with high bias values
- Fix SSAO resolution scale parameter handling
- Improve outlines, visually more stable at different view distances

## [v3.29.0] - 2023-01-15

- `meshes` extension: Fixed a bug in mesh visualization (show backfaces when opacity < 1)
- Add color quick select control to Volume controls
- Fix `dropFiles` bug
- Fix some cyclic imports and reduce the use of const enums. This should make it easier to use the library with the `isolatedModules: true` TS config.
- Fix `dropFiles` bug (#679)
- Add `input type='color'` picker to `CombinedColorControl`
- Set `ParameterMappingControl` disabled when state is updating
- Performance tweaks
    - Update clip `defines` only when changed
    - Check for identity in structure/unit areEqual methods
    - Avoid cloning of structure representation parameters
    - Make SymmetryOperator.createMapping monomorphic
    - Improve bonding-sphere calculation
    - Defer Scene properties calculation (markerAverage, opacityAverage, hasOpaque)
    - Improve checks in in UnitsRepresentation setVisualState
- Add StructureElement.Loci.forEachLocation
- Add RepresentationRegistry.clear and ThemeRegistry.clear
- Add generic Loci support for overpaint, substance, clipping themes
- Add `.getCenter` and `.center` to `Camera`
- Add support to dim unmarked groups
- Add support for marker edge strength

## [v3.28.0] - 2022-12-20

- Show histogram in direct volume control point settings
- Add `solidInterior` parameter to sphere/cylinder impostors
- [Breaking] Tweak `ignoreHydrogens` non-polar handling (introduced in 3.27.0)
- Add `meshes` and `volumes-and-segmentations` extensions
    - See https://molstarvolseg.ncbr.muni.cz/ for more info
- Fix missing support for info in `ParamDefinition.Converted`
- Add support for multi-visual volume representations
- Improve volume isosurface bounding-sphere
- Add basic volume segmentation support to core
    - Add `Volume.Segment` model
    - Add `Segmentation` custom volume property
    - Add `SegmentRepresentation` representation
    - Add `volume-segment` color theme
- Fix GPU marching cubes failing for large meshes with webgl2 (due to use of float16)

## [v3.27.0] - 2022-12-15

- Add an `includeTransparent` parameter to hide/show outlines of components that are transparent
- Fix 'once' for animations of systems with many frames
- Better guard against issue (black fringes) with bumpiness in impostors
- Improve impostor shaders
    - Fix sphere near-clipping with orthographic projection
    - Fix cylinder near-clipping
    - Add interior cylinder caps
    - Add per-pixel object clipping
- Fix `QualityAssessment` assignment bug for structures with different auth vs label sequence numbering
- Refresh `ApplyActionControl`'s param definition when toggling expanded state
- Fix `struct_conn` bond assignment for ions
- Ability to show only polar hydrogens

## [v3.26.0] - 2022-12-04

- Support for ``powerPreference`` webgl attribute. Add ``PluginConfig.General.PowerPreference`` and ``power-preference`` Viewer GET param.
- Excluded common protein caps `NME` and `ACE` from the ligand selection query
- Add screen-space shadow post-processing effect
- Add "Structure Molecular Surface" visual
- Add `external-volume` theme (coloring of arbitrary geometries by user-selected volume)

## [v3.25.1] - 2022-11-20

- Fix edge-case in `Structure.eachUnitPair` with single-element units
- Fix 'auto' structure-quality for coarse models

## [v3.25.0] - 2022-11-16

- Fix handling of gzipped assets (reverts #615)

## [v3.24.0] - 2022-11-13

- Make `PluginContext.initContainer` checkered canvas background optional
- Store URL of downloaded assets to detect zip/gzip based on extension (#615)
- Add optional `operator.key`; can be referenced in `IndexPairBonds`
- Add overpaint/transparency/substance theme strength to representations
- Fix viewport color for transparent background

## [v3.23.0] - 2022-10-19

- Add `PluginContext.initContainer/mount/unmount` methods; these should make it easier to reuse a plugin context with both custom and built-in UI
- Add `PluginContext.canvas3dInitialized`
- `createPluginUI` now resolves after the 3d canvas has been initialized
- Change EM Volume Streaming default from `Whole Structure` to `Auto`

## [v3.22.0] - 2022-10-17

- Replace `VolumeIsosurfaceParams.pickingGranularity` param with `Volume.PickingGranuality`

## [v3.21.0] - 2022-10-17

- Add `VolumeIsosurfaceParams.pickingGranularity` param
- Prevent component controls collapsing when option is selected

## [v3.20.0] - 2022-10-16

- [Breaking] Rename the ``model-index`` color theme to ``trajectory-index``
- Add a new ``model-index`` color theme that uniquely colors each loaded model
- Add the new ``model-index`` and ``structure-index`` color themes as an option for the carbon color in the ``element-symbol`` and ``ilustrative`` color themes
- Add ``structure-index`` color theme that uniquely colors each root structure
- Add ``nearest`` method to ``Lookup3D``
- Add mipmap-based blur for skybox backgrounds

## [v3.19.0] - 2022-10-01

- Fix "empty textures" error on empty canvas
- Optimize BinaryCIF integer packing encoder
- Fix dual depth peeling when post-processing is off or when rendering direct-volumes
- Add ``cameraClipping.minNear`` parameter
- Fix black artifacts on specular highlights with transparent background

## [v3.18.0] - 2022-09-17

- Integration of Dual depth peeling - OIT method
- Stereo camera improvements
    - Fix param updates not applied
    - Better param ranges and description
    - Add timer.mark for left/right camera

## [v3.17.0] - 2022-09-11

- [Fix] Clone ``Canvas3DParams`` when creating a ``Canvas3D`` instance to prevent shared state between multiple instances
- Add ``includeResidueTest`` option to ``alignAndSuperposeWithSIFTSMapping``
- Add ``parentDisplay`` param for interactions representation.
- [Experimental] Add support for PyMOL, VMD, and Jmol atom expressions in selection scripts
- Support for ``failIfMajorPerformanceCaveat`` webgl attribute. Add ``PluginConfig.General.AllowMajorPerformanceCaveat`` and ``allow-major-performance-caveat`` Viewer GET param.
- Fix handling of PDB TER records (#549)
- Add support for getting multiple loci from a representation (``.getAllLoci()``)
- Add ``key`` property to intra- and inter-bonds for referencing source data
- Fix click event triggered after move

## [v3.16.0] - 2022-08-25

- Support ``globalColorParams`` and ``globalSymmetryParams`` in common representation params
- Support ``label`` parameter in ``Viewer.loadStructureFromUrl``
- Fix ``ViewportHelpContent`` Mouse Controls section

## [v3.15.0] - 2022-08-23

- Fix wboit in Safari >=15 (add missing depth renderbuffer to wboit pass)
- Add 'Around Camera' option to Volume streaming
- Avoid queuing more than one update in Volume streaming

## [v3.14.0] - 2022-08-20

- Expose inter-bonds compute params in structure
- Improve performance of inter/intra-bonds compute
- Fix defaultAttribs handling in Canvas3DContext.fromCanvas
- Confal pyramids extension improvements
    - Add custom labels to Confal pyramids
    - Improve naming of some internal types in Confal pyramids extension coordinate
    - Add example mmCIF file with categories necessary to display Confal pyramids
    - Change the lookup logic of NtC steps from residues
- Add support for download of gzipped files
- Don't filter IndexPairBonds by element-based rules in MOL/SDF and MOL2 (without symmetry) models
- Fix Glycam Saccharide Names used by default
- Fix GPU surfaces rendering in Safari with WebGL2
- Add ``fov`` (Field of View) Canvas3D parameter
- Add ``sceneRadiusFactor`` Canvas3D parameter
- Add background pass (skybox, image, horizontal/radial gradient)
    - Set simple-settings presets via ``PluginConfig.Background.Styles``
    - Example presets in new backgrounds extension
    - Load skybox/image from URL or File (saved in session)
    - Opacity, saturation, lightness controls for skybox/image
    - Coverage (viewport or canvas) controls for image/gradient
- [Breaking] ``AssetManager`` needs to be passed to various graphics related classes
- Fix SSAO renderable initialization
- Reduce number of webgl state changes
    - Add ``viewport`` and ``scissor`` to state object
    - Add ``hasOpaque`` to scene object
- Handle edge cases where some renderables would not get (correctly) rendered
    - Fix text background rendering for opaque text
    - Fix helper scenes not shown when rendering directly to draw target
- Fix ``CustomElementProperty`` coloring not working

## [v3.13.0] - 2022-07-24

- Fix: only update camera state if manualReset is off (#494)
- Improve handling principal axes of points in a plane
- Add 'material' annotation support for textures
- More effort to avoid using ``flat`` qualifier in shaders: add ``dVaryingGroup``
- Enable ``immediateUpdate`` for iso level in isosurface and volume streaming controls
- Add support to download CCD from configurable URL

## [v3.12.1] - 2022-07-20

- Fix plugin behavior dispose logic to correctly unsubscribe observables.

## [v3.12.0] - 2022-07-17

- Add ``colorMarker`` option to Renderer. This disables the highlight and select marker at a shader level for faster rendering of large scenes in some cases.
- Bind shared textures only once per pass, not for each render item
- Fix missing 'material' annotation for some uniforms, causing unnecessary uniform updates
- Remove use of ``isnan`` in impostor shaders, not needed and causing slowdown
- Avoid using ``flat`` qualifier in shaders, causing slowdown
- Improve CellPack's ``adjustStyle`` option (disable ``colorMarker``, set component options, enable marking w/o ghost)
- Scan all entities when looking for ``struct_conn`` entries (fixes issue when the same ``label_asym_id`` is used in more than one entity)

## [v3.11.0] - 2022-07-04

- Add ``instanceGranularity`` option for marker, transparency, clipping, overpaint, substance data to save memory
- CellPack extension tweaks
    - Use instancing to create DNA/RNA curves to save memory
    - Enable ``instanceGranularity`` by default
    - Add ``adjustStyle`` option to LoadCellPackModel action (stylized, no multi-sample, no far clipping, chain picking)
- Structure Superposition now respects pivot's coordinate system

## [v3.10.2] - 2022-06-26

- Fix superfluous shader varying
- Improve use of gl_VertexID when possible

## [v3.10.1] - 2022-06-26

- Fix groupCount when updating TextureMesh-based visuals

## [v3.10.0] - 2022-06-24

- Add support for Glycam saccharide names
- Add ``PluginConfig.Viewport.ShowTrajectoryControls`` config option

## [v3.9.1] - 2022-06-19

- Fix missing ``super.componentWillUnmount()`` calls (@simeonborko)
- Fix missing ``uGroupCount`` update for visuals
- Fix missing aromatic bond display

## [v3.9.0] - 2022-05-30

- Improve picking by using drawbuffers (when available) to reduce number of drawcalls
- GPU timing support
    - Add ``timing-mode`` Viewer GET param
    - Add support for webgl timer queries
    - Add timer marks around GPU render & compute operations
- Volume Server CIF: Add check that a data block contains volume data before parsing
- Fix ``Scene.clear`` not clearing primitives & volumes arrays (@JonStargaryen)
- Fix rendering volumes when wboit is switched off and postprocessing is enabled

## [v3.8.2] - 2022-05-22

- Fix ``Scene.opacityAverage`` not taking xray shaded into account

## [v3.8.1] - 2022-05-14

- Fix issues with marking camera/handle helper (#433)
- Fix issues with array uniforms when running with headless-gl
- Fix Polymer Chain Instance coloring
- Improve performance of scene marker/opacity average calculation

## [v3.8.0] - 2022-04-30

- Add support for outlines around transparent objects
- Improve per-group transparency when wboit is switched off
- Improve ``ColorTheme`` typing with ``ColorType`` generic.
    - Defaults to ``ColorTypeLocation``
    - Set when using ``ColorTypeDirect`` or ``ColorTypeGrid``
- Fix case handling of ``struct_conf`` mmCIF enumeration field (#425)
- Fix ``allowTransparentBackfaces`` for per-group transparency
- Fix ``FormatRegistry.isApplicable`` returning true for unregistered formats
- Fix: handle building of ``GridLookup3D`` with zero cell size
- Fix ``ignoreLight`` for direct-volume rendering with webgl1
- Fix (non-black) outlines when using transparent background

## [v3.7.0] - 2022-04-13

- Fix ``xrayShaded`` for texture-mesh geometries
- [Breaking] Change ``allowTransparentBackfaces`` to ``transparentBackfaces`` with options ``off``, ``on``, ``opaque``. This was only added in 3.6.0, so allowing a breaking change here.
    - ``off``: don't show (default)
    - ``on``: show with transparency
    - ``opaque``: show fully opaque
- Add option to disable file drop overlay.

## [v3.6.2] - 2022-04-05

- ModelServer ligand queries: fixes for alternate locations, additional atoms & UNL ligand
- React 18 friendly ``useBehavior`` hook.

## [v3.6.1] - 2022-04-03

- Fix React18 related UI regressions.

## [v3.6.0] - 2022-04-03

- Check that model and coordinates have same element count when creating a trajectory
- Fix aromatic rings assignment: do not mix flags and planarity test
- Improve bonds assignment of coarse grained models: check for IndexPairBonds and exhaustive StructConn
- Fix unit mapping in bondedAtomicPairs MolScript query
- Improve pdb parsing: handle non unique atom and chain names (fixes #156)
- Fix volume streaming for entries with multiple contour lists
- Add ``allowTransparentBackfaces`` parameter to support double-sided rendering of transparent geometries
- Fix handling of case insensitive mmCIF enumeration fields (including entity.type)
- Fix ``disable-wboit`` Viewer GET param
- Add support for React 18.
    - Used by importing ``createPluginUI`` from ``mol-plugin-ui/react18``;
    - In Mol* 4.0, React 18 will become the default option.

## [v3.5.0] - 2022-03-25

- Fix issues with bounding-sphere & color-smoothing (mostly for small geometries)
- Support BCIF => CIF conversion in ``cif2bcif`` CLI tool

## [v3.4.0] - 2022-03-13

- Fix handling of mmcif with empty ``label_*`` fields
- Improve saccharide detection (compare against list from CCD)
- Fix legend label of hydrophobicity color theme
- Add ``LoadTrajectory`` action
- Add ``CustomImportControls`` to left panel
- Add Zenodo import extension (load structures, trajectories, volumes, and zip files)
- Fix loading of some compressed files within sessions
- Fix wrong element assignment for atoms with Charmm ion names
- Fix handling of empty symmetry cell data
- Add support for ``trr`` and ``nctraj`` coordinates files
- Add support for ``prmtop`` and ``top`` topology files

## [v3.3.1] - 2022-02-27

- Fix issue with unit boundary reuse (do at visual level instead)
- Add option to ignore ions for inter-unit bond computation

## [v3.3.0] - 2022-02-27

- Fix parsing contour-level from emdb v3 header files
- Fix invalid CSS (#376)
- Fix "texture not renderable" & "texture not bound" warnings (#319)
- Fix visual for bonds between two aromatic rings
- Fix visual for delocalized bonds (parsed from mmcif and mol2)
- Fix ring computation algorithm
- Add ``UnitResonance`` property with info about delocalized triplets
- Resolve marking in main renderer loop to improve overall performance
- Use ``throttleTime`` instead of ``debounceTime`` in sequence viewer for better responsiveness
- Change line geometry default ``scaleFactor`` to 2 (3 is too big after fixing line rendering)
- Trajectory animation performance improvements
    - Reuse ``Model.CoarseGrained`` for coordinate trajectories
    - Avoid calculating ``InterUnitBonds`` when ``Structure.parent`` ones are empty
    - Reuse unit boundary if sphere has not changed too much
    - Don't show 'inter-bond' and 'element-cross' visuals in line representations of polymerAndLigand preset
- Fix additional mononucleotides detected as polymer components
- Fix and improve ``canRemap`` handling in ``IntraUnitBonds``
- Reuse occlusion for secondary passes during multi-sampling
- Check if marking passes are needed before doing them
- Add ``resolutionScale`` parameter to allow trading quality of occlusion for performance

## [v3.2.0] - 2022-02-17

- Rename "best database mapping" to "SIFTS Mapping"
- Add schema and export support for ``atom_site.pdbx_sifts_xref_*`` fields
- Add schema export support for ``atom_site.pdbx_label_index`` field
- Add `traceOnly` parameter to chain/UniProt-based structure alignment
- Store ``IndexPairBonds`` as a dynamic property.

## [v3.1.0] - 2022-02-06

- Fix ``xrayShaded`` & ``ignoreLight`` params not working at the same time
- Add ``ignoreLight`` to component params
- Tweaks for cleaner default representation style
    - Cartoon: use ``nucleotide-ring`` instead of ``nucleotide-block``
    - Focus: use ``xrayShaded`` instead of opacity; adjust target size; don't show non-covalent interactions twice
- Fix representation preset side effects (changing post-processing parameters, see #363)
- Add Quick Styles panel (default, illustrative, stylized)
- Fix exported structure missing secondary-structure categories (#364)
- Fix volume streaming error message: distinguish between missing data and server error (#364)

## [v3.0.2] - 2022-01-30

- Fix color smoothing of elongated structures (by fixing ``Sphere.expand`` for spheres with highly directional extrema)
- Fix entity label not displayed when multiple instances of the same entity are highlighted
- Fix empty elements created in ``StructureElement.Loci.extendToAllInstances``
- Measurement options tweaks (allow larger ``textSize``; make ``customText`` essential)
- Fix visual visibility sync edge case when changing state snapshots

## [v3.0.1] - 2022-01-27

- Fix marking pass not working with ``transparentBackground``
- Fix pdbe xray maps url not https
- Fix entity-id color theme broken for non-IHM models
- Improve/fix marking of ``InteractionsInterUnitVisual`` (mark when all contact-feature members are given)
- Add missing "entity-id" and "enity-source" options for carbon coloring to "element-symbol" color theme
- Fix VolumeServer/query CLI
- Support automatic iso-value adjustment for VolumeServer data in ``Viewer.loadVolumeFromUrl``
- Emit drag event whenever started within viewport (not only for non-empty loci)

## [v3.0.0] - 2022-01-23

- Assembly handling tweaks:
    - Do not include suffix for "identity assembly operators"
    - Do not include assembly-related categories to export if the structure was composed from an assembly
    - Special case for ``structAsymMap`` if Mol* asym id operator mapping is present
- Support for opening ZIP files with multiple entries
- Add Model Export extension
- Bugfix: Automatically treat empty string as "non-present" value in BinaryCIF writer.
- Fix coarse model support in entity-id color theme
- Fix marking of carbohydrate visuals (whole chain could get marked instead of single residue)
- Add custom colors to "element-symbol", "molecule-type", "residue-name", and "secondary-structure" themes
- Support/bugfixes for ``atom_site.pdbx_sifts_xref`` categories
- Improve/fix marking of ``InteractionsIntraUnitVisual`` (mark when all contact-feature members are given)

## [v3.0.0-dev.10] - 2022-01-17

- Fix ``getOperatorsForIndex``
- Pass animation info (current frame & count) to state animations
    - Fix camera stutter for "camera spin" animation
- Add formal charge parsing support for MOL/SDF files (thanks @ptourlas)
- [Breaking] Cleaner looking ``MembraneOrientationVisuals`` defaults
- [Breaking] Add rock animation to trackball controls
    - Add ``animate`` to ``TrackballControlsParams``, remove ``spin`` and ``spinSpeed``
    - Add ``animate`` to ``SimpleSettingsParams``, remove ``spin``
- Add "camera rock" state animation
- Add support for custom colors to "molecule-type" theme
- [Breaking] Add style parameter to "illustrative" color theme
    - Defaults to "entity-id" style instead of "chain-id"
- Add "illustrative" representation preset

## [v3.0.0-dev.9] - 2022-01-09

- Add PDBj as a ``pdb-provider`` option
- Move Viewer APP to a separate file to allow use without importing light theme & index.html
- Add symmetry support for mol2 files (only spacegroup setting 1)
- Fix mol2 files element symbol assignment
- Improve bond assignment from ``IndexPairBonds``
    - Add ``key`` field for mapping to source data
    - Fix assignment of bonds with unphysical length
- Fix label/stats of single atom selection in multi-chain units

## [v3.0.0-dev.8] - 2021-12-31

- Add ``PluginFeatureDetection`` and disable WBOIT in Safari 15.
- Add ``disable-wboit`` Viewer GET param
- Add ``prefer-webgl1`` Viewer GET param
- [Breaking] Refactor direct-volume rendering
    - Remove isosurface render-mode (use GPU MC instead)
    - Move coloring into theme (like for other geometries/renderables)
        - Add ``direct`` color type
        - Remove color from transfer-function (now only alpha)
        - Add direct-volume color theme support
        - Add volume-value color theme
- [Breaking] Use size theme in molecular/gaussian surface & label representations
    - This is breaking because it was hardcoded to ``physical`` internally but the repr size theme default was ``uniform`` (now ``physical``)

## [v3.0.0-dev.7] - 2021-12-20

- Reduce number of created programs/shaders
    - Support specifying variants when creating graphics render-items
    - Change double-side shader param from define to uniform
    - Remove dMarkerType shader define (use uMarker as needed)
    - Support to ignore defines depending on the shader variant
    - Combine pickObject/pickInstance/pickGroup shader variants into one
    - Combine markingDepth/markingMask shader variants into one
    - Correctly set shader define flags for overpaint, transparency, substance, clipping
- [Breaking] Add per-object clip rendering properties (variant/objects)
    - ``SimpleSettingsParams.clipping.variant/objects`` and ``RendererParams.clip`` were removed

## [v3.0.0-dev.6] - 2021-12-19

- Enable temporal multi-sampling by default
    - Fix flickering during marking with camera at rest
- Enable ``aromaticBonds`` in structure representations by default
- Add ``PluginConfig.Structure.DefaultRepresentationPreset``
- Add ModelArchive support
    - schema extensions (e.g., AlphaFold uses it for the pLDDT score)
    - ModelArchive option in DownloadStructure action
    - ``model-archive`` GET parameter for Viewer app
    - ``Viewer.loadModelArchive`` method
- Improve support for loading AlphaFold structures
    - Automatic coloring by pLDDT
    - AlphaFold DB option in DownloadStructure action
    - ``afdb`` GET parameter for Viewer app
    - ``Viewer.loadAlphaFoldDb`` method
- Add QualityAssessment extension (using data from ma_qa_metric_local mmcif category)
    - pLDDT & qmean score: coloring, repr presets, molql symbol, loci labels (including avg for mutli-residue selections)
    - pLDDT: selection query
- Warn about erroneous symmetry operator matrix (instead of throwing an error)
- Added ``createPluginUI`` to ``mol-plugin-ui``
    - Support ``onBeforeUIRender`` to make sure initial UI works with custom presets and similar features.
- [Breaking] Removed ``createPlugin`` and ``createPluginAsync`` from ``mol-plugin-ui``
    - Please use ``createPluginUI`` instead
- Improve aromatic bonds handling
    - Don't detect aromatic bonds for rings < 5 atoms based on planarity
    - Prefer atoms in aromatic rings as bond reference positions

## [v3.0.0-dev.5] - 2021-12-16

- Fix initial camera reset not triggering for some entries.

## [v3.0.0-dev.4] - 2021-12-14

- Add ``bumpiness`` (per-object and per-group), ``bumpFrequency`` & ``bumpAmplitude`` (per-object) render parameters (#299)
- Change ``label`` representation defaults: Use text border instead of rectangle background
- Add outline color option to renderer
- Fix false positives in Model.isFromPdbArchive
- Add drag and drop support for loading any file, including multiple at once
    - If there are session files (.molx or .molj) among the dropped files, only the first session will be loaded
- Add drag and drop overlay
- Safari 15.1 - 15.3 WebGL 2 support workaround
- [Breaking] Move ``react`` and ``react-dom`` to ``peerDependencies``. This might break some builds.

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
