# v6 compatibility changes

Record public API and behavior changes introduced by the refactor here. This is
an initial ledger, not a complete API audit. Remaining audit work is tracked in
[checklist.md](checklist.md).

## Accepted: classic `molstar.lib.shape`

The classic bundle still exposes `molstar.lib`, but its shape objects now expose
the model-side API following the separation of model data and graphics helpers.

- `Shape.createRenderObject`, `Shape.createTransform`, `Shape.getTheme`, and
  `Shape.groupIterator` are no longer on the classic model-side `Shape` object.
- `ShapeGroup.getBoundingSphere` is no longer on its model-side object.
- Model-side `Shape.create` requires `groupCount` as argument 7, with optional
  `transforms` as argument 8. The old seventh-argument transforms call must change.
  JavaScript callers passing the old signature can silently store the transforms
  array as `groupCount`.

Graphics integrations belong to `@molstar/graphics/geo/shape/shape`. These classic API
changes are accepted for v6; do not restore them through a compatibility shim.

## Temporary: standalone MolQL validation

The standalone builder and `mvs-validate` check MolQL expression structure but do
not compile expressions. An unknown symbol can therefore pass CLI validation.
The MVS runtime still performs compiler validation and rejects it when loading.

Full CLI validation will return after designing the mol-script import/integration.
See the deferred validation item in [checklist.md](checklist.md).

## Headless MP4 integration

`HeadlessPluginContext` now owns image rendering only. Its former `getAnimation`
and `saveAnimation` methods live on `Mp4HeadlessPluginContext`, imported from
`@molstar/mp4-export-extension/headless`. Register `Mp4Export` in the plugin spec
when using those methods. Native GL remains opt-in rather than auto-installed.

## Shape bounds and graphics helpers

Model-side shapes use a structural `ShapeGeometry` contract. Geometry-aware shape
creation and group bounds belong to `@molstar/graphics/geo/shape/shape`.
`Loci.getBoundingSphere` for shape groups delegates to the optional
`shape.getGroupBoundingSphere` callback; the graphics factory supplies it.
Shapes created directly through the model factory need that callback for group bounds.

Streamline picking and location-iterator helpers moved from model props to
`@molstar/graphics/geo/streamlines`: `getStreamlinesVisualLoci`, `getStreamlinesLoci`,
`eachStreamlines`, and `createStreamlinesLocationIterator`.

`PositionLocation` and `isPositionLocation` are owned by the model package and
re-exported by the graphics location iterator. `Viewport` is owned by core and
re-exported by the graphics camera utility. Those re-exports preserve existing
behavior at the mapped module paths.

## Color theme registration

`ColorTheme.createRegistry()` in graphics no longer includes `external-structure`
or `external-volume`. `PluginContext` adds those plugin-owned providers to its
registries, preserving default plugin behavior. Consumers constructing a bare
graphics registry must register those providers explicitly if needed.

## Declaration contracts

`ExternalModules['jpeg-js']` exposes the injected codec's `encode` contract instead
of the entire codec module type. `JpegBufferRet` is owned by headless and describes
the returned `width`, `height`, and `Buffer` data.

Standalone MVS `validateTree` accepts the logger contract
`{ log: { error(message: string): void } }` rather than requiring `PluginContext`.
Existing plugin arguments still satisfy that contract.

## Distribution and import paths

Compiled packages use ESM and new package subpaths. Old monolithic CommonJS and
source paths require migration; use `migration-map.json` for ownership mappings.
Classic Viewer/MVS Stories asset paths and globals remain available, subject to
the shape API changes above. Package exports exclude tests and build caches.

## MVSX ZIP options

`createMVSX(data, assets, options?: { zip?: ZipOptions })` accepts fflate ZIP options
under `options.zip`, including `mtime` and compression level. The default uses
the current time rather than v5's fixed
ZIP timestamps. Callers needing reproducible archives should provide `zip.mtime`;
use a local-calendar date when identical ZIP date fields across time zones matter.
