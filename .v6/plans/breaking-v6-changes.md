# v6 compatibility changes

Record public API and behavior changes introduced by the refactor here. The
current prototype audit is complete for the scope below. Update this ledger as
implementation continues; release/downstream validation remains in [checklist.md](checklist.md).

## Audit checkpoint

Compared v5 master `a1bf00ba55d737461e7e149ccc1996d979523b34` with prototype source
at `094c412a2f40692440434b5cd6eea8f1d0b77821` on 2026-10-04. Using the migration
map, compared 1,547 production TypeScript/TSX modules across packages, extensions,
apps, CLI tools, and servers. No mapped targets or upstream production source
mappings were missing. Tests/performance fixtures were outside this source scan.

A TypeScript syntax/token scan removed top-level imports, normalized export module
paths, and ignored comments/formatting. It identified 63 remaining source
differences, reviewed for exported names, namespace/class members, signatures,
and behavior. Reviewed the 17 new/extracted source modules separately. Internal
shape-factory call changes, ESM inline import suffixes, alias-only type changes,
and CLI path rewrites account for much of the remaining diff.

Checked compiled public modules and actual classic Viewer globals in local Chrome.
Node probes confirmed the MVS cache/file-URL changes and the standalone/runtime
MolQL validation distinction. This checkpoint records the API/behavior diff; it
does not certify arbitrary downstream applications, native GL, or browser platforms.
Changes already present in the merged v5 master, including rendering/marking and
BinaryCIF/PDB fixes, are not attributed to this refactor.

## Preserved classic globals

The classic Viewer still exposes `molstar.Viewer`, `ViewerAutoPreset`,
`PluginExtensions`, `ExtensionMap`, `lib`, `version`, `consoleStats`,
`isDebugMode`, `isProductionMode`, `isTimingMode`, `setDebugMode`,
`setProductionMode`, and `setTimingMode`.

`molstar.lib` retains `structure`, `volume`, `shape`, `loci`, `math`, `plugin`, and
`extensions`. The Viewer and MVS Stories entrypoint implementations retain their
public exports. MVS Stories still registers its custom elements and exports its
context/load/download functions and `MVSData`. Retaining these namespaces does
not preserve every nested API, particularly the shape changes below.

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

The graphics module exports functions directly, rather than a replacement
`Shape`/`ShapeGroup` namespace:

| v5 call | v6 graphics call |
| --- | --- |
| `Shape.create(..., transforms?, groupCount?)` | `create(..., transforms?, groupCount?)` |
| `Shape.createRenderObject(shape, props)` | `createRenderObject(shape, props)` |
| `Shape.createTransform(...)` | `createTransform(...)` |
| `Shape.getTheme(shape)` | `getTheme(shape)` |
| `Shape.groupIterator(shape)` | `groupIterator(shape)` |
| `ShapeGroup.getBoundingSphere(loci, out?)` | `getGroupBoundingSphere(loci, out?)` |

For rendered shapes, import that graphics module as a namespace and use its
`create` function; it preserves the old optional-argument order and computes the
initial group count from geometry when omitted. Model-side `Shape` defaults its
generic geometry type to `ShapeGeometry`, not the full graphics `Geometry` union.
Supply a concrete generic type when accessing geometry-specific fields.

`shape.groupCount` is now a stored value rather than a getter deriving the count
from geometry on each access. Consumers changing geometry group counts after
creation must recreate the shape with the correct count instead of relying on
the old getter. The model contract still marks `groupCount` readonly.

## Temporary: standalone MolQL validation

The standalone builder and `mvs-validate` check MolQL expression structure but do
not compile expressions. An unknown symbol can therefore pass CLI validation.
The MVS runtime still performs compiler validation and rejects it when loading.

Full CLI validation will return after designing the mol-script import/integration.
See the deferred validation item in [checklist.md](checklist.md).

This also affects `MVSData.validationIssues`/`isValid` and schema decoding of MolQL
selectors in the standalone builder. Runtime sanity checks invoke the compiler
validation for both scene and animation trees, and selector loading still compiles
expressions. A syntactically valid unknown symbol can pass standalone checks.

The standalone `@molstar/mvs-builder/expression` module exposes JSON expression
construction/shape checks only. Use `@molstar/model/script/...` for MolScript
compiler/runtime functionality. New `@molstar/mvs/behavior-id` constants keep
the existing `molviewspec` name and `ms-plugin.molviewspec` transformer ID; the
loader now checks that stable ID without importing the behavior module. Serialized
snapshot identity is unchanged.

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

## Other helper moves requiring new imports

These are symbol-level moves within formerly larger modules; the whole-file
migration map alone is insufficient for these imports.

| Public API | v5 module | v6 import |
| --- | --- | --- |
| `CollapsableControls`, `CollapsableProps`, `CollapsableState` | `mol-plugin-ui/base` | `@molstar/plugin-ui/controls/collapsable` |
| `GaussianDensityTextureData`, `computeGaussianDensityTexture`, `computeGaussianDensityTexture2d`, `computeGaussianDensityTexture3d` | `mol-math/geometry/gaussian-density` | `@molstar/graphics/geo/gaussian-density` |
| `fillIdentityTransform` | `mol-geo/geometry/transform-data` | `@molstar/core/math/linear-algebra/3d/transform-array` |

Collapsible controls are no longer exported from `@molstar/plugin-ui/base`.
Gaussian-density CPU APIs remain in `@molstar/core/math/geometry/gaussian-density`;
their GPU functions/data types moved to graphics with the existing signatures.
The GPU helper behavior is carried over from v5; this audit does not redesign it.
`fillIdentityTransform` is no longer exported by the graphics transform-data module.

`Tokens` now lives in `@molstar/core/util/token-ranges`, with the same `data`,
`count`, and `indices` fields. The IO tokenizer retains its type re-export, so
the mapped tokenizer import remains valid. `openRead` now lives in
`@molstar/common-server/open-read`; the volume server's mapped common/file module
retains its re-export and behavior.

## Temporary: external color themes unavailable by default

`ColorTheme.createRegistry()` in graphics no longer includes `external-structure`
or `external-volume`. They are also intentionally absent from the default plugin
registries for this prototype. The special registration helper in `PluginContext`
has been removed; it uses the ordinary graphics registry factory.

Implementations remain at `@molstar/plugin/themes/external-structure` and
`@molstar/plugin/themes/external-volume` because they depend on plugin state
objects and selection queries. Bring them back through the planned registry
composition work, as tracked in the checklist. Existing presets/snapshots that
request these themes cannot rely on default registration until that work is done.

## Declaration contracts

`ExternalModules['jpeg-js']` exposes the injected codec's `encode` contract instead
of the entire codec module type. `JpegBufferRet` is owned by headless and describes
the returned `width`, `height`, and `Buffer` data.

Standalone MVS `validateTree` accepts the logger contract
`{ log: { error(message: string): void } }` rather than requiring `PluginContext`.
Existing plugin arguments still satisfy that contract.

`Clip.Rotation` is a new exported alias for `{ axis: Vec3, angle: number }`.
Explicit parameter-schema annotations retain the existing rotation value shape.
The background image asset declaration changes from TypeScript `export =` to a
default export to match ESM asset imports.

Standalone MVS also exposes package-local `object`, `json`, and `color-names`
helpers. The builder's `ColorNames` values are membership flags (`true`), not
the numeric colors in core's color-name table; use core when a color value is needed.
The MVS JSON expression shape remains structurally compatible with MolScript
expressions, but these builder helpers do not provide compiler semantics.

## Distribution and import paths

Compiled packages use ESM and new package subpaths. Old monolithic CommonJS and
source paths require migration; use `migration-map.json` for ownership mappings.
Classic Viewer/MVS Stories asset paths and globals remain available, subject to
the shape API changes above. Package exports exclude tests and build caches.

`molstar` now packages only browser distributions. Library/server/CLI consumers
must install the owning `@molstar/*` packages instead of importing
`molstar/lib/...` or `molstar/lib/commonjs/...`. There is no CommonJS build or
`require` export condition. Published library modules and the browser ESM
distribution target ES2022. All 43 public packages share one release version.

Native `gl`/`canvas` remain injected by headless callers and are optional peers.
`@molstar/plugin-headless/native` adds lazy `loadNativeModule` resolution for the
MVS rendering CLI; importing it or requesting `--help` no longer eagerly loads native
modules. A rendering command still needs the caller to install them explicitly.

## Version metadata

`PLUGIN_VERSION` is generated from `version.json` rather than falling back to
`'(development)'` outside a bundle. `PLUGIN_VERSION_DATE` comes from version
synchronization and is reused for the same version unless
`MOLSTAR_BUILD_TIMESTAMP` is explicitly provided to `version:sync`/`build:lib`.
It is no longer automatically the time of every application bundle build.
The old `__MOLSTAR_PLUGIN_VERSION__`/`__MOLSTAR_BUILD_TIMESTAMP__` definitions
do not control the generated module's values.

## CLI generated imports and entrypoints

`cifschema --moldataImportPath` now defaults to `@molstar/core/data`, replacing
`molstar/lib/mol-data` in generated TypeScript. Pass the option explicitly to
generate another import path. Its field-name CSV presets read canonical root
data in the workspace and staged package assets when installed.

CLI/server entrypoints live in their owning packages and use checked-in
`bin/*.mjs` launchers. Their existing command argument parsing is retained;
repository maintenance generators such as `syminfo` target the new package/test
paths. Model-server TAR helpers now expose `decodeLongPath`, `encodePax`, and
`decodePax` as named ESM/TypeScript exports instead of CommonJS `exports` assignments.

## MVSX ZIP options

`createMVSX(data, assets, options?: { zip?: ZipOptions })` accepts fflate ZIP options
under `options.zip`, including `mtime` and compression level. The default uses
the current time rather than v5's fixed
ZIP timestamps. Callers needing reproducible archives should provide `zip.mtime`;
use a local-calendar date when identical ZIP date fields across time zones matter.

The return value remains an asynchronous `Uint8Array` archive with the same
document/assets layout. Compression implementation and timestamp defaults changed,
so matching archive contents does not imply byte identity with v5.
`MVSData.toMVSX` currently does not forward ZIP options; the direct `createMVSX`
helper exposes this control.

## Follow-up: automatic MVS asset loading

`MVSData.toMVSX` now uses platform `fetch` instead of core's `ajaxGet` task.
Explicit `options.assets` still bypasses automatic fetching.

The options argument is optional and now accepts `fetch?: typeof globalThis.fetch`.
The callback receives the resolved URI and returns a platform-compatible `Response`.
Omitting it uses platform `fetch`. Cached assets are consulted before fetching,
including empty content, so a shared cache avoids repeated requests and permits
exports when the source is unavailable. Failed responses are not cached.

Under Node, platform `fetch` does not load `file://` assets. Callers can provide a
file-aware `options.fetch` adapter (for example, returning a `Response` containing
bytes read with Node's `fs` APIs), or supply `options.assets` directly. The old
`setFSModule(fs)` core IO hook no longer controls standalone builder fetching.
This explicit adapter contract keeps the builder independent of plugin/rendering
code. Tests cover default and custom fetching, resolved file URIs, cache reuse,
explicit assets, skipped external URIs, and HTTP errors.
