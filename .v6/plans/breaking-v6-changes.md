# v6 compatibility changes

Record public API and behavior changes introduced by the refactor here. The current prototype audit is complete for the
scope below. Update this ledger as implementation continues; release/downstream validation remains in
[checklist.md](checklist.md).

## Audit checkpoint

Compared v5 master `a1bf00ba55d737461e7e149ccc1996d979523b34` with prototype source at
`094c412a2f40692440434b5cd6eea8f1d0b77821` on 2026-10-04. Using the migration map, compared 1,547 production
TypeScript/TSX modules across packages, extensions, apps, CLI tools, and servers. No mapped targets or upstream
production source mappings were missing. Tests/performance fixtures were outside this source scan.

A TypeScript syntax/token scan removed top-level imports, normalized export module paths, and ignored
comments/formatting. It identified 63 remaining source differences, reviewed for exported names, namespace/class
members, signatures, and behavior. Reviewed the 17 new/extracted source modules separately. Internal shape-factory call
changes, ESM inline import suffixes, alias-only type changes, and CLI path rewrites account for much of the remaining
diff.

Checked compiled public modules and actual classic Viewer globals in local Chrome. Node probes confirmed the MVS
cache/file-URL changes and the standalone/runtime MolQL validation distinction. This checkpoint records the API/behavior
diff; it does not certify arbitrary downstream applications, native GL, or browser platforms. Changes already present in
the merged v5 master, including rendering/marking and BinaryCIF/PDB fixes, are not attributed to this refactor.

## Preserved classic globals

The classic Viewer still exposes `molstar.Viewer`, `ViewerAutoPreset`, `PluginExtensions`, `ExtensionMap`, `lib`,
`version`, `consoleStats`, `isDebugMode`, `isProductionMode`, `isTimingMode`, `setDebugMode`, `setProductionMode`, and
`setTimingMode`.

`molstar.lib` retains `structure`, `volume`, `shape`, `loci`, `math`, `plugin`, and `extensions`. The Viewer and MVS
Stories entrypoint implementations retain their public exports. MVS Stories still registers its custom elements and
exports its context/load/download functions and `MVSData`. Retaining these namespaces does not preserve every nested
API, particularly the shape changes below.

## Accepted: classic `molstar.lib.shape`

The classic bundle still exposes `molstar.lib`, but its shape objects now expose the model-side API following the
separation of model data and graphics helpers.

- `Shape.createRenderObject`, `Shape.createTransform`, `Shape.getTheme`, and `Shape.groupIterator` are no longer on the
  classic model-side `Shape` object.
- `ShapeGroup.getBoundingSphere` is no longer on its model-side object.
- Model-side `Shape.create` requires `groupCount` as argument 7, with optional `transforms` as argument 8. The old
  seventh-argument transforms call must change. JavaScript callers passing the old signature can silently store the
  transforms array as `groupCount`.

Graphics integrations belong to `@molstar/graphics/geo/shape/shape`. These classic API changes are accepted for v6; do
not restore them through a compatibility shim.

The graphics module exports functions directly, rather than a replacement `Shape`/`ShapeGroup` namespace:

| v5 call                                       | v6 graphics call                        |
| --------------------------------------------- | --------------------------------------- |
| `Shape.create(..., transforms?, groupCount?)` | `create(..., transforms?, groupCount?)` |
| `Shape.createRenderObject(shape, props)`      | `createRenderObject(shape, props)`      |
| `Shape.createTransform(...)`                  | `createTransform(...)`                  |
| `Shape.getTheme(shape)`                       | `getTheme(shape)`                       |
| `Shape.groupIterator(shape)`                  | `groupIterator(shape)`                  |
| `ShapeGroup.getBoundingSphere(loci, out?)`    | `getGroupBoundingSphere(loci, out?)`    |

For rendered shapes, import that graphics module as a namespace and use its `create` function; it preserves the old
optional-argument order and computes the initial group count from geometry when omitted. Model-side `Shape` defaults its
generic geometry type to `ShapeGeometry`, not the full graphics `Geometry` union. Supply a concrete generic type when
accessing geometry-specific fields.

`shape.groupCount` is now a stored value rather than a getter deriving the count from geometry on each access. Consumers
changing geometry group counts after creation must recreate the shape with the correct count instead of relying on the
old getter. The model contract still marks `groupCount` readonly.

## Temporary: standalone MolQL validation

The standalone builder and `mvs-validate` check MolQL expression structure but do not compile expressions. An unknown
symbol can therefore pass CLI validation. The MVS runtime still performs compiler validation and rejects it when
loading.

Full CLI validation will return after designing the mol-script import/integration. See the deferred validation item in
[checklist.md](checklist.md).

This also affects `MVSData.validationIssues`/`isValid` and schema decoding of MolQL selectors in the standalone builder.
Runtime sanity checks invoke the compiler validation for both scene and animation trees, and selector loading still
compiles expressions. A syntactically valid unknown symbol can pass standalone checks.

The standalone `@molstar/mvs-builder/expression` module exposes JSON expression construction/shape checks only. Use
`@molstar/model/script/...` for MolScript compiler/runtime functionality. New `@molstar/mvs/behavior-id` constants keep
the existing `molviewspec` name and `ms-plugin.molviewspec` transformer ID; the loader now checks that stable ID without
importing the behavior module. Serialized snapshot identity is unchanged.

## Headless MP4 integration

`HeadlessPluginContext` now owns image rendering only. Its former `getAnimation` and `saveAnimation` methods live on
`Mp4HeadlessPluginContext`, imported from `@molstar/mp4-export-extension/headless`. Register `Mp4Export` in the plugin
spec when using those methods. Native GL remains opt-in rather than auto-installed.

## Shape bounds and graphics helpers

Model-side shapes use a structural `ShapeGeometry` contract. Geometry-aware shape creation and group bounds belong to
`@molstar/graphics/geo/shape/shape`. `Loci.getBoundingSphere` for shape groups delegates to the optional
`shape.getGroupBoundingSphere` callback; the graphics factory supplies it. Shapes created directly through the model
factory need that callback for group bounds.

Streamline picking and location-iterator helpers moved from model props to `@molstar/graphics/geo/streamlines`:
`getStreamlinesVisualLoci`, `getStreamlinesLoci`, `eachStreamlines`, and `createStreamlinesLocationIterator`.

`PositionLocation` and `isPositionLocation` are owned by the model package and re-exported by the graphics location
iterator. `Viewport` is owned by core and re-exported by the graphics camera utility. Those re-exports preserve existing
behavior at the mapped module paths.

## Other helper moves requiring new imports

These are symbol-level moves within formerly larger modules; the whole-file migration map alone is insufficient for
these imports. Symbol moves from the plugin composition split are recorded in `migration-symbols.json` and described in
_Plugin composition step 1_ below.

| Public API                                                                                                                          | v5 module                            | v6 import                                              |
| ----------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------ | ------------------------------------------------------ |
| `CollapsableControls`, `CollapsableProps`, `CollapsableState`                                                                       | `mol-plugin-ui/base`                 | `@molstar/plugin-ui/controls/collapsable`              |
| `GaussianDensityTextureData`, `computeGaussianDensityTexture`, `computeGaussianDensityTexture2d`, `computeGaussianDensityTexture3d` | `mol-math/geometry/gaussian-density` | `@molstar/graphics/geo/gaussian-density`               |
| `fillIdentityTransform`                                                                                                             | `mol-geo/geometry/transform-data`    | `@molstar/core/math/linear-algebra/3d/transform-array` |

Collapsible controls are no longer exported from `@molstar/plugin-ui/base`. Gaussian-density CPU APIs remain in
`@molstar/core/math/geometry/gaussian-density`; their GPU functions/data types moved to graphics with the existing
signatures. The GPU helper behavior is carried over from v5; this audit does not redesign it. `fillIdentityTransform` is
no longer exported by the graphics transform-data module.

`Tokens` now lives in `@molstar/core/util/token-ranges`, with the same `data`, `count`, and `indices` fields. The IO
tokenizer retains its type re-export, so the mapped tokenizer import remains valid. `openRead` now lives in
`@molstar/common-server/open-read`; the volume server's mapped common/file module retains its re-export and behavior.

## Temporary: external color themes unavailable by default

`ColorTheme.createRegistry()` in graphics no longer includes `external-structure` or `external-volume`. They are also
intentionally absent from the default plugin registries for this prototype. The special registration helper in
`PluginContext` has been removed; it uses the ordinary graphics registry factory.

Implementations remain at `@molstar/plugin/themes/external-structure` and `@molstar/plugin/themes/external-volume`
because they depend on plugin state objects and selection queries. Bring them back through the planned registry
composition work, as tracked in the checklist. Existing presets/snapshots that request these themes cannot rely on
default registration until that work is done.

## Plugin composition step 1

Step 1 of [plugin-composition.md](plugin-composition.md) split facades and large modules so that applications compose
the plugin from catalogs instead of importing everything. Whole-file ownership stays in `migration-map.json`; where a v5
module was split, its entry maps to an array of the files that received its contents, or to `null` when nothing succeeds
it. [`migration-symbols.json`](migration-symbols.json) maps each exported symbol of an affected v5 module to its current
defining module (`null` when removed). `scripts/workspace/check.mjs` verifies that every target and symbol in both files
still exists.

### `StateTransforms` removed

`StateTransforms` and `@molstar/plugin/state/transforms` no longer exist. Import each transformer from its defining
module; the `StateTransforms.<Namespace>.<Name>` keys in `migration-symbols.json` list all 120 members.

| v5                                                         | v6 import                                                   |
| ---------------------------------------------------------- | ----------------------------------------------------------- |
| `StateTransforms.Data.Download`                            | `@molstar/plugin/state/transforms/data/fetch`               |
| `StateTransforms.Data.ParseCif`                            | `@molstar/plugin/state/formats/cif`                         |
| `StateTransforms.Model.TrajectoryFromMmCif`                | `@molstar/plugin/state/formats/trajectory/mmcif`            |
| `StateTransforms.Model.StructureFromModel`                 | `@molstar/plugin/state/transforms/structure/hierarchy`      |
| `StateTransforms.Model.StructureComponent`                 | `@molstar/plugin/state/transforms/structure/selection`      |
| `StateTransforms.Representation.StructureRepresentation3D` | `@molstar/plugin/state/transforms/structure/representation` |
| `StateTransforms.Volume.VolumeFromCcp4`                    | `@molstar/plugin/state/formats/volume/ccp4`                 |

The classic Viewer global `molstar.lib.plugin.StateTransforms` keeps its shape (the same seven namespaces and 120
members). It is now an app-level object in `@molstar/viewer/state-transforms` for classic-script pages that cannot
import leaf modules. The `StateTransforms` type alias (`typeof` the facade) is gone with the facade. Neither
`@molstar/viewer/lib` nor `@molstar/viewer/state-transforms` exports a `StateTransforms` type by name; where the type is
needed, use `typeof StateTransforms` imported from `@molstar/viewer/state-transforms`.

### Transform and format modules split

The seven transform modules and the six format modules are split by functionality:
`state/transforms/{data,misc, particles,shape,structure,volume}/...` (structure effects live in `structure/effects/`)
and `state/formats/{trajectory,coordinates,topology,volume,shape,particles}/<format>.ts`. Format-specific transformers
(`ParseCcp4`, `VolumeFromCcp4`, `TrajectoryFromPDB`, ...) live next to their provider in the format module. Each format
family also has `provider.ts`, `category.ts`, and `catalog.ts` modules for shared helpers, the category, and the
built-in provider list.

`@molstar/plugin/state/transforms/catalog` imports every built-in transform and format module for its registration side
effects and exports nothing. The default spec imports it, so a plugin built from `DefaultPluginSpec` registers all
built-in transformer ids. A custom spec must import the modules it uses, or the catalog, before it restores snapshots
(see _Snapshot loading_).

### Default specs

`DefaultPluginSpec` moved to `@molstar/plugin/default-spec` and `DefaultPluginUISpec` to
`@molstar/plugin-ui/default-spec`. `@molstar/plugin/spec` and `@molstar/plugin-ui/spec` keep only types and helpers
(`PluginSpec`, `PluginUISpec`). The classic Viewer globals `molstar.lib.plugin.DefaultPluginSpec` and
`DefaultPluginUISpec` are unchanged.

Both default-spec modules load `DefaultRegistry` and every built-in catalog. Apps that compose their own registry and
only want the default UI parts import `DefaultPluginUIComponents()` and `DefaultPluginUICustomParamEditors()` from
`@molstar/plugin-ui/default-ui` instead; `DefaultPluginUISpec()` is composed from them.

### Script languages

PyMOL, VMD, and Jmol are no longer enabled implicitly. Import `@molstar/model/script/transpilers/pymol`, `.../vmd`,
`.../jmol`, or `@molstar/model/script/transpilers/all` to enable a language; the default plugin spec and the Viewer
import `all`. Using a language that was not imported throws `Script language '<x>' is not available in this build`.
`mol-script` is always available. `Script.getAvailableLanguages()` returns `mol-script` plus the registered languages,
and the script-language select in the plugin UI lists only those. The `_transpiler` object exported by `transpilers/all`
is removed; after enabling the language, call `parse(language, text)` from `@molstar/model/script/transpile` (or
`Script.toExpression`), or import `transpiler` from `.../transpilers/<lang>/parser` directly.

### Data format providers

`DataFormatProvider` requires a `name` and gains an `Id extends string` type parameter so the name keeps its literal
type; the `DataFormatProvider(...)` helper takes a `const` type parameter for the same reason. Custom providers must add
`name`. `BuiltInTrajectoryFormats`, `BuiltInVolumeFormats`, `BuiltInShapeFormats`, `BuiltInTopologyFormats`,
`BuiltInCoordinatesFormats`, and `BuiltInParticlesFormats` are now `as const` provider arrays instead of
`[name, provider]` tuples; read `provider.name` instead of the first tuple element. `PluginSpec.customFormats` and
`DataFormatRegistry.add(name, provider)` keep their `[name, provider]` shape in this step.

The format name types (`BuiltInTrajectoryFormat`, `BuiltInVolumeFormat`, ...) are derived from the catalogs and live in
`@molstar/plugin/state/formats/<family>/catalog`. `BuiltInVolumeFormat` and `BuiltInShapeFormat` are new; the misspelled
`BuildInVolumeFormat` and `BuildInShapeFormat` remain as deprecated aliases. `MmcifProvider.parse` and
`CifCoreProvider.parse` take an optional `params` argument. The volume formats' parameter type, formerly an unexported
`Params`, is exported as `VolumeFormatParams` from `@molstar/plugin/state/formats/volume/provider`.

### Graphics catalogs

The `BuiltIn` values were removed from the `ColorTheme` and `SizeTheme` namespaces and from
`StructureRepresentationRegistry`, `VolumeRepresentationRegistry`, and `ParticleRepresentationRegistry`. The `BuiltIn`
types, `BuiltInParams`, and `createRegistry()` keep their names; `ColorTheme.createRegistry()` and
`SizeTheme.createRegistry()` return empty registries
([Registries start empty](#plugin-composition-step-3-registries-start-empty)). Read the values from the catalogs:

| v5 value                                  | v6 catalog export                                                                 |
| ----------------------------------------- | --------------------------------------------------------------------------------- |
| `ColorTheme.BuiltIn`                      | `BuiltInColorThemes` from `@molstar/graphics/theme/color/catalog`                 |
| `SizeTheme.BuiltIn`                       | `BuiltInSizeThemes` from `@molstar/graphics/theme/size/catalog`                   |
| `StructureRepresentationRegistry.BuiltIn` | `BuiltInStructureRepresentations` from `@molstar/graphics/repr/structure/catalog` |
| `VolumeRepresentationRegistry.BuiltIn`    | `BuiltInVolumeRepresentations` from `@molstar/graphics/repr/volume/catalog`       |
| `ParticleRepresentationRegistry.BuiltIn`  | `BuiltInParticleRepresentations` from `@molstar/graphics/repr/particles/catalog`  |

Catalogs are built with the new `namedCatalog` helper (`@molstar/graphics/util/named-catalog`), which checks at compile
time that each key equals its provider's `name`; the runtime key/name check in the registry constructors is removed.

### Selection queries

`@molstar/plugin/state/helpers/structure-selection-query` is split into `@molstar/plugin/state/queries/structure/*`:
`catalog` (`StructureSelectionQueries`), `registry` (`StructureSelectionQueryRegistry`), `query`
(`StructureSelectionQuery`, `StructureSelectionCategory`), `dynamic` (`ResidueQuery`, `ElementSymbolQuery`,
`EntityDescriptionQuery`, and the `get*Queries` helpers), and the category modules `basic`, `bond`, `common`,
`manipulate`, `residue`, `structure`, and `type`. `applyBuiltInSelection` was unused and is removed; apply
`StructureSelectionFromExpression` with `StructureSelectionQueries[name].expression` instead.

### Structure presets

Representation presets moved from `state/builder/structure/representation-preset` to
`state/builder/structure/representation-presets/`, and hierarchy presets from `.../hierarchy-preset` to
`.../hierarchy-presets/`: one module per preset, plus `types` (`StructureRepresentationPresetProvider`,
`TrajectoryHierarchyPresetProvider`, `presetStaticComponent`, `presetSelectionComponent`) and `catalog`
(`PresetStructureRepresentations`, `PresetTrajectoryHierarchy`). Presets gain an optional `alias` (for example `auto`,
`default`, `all-models`) next to `id`; `PresetProvider` and both provider helpers take `Id` and `Alias` type parameters.
The catalogs export the name types `BuiltInStructureRepresentationPresetId`,
`BuiltInStructureRepresentationPresetAlias`, `BuiltInTrajectoryHierarchyPresetId`, and
`BuiltInTrajectoryHierarchyPresetAlias`. The `auto` preset is its own module, `representation-presets/auto`.

### Behaviors and markdown extensions

`BuiltInPluginBehaviors` moved to `@molstar/plugin/behavior/built-in`; `PluginBehaviors` stays in
`@molstar/plugin/behavior`. `@molstar/plugin/behavior` no longer re-exports `PluginBehavior`; import it from
`@molstar/plugin/behavior/behavior`. The built-in markdown extensions moved out of
`@molstar/plugin/state/manager/markdown-extensions` into `@molstar/plugin/state/markdown/*`: `BuiltInMarkdownExtension`
is exported from `.../markdown/catalog`, and the individual extensions are in the sibling modules (`audio`, `camera`,
`highlight`, `query`, `snapshots`). `MarkdownExtension`, `MarkdownExtensionManager`, and
`defaultParseMarkdownCommandArgs` stay in the manager module.

### Snapshot loading

Restoring a snapshot (`PluginState.setSnapshot`, and every entry of `PluginStateSnapshotManager.setStateSnapshot`) first
checks that each transformer id in the behavior tree, data tree, and transition frames is registered. A snapshot that
names unregistered transformers throws one error listing all missing ids
(`Snapshot uses transformers that are not available in this plugin: <ids>. Import the modules that define them.`) before
the plugin or the snapshot manager changes. `StateTransformer.has(id)` is a non-throwing registration check;
`StateTransformer.get` still throws for unknown ids.

## Plugin composition step 2

Step 2 of [plugin-composition.md](plugin-composition.md) replaces three spec fields with registry entries. Step 3
([Registries start empty](#plugin-composition-step-3-registries-start-empty)) then removes the constructor preloads
(formats, representations, themes, presets, selection queries, markdown extensions).

### `PluginSpec.actions`, `animations`, and `customFormats` removed

`PluginSpec` no longer has `actions`, `animations`, or `customFormats`. The `PluginContext` constructor throws
`PluginSpec.<key> was removed in 6.0; use registry entries (see the migration guide)` when any of the three keys has a
value other than `undefined`; a key present with `undefined` (for example `customFormats: o?.customFormats`) is ignored.
This is a check, not an alias. The replacements are entries of `PluginSpec.registry`:

| 5.x                                 | 6.0                                                                                                                                     |
| ----------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------- |
| `actions: [PluginSpec.Action(a)]`   | `registry: [..., { actions: [a] }]`                                                                                                     |
| `animations: [...]`                 | `registry: [..., { animations: [...] }]`                                                                                                |
| `customFormats: [[name, provider]]` | `registry: [..., { formats: [DataFormatProvider.withName(provider, name)] }]`                                                           |
| `actions: defaultSpec.actions`      | `registry: DefaultRegistry` (or `defaultSpec.registry`), or the `DefaultActions` entry                                                  |
| `animations: [AnimateModelIndex]`   | replace the default entry: `registry: [...DefaultRegistry.filter((e) => e !== DefaultAnimations), { animations: [AnimateModelIndex] }]` |

`DefaultPluginSpec()` and `DefaultPluginUISpec()` return `{ registry: DefaultRegistry, behaviors, ... }` and no longer
carry `actions` or `animations`, so code that read `defaultSpec.actions`/`defaultSpec.animations` (or
`plugin.spec.actions`) must read the registry entries (`DefaultActions.actions`, `DefaultAnimations.animations`) or the
plugin (`plugin.state.data.actions`, `plugin.managers.animation.animations`). A literal spec that did not spread the
default spec gets no actions and no animations until it lists them (for example `registry: DefaultRegistry`); a spec
that listed `actions: defaultSpec.actions` and `animations: defaultSpec.animations` becomes
`registry: defaultSpec.registry`.

What changes now: actions and animations come only from `registry` (and behaviors); `PluginContext` no longer registers
`AnimateStateSnapshotTransition` when no animations are listed, so `plugin.managers.animation.current` is `undefined`
(it is typed `Current | undefined`) until an animation is registered. `AnimateStateSnapshotTransition` stays base
functionality: snapshot transitions and the UI snapshot controls play it directly and `play()` registers it on demand.
`DefaultAnimations` keeps it in its current position, so the default animation remains `AnimateModelIndex`.
`PluginAnimationManager.setSnapshot` skips `snapshot.current` and logs a warning when the snapshot names an unregistered
animation, instead of writing `paramValues` into the previously selected animation. `PluginContext.initCustomFormats`,
`initDataActions`, and `initAnimations` are removed (they were private).

What step 3 changes: see [Registries start empty](#plugin-composition-step-3-registries-start-empty); a literal spec
gets no formats, representations, themes, presets, selection queries, or markdown extensions either, and must list
`DefaultRegistry` (or its own entries).

### Viewer and mesoscale-explorer `customFormats`

`Viewer.create` options are unchanged, including `customFormats: [name, provider][]` with `DataFormatProvider.Unnamed`
providers. The Viewer builds
`registry: [...DefaultRegistry, { formats: customFormats.map(([n, p]) => withName(p, n)) }]` with the custom-formats
entry last. A custom name that equals a built-in format name (for example `['pdb', MyPdb]`) still overrides it:
`DefaultFormats` is replaced by a copy without that provider, so `dataFormats.get(name)` returns the custom provider and
no name is registered twice. The overriding provider is listed after the built-ins, so for that format `auto()`
tie-breaking by registration order sees it last, as in 5.x. mesoscale-explorer keeps its own `customFormats` option
shape and applies the same rule.

## Plugin composition step 3: registries start empty

`new PluginContext(spec)` (and `PluginUIContext`, `HeadlessPluginContext`) starts with empty registries. After `init()`
a plugin whose spec lists nothing has no representations, color or size themes (in the structure, volume, and particles
scopes), data formats, trajectory-hierarchy or representation presets, selection queries, markdown extensions,
drag-and-drop handlers, actions, or animations; it has only what `spec.registry` and the behaviors register. A literal
spec that relied on the implicit catalog adds `DefaultRegistry` (or its own entries):

```ts
// 5.x: a literal spec received every built-in provider
const spec = { behaviors: [...] };
// 6.0
import { DefaultRegistry } from '@molstar/plugin/default-registry';
const spec = { registry: DefaultRegistry, behaviors: [...] };
```

Script-tag code using the classic Viewer global does the same with `molstar.lib.plugin.DefaultRegistry`, which the
Viewer now exposes next to `DefaultPluginSpec` and `DefaultPluginUISpec`.

Specs built from `DefaultPluginSpec()` or `DefaultPluginUISpec()` list `DefaultRegistry` already and are unchanged.

- `ColorTheme.createRegistry()` and `SizeTheme.createRegistry()` return empty registries. Code that creates a
  `RepresentationContext` by hand (for example `colorThemeRegistry: ColorTheme.createRegistry()`) must add the themes it
  uses (`BuiltInColorThemes` from `@molstar/graphics/theme/color/catalog`, or single providers).
- `StructureRepresentationRegistry`, `VolumeRepresentationRegistry`, and `ParticleRepresentationRegistry` constructors
  add nothing; `DataFormatRegistry`, `TrajectoryHierarchyBuilder`, `StructureRepresentationBuilder`,
  `StructureSelectionQueryRegistry`, and `MarkdownExtensionManager` constructors add nothing either. The modules that
  define these classes no longer import any catalog.
- `PluginContext.dropOverriddenPreloadedFormats` (the transitional handling that let a `customFormats` entry replace a
  preloaded built-in format) is removed. The Viewer and mesoscale-explorer replace the matching default entry as
  described above, so a custom format with a built-in name still overrides it; a different provider under an existing
  format name in any other spec is an error, as for every registry.
- mesoscale-explorer no longer clears registries after `init()`; it lists spacefill with its themes, the `uniform` and
  `illustrative` color themes, the mmCIF format, its three animations, and its custom formats.
- The import-graph check requires `@molstar/plugin/context`, `@molstar/plugin/spec`, `@molstar/plugin-ui`, and the base
  UI entry points to reach no catalog module by value.

## Plugin composition step 3: UI

### `createPluginUI` requires a spec

`createPluginUI({ target, render, spec })` no longer falls back to `DefaultPluginUISpec()`; `spec` is a required option
and `@molstar/plugin-ui` (`index.ts`) no longer imports `@molstar/plugin-ui/default-spec`. Pass
`spec: DefaultPluginUISpec()` (import it from `@molstar/plugin-ui/default-spec`) for the previous behavior.

### Base UI renders minimal structure tools

The base `Plugin` layout (`ControlsWrapper`) no longer falls back to `DefaultStructureTools`. Without
`spec.components.structureTools` it renders `MinimalStructureTools`: structure source, measurements, components, and the
behavior-registered custom structure controls. The full tools component (also adding superposition, quick styles,
procedural animation, volume streaming, volume source, and particle source) moved from `@molstar/plugin-ui/controls` to
`@molstar/plugin-ui/default-structure-tools` and is set as `DefaultPluginUISpec().components.structureTools`.

A spec that spreads `DefaultPluginUISpec()` but replaces `components` wholesale must spread the default components to
keep the full tools, or it silently gets the minimal list:

```ts
const defaultSpec = DefaultPluginUISpec();
const spec = { ...defaultSpec, components: { ...defaultSpec.components, remoteState: 'none' } };
```

### Quick styles and volume controls

Quick Styles resolves its presets by id through the representation preset registry (`DefaultRepresentationPreset`
config, then `preset-structure-representation-auto` for Default; `-polymer-and-ligand`, `-illustrative`, and
`-molecular-surface` for Cartoon, Spacefill, and Surface) and hides the buttons whose preset is not registered.
Lazy-volume loading in the volume source controls uses the isosurface representation and uniform color theme only when
they are registered, and otherwise falls back to the registry default.

## Plugin composition step 3

Step 3 of [plugin-composition.md](plugin-composition.md) removes implicit defaults from the plugin context.

### View models and hooks require a spec

`PluginViewModel` and `PluginUIViewModel` (`@molstar/plugin-extension`) take `{ spec }` as a required constructor option
and no longer fall back to `DefaultPluginSpec()` / `DefaultPluginUISpec()`. `useCreatePluginViewModel` and
`useCreatePluginUIViewModel` require `options.spec` (a spec object) and no longer accept the
`spec: (defaultSpec) => spec` callback form. Callers import the default spec themselves. The classic Viewer global
(`molstar.PluginExtensions.plugin.models`) exposes app-level subclasses that keep the 5.x default when no spec is given
(`apps/viewer/src/view-models.ts`).

```ts
// 5.x
useCreatePluginUIViewModel({ spec: (spec) => ({ ...spec, behaviors: [...spec.behaviors, MyBehavior] }) });
// 6.0
const spec = DefaultPluginUISpec();
useCreatePluginUIViewModel({ spec: { ...spec, behaviors: [...spec.behaviors, MyBehavior] } });
```

### Drag-and-drop open-files fallback is a registry entry

`DragAndDropManager.handle()` no longer opens unrecognized files itself. The open-anything handler is the
`{ name: 'open-files', handle, fallback: true }` entry of `DefaultDragAndDrop` (part of `DefaultRegistry`), so it still
runs after the handlers and the built-in session handling (`.molx`/`.molj`). A plugin without that entry ignores drops
it does not recognize; add `DefaultDragAndDrop` (or your own fallback handler) to open arbitrary files.

### Snapshot name warnings

`PluginState.setSnapshot` calls the new `PluginState.reportUnregisteredNames(snapshot)` after the behavior tree is
applied and before the data tree. It logs a warning for each representation type, color theme, and size theme name of
the structure, volume, and particles representation transformers (data tree and transition frames) that is not
registered in the matching scope; the registry default replaces it when the data tree is normalized. It never throws.

## Plugin composition step 3: presets

### Mixed name and provider props

`createStructureRepresentationParams`, `createVolumeRepresentationParams`, and the structure representation builder
(`addRepresentation`, `buildRepresentation`) accept a name or a provider in each of `type`, `color`, and `size` and
resolve each field on its own, so `{ type: CartoonRepresentationProvider, color: 'element-symbol' }` works. The "any
string field selects the by-name path" dispatch is gone. The params types follow the field: a provider gives its params,
a built-in name gives `BuiltInParams`, and any other string gives `{}`. The prop types are now
`StructureRepresentationProps<R, C, S>` and `VolumeRepresentationProps<R, C, S>` with `R`, `C`, and `S` constrained to a
provider or a string (the new `*RepresentationRef`, `*ColorThemeRef`, and `*SizeThemeRef` types);
`*RepresentationBuiltInProps` is kept as the name-only alias. The builder methods are generic in `R`, `C`, and `S`
instead of in the whole props object, so the params of each field are checked against that field; a call that previously
compiled only because the props were checked against the union of every built-in's params can now fail to compile (for
example `colorParams: { palette }` on `illustrative`). Names that are not registered still warn and use the registry
default, and an empty registry still throws.

### Presets import what they run and export entries

Each representation preset passes provider objects (and the imported theme providers for its fixed color themes) to the
builder instead of names, and its module exports an entry next to the preset: `EmptyPresetEntry`, `AutoPresetEntry`,
`AtomicDetailPresetEntry`, `PolymerCartoonPresetEntry`, `PolymerAndLigandPresetEntry`, `ProteinAndNucleicPresetEntry`,
`CoarseSurfacePresetEntry`, `IllustrativePresetEntry`, `MolecularSurfacePresetEntry`, `AutoLodPresetEntry`, and
`MesoscalePresetEntry`. An entry lists its own preset first, followed by the representation entries (and theme
providers) it builds; `AutoPresetEntry` also includes the entries of the presets `auto` chooses between. The hierarchy
presets export `DefaultHierarchyPresetEntry`, `AllModelsHierarchyPresetEntry`, `UnitcellHierarchyPresetEntry`,
`SupercellHierarchyPresetEntry`, and `CrystalContactsHierarchyPresetEntry` with the color themes they apply (and, for
all-models, the default preset it delegates to). Hierarchy preset entries do not include a representation preset; that
is the configured one. `@molstar/plugin/registry/merge` exports `mergeRegistryEntries(...entries)` for building such
entries. `DefaultPresets` is built from these entries and lists the same 5 hierarchy and 11 representation presets in
the same order. The preset providers keep their names.

### Hierarchy preset `representationPreset`

The `representationPreset` param of every hierarchy preset is typed with the built-in id and alias unions (`import type`
from the catalog) and no longer defaults to `'auto'`: when it is absent or empty the preset applies
`PluginConfig.Structure.DefaultRepresentationPreset`, and the hierarchy presets no longer import the `auto` preset. The
configured preset is resolved through the registry, so a hierarchy preset applied to a plugin that has not registered it
fails with `Preset '<id>' is not registered in this plugin`. Calls that relied on the old `'auto'` param default and a
config that names another preset now get the configured preset.

### `presetSelectionComponent`

`presetSelectionComponent(plugin, structure, query, tag, params?)` takes a `StructureSelectionQuery` object and a tag
instead of a `StructureSelectionQueries` key; the component key stays `selection-<tag>`. Import the query from its
module (for example `protein` and `nucleic` from `@molstar/plugin/state/queries/structure/type`). The module no longer
imports the query catalog. Replace `presetSelectionComponent(plugin, s, 'protein')` with
`presetSelectionComponent(plugin, s, protein, 'protein')`.

### Delegating presets

The Viewer's `ViewerAutoPreset` looks up the model-archive quality-assessment and SB-NCBR partial-charges presets by id
through `plugin.builders.structure.representation.resolveProvider` and skips them when they are not registered, instead
of importing them; it still falls back to the imported `AutoPreset`.

## Plugin composition step 3: formats and post-load presets

### Format entries

Every built-in data format module exports a registry entry next to its provider: the provider name without `Provider`
(`SdfProvider` and `Sdf`, `MmcifProvider` and `Mmcif`, `Ccp4Provider` and `Ccp4`, `RelionStarParticlesProvider` and
`RelionStarParticles`). An entry holds `formats: [provider]` and the actions `DefaultActions` lists for the format
(CCP4: `ParseCcp4`, `VolumeFromCcp4`; mmCIF: `ParseCif`, `TrajectoryFromMmCif`; PDB, PDBQT, and PQR:
`TrajectoryFromPDB`; SDF: none), so a slim plugin registers `registry: [Sdf, Ccp4, ...]`. Volume and particle format
entries also include the representations and themes their `visuals` step builds: `Isosurface` for CCP4, DSN6, DX, Cube,
density-server CIF, structure-factor CIF, and MTZ, `Segment` for segmentation CIF, and the particle spacefill
representation or the spacefill, fibers, and target representations with their color themes for the particle formats.
Shape, topology, coordinates, and trajectory formats carry no representations. The `visuals` steps pass provider objects
(`IsosurfaceRepresentationProvider`, `UniformColorThemeProvider`) instead of names, and
`VolumeRepresentation3DHelpers.getDefaultParams` and `getDefaultParamsStatic` accept a provider object or a name for the
representation and the color and size themes. `DefaultFormats` and `DefaultActions` are unchanged, and the built-in
catalogs (`BuiltIn*Formats`) still list the providers.

### Post-load presets come from config

Trajectory format `visuals` and the `DownloadStructure` action apply `PluginConfig.Structure.DefaultHierarchyPreset`
(default `preset-trajectory-default`) instead of `'default'`. `LoadTrajectory`, `AddTrajectory`,
`StructureHierarchyManager.updateStructure`, and the Cube format apply
`PluginConfig.Structure.DefaultRepresentationPreset` (default `preset-structure-representation-auto`) instead of
`'auto'`, and `DownloadStructure` takes both the default of its representation select and the empty-preset check from
config. A config value that is not registered fails the load with `Preset '<id>' is not registered in this plugin` (the
actions log it and revert, as for any other error in their transaction); there is no fallback to the default preset. A
plugin that does not register the default presets sets the two config items.

### `DownloadStructure` depends on the registered formats

`DownloadStructure` offers only the sources whose format is registered (`plugin.dataFormats.has`): PDB, PDB-IHM,
AlphaFold DB, and Model Archive need `mmcif`, SWISS-MODEL needs `pdb`, PubChem needs `mol`, and URL needs at least one
trajectory format. The URL format list is the registered formats of the trajectory category in registration order
(formats that extensions register in that category, such as `g3d` in the Viewer, now appear in it;
`BuiltInTrajectoryFormat` no longer types the `format` param, which is a `string`), and its default is `mmcif` when
registered. The "Load all entries into a single trajectory" option is hidden unless `mmcif` is registered. With no
usable source the action's source select has a single empty option.

### `builders.structure.parseTrajectory(blob)`

The blob overload no longer imports the mmCIF transformers. It calls `parseBlob` of the registered `mmcif` format and
fails with `parseTrajectory(blob) requires the 'mmcif' data format to be registered in this plugin.` when mmCIF is not
registered. `parseBlob` is a new optional method of `TrajectoryFormatProvider` implemented by `MmcifProvider`, and
`parseMmcifBlob(plugin, blob, params)` in `@molstar/plugin/state/formats/trajectory/mmcif` is the same function.

## Declaration contracts

`ExternalModules['jpeg-js']` exposes the injected codec's `encode` contract instead of the entire codec module type.
`JpegBufferRet` is owned by headless and describes the returned `width`, `height`, and `Buffer` data.

Standalone MVS `validateTree` accepts the logger contract `{ log: { error(message: string): void } }` rather than
requiring `PluginContext`. Existing plugin arguments still satisfy that contract.

`Clip.Rotation` is a new exported alias for `{ axis: Vec3, angle: number }`. Explicit parameter-schema annotations
retain the existing rotation value shape. The background image asset declaration changes from TypeScript `export =` to a
default export to match ESM asset imports.

Standalone MVS also exposes package-local `object`, `json`, and `color-names` helpers. The builder's `ColorNames` values
are membership flags (`true`), not the numeric colors in core's color-name table; use core when a color value is needed.
The MVS JSON expression shape remains structurally compatible with MolScript expressions, but these builder helpers do
not provide compiler semantics.

## Distribution and import paths

Compiled packages use ESM and new package subpaths. Old monolithic CommonJS and source paths require migration; use
`migration-map.json` for ownership mappings and `migration-symbols.json` for symbols whose module was split or removed.
Classic Viewer/MVS Stories asset paths and globals remain available, subject to the shape API changes above. Package
exports exclude tests and build caches.

`molstar` now packages only browser distributions. Library/server/CLI consumers must install the owning `@molstar/*`
packages instead of importing `molstar/lib/...` or `molstar/lib/commonjs/...`. There is no CommonJS build or `require`
export condition. Published library modules and the browser ESM distribution target ES2022. All 43 public packages share
one release version.

Native `gl`/`canvas` remain injected by headless callers and are optional peers. `@molstar/plugin-headless/native` adds
lazy `loadNativeModule` resolution for the MVS rendering CLI; importing it or requesting `--help` no longer eagerly
loads native modules. A rendering command still needs the caller to install them explicitly.

## Version metadata

`PLUGIN_VERSION` is generated from `version.json` rather than falling back to `'(development)'` outside a bundle. The
plugin startup log displays only this version. `PLUGIN_VERSION_DATE` and the build timestamp override have been removed.
Additional build identification metadata is deferred in the v6 checklist. The old `__MOLSTAR_PLUGIN_VERSION__`
definition does not control the generated module's value.

## CLI generated imports and entrypoints

`cifschema --moldataImportPath` now defaults to `@molstar/core/data`, replacing `molstar/lib/mol-data` in generated
TypeScript. Pass the option explicitly to generate another import path. Its field-name CSV presets read canonical root
data in the workspace and staged package assets when installed.

CLI/server entrypoints live in their owning packages and use checked-in `bin/*.mjs` launchers. Their existing command
argument parsing is retained; repository maintenance generators such as `syminfo` target the new package/test paths.
Model-server TAR helpers now expose `decodeLongPath`, `encodePax`, and `decodePax` as named ESM/TypeScript exports
instead of CommonJS `exports` assignments.

## MVSX ZIP options

`createMVSX(data, assets, options?: { zip?: ZipOptions })` accepts fflate ZIP options under `options.zip`, including
`mtime` and compression level. The default uses the current time rather than v5's fixed ZIP timestamps. Callers needing
reproducible archives should provide `zip.mtime`; use a local-calendar date when identical ZIP date fields across time
zones matter.

The return value remains an asynchronous `Uint8Array` archive with the same document/assets layout. Compression
implementation and timestamp defaults changed, so matching archive contents does not imply byte identity with v5.
`MVSData.toMVSX` currently does not forward ZIP options; the direct `createMVSX` helper exposes this control.

## Follow-up: automatic MVS asset loading

`MVSData.toMVSX` now uses platform `fetch` instead of core's `ajaxGet` task. Explicit `options.assets` still bypasses
automatic fetching.

The options argument is optional and now accepts `fetch?: typeof globalThis.fetch`. The callback receives the resolved
URI and returns a platform-compatible `Response`. Omitting it uses platform `fetch`. Cached assets are consulted before
fetching, including empty content, so a shared cache avoids repeated requests and permits exports when the source is
unavailable. Failed responses are not cached.

Under Node, platform `fetch` does not load `file://` assets. Callers can provide a file-aware `options.fetch` adapter
(for example, returning a `Response` containing bytes read with Node's `fs` APIs), or supply `options.assets` directly.
The old `setFSModule(fs)` core IO hook no longer controls standalone builder fetching. This explicit adapter contract
keeps the builder independent of plugin/rendering code. Tests cover default and custom fetching, resolved file URIs,
cache reuse, explicit assets, skipped external URIs, and HTTP errors.
