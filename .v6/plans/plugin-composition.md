# Plugin composition implementation plan

Status: planned, not started (2026-10-06). Line numbers below refer to commit `477703ca5`.

This plan implements the [plugin-composition design](../designs/plugin-composition.md) ("the spec"; `spec §N` refers to
its sections). The spec holds the contract, rules, and decisions; this plan holds the ordered steps, call-site
inventories, per-consumer migration, tooling, and the acceptance checklist. Remaining v6 work is tracked in
[checklist.md](checklist.md); accepted API changes go to [breaking-v6-changes.md](breaking-v6-changes.md) as listed in
spec §13. Paths without a package prefix are relative to `packages/plugin/core/src/`.

## 1. Inventories

### 1.1 Current coupling

`new PluginContext(spec)` loads the full built-in catalog through these edges. Step 3 removes each of them.

| Edge                                                            | Mechanism                                                                                        | What it loads                                                                        |
| --------------------------------------------------------------- | ------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------ |
| `context.ts` → graphics representation registries               | Constructors add every `BuiltIn` entry; class and catalog share a module                         | 17 structure, 5 volume, 4 particle representations                                   |
| `context.ts` → `ColorTheme`/`SizeTheme.createRegistry()`        | `ThemeRegistry` requires a built-in map; created for three scopes                                | 37 color and 6 size themes, three times                                              |
| `context.ts` → `DataFormatRegistry`                             | `state/formats/registry.ts` value-imports all six format lists; format modules import the facade | 37 formats, the 116 transformers of the seven facade modules                         |
| `context.ts` → `StructureBuilder` → preset builders             | Constructors add all presets; `resolveProvider` checks the static map first                      | 5 hierarchy and 11 representation presets                                            |
| `hierarchy-preset.ts` → `PresetStructureRepresentations`        | Value import used as `PresetStructureRepresentations.auto.id` fallback                           | Every representation preset                                                          |
| `context.ts` → `StructureSelectionQueryRegistry`                | Constructor preload; registry, query type, and catalog share a module that imports the facade    | 65 queries (35 catalog queries, 30 residue queries)                                  |
| `context.ts`, `state.ts`, focus representation → `behavior.ts`  | `export *` of `behavior/behavior` plus the `PluginBehaviors` map of dynamic behaviors            | All dynamic behaviors, including eight custom-property ones                          |
| `context.ts` → `StructureComponentManager`                      | Direct value imports                                                                             | `InteractionsProvider`, focus representation behavior, the query catalog, transforms |
| `context.ts` → `MarkdownExtensionManager`                       | Constructor preload; `BuiltInMarkdownExtension` lives in the manager module                      | 14 extensions, `Script` and every transpiler, the snapshot transition animation      |
| `context.ts` → `DragAndDropManager`                             | Hard-coded fallback called from `handle()`                                                       | `OpenFiles` action, `PluginCommands` (session open)                                  |
| `state.ts`, markdown manager, UI controls → snapshot transition | Value import of `AnimateStateSnapshotTransition`, played directly                                | One animation (small, no catalogs)                                                   |
| `spec.ts` (`@molstar/plugin/spec`)                              | `DefaultPluginSpec` shares the module with the `PluginSpec` type                                 | The facade, `StateActions`, all behaviors, all animations                            |
| `plugin-ui/spec.ts`, `createPluginUI`                           | `DefaultPluginUISpec` shares the module with `PluginUISpec`; `createPluginUI` falls back to it   | The full default spec and volume-streaming UI                                        |
| `plugin-ui/plugin.tsx` → `DefaultStructureTools`                | Fallback when `components.structureTools` is not set                                             | Volume streaming, quick styles (preset map), volume/particle UIs                     |
| Facade cycle: `transforms/model` → selection queries → facade   | `StateTransforms` lazy getters work around it (issue #1791)                                      | —                                                                                    |

### 1.2 `registry.default` call sites

All 17 are on structure, volume, or particle representation registries: nine in transformer param definitions (six
`PD.Mapped` type defaults, three theme-default lookups) and eight in the representation-param helpers. Spec §4.3 gives
the empty-registry rules for both groups.

| Group                                    | Sites                                                                                                                                                                                   |
| ---------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Param definitions, `PD.Mapped` defaults  | `state/transforms/representation.ts:103, 129, 1274, 1288`; `state/transforms/particles.ts:533, 559`                                                                                     |
| Param definitions, theme-default lookups | `registry.get(registry.default.name)` at `state/transforms/representation.ts:90, 1270`; `state/transforms/particles.ts:520` (seed `defaultColorTheme.name` and `defaultSizeTheme.name`) |
| Helpers                                  | `state/helpers/structure-representation-params.ts:108, 137, 151, 177`; `state/helpers/volume-representation-params.ts:103, 132, 146, 172`                                               |

### 1.3 Preset short keys

Short keys become the `alias` of each built-in preset (spec §4.3); the call sites below keep working unchanged and are
listed for verification. Keys do not overlap between the two builders.

| Builder        | Alias (short key) → id                                                                                                                                                                                                                                                                   |
| -------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Hierarchy      | `default` → `preset-trajectory-default`, `all-models` → `preset-trajectory-all-models`, `unitcell` → `preset-trajectory-unitcell`, `supercell` → `preset-trajectory-supercell`, `crystalContacts` → `preset-trajectory-crystal-contacts`                                                 |
| Representation | `empty`, `auto`, `atomic-detail`, `polymer-cartoon`, `polymer-and-ligand`, `protein-and-nucleic`, `coarse-surface`, `illustrative`, `molecular-surface`, `auto-lod`, `mesoscale` → `preset-structure-representation-<key>` (for example `auto` → `preset-structure-representation-auto`) |

`applyPreset(x, '<short>')` calls in code (19), plus one computed default:

- `state/formats/trajectory.ts:28` (`'default'`), `state/formats/volume.ts:239` (`'auto'`, Cube path),
  `state/manager/structure/hierarchy.ts:229` (`'auto'`, `updateStructure`).
- `state/actions/structure.ts:308, 321` (`'default'`, `DownloadStructure`), `:501` (`'auto'`, `AddTrajectory`), `:654`
  (`'auto'`, `LoadTrajectory`).
- `extensions/plugin/src/loaders.ts:80` (`'all-models'`), `:95` (`'default'`), `:390` (`params.preset ?? 'default'`).
- `extensions/volumes-and-segmentations/src/entry-models.ts:31` (`'default'`).
- `examples/interactions/src/index.tsx:155, 161, 250, 254, 309, 313` (`'default'`).
- `examples/basic-wrapper/src/index.ts:82` (`'default'`).
- `smoke/headless/capture.mjs:54` and `smoke/browser/library/index.html:15` (`'default'`).

Documentation examples with `applyPreset(trajectory, 'default')` (unchanged; resolved through the alias):
`docs/docs/plugin/custom-library.md:54, 144`, `docs/docs/plugin/instance.md:151, 201-204, 311`,
`docs/docs/plugin/transforms/custom-conformation.md:12`, and `docs/docs/plugin/transforms/custom-trajectory.md:47, 73`.

Short keys in `representationPreset` values (6): `extensions/volumes-and-segmentations/src/entry-models.ts:32`
(`'polymer-cartoon'`), `examples/interactions/src/index.tsx:156` (`'polymer-cartoon'`), `:162, 255, 314`
(`'atomic-detail'`), `examples/basic-wrapper/src/index.ts:93` (`'auto'`).

Short-key types, defaults, and lookups:

- `state/builder/structure/hierarchy.ts`: the `keyof PresetTrajectoryHierarchy` union (`:22`), `defaultProvider`
  (`:30`), the static-map branch of `resolveProvider` (`:34`), the `applyPreset<K extends keyof ...>` overload (`:95`),
  and `getPresetSelect`'s `PD.Select('auto', …)` default (`:66`), whose options are ids.
- `state/builder/structure/representation.ts`: the same at `:28`, `:39`, `:43`, `:104`, and `getPresetSelect` (`:69`).
- `state/builder/structure/hierarchy-preset.ts:40`:
  `representationPreset: PD.Text<keyof PresetStructureRepresentations>('auto')`; the
  `PresetStructureRepresentations.auto.id` fallbacks at `:80-82`, `:143-145`, `:182-184`, `:275-277`.
- `extensions/plugin/src/loaders.ts:513`: `LoadTrajectoryParams.preset?: keyof PresetTrajectoryHierarchy`, and the
  `preset: 'all-models' // or 'default'` doc comments at `loaders.ts:330` and `apps/viewer/src/app.ts:204` (Viewer API,
  §3.1).
- `state/actions/structure.ts:41, 290-291`: `PresetStructureRepresentations.auto.id`/`.empty.id`.
- `packages/plugin/ui/src/structure/quick-styles.tsx:140-156`: `PresetStructureRepresentations.auto.id`,
  `.illustrative`, `['polymer-and-ligand']`, `['molecular-surface']`.

Value calls into the representation-preset catalog outside it (`PresetStructureRepresentations.auto.apply`), which
import the auto preset's defining module instead; `auto` pulls in the presets it composes (spec §5.1):
`extensions/dnatco/src/confal-pyramids/behavior.ts:45`, `extensions/dnatco/src/ntc-tube/behavior.ts:45`,
`extensions/sb-ncbr/src/tunnels/behavior.ts:96`, `extensions/sb-ncbr/src/partial-charges/preset.ts:26`,
`extensions/model-archive/src/quality-assessment/behavior.ts:229, 254`,
`extensions/assembly-symmetry/src/behavior.ts:222`, `extensions/anvil/src/behavior.ts:198`,
`extensions/rcsb/src/validation-report/behavior.ts:335, 405, 442`, and `apps/viewer/src/presets.ts:47` (an app, allowed
to keep the catalog, but switched for consistency).

### 1.4 Data-format `get` checks

`dataFormats.get` throws for an unknown name, so these checks are dead and become `has(name)` checks before `get`:

- `state/actions/structure.ts:318` (`DownloadStructure`) and `:591` (`LoadTrajectory`).
- `state/actions/file.ts:134` (`DownloadFile`).
- `state/builder/structure.ts:40` (`parseTrajectory`).
- `extensions/plugin/src/loaders.ts:244` (`loadVolumeFromUrl`).

These stay, because they follow a possible `auto()` result: the `'auto'`-or-`get` helper at the top of
`state/actions/file.ts` (`:20-22`) and the `auto` branch of `state/helpers/particle-targets.ts:180-183`.

### 1.5 Script transpiler reachers

Modules that reach `transpilers/all.ts` through `@molstar/model/script/script` today: the graphics
`theme/{overpaint,transparency,clipping,substance,emissive,wiggle}` modules (which every representation imports through
`repr/representation.ts`), `state/helpers/structure-component.ts`, `state/helpers/structure-query.ts`,
`state/transforms/model.ts`, `state/transforms/representation.ts`, `state/manager/markdown-extensions.ts` (the `query`
extension), and the UI parameter controls. After the script split (step 1) none of them reaches a transpiler.

Direct `parse` users of `@molstar/model/script/transpile`, which reaches `transpilers/all.ts` today:
`packages/model/src/script/script.ts`, `examples/mvs-stories/src/stories/molql.ts` (`parse('pymol', ...)` at module
level), and `packages/mvs/runtime/src/_test/molql.test.ts` (`parse('pymol', ...)`). After the split the last two import
`@molstar/model/script/transpilers/pymol` (or `all`) themselves.

### 1.6 Hard-coded policy and provider names

Sites that name a specific preset or provider and change in step 3 (spec §5.2):

- Hierarchy preset after loading: trajectory format `visuals` (`state/formats/trajectory.ts`) and `DownloadStructure`
  read `DefaultHierarchyPreset` instead of `'default'`.
- Representation preset after loading: `LoadTrajectory` and `AddTrajectory` (`state/actions/structure.ts`),
  `StructureHierarchyManager.updateStructure` (`state/manager/structure/hierarchy.ts`), the Cube path in
  `state/formats/volume.ts`, and the hierarchy presets read `DefaultRepresentationPreset` instead of `'auto'` or
  `.auto.id`. `state/actions/structure.ts` also replaces `.auto.id`/`.empty.id` with the config value and a string id
  constant.
- Focus representation behavior (`behavior/dynamic/selection/structure-focus-representation.ts:37-63`):
  `ball-and-stick`, `physical`, and the `SizeTheme.BuiltIn.uniform` catalog read.
- Volume streaming (`behavior/dynamic/volume-streaming/transformers.ts:347`):
  `VolumeRepresentationRegistry.BuiltIn.isosurface`.
- `updateFocusRepr` (representation-preset helpers): `ball-and-stick` and `element-symbol`, and a value import of the
  focus behavior.
- Volume and particle format `visuals`: representation and theme names, replaced by provider objects.
- `builders.structure.parseTrajectory(blob)`: mmCIF only. `DownloadStructure`: URL format list from
  `BuiltInTrajectoryFormats` and per-source hard-coded formats. Cube format: builds a structure through the model
  transforms.

### 1.7 Files touched

- `packages/plugin/core/src`: `spec.ts`, new `default-spec.ts` and `default-registry.ts`, `context.ts`, `state.ts`,
  `config.ts`, `behavior.ts` and `behavior/behavior.ts`, new built-in behaviors module, the focus representation
  behavior and its id module, volume streaming, `state/transforms/**` (split, catalog), `state/formats/**` (split,
  catalogs), `state/actions/**`, `state/builder/structure/**` (preset split, builders), `state/helpers/*-params.ts`,
  `state/helpers/structure-component.ts`, `state/helpers/particle-targets.ts`, `state/manager/**` (component, hierarchy,
  markdown, drag-and-drop, animation, loci labels, snapshots), `themes/external-*`, the selection-query modules.
- `packages/core/src/state`: `StateActionManager` counting and `toAction()` caching only (no normalization change).
- `packages/graphics/src`: representation and theme registries (`has`, counting, no preloads, typed `default`), catalog
  modules for representations and themes, namespace type imports.
- `packages/model/src/script`: `script.ts` split, transpiler table in `transpile.ts`,
  `transpilers/{pymol,vmd,jmol,all}.ts`; MolScript runtime-table counting.
- Direct `parse` users: `examples/mvs-stories/src/stories/molql.ts`, `packages/mvs/runtime/src/_test/molql.test.ts`.
- `docs/docs/plugin/**`: the preset short keys in §1.3.
- `packages/plugin/ui/src`: `spec.ts`, new `default-spec.tsx`, `index.ts`, `plugin.tsx`/`ControlsWrapper`, minimal
  structure tools, quick styles, volume controls, script-param control, snapshot controls.
- `packages/mvs/runtime/src`: `behavior.ts`, `MVSRuntimeRegistry` catalog module, MVSJ/MVSX providers.
- `extensions/*`: behaviors listed in §1.3, §3.8; `extensions/plugin` loaders and view models; g3d provider.
- `apps/*`, `examples/*`, `cli/*`, `smoke/*`: §3.
- `scripts/workspace/`: `registry-dump.mjs`, `import-graph.mjs`, `catalogs.json`; `.v6/baselines/*.json`.

### 1.8 Duplicate handling today

What each registry does today with a duplicate registration; step 2 replaces it with the counting rules of spec §4.2.
Removing an unknown provider from the representation, theme, or format registries today drops the last entry
(`splice(-1, 1)`), which behaviors that call `registry.remove(p)` in `unregister()`, such as MVS, can hit.

| Entry field                                    | Today on duplicate                 |
| ---------------------------------------------- | ---------------------------------- |
| `<scope>.themes.color/size`                    | Throws, even for the same object   |
| `<scope>.representations`                      | Throws, even for the same object   |
| `formats`                                      | Appends a second entry             |
| `structure.presets.hierarchy`/`representation` | Throws                             |
| `structure.selectionQueries`                   | Appends                            |
| `lociLabels`                                   | Appends                            |
| `markdownExtensions`                           | Replaces                           |
| `dragAndDrop`                                  | Replaces                           |
| `actions`                                      | Ignored; `remove` is unconditional |
| `animations`                                   | Logs an error; no unregister       |

## 2. Steps

Each step keeps the build, in-repo apps, and the full default Viewer working. Step 1 corresponds to
[architecture phase 2](../designs/architecture.md#102-technical-phases), steps 2–4 to phase 3.

### Step 0: baselines

- [x] Before step 1 changes which modules register transformers, run `scripts/workspace/registry-dump.mjs` (§4.1) on the
      v6 prototype head for both `DefaultPluginSpec()` and the Viewer with default options.
- [x] Commit `.v6/baselines/registry.json` and `.v6/baselines/transformer-ids.json`, generated at `f088fae56`. The
      prototype is expected to match 5.x apart from the external color themes, which §6.1 lists; a 5.x dump was not
      generated.

### Step 1: module splits

- [x] Remove the `StateTransforms` facade and its lazy getters (and with them the #1791 cycle). Split the seven
      transform modules into this layout (exact file names are settled during the split):

      ```text
      state/transforms/
        catalog.ts   # imports every transform and format module (spec §6.1)
        data/        fetch.ts, json.ts, ...
        structure/   hierarchy.ts, selection.ts, representation.ts, animation.ts,
                     measurement.ts, unitcell.ts, bounding-box.ts,
                     effects/{overpaint,transparency,emissive,substance,clipping,wiggle,theme-strength}.ts
        volume/      ops.ts, representation.ts
        particles/   ops.ts, representation.ts, unitcell.ts
        shape/       representation.ts, box.ts
        misc/        group.ts
      state/formats/
        cif.ts                       # shared ParseCif, ParseBlob
        trajectory/{catalog,mmcif,cif-core,pdb,gro,xyz,lammps,mol,sdf,mol2}.ts
        coordinates/{catalog,dcd,xtc,trr,nctraj,lammps}.ts
        topology/{catalog,psf,prmtop,top}.ts
        volume/{catalog,ccp4,dsn6,mtz,dx,cube,density-server,segmentation,structure-factors}.ts
        shape/{catalog,ply,obj,vtp}.ts
        particles/{catalog,star,tbl,ndjson,em,simularium,mmcif-assembly}.ts
      ```

- [x] Move format-specific parse and conversion transformers next to their providers (`formats/trajectory/sdf.ts` holds
      `TrajectoryFromSDF`, `SdfProvider`, and the `Sdf` entry). Move `VolumeRepresentation3DHelpers`, `getTrajectory`,
      `getBoxMesh`, and the measurement-data helpers out of the large modules. Move format category constants to small
      modules.
- [x] Add `state/transforms/catalog.ts`, which value-imports every built-in transform and format module and exports
      nothing (spec §6.1). Several transformers (`ImportString`, `ImportJson`, `ParseJson`) are imported by nothing
      else.
- [x] Format catalogs become provider arrays and derive the format name types (spec §8); add `BuiltInVolumeFormat` and
      `BuiltInShapeFormat`, keep `BuildIn*` as deprecated aliases. Add the `Id` parameter and `const` helper to
      `DataFormatProvider` and forward it from `TrajectoryFormatProvider` and its siblings.
- [x] Split the preset modules (spec §7): types-and-helpers module plus one module per preset; the auto preset module
      imports the presets it composes. Preset catalogs derive the id and alias unions
      (`BuiltInTrajectoryHierarchyPresetId`/`Alias`, `BuiltInStructureRepresentationPresetId`/`Alias`).
- [x] Split `state/helpers/structure-selection-query.ts` into the `state/queries/structure/` group modules of spec §7,
      each exporting its queries and a group entry, plus `catalog.ts`. Create the residue queries once at module scope
      in `residue.ts`. Move the per-structure `get*Queries` helpers to `dynamic.ts`. Remove the unused
      `applyBuiltInSelection`. `helpers/structure-component.ts`, the `StructureComponent` transform, the component
      manager, and the UI import the group modules they use; the external-structure theme imports `backbone` from
      `structure.ts`. Build `DefaultSelectionQueries` from the individual queries in today's order and check it against
      the step-0 baseline.
- [x] Split behaviors: `BuiltInPluginBehaviors` to its own module; `behavior.ts` becomes only the `PluginBehaviors`
      catalog; base modules import `PluginBehavior` from `@molstar/plugin/behavior/behavior`; library code and slim apps
      import behaviors from their defining modules. Add the focus representation id module
      `behavior/dynamic/selection/structure-focus-representation/id.ts`.
- [x] Make `StateActions` a catalog module.
- [x] Move `BuiltInMarkdownExtension` into catalog entries; the `query` extension is its own entry and imports no
      transpiler.
- [x] Split `@molstar/model/script/script` (spec §7): the core module handles `mol-script` only; `transpile.ts` holds
      the module-global transpiler table and no longer imports `transpilers/all.ts`; add
      `transpilers/{pymol,vmd,jmol}.ts` registering modules (the parsers stay in the `<lang>/` directories), and make
      `transpilers/all.ts` import all three instead of exporting `_transpiler`. `parse` throws the spec §7
      script-language error for an unregistered language; the script-param UI lists registered languages only. The
      direct `parse` users in §1.5 import `transpilers/pymol`. Add the `transpilers/all` import to the default spec
      modules now (the default registry takes it over in step 2) and to the Viewer, so no consumer loses a language.
- [x] Move the graphics built-in maps to catalog modules (`@molstar/graphics/repr/structure/catalog`,
      `@molstar/graphics/theme/color/catalog`, ...) wrapped in `namedCatalog`; namespace types type-import them;
      namespace values are removed. Constructors kept preloading until step 3 removed the preloads.
- [x] Move `DefaultPluginSpec` to `@molstar/plugin/default-spec` and `DefaultPluginUISpec` to
      `@molstar/plugin-ui/default-spec`; the `spec` modules keep only types and helpers.
- [x] Land transformer-id snapshot validation (`validateSnapshotTransformers`, spec §6.2) in this step, since the split
      changes which modules register transformers.
- [x] Migrate every consumer of the facade, the moved default specs, and moved modules, including `cli/mvs-render`, the
      headless examples, and the smoke fixtures (§3.6, §3.7).
- [x] Switch the extension and app calls of `PresetStructureRepresentations.auto.apply` in §1.3 to the auto preset's
      defining module.
- [x] Build the import-graph check with the catalog manifest (§4.2, §4.3), run from `check:workspace`; later steps
      extend its rules. First version: rules a (catalog value imports), b (base entry points reach default-composition
      modules) and d (type imports to a higher package), with the planned violations in the manifest allowlist (each
      with a reason and the removing step; stale entries fail).

### Step 2: registries

- [x] Reference counting in every registry and manager with the per-registry keys and duplicate rules of spec §4.2;
      `remove` of an unknown provider is a no-op. `clear()` drops providers and counts.
- [x] `has(nameOrProvider)` on representation and theme registries; `RepresentationRegistry.default` typed as possibly
      `undefined`; `ThemeRegistry` without a built-in map.
- [x] `DataFormatProvider.withName` (copies memoized in a module-level `WeakMap<provider, Map<name, provider>>`) and
      `DataFormatProvider.Unnamed`, and `DataFormatRegistry` as in spec §4.3: `add(provider)` uses `provider.name`;
      `add(name, provider)` registers `provider` when the names match, otherwise the named copy with a deprecation
      warning; `remove(provider)` decrements by identity; `remove(name)` resolves the provider, then decrements, and is
      a no-op for an unknown name; `list` keeps `{ name, provider }` items; `has(name)` is added. Apply the dead-check
      rewrites in §1.4. (`DataFormatProvider.name`, the `Id` parameter, and names on every in-repo provider landed in
      step 1.)
- [x] `PluginAnimationManager.unregister`; `PluginDragAndDropEntry`, `fallback` ordering, and `addHandler` options.
- [x] Cached `toAction()` and id-counted `StateActionManager` in core.
- [x] `plugin.register` with atomic conflict checks and idempotent undo; `spec.registry` registered in `init()`.
- [x] Preset builders resolve through the registry by `id`, then `alias`: index `alias` in the builders with the §4.2
      conflict rule, remove `defaultProvider` and the static-map lookup, add the id/alias-union and `string` overloads
      (`TrajectoryHierarchyBuilder` gains `string`), throw for an unresolved string, and fix `getPresetSelect` defaults.
      Verify that every §1.3 call site, smoke fixture, and documentation example still resolves. (Step 1 already added
      the optional `alias` to `PresetProvider`, set it on every built-in preset, derived the id/alias unions in the
      catalogs under `state/builder/structure/{representation,hierarchy}-presets/`, and typed
      `LoadTrajectoryParams.preset` with the hierarchy alias union.)
- [x] Representation-name handling (spec §4.3): the nine param-definition sites fall back to empty mapped params for
      empty registries; the eight helper sites and `buildRepresentation` check names with `has` and warn with provider
      and scope; they throw for empty registries when given data, and return empty mapped params when called without
      data (as `StructureFocusRepresentation` and `extensions/meshes/src/examples.ts:148` do); the three representation
      transformers throw in `apply`/`update` for empty registries. No core-state change.
- [x] `PluginConfig.Structure.DefaultHierarchyPreset`.
- [x] Add catalogs and `DefaultRegistry` (spec §5.3) and register it alongside the existing constructor preloads;
      identity counting makes the duplicates harmless, provided the `StructureSelectionQueryRegistry` constructor pushes
      the catalog's module-scope residue queries instead of creating new ones (otherwise each residue query would be
      listed twice until step 3). `DefaultThemes` includes the `ExternalColorThemes` entry (spec §5.4), exported from
      `themes/external-*`.
- [x] Remove the three spec fields with the constructor check (spec §3.1); convert the Viewer and mesoscale
      `customFormats` (§3.1, §3.2); migrate every in-repo spec (Viewer, docking-viewer, mesoscale-explorer, mvs-stories,
      examples, smoke fixtures, `cli/state-docs`, `cli/mvs-render`).
- [x] In the same change that removes `spec.animations`, delete `initAnimations` and its implicit snapshot-transition
      registration (spec §4.4) and type `current` as possibly undefined; otherwise the implicit registration would
      always fire first and make the snapshot transition the default animation instead of `AnimateModelIndex`.
- [x] `setSnapshot` skips `snapshot.current` with a warning for an unregistered animation name.
- [x] Separate fix: count `addCustomProp`/`removeCustomProp` and `addSymbol`/`removeSymbol` in
      `DefaultQueryRuntimeTable`.

### Step 3: context decoupling

- [x] Check that every in-repo spec lists `DefaultRegistry` or its own entries, then remove the constructor preloads and
      implicit imports of §1.1 (spec §9). Also delete `PluginContext.dropOverriddenPreloadedFormats` (the transitional
      handling that lets a `customFormats` override replace a preloaded built-in format), and have mesoscale-explorer
      list only the formats, representations, and themes it needs. Done: the registries, builders, selection-query
      registry, markdown manager, and `ColorTheme`/`SizeTheme.createRegistry()` start empty and import no catalog; the
      import-graph allowlist is empty, and the base entry points may reach no catalog (rule c).
- [x] Remove the volume-streaming behavior from the `PluginContext` closure. After step 1 it is still reached through
      `state/manager/structure/hierarchy.ts` → `hierarchy-state.ts` → `behavior/dynamic/volume-streaming/behavior.ts`,
      and `behavior/dynamic/volume-streaming/util.ts` through `formats/registry.ts` → the volume catalog →
      `volume/ccp4.ts` → `volume/provider.ts`. Move the shared helpers and type checks to modules that do not import the
      behavior.
- [x] `state/builder/structure/representation-presets/types.ts` value-imports the selection-query catalog for
      `presetSelectionComponent` (keyed lookup by catalog key); import the individual queries or take a query object, so
      preset helpers do not load every query (found by the step 1 import-graph check).
- [x] Presets import what they run and list it in their entries; add representation entries with their default themes in
      plugin-layer modules; mixed name/provider props in the helpers (spec §5.1); development-mode default-theme check.
- [x] Format entries with actions and, for volume, particle, and shape formats, the representation entries and themes
      their `visuals` apply as provider objects (CCP4 includes `Isosurface`). List only actions in today's default
      action list, for example CCP4: `ParseCcp4`, `VolumeFromCcp4`; mmCIF: `ParseCif`, `TrajectoryFromMmCif`; SDF: none.
- [x] Policy through config for every site in §1.6. When the configured preset is not registered, loading fails with the
      preset error.
- [x] Focus representation behavior: its `nciParams` group passes `InteractionsRepresentationProvider` and
      `InteractionTypeColorThemeProvider`, which only the toggleable `CustomProps.Interactions` behavior registers (and
      `DefaultPluginSpec` lists after the focus behavior). Its param defaults are data-less and do not assert
      registration. Before applying the `interactions` focus component it checks
      `registry.has(InteractionsRepresentationProvider)` and skips the component when it is absent, instead of
      substituting another representation. Replace the `SizeTheme.BuiltIn.uniform` read with an import.
- [x] `updateFocusRepr` (spec §7): stop naming `ball-and-stick` and `element-symbol` (which would warn and substitute on
      every preset in a slim plugin) and value-importing the focus behavior; read the behavior's current representation
      type from its params, skip when that type or the requested theme is not registered (`has`), and reach the behavior
      through `plugin.state.hasBehavior`/`updateBehavior` with the transformer id from the id module and an
      `import type` for its params.
- [x] `StructureComponentManager`: move the interactions code out. Keep the `interactions` option slot so
      `structureComponentManager.options.interactions` round-trips through snapshot JSON; take its param definition from
      a params-only module instead of `InteractionsProvider`; apply it through a handler the Interactions
      custom-property behavior registers with `managers.structure.component.registerOptionHandler('interactions', ...)`,
      which returns an undo. Without the behavior the option is preserved but ignored. The focus representation update
      in `setOptions` goes through the same id module and runs only when `plugin.state.hasBehavior(id)`; today
      `updateBehavior` inserts the focus behavior when it is absent. Pick the default query by identity (`current` if
      registered, else the first option) instead of `options[1][0]`, and return an empty select for an empty registry.
- [x] `ViewerAutoPreset` and other delegating presets look up optional presets by id and skip them when absent.
- [x] Move the open-files drag-and-drop fallback into `DefaultDragAndDrop`.
- [x] `parseTrajectory(blob)` moves into the mmCIF entry's module or becomes format-neutral; the Cube structure path
      imports its transforms. `DownloadStructure` (spec §9), when building params, offers only the sources whose formats
      are registered (`dataFormats.has(name)`), builds the URL format list from the registry, and hides the multi-source
      blob option unless mmCIF is registered.
- [x] `PluginViewModel`, `PluginUIViewModel`, and their hooks require an explicit spec.
- [x] UI (spec §10): `createPluginUI` requires a spec; base layout renders minimal structure tools; full tools move to
      `DefaultPluginUISpec().components.structureTools`; quick styles resolve presets by id; volume controls use the
      registry or config; apply the `components` fixes in §3.5.
- [x] Add `reportUnregisteredNames` to `setSnapshot` (spec §6.2).

### Step 4: slim acceptance and boundary checks

- [x] Write the new `BallAndStickPreset` (`preset-structure-representation-ball-and-stick`) and the slim example of spec
      §12. Done: `representation-presets/ball-and-stick.ts` (`BallAndStickPreset`, `BallAndStickPresetEntry`, not in
      `DefaultPresets`) and `examples/slim-plugin`.
- [x] Pin the excluded-module list in the import-graph check; add the esbuild metafile checks for a split and a
      single-file build. Done: `scripts/workspace/slim-exclusions.json` (shared by `import-graph.mjs` rule e and
      `slim-bundle.mjs`, which `check:workspace` runs); the import-graph check now follows `verbatimModuleSyntax`
      (`import { type A }` is a value import).
- [x] Rendering smoke test. Done: `node smoke/run.mjs slim` (`smoke/slim/`, included in `all`, `smoke:slim`) installs
      the packed packages the example imports, bundles a copy of the example sources with esbuild and the `molstar-src`
      condition, and drives the page with Playwright: no page or console errors, a visible ball-and-stick representation
      with render objects, exactly the slim providers, the unregistered-`cartoon` snapshot warning and default, and the
      PyMOL script error.
- [x] Add an import-graph assertion that `themes/external-structure`, `themes/external-volume`, and
      `state/queries/structure/*` reach no transform module, `PluginContext`, or catalog by value (spec §5.4). Done:
      rule f, from the `boundary` section of `slim-exclusions.json`.
- [x] Fix the leaks the checks found: `MmcifFormat` moved to `formats/structure/mmcif-format.ts` (the mmCIF parser
      module `mmcif.ts` re-exports it and keeps the property-provider registrations);
      `model/structure/export/categories/utils.ts` imports `getCifFieldType` from `data-model` instead of `reader/cif`;
      `state/builder/structure/{representation,hierarchy}.ts` used `import { type ... }` of the preset catalogs, which
      `verbatimModuleSyntax` keeps as a side-effect import (now `import type`).
- [x] Complete the acceptance checklist (§6): audited in step 4; see the evidence under §6.

### Step 5: extensions and MVS

- [x] Move extensions to `plugin.register` where it simplifies them (§3.8).
- [ ] Move MVS to `plugin.register` (§3.9).

## 3. Per-consumer migration

`check:publish` runs or builds `cli/mvs-render`, `examples/image-renderer`, `examples/glb-export`, the
`extensions/plugin` view models, and the smoke fixtures, so step 1 and step 2 must migrate them (§3.5–§3.8).

### 3.1 Viewer (`apps/viewer`)

- [ ] Spec: `registry: [...DefaultRegistry, ViewerEntry, customFormatsEntry]`, with the custom-formats entry last so
      built-in format order and `auto()` tie-breaking are unchanged. `ViewerAutoPreset` may move into `ViewerEntry` or
      stay registered in `onBeforeUIRender`; its entry includes only the auto preset it falls back to.
- [x] Convert `customFormats` as in spec §11: type it `[name: string, provider: DataFormatProvider.Unnamed][]` and
      register `DataFormatProvider.withName(provider, name)` for each tuple. `G3dProvider` gains `name: 'g3d'`, so the
      default tuple passes through unchanged.
- [x] Overriding a built-in format: in 5.x a tuple whose name matches a built-in format (`['pdb', MyPdbProvider]`)
      overrides it for `dataFormats.get(name)`, because the custom formats are added after the built-ins and the map
      keeps the last one. In 6.0 that would be a different provider under an existing key, which rejects `init()` and
      would break `Viewer.create`. When a `customFormats` name equals a format in `DefaultFormats`, the Viewer registers
      `{ ...DefaultFormats, formats: DefaultFormats.formats.filter((p) => !customNames.has(p.name)) }` in its place,
      then appends the custom-formats entry. `get(name)` returns the custom provider as in 5.x.
- [ ] `loadTrajectory({ preset })` keeps accepting the hierarchy short keys through the preset aliases, and also accepts
      ids.
- [ ] Import `@molstar/model/script/transpilers/all` at the top of the Viewer entry.
- [ ] `molstar.lib.plugin` (`apps/viewer/src/lib.ts`): keep `StateTransforms` as an app-level object literal with the
      same `Data`/`Misc`/`Model`/`Particles`/`Volume`/`Representation`/`Shape` keys and member names, assembled from the
      split modules; no library module imports it. Keep `StateActions`, `DefaultPluginSpec`, and `DefaultPluginUISpec`
      as values from their new modules, and add `DefaultRegistry` so script-tag code can migrate a literal spec with
      `registry: molstar.lib.plugin.DefaultRegistry`. Record the carve-out in the ledger.

### 3.2 mesoscale-explorer

- [x] Keep the `customFormats` option and convert it like the Viewer, including the override handling.
- [x] Replace the literal spec (`actions: defaultSpec.actions`, three animations) and the post-`init()`
      `registry.clear()`/`lociLabels.clearProviders()` with an explicit registry: spacefill with its themes, the
      `uniform` and `illustrative` structure color themes it applies by name (`data/state.ts`, `ui/entities.tsx`,
      `ui/states.tsx`), any other theme it names, its formats and actions, and its three animations. Today it clears
      only the structure representation registry and keeps every theme. It sets themes through direct state updates and
      `StructureRepresentation3D` applies (`data/*/preset.ts`), not the builders or helpers, so a theme it names but
      does not register would be substituted with the registry default silently (spec §4.3); the registry must list
      every theme it names. Done: it lists `Spacefill` (element-symbol and `physical` themes), `uniform` and
      `illustrative`, the mmCIF entry (`Mmcif`, with its actions), its three animations, and the custom formats; the
      loci-label clearing was a no-op (no listed behavior adds a provider) and is gone, and it lists no other actions
      (its own actions are applied directly, not looked up in the registry).

### 3.3 docking-viewer

- [x] Add `DefaultRegistry` (or a narrower list) to its literal spec.

### 3.4 mvs-stories

- [x] It spreads `DefaultPluginUISpec()`, so it inherits `registry` and needs no registry change. Move the
      `DefaultPluginUISpec` import to `default-spec` and fix `components` in `elements/viewer.tsx` (§3.5).

### 3.5 Examples

- [x] `components` overrides: `examples/proteopedia-wrapper`, `examples/basic-wrapper`, `examples/interactions`,
      `examples/ihm-restraints`, `examples/ligand-editor`, and `apps/mvs-stories` (`elements/viewer.tsx`) spread the
      default spec but replace `components` wholesale; each spreads `...defaultSpec.components` into its `components`,
      or it would silently get the minimal tools.
- [x] `proteopedia-wrapper` replaces `DefaultAnimations` instead of overriding `animations`.
- [ ] `examples/interactions` and `examples/basic-wrapper` use preset ids (§1.3).
- [ ] `examples/image-renderer` and `examples/glb-export`: default specs from `default-spec`, leaf transformer imports.

### 3.6 Smoke fixtures

- [x] `smoke/headless/capture.mjs`, `smoke/browser/viewer/index.html`, `smoke/browser/library/index.html`,
      `smoke/node/runtime.mjs`, and `smoke/types/consumer.tsx` import the default specs from `default-spec`, replace
      `spec.actions` checks with `plugin.state.data.actions` lookups, and replace `StateTransforms` with leaf imports.
- [ ] `smoke/headless/capture.mjs` and `smoke/browser/library/index.html` rewrite `applyPreset(trajectory, 'default')`
      to `'preset-trajectory-default'` (§1.3).

### 3.7 CLI tools

- [x] `cli/state-docs` imports the transformer catalog instead of the facade, builds its context from
      `DefaultPluginSpec()`, and awaits `init()`, so the generated docs keep listing every built-in representation and
      theme.
- [x] `cli/mvs-render`: same treatment as the smoke fixtures.

### 3.8 Extensions

- [x] Extensions that add an action in `register()` and remove it in `unregister()` (assembly symmetry, SB-NCBR tunnels,
      g3d, volumes and segmentations; MVS is §3.9) need no change once actions are counted; verified. The four now
      register the action inside their entry (below), so the count is one per registration and the spec's own entry
      keeps the action when the behavior is toggled off. `apps/viewer/src/_test/extension-registry.test.ts` pins this
      for assembly symmetry, tunnels and g3d (count 1 with the behavior, 2 with an app entry listing the same action, 1
      after the behavior is removed, 0 after the entry's undo);
      `extensions/volumes-and-segmentations/src/_test/behavior.test.ts` does the same for `LoadVolseg`.
- [x] `extensions/plugin` view models and hooks require an explicit spec; the hooks' `spec: (defaultSpec) => spec` form
      is removed, and callers import the default themselves.
- [x] Replace catalog value imports with defining-module imports (§1.3): nothing remained in `extensions/` (the
      `AutoPreset` imports were switched in step 1; `extensions/plugin` only has `import type` of catalog modules, which
      the import-graph check allows).
- [x] Replace hand-written provider registration with one `plugin.register(entry)` call and its undo, built in
      `register()` (so it can close over the behavior, as the loci label providers do) and undone in `unregister()`:
      anvil (representation, selection query, preset), assembly symmetry (cluster color theme, preset, action), dnatco
      (color themes, representations, presets), g3d (loci label, action), kinemage (`KIN` format, drag-and-drop
      handler), model-archive quality assessment (color themes, query, presets, loci label), PDBe structure quality
      (color theme, loci label), RCSB validation report (color themes, representation, query, presets, loci label),
      SB-NCBR partial charges (color theme, preset, loci label) and tunnels (preset, action), volume tools mask and
      segmentor (volume color themes), volumes and segmentations (action), wwPDB CCD (hierarchy and representation
      presets). Each entry lists its providers in the order the hand-written code registered them, so registry order is
      unchanged (the registry comparison in `check:workspace` passes). Public exports are unchanged.
- [x] Left imperative, as the spec (§3.3) says or because an entry would not simplify them: custom model, structure and
      volume properties with their `autoAttach` retuning in `update()`, `DefaultQueryRuntimeTable` custom properties and
      symbols (RCSB, PDBe, MA, g3d, anvil, kinemage), UI maps (`customStructureControls`,
      `genericRepresentationControls`, `customImportControls`: assembly symmetry, anvil, kinemage, MA pairwise plot,
      geo/model/mp4 export, zenodo import, volseg), the backgrounds config, debug helpers (canvas debug registry), the
      mesh streaming behavior, and the volseg entry's loci label provider (`entry-root.ts`), which is created and
      removed with each loaded entry rather than by the behavior. The geo, model and mp4 export and the Zenodo import
      behaviors only register UI map entries, so they have nothing to put in an entry.
- [x] Jest loads `.tsx` modules (`tsx` module extension, the `react-jsx` transform, an image stub) and maps
      `@molstar/volseg-api-extension`, so a test can load the extension behaviors that import their UI.

### 3.9 MVS (`packages/mvs/runtime`)

- [ ] `MolViewSpec.register()` builds one entry from the provider part of `Registrables`
      (`mvs/runtime/src/behavior.ts`): representations, color themes (including the multilayer theme that
      `makeMultilayerColorThemeProvider` builds from this plugin's registry), loci labels, the drag-and-drop handler
      (already `{ name, handle }`), the `MVSJ`/`MVSX` formats (providers gain `name`; the `{ name, provider }` wrapper
      goes), and its actions. It passes the entry to `plugin.register` and calls the undo in `unregister()`.
- [ ] Custom model and structure properties (registered with the `autoAttach` param and retuned in `update()`) and the
      state and markdown ref/URI resolvers stay imperative.
- [ ] Add `MVSRuntimeRegistry` in a catalog module listed in the manifest: every representation and color/size theme
      provider MVS can name. MVS names them in `load-helpers.ts` and in `load-extensions/non-covalent-interactions.ts`
      (`interactions` and `interaction-type` from the non-covalent-interactions extension).

## 4. Tooling

### 4.1 Registry baseline dump

`scripts/workspace/registry-dump.mjs` (run with `node`, after `pnpm build:lib`) writes `.v6/baselines/registry.json` and
`.v6/baselines/transformer-ids.json` for two targets: `default` (`new PluginContext(DefaultPluginSpec())`) and `viewer`
(`new PluginUIContext(createViewerSpec({}))` plus the `ViewerAutoPreset` registration that `Viewer.create` does). It
runs once in step 0 and again for the acceptance comparison; `--out <dir>` writes elsewhere.

- Each target runs in its own Node process after `init()`, so the dump includes what behaviors register and the global
  transformer registry holds only what that target imports. No UI is rendered and no canvas is created.
- Node needs two shims for browser-oriented modules: `globalThis.window = globalThis`, and a module load hook that turns
  image and style imports into URL strings.
- `registry.json` records, per target: representations and color/size themes per scope (registry order), format names,
  hierarchy and representation preset ids, selection queries (`category: label`), the number of loci label providers
  (they have no names), markdown extension, drag-and-drop, and animation names, data and behavior state actions (all,
  and per `from` type), spec behavior ids, and custom model/structure/volume property names. Actions have random UUIDs,
  so a transformer-derived action is keyed by its transformer id and any other action by `action:<display name>`.
- `transformer-ids.json` records the sorted ids from `StateTransformer.getAll()`.
- Both files include the commit they were generated from; the rest of the output is deterministic.
- `scripts/workspace/registry-compare.mjs` (run by `check:workspace`) dumps both targets again and compares them with
  the baseline, allowing only the §6.1 differences.

### 4.2 Catalog manifest

`scripts/workspace/catalogs.json` lists every catalog module (`behavior.ts`, `state/actions.ts`, the preset and query
catalogs, the `catalog.ts` files, the graphics catalogs, `MVSRuntimeRegistry`'s module) and the three
default-composition modules; apps and base modules are classified as in spec §8. A new catalog module is added to the
manifest in the same change.

### 4.3 Import-graph check

`scripts/workspace/import-graph.mjs`, run from `check:workspace`:

- [ ] Rejects value imports of catalog modules outside catalog modules, default-composition modules, and apps.
- [ ] Rejects type imports that point to a higher package.
- [x] Rejects value imports of default specs or catalogs from the base entry points (`@molstar/plugin/context`,
      `@molstar/plugin/spec`, `@molstar/plugin-ui`, `@molstar/plugin-ui/spec`), including the modules split in step 1.
- [ ] Includes extensions, servers, and CLI packages.
- [x] Checks the slim example's graph against the excluded-module list (spec §12), including UI modules (rule e), and
      the spec §5.4 boundary (rule f).

## 5. Migration-tool additions

Add to `@molstar/migrate-6-cli` ([architecture §9.2](../designs/architecture.md#92-migration-tool)):

- [ ] Rewrite literal `actions`/`animations`/`customFormats` spec fields into `registry` entries.
- [ ] Rewrite literal specs that relied on the implicit catalog to include `DefaultRegistry` (spec §3.1).
- [ ] Rewrite spread-and-override specs to filter the matching default entry.
- [ ] Report spec objects it cannot classify.
- [ ] Split combined `{ DefaultPluginUISpec, type PluginUISpec }` imports.
- [ ] Map the moved name-type, category-constant, and provider paths.
- [ ] When a spec spreads `DefaultPluginUISpec()` and sets a `components` literal without spreading the default
      components, insert `structureTools: DefaultPluginUISpec().components!.structureTools`, or report the spec when it
      cannot.
- [ ] Rewrite short-key `PresetTrajectoryHierarchy`/`PresetStructureRepresentations` type uses to the id and alias
      unions. Short-key string values need no rewrite; they resolve through aliases.
- [ ] Rewrite `PresetStructureRepresentations.<key>` and other catalog value uses in library code to defining-module
      imports, and report what it cannot resolve.
- [ ] Report PyMOL, VMD, or Jmol script use without a transpiler import.
- [ ] Emit the `DefaultRegistry` export names fixed in spec §5.3.

## 6. Acceptance checklist

- [x] The slim example renders an SDF ligand, and the import graph and the split and single-file bundles exclude every
      module in the spec §12 exclusion table (`pnpm check:workspace`; rendering checked by `node smoke/run.mjs slim`).
- [x] The providers registered by `DefaultPluginSpec`, and by the Viewer with default options, match the step-0 baseline
      of each registry (names and order), apart from the differences in §6.1 and the Viewer's own entries.
- [x] The `StateTransformer` ids registered after importing `DefaultPluginSpec` contain the baseline id set.
- [x] Existing snapshots restore in the Viewer. A snapshot needing an unimported transformer fails before changing
      state.
- [x] Built-in preset short keys resolve through their aliases, including `loadTrajectory({ preset: 'all-models' })`.
      Loading an unregistered format or an unknown preset id or alias produces a clear error.
- [x] Building a representation through the builders or helpers with an unregistered type, color, or size name logs a
      warning naming the provider and scope and renders the registry default. Restoring a snapshot that names an
      unregistered representation or theme reports it the same way.
- [x] With an empty representation, color theme, or size theme registry in a scope, the helpers and a direct `apply` of
      `StructureRepresentation3D`, `VolumeRepresentation3D`, or `ParticlesRepresentation3D` fail with a clear error.
- [x] In the slim app, loading a snapshot with a cartoon representation reports `cartoon` as unregistered and renders
      the registry default (`node smoke/run.mjs slim`).
- [x] A PyMOL script in the slim app fails with the script-language error (`node smoke/run.mjs slim`); the default
      plugin and the Viewer evaluate it (`script-languages.test.ts`, `composition-check.mjs`).
- [x] `register` with a conflicting entry throws and changes nothing; the undo is idempotent; a behavior's
      `unregister()` does not remove a provider the spec also registered.
- [x] `Viewer.create` with a `customFormats` tuple that overrides a built-in name (for example `['pdb', MyPdbProvider]`)
      loads that format with the custom provider.
- [x] `viewer.loadTrajectory({ ..., preset: 'all-models' })` still works.
- [x] mesoscale-explorer, docking-viewer, mvs-stories, proteopedia-wrapper, the examples and smoke fixtures in §3,
      `cli/mvs-render`, and `cli/state-docs` build and behave as before.

Left open. The slim-app items (the cartoon snapshot and the PyMOL script) are covered by the slim browser smoke test,
not by this audit; the default plugin and Viewer half of the PyMOL item is verified below. The last item is verified for
the builds, the type checks, the smoke fixtures, `cli/mvs-render` and `cli/state-docs`, but nothing exercises the
mesoscale-explorer, docking-viewer, proteopedia-wrapper, or the other examples at run time: their specs are built inside
`create`/`main` functions that need a browser, and the native headless smoke (`smoke/headless`, needs `gl`) was not run.

#### Evidence

Run with `pnpm build:lib` first. Jest tests are named by file and test; `pnpm check:workspace` runs the two new
workspace checks (`registry-compare.mjs`, `composition-check.mjs`) after the existing ones.

1. Slim example: `scripts/workspace/import-graph.mjs` rule e/f, `scripts/workspace/slim-bundle.mjs`, and the browser
   render, as recorded under step 4.
2. Registry baseline: `scripts/workspace/registry-compare.mjs` dumps both targets (`registry-dump.mjs`: default spec,
   and the Viewer spec plus `ViewerAutoPreset`) and compares them with `.v6/baselines/registry.json`. The only
   differences are the two in §6.1: `external-structure` and `external-volume` in the color themes of all three scopes
   (both targets) and `open-files` in drag and drop (both targets); the Viewer adds nothing of its own. Complements:
   `_test/default-registry.test.ts` ("registers into a plugin with the default spec", "DefaultThemes includes the
   external color themes in all three scopes") and `_test/empty-registries.test.ts` ("the default spec registers the 5.x
   provider sets in the 5.x order").
3. Transformer ids: the same script checks that every id of `.v6/baselines/transformer-ids.json` is still registered
   (the sets are identical today).
4. Snapshots: `packages/plugin/core/src/state/_test/snapshot-restore.test.ts` ("rebuilds the data tree and the
   representations", "restores through the snapshot manager entry format as well", "fails before changing anything when
   a transformer is not imported") builds a snapshot from a real default plugin (crambin, JSON round trip) and restores
   it into a fresh one. `scripts/workspace/composition-check.mjs` restores the same snapshot in a plugin with the real
   `createViewerSpec` ("a default-plugin snapshot restores in a plugin with the Viewer spec", "a snapshot needing an
   unregistered transformer fails before changing the Viewer plugin"). Earlier unit tests:
   `_test/snapshot-validation.test.ts` ("fails before any side effect when a transformer is not registered") and
   `manager/_test/snapshots.test.ts` ("leaves the manager unchanged when any entry names an unregistered transformer").
   The repository has no `.molj` from before the refactor, so the snapshot comes from this tree; the transformer id set
   (item 3) and the unchanged snapshot format are what make older snapshots restore.
5. Aliases and errors: `state/builder/structure/_test/preset-aliases.test.ts` (every catalog key and id resolves in the
   default plugin; "give a clear error for an unknown id or alias"), `preset-names.test.ts`, `preset-registry.test.ts`
   ("resolves by id, then alias", "throws for an unresolved string"); `extensions/plugin/src/_test/loaders.test.ts`
   ("applies the 'all-models' preset through its alias, with one structure per frame", "fails with a clear error for an
   unknown preset", "fails with a clear error naming the unregistered format"); `_test/empty-registries.test.ts` ("does
   not use a built-in provider when a name is requested").
6. Unregistered names: `state/builder/structure/_test/unregistered-names.test.ts` (a real plugin, `addRepresentation`
   with an unregistered type, color and size warns with scope and kind and stores the registry defaults),
   `state/helpers/_test/representation-registry.test.ts` ("warn for an unregistered structure representation...", "warn
   in the volume helpers"), `state/_test/report-unregistered-names.test.ts`, and the cartoon-less restore in
   `snapshot-restore.test.ts` ("reports a representation that is not registered and renders the registry default").
7. Empty registries: `state/helpers/_test/representation-registry.test.ts` ("representation param helpers with empty
   registries") and `state/transforms/_test/representation-registry.test.ts` ("apply throws without representations",
   "apply and update throw without color themes", "... without size themes", for the structure, volume and particles
   transformers).
8. Scripts (default plugin and Viewer): `_test/script-languages.test.ts` (the default plugin evaluates PyMOL, VMD and
   Jmol scripts), `packages/model/src/script/_test/languages.test.ts` (without a transpiler import only `mol-script` is
   available and the others throw "is not available in this build"; `transpilers/all` enables all), and
   `composition-check.mjs` ("the Viewer spec evaluates a PyMOL script"). The slim half stays with the browser smoke.
9. `register`: `_test/register.test.ts` ("throws one error listing every conflict and changes nothing", "does not
   register anything before a conflict later in the input", "decrements exactly the registrations it made and is
   idempotent", "does not let a behavior-style removal remove a provider the spec also registered").
10. Viewer `customFormats`: `apps/viewer/src/_test/registry.test.ts` ("lets a custom format override a built-in one, as
    in 5.x", "loads data with the custom provider of an overridden built-in format") and `composition-check.mjs`
    ("Viewer customFormats overrides a built-in format and loads the data with the custom provider", with
    `createViewerSpec`).
11. `loadTrajectory`: `loaders.test.ts` as in item 5, and `composition-check.mjs` ("Viewer loadTrajectory({ preset:
    'all-models' }) works").
12. Consumers: `pnpm build:apps` builds the Viewer, mesoscale-explorer, docking-viewer, mvs-stories and every example
    (including proteopedia-wrapper, basic-wrapper, interactions, ihm-restraints, ligand-editor, and slim-plugin);
    `pnpm check:types:full` covers the rest of the workspace; `smoke/run.mjs node`, `types`, `source`, `cli` and
    `browser` (Viewer and library pages, classic Viewer and MVS Stories) pass on the packed artifacts;
    `composition-check.mjs` runs `cli/mvs-render --help` and `cli/state-docs`, whose output must list every registered
    representation and theme (`cli/state-docs/src/_test/state-docs.test.ts` does the same in jest).

Failures the audit found and fixed: `cli/state-docs` threw with the default plugin because three param getters
dereferenced their data when called without it (`getOrientationParticlesParams`, `getParticleTargetParams`,
`getOperatorHklColorThemeParams` with `Structure.Empty`); they now accept the missing data.

- Item 13, runtime: on the step 4 tree, the built Viewer (1CBS), mesoscale-explorer (local 1CRN mmCIF via `url`),
  docking-viewer (local `ace2.pdbqt` + `ace2-hit.mol2`), mvs-stories (`examples/mvs/kinase-story.mvsj`), and the
  proteopedia-wrapper, basic-wrapper, interactions, ligand-editor, ihm-restraints, lighting, alpha-orbitals, react,
  alphafolddb-pae, volume-tools, and slim-plugin examples were opened in a browser from a static server. All loaded and
  rendered without plugin warnings or errors. Unrelated: the alphafolddb-pae example's PAE request
  (`..._predicted_aligned_error_v4.json`) now returns 404 from AlphaFold DB; mesoscale-explorer's deploy-time extras
  (`../extras/driver.*`, `../examples/list.json`) are absent from a plain build.

### 6.1 Intentional registry differences

The registry dump comparison allows only these differences. A step that introduces another one adds it here.

| Registry      | Difference from the baseline                                                                                      |
| ------------- | ----------------------------------------------------------------------------------------------------------------- |
| Color themes  | `external-structure` and `external-volume` are present in all three scopes, as in 5.x; the prototype dropped them |
| Drag and drop | The open-anything handler is listed as `open-files` (`DefaultDragAndDrop`, fallback) instead of being hard-coded  |
