# Mol\* 6.0: plugin composition

Design specification, 2026-10-06. [Architecture §6](architecture.md#6-plugin-composition) and the
[summary](summary.md#plugin-composition) summarize this design and link here. The
[implementation plan](../plans/plugin-composition.md) holds the ordered steps, call-site inventories, per-consumer
migration, tooling, and the acceptance checklist. This spec describes the target; none of it is implemented yet. Paths
without a package prefix are relative to `packages/plugin/core/src/`.

## 1. Goals and settled decisions

Goals:

- **Include only what is imported and listed.** An application's plugin contains the formats, representations, themes,
  presets, and other providers that its spec lists. Those values come from static imports, so the bundle contains them
  and nothing else.
- **Cheap plugin creation.** Given a spec, `PluginContext` registers the listed entries in order. It does not resolve
  dependencies, look up providers by string name, or load catalogs.
- **Unchanged defaults.** `DefaultPluginSpec`, `DefaultPluginUISpec`, and the Viewer keep today's providers, order, and
  defaults.

Settled decisions. Later sections implement them; where a section adds a mitigation, the decision itself is unchanged.

1. `PluginSpec` gets `registry?: PluginRegistryEntry[]`. Entries are pure data: providers only, with no behaviors,
   config, or functions run at registration. There is no `requires` and no dependency resolution; what is included is
   what the spec file imports and lists. Any dynamic resolution would live above `PluginContext`, and none is planned.
2. `PluginSpec` keeps only `behaviors`, `config`, `canvas3d`, and `layout` besides `registry`. Spec-level `actions`,
   `animations`, and `customFormats` move into entries and are removed with no aliases. `DataFormatProvider` gains a
   required `name`.
3. The public `Viewer.create` API is not broken. Its `customFormats` option keeps the `[name, provider]` shape, and the
   Viewer converts it into a registry entry.
4. Behaviors stay as they are (state tree, snapshots, toggling). `plugin.register(entry)` returns an undo function so
   behaviors can reuse the entry format; migrating to it is optional.
5. Registries start empty. Registration is by provider identity with reference counting. A different provider under an
   existing key is an error.
6. Transformers stay in the global `StateTransformer` registry and register when their module is imported. They are not
   listed in the spec. Snapshot loading checks all transformer ids before changing state.
7. Presets statically import the providers they run and include them in their entry. Policy choices, such as the preset
   a format applies after loading, go through config and the registries by id.
8. The `StateTransforms` facade is removed. Transform modules are grouped by functionality, with subgroups;
   format-specific transformers live next to their format providers.
9. `BuiltInTrajectoryFormat` and sibling name types are kept. They are derived in catalog modules and imported elsewhere
   only with `import type`, which `verbatimModuleSyntax` erases.
10. Static single-file builds (one IIFE or ESM file) work: static imports only, synchronous registration, and no runtime
    selection by string name except where an app deliberately imports a full catalog.
11. Unregistered representation and theme names are lenient: parameter normalization substitutes the registry default as
    today, for direct applies and snapshot restore alike, and the plugin warns, naming the missing provider and its
    scope, where it sees the name before normalization. An empty registry is an error. Core state
    (`packages/core/src/state`) does not change for this (§4.3).
12. Presets are resolved by `id` or by an optional `alias` on the preset definition. Built-in presets use their current
    short keys (`'default'`, `'auto'`, `'all-models'`, ...) as aliases, so existing calls keep working; an unresolved
    string throws (§4.3).
13. Script languages other than MolScript are enabled by a top-level import of their transpiler module, the same model
    as transformer registration. Nothing about them goes in the spec (§7).

Transformer ids, snapshot JSON, provider names, `PluginSpec.Action`/`Behavior`, and the `createPluginUI`/`Viewer.create`
entry points remain ([architecture §9.3](architecture.md#93-compatibility-contract)).

### 1.1 Terms

| Term                       | Meaning                                                                                                                                                                   |
| -------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Provider                   | A value a registry holds: a representation, theme, format, preset, selection query, loci label provider, markdown extension, drag-and-drop handler, action, or animation. |
| Entry                      | A `PluginRegistryEntry`: a plain record of providers grouped by registry.                                                                                                 |
| Registry key               | What a registry uses to detect conflicts: a name, a preset id, an action id, or identity only (§4.2).                                                                     |
| Catalog module             | A module that value-imports a complete provider set, such as all trajectory formats or all color themes.                                                                  |
| Default-composition module | A module that assembles the default plugin from catalogs: `@molstar/plugin/default-registry`, `@molstar/plugin/default-spec`, and `@molstar/plugin-ui/default-spec`.      |
| Base module                | Any other library module, including `PluginContext`, managers, builders, transformers, providers, and generic UI.                                                         |

Naming convention: a provider is `<Thing>Provider` (`SdfProvider`) and the entry exported next to it is the bare name
(`Sdf`). A preset entry is `<Thing>Preset` (`BallAndStickPreset`) and its provider `<Thing>PresetProvider`. The
convention applies to new built-in exports. Existing exported preset providers that already use the `<Thing>Preset` name
(`QualityAssessmentPLDDTPreset`, `AssemblySymmetryPreset`, `ViewerAutoPreset`, and others in extensions and apps) keep
their names; an entry exported next to one is named `<Thing>PresetEntry`.

## 2. Current state

Only facts that motivate the design are listed here. The coupling edges with their counts and the call-site inventories
are in the [plan](../plans/plugin-composition.md#1-inventories).

- **Implicit catalog.** `new PluginContext(spec)` loads the full built-in catalog even with an empty spec: registry and
  manager constructors preload every built-in representation, theme, format (with the facade's transformers), preset,
  selection query, and markdown extension, and `context.ts` imports every dynamic behavior. `DefaultPluginSpec` adds
  only actions, behaviors, and animations. No spec field can remove a built-in provider.
- **Preset keys.** `PresetTrajectoryHierarchy` and `PresetStructureRepresentations` are keyed by short names (`default`,
  `all-models`, `auto`, ...); the builders' registries are keyed by `id` (`preset-trajectory-default`).
  `resolveProvider` checks the static map, then the registry; the typed `applyPreset` overloads use the short keys; an
  unresolved string returns `undefined` silently. `PluginConfig.Structure.DefaultRepresentationPreset` (default
  `'preset-structure-representation-auto'`) exists, but there is no hierarchy-preset item and format visuals hard-code
  `'default'`.
- **`registry.default`.** Transformer param definitions and the representation-param helpers read
  `registry.default.name` or `.provider` on representation registries (plan §1.2). It is the first registered entry
  (`direct-volume` for volumes). The getter already returns `undefined` for an empty registry, but its type does not say
  so. `ThemeRegistry.default` (first after sorting by category and label) has no callers.
- **Lenient names.** `RepresentationRegistry.get(name)` returns `EmptyRepresentationProvider` and `ThemeRegistry.get`
  returns the empty theme for an unknown name. Core normalizes the params of every new or changed transform before
  `apply`, not only on snapshot restore: `resolveParams` in `core/src/state/state.ts` calls
  `PD.normalizeParams(definition, params, 'all')` whenever the transform version changes, and `PD.Mapped` replaces a
  name that is not an option with the default. A direct
  `.apply(StructureRepresentation3D, { type: { name: 'cartoon' } })` in a plugin without cartoon therefore renders the
  default representation. Representation providers name their default themes as strings (`defaultColorTheme: { name }`).
  `dataFormats.get(name)` already throws `unknown data format name`.
- **Format and provider knowledge in generic code.** `parseTrajectory(blob)` handles mmCIF only, `DownloadStructure`
  uses static format lists, volume and particle `visuals` apply names, and `@molstar/model/script/script` value-imports
  every script transpiler, which every representation reaches (plan §1.5, §1.6).

## 3. The spec and entry contract

`PluginSpec` and `PluginRegistryEntry` live in `@molstar/plugin/spec`. That module holds only types and the
`PluginSpec.Action`/`Behavior` helpers; `DefaultPluginSpec` moves to `@molstar/plugin/default-spec`
([architecture §3.3, §9.1](architecture.md#33-plugin-initialization-cycles)). Graphics and model types are referenced
with `import type`; plugin already depends on those packages, so no package cycle arises.

```ts
interface PluginSpec {
  registry?: readonly PluginRegistryEntry[];
  behaviors: PluginSpec.Behavior[];
  config?: [PluginConfigItem, unknown][];
  canvas3d?: PartialCanvas3DProps;
  layout?: { initial?: Partial<PluginLayoutStateProps> };
}

interface PluginRegistryEntry {
  readonly structure?: {
    readonly themes?: PluginRegistryEntry.Themes;
    readonly representations?: readonly StructureRepresentationProvider<any>[];
    readonly presets?: {
      readonly hierarchy?: readonly TrajectoryHierarchyPresetProvider<any, any>[];
      readonly representation?: readonly StructureRepresentationPresetProvider<any, any>[];
    };
    readonly selectionQueries?: readonly StructureSelectionQuery[];
  };
  readonly volume?: {
    readonly themes?: PluginRegistryEntry.Themes;
    readonly representations?: readonly VolumeRepresentationProvider<any>[];
  };
  readonly particles?: {
    readonly themes?: PluginRegistryEntry.Themes;
    readonly representations?: readonly ParticleRepresentationProvider<any>[];
  };
  readonly formats?: readonly DataFormatProvider[];
  readonly lociLabels?: readonly LociLabelProvider[];
  readonly markdownExtensions?: readonly MarkdownExtension[];
  readonly dragAndDrop?: readonly PluginDragAndDropEntry[];
  readonly actions?: readonly (StateAction | StateTransformer | PluginSpec.Action)[];
  readonly animations?: readonly PluginStateAnimation[];
}

namespace PluginRegistryEntry {
  export interface Themes {
    readonly color?: readonly ColorTheme.Provider<any, any>[];
    readonly size?: readonly SizeTheme.Provider<any, any>[];
  }
}

// state/manager/drag-and-drop.ts; MVS already uses this shape
interface PluginDragAndDropEntry {
  readonly name: string;
  readonly handle: PluginDragAndDropHandler; // unchanged: (files, plugin) => Promise<boolean> | boolean
  readonly fallback?: boolean; // §4.5
}
```

Representation provider types have no default for their first type parameter, so entry fields use `<any>`. Fields are
`readonly` so `as const` catalogs and default entries are assignable.

### 3.1 Spec field changes

| 5.x spec field               | 6.0                                                             |
| ---------------------------- | --------------------------------------------------------------- |
| `actions: [...]`             | `registry: [..., { actions: [...] }]`                           |
| `animations: [...]`          | `registry: [..., { animations: [...] }]`                        |
| `customFormats: [[name, p]]` | `registry: [..., { formats: [p] }]`, where `p.name` is the name |

The fields are removed rather than kept as aliases. TypeScript rejects them; for plain JavaScript and script-tag users,
the `PluginContext` constructor throws when any of the three keys has a value other than `undefined`, with the message
`PluginSpec.<key> was removed in 6.0; use registry entries (see the migration guide)`. A key present with an `undefined`
value, as in `customFormats: o?.customFormats`, is ignored. This is a check, not an alias.

A 5.x spec that does not spread `DefaultPluginSpec()`/`DefaultPluginUISpec()` received the full built-in catalog from
the `PluginContext` constructor. In 6.0 the same literal gets an empty plugin. Two migration forms keep 5.x behavior:

- A literal spec, with or without `actions: defaultSpec.actions`/`animations: defaultSpec.animations`, becomes
  `registry: DefaultRegistry` or `registry: defaultSpec.registry`, plus entries for its own providers.
- A spec that spreads the default spec and replaces one list replaces the matching default entry instead, for example
  `registry: [...DefaultRegistry.filter((e) => e !== DefaultAnimations), { animations: [AnimateModelIndex] }]`.

The migration tool rewrites both forms and reports any spec it cannot classify.

`DataFormatProvider` gains a required `name`, like representation and theme providers. Each provider's `name` equals its
current registration key: the built-in list keys (`mmcif`, `sdf`, `cifCore`, ...) and the extension keys (`MVSJ`,
`MVSX`, `KIN`, `g3d`, ...). These names appear in format params and `DownloadFile` options and do not change. To keep
names literal, the interface gets an id parameter and the helper becomes `const`-generic:

```ts
interface DataFormatProvider<P = any, R = any, V = any, D = any, Id extends string = string> {
  readonly name: Id;
  // ... existing fields
}
function DataFormatProvider<const T extends DataFormatProvider>(p: T): T;
```

`TrajectoryFormatProvider` and its siblings forward `Id`. Built-in providers use the helper instead of a widening
`: TrajectoryFormatProvider` annotation. `DataFormatProvider.withName(p, name)` returns `p` when `p.name === name`,
otherwise a memoized `{ ...p, name }` copy, which is a separate identity. `DataFormatProvider.Unnamed`
(`Omit<DataFormatProvider, 'name'> & { name?: string }`) types inputs that may lack a name, such as the Viewer's
`customFormats`, so TypeScript callers that pass providers without `name` still compile.

### 3.2 Rules

- **Pure data.** An entry contains provider values only: no behaviors, no config, no functions run at registration. A
  preset entry cannot suggest itself as the default; config stays in the spec.
- **No resolution.** Entries have no `requires`. Each entry lists what it directly needs (§5). What is registered is
  exactly what the spec lists.
- **Not serialized.** Entries are not part of the state tree or snapshots and are never toggled.
- **Runtime-built entries are allowed.** An entry can be built at run time, for example in a behavior's `register()`
  from plugin state, and passed to `plugin.register`. It is still plain data when registered.

### 3.3 Not covered by entries

These stay behavior-managed (or global), because they take behavior params, are global, or belong to the UI:

- Custom model, structure, and volume properties. `CustomProperty.Registry.register(provider, defaultAutoAttach)` takes
  a behavior param that `update()` changes. Loci label providers that close over behavior params (for example the
  accessible-surface-area tooltip) also stay in the behavior.
- Markdown ref and URI resolvers, and `State.registerRefResolver`.
- MolScript symbols in the module-global `DefaultQueryRuntimeTable` and script languages (§7). Like transformers, these
  tables are shared by every plugin on the page. Separate fix: `addCustomProp`/`removeCustomProp` and
  `addSymbol`/`removeSymbol` count registrations, so one plugin's behavior unregistering does not remove symbols another
  plugin still uses.
- UI maps (§10).

## 4. Registration semantics

### 4.1 Lifecycle

`init()` creates `interactivity`, `lociLabels`, and `builders.structure` as today (they are `void 0` before), then
registers `spec.registry`, then initializes behaviors. Within one entry the order is fixed: for each scope (structure,
volume, particles) themes then representations; then formats; then structure presets and selection queries; then loci
labels, markdown extensions, drag-and-drop handlers, actions, and animations. Entries register in spec order. Behaviors
can rely on registered providers.

```ts
class PluginContext {
  register(entry: PluginRegistryEntry | readonly PluginRegistryEntry[]): () => void;
}
```

- **Atomic.** `register` checks every provider in its input against existing registrations and against each other before
  changing any registry. On a conflict it throws one error listing every conflict and leaves the plugin unchanged.
  `init()` registers the spec's entries through the same path; a conflict rejects `init()`.
- **Undo.** The returned function decrements exactly the registrations it made. It is idempotent.
- **Timing.** Calling `register` before `init()` has created the managers and builders throws
  `PluginContext.register called before init()`. Behaviors run after this point.
- **Ownership.** Reference counting lives in each registry and manager, not only in `register`, so imperative
  `add`/`remove` calls take part in it. A behavior that calls `registry.remove(p)` does not remove a provider the spec
  also registered.

### 4.2 Keys, identity, and reference counting

Registering the same object again increments a count. Unregistering decrements it and removes the provider at zero.
Registering a different object under an existing key is an error. Removing an unknown provider is a no-op. Today's
duplicate and unknown-removal behavior of each registry is listed in the
[plan](../plans/plugin-composition.md#18-duplicate-handling-today).

| Entry field                                    | Registry or manager                                       | Key                            |
| ---------------------------------------------- | --------------------------------------------------------- | ------------------------------ |
| `<scope>.themes.color/size`                    | `representation.<scope>.themes.{color,size}ThemeRegistry` | `name`                         |
| `<scope>.representations`                      | `representation.<scope>.registry`                         | `name`                         |
| `formats`                                      | `dataFormats`                                             | `name`                         |
| `structure.presets.hierarchy`/`representation` | `builders.structure.hierarchy`/`.representation`          | `id`, and `alias` when present |
| `structure.selectionQueries`                   | `query.structure.registry`                                | Identity only                  |
| `lociLabels`                                   | `managers.lociLabels`                                     | Identity only                  |
| `markdownExtensions`                           | `managers.markdownExtensions`                             | `name`                         |
| `dragAndDrop`                                  | `managers.dragAndDrop`                                    | `name`                         |
| `actions`                                      | `state.data.actions` (core `StateActionManager`)          | Id of the cached action (§6.4) |
| `animations`                                   | `managers.animation`                                      | `name`                         |

- Actions register on the data state only; no code registers actions on the behavior state. A `PluginSpec.Action`
  wrapper is compared by the action it carries (`PluginSpec.Action(X)` allocates a new wrapper). `init()` reads only
  `a.action`, so `customControl` and `autoUpdate` have no reader; they stay on the type, unused.
- Selection queries and loci labels have no name, so they have no conflict check; identity counting deduplicates them
  once each query is created once (plan step 1).
- A drag-and-drop provider's identity is its `handle` function: entries with the same `name` and `handle` are the same
  provider, whatever the wrapper. `addHandler(name, fn, options?: { fallback?: boolean })` registers
  `{ name, handle: fn, ...options }` with counting; a different function under an existing name is a conflict.
  `removeHandler(name)` decrements the provider under `name` and is a no-op for an unknown name.
- `PluginAnimationManager` gains `unregister`, with counting. Removing the current animation selects the first remaining
  one, or none, and resets the cached params.
- `clear()` (representation and theme registries, `lociLabels.clearProviders()`) drops all providers and counts; later
  decrements of a cleared provider are no-ops. Apps should prefer registering a narrower set over clearing.

### 4.3 Defaults and lookup

**Representation defaults.** `RepresentationRegistry.default` returns the first registered `{ name, provider }` entry,
or `undefined` when the registry is empty; its type now includes `undefined`. There is no config item and no
"applicable" test; `default` takes no data. `DefaultRegistry` registers each representation catalog in today's order, so
the defaults stay `cartoon`, `direct-volume`, and particle `spacefill`. `ThemeRegistry` no longer takes a built-in map;
its `default` keeps the sorted order. Helpers use the requested name, and `default` when none is given.

**Unregistered names are lenient.** Graphics registries keep the lenient `get(name)`. `ThemeRegistry.has(provider)`
compares `provider.name` today; both registry kinds get `has(nameOrProvider: string | Provider)`. There is no checked
`getOrThrow`, no transformer hook, and no core-state option:

- **Normalization is unchanged.** `PD.Mapped` substitutes the registry default for a type, color, or size name that is
  not registered, for direct applies and snapshot restore alike.
- **Warnings where the name is seen.** Two places receive names before normalization and warn through `plugin.log.warn`,
  naming the missing provider and its scope (for example
  `Structure representation 'cartoon' is not registered in this plugin; the registry default is used`): snapshot restore
  (`reportUnregisteredNames`, §6.2), and the representation builders and param helpers (`buildRepresentation`,
  `createStructureRepresentationParams`, `createStructure{Color,Size}ThemeParams`, and the volume equivalents). The
  helpers check each requested name and the representation's default theme names with `has` in the scope, warn, and
  otherwise build params as today, so normalization substitutes.
- **Direct applies are silent.** A transformer applied with raw params, without the helpers, is substituted without a
  warning, as in 5.x.
- **Empty registries are errors.** When the scope's representation, color theme, or size theme registry is empty, there
  is nothing to substitute. `StructureRepresentation3D`, `VolumeRepresentation3D`, and `ParticlesRepresentation3D` throw
  in `apply`/`update` (for example `No structure representations are registered in this plugin`) instead of rendering
  `EmptyRepresentationProvider` or the empty theme, and the builders and helper calls given data (`structure` or
  `volume`) throw the same error.
- **Param definitions never assert.** Transformer param definitions and helper calls without data must not throw,
  because `StructureFocusRepresentation` (whose `customDefault` params call the helper with no structure) and
  `cli/state-docs` build them outside any load. With an empty registry, the `PD.Mapped` type defaults fall back to
  `PD.MappedStatic('', { '': PD.EmptyGroup() })` and the `registry.get(registry.default.name)` theme lookups are
  skipped, so the theme params seeded from them become empty mapped params too. A data-less helper call returns such
  params; applying them fails with the transformer error above.
- **Provider objects** passed to a builder or helper are used directly; `getName` keeps throwing for an unregistered
  provider, because an imported but unlisted provider is a composition bug, not a degraded snapshot.

**Preset lookup.** Preset builders resolve strings through their registry only; the static-map lookup and the fixed
`defaultProvider` are removed. `PresetProvider` gains an optional `alias`:

```ts
interface PresetProvider<O extends StateObject = StateObject, P = any, S = {}, Id extends string = string> {
  id: Id;
  alias?: string; // e.g. 'default', 'auto', 'all-models'; indexed by the registry under the same conflict rule as id
  // ... existing fields
}
```

Built-in presets set `alias` to their current static-map key (plan §1.3), so short-key calls such as
`applyPreset(trajectory, 'default')`, `representationPreset: 'auto'`, and the Viewer's `loadTrajectory({ preset })` keep
working without a static map. A string resolves by `id`, then by `alias`. An `alias` that equals another registered
preset's `id` or `alias` is a conflict (§4.2). The `TrajectoryHierarchyPresetProvider` and
`StructureRepresentationPresetProvider` helpers become `const`-generic so ids and aliases stay literal, and the preset
catalogs derive the id and alias unions (`BuiltInTrajectoryHierarchyPresetId`/`Alias`,
`BuiltInStructureRepresentationPresetId`/`Alias`). Both builders have three public `applyPreset` overloads, typed
through `import type`: a built-in id or alias (typed params), a provider object, and `string` (any params).
`TrajectoryHierarchyBuilder` gains the `string` overload, so calls that pass a config value compile. An unresolved
string throws `Preset '<x>' is not registered in this plugin` instead of returning `undefined`. `getPresetSelect`
defaults to the configured preset id when registered, otherwise the first option. Config items hold ids.

**Format lookup.** `DataFormatRegistry` is keyed by `provider.name` and keeps priority-then-registration-order
resolution for `auto()`. The two-argument `add(name, provider)` is kept; with a name other than `provider.name` it
registers `DataFormatProvider.withName(provider, name)` (§3.1) and logs a deprecation warning. Registering an original
and a renamed copy under the same name is a conflict. `remove` takes a provider or a name, `has(name)` is added, and
`get(name)` keeps throwing, saying the format is not registered in this plugin, so `if (!get(...))` checks are dead and
become `has(name)` checks (plan §1.4).

### 4.4 Animations

- `initAnimations` no longer registers `AnimateStateSnapshotTransition` when the spec lists no animations (today it
  does, only so that `current`, typed `this._current!`, is defined). With none registered, `current` is `undefined`,
  typed so, and readers guard it as `isAnimatingStateTransition` already does.
- `AnimateStateSnapshotTransition` is base functionality: `state.ts`, the markdown-extension manager, and the UI
  snapshot controls keep their direct import of this small module and play it.
- `play(animation)` keeps registering an unregistered animation on demand, outside reference counting. A later counted
  `register` of the same object adopts it.
- `DefaultAnimations` lists it in its current position (sixth), so the default animation stays `AnimateModelIndex`.
- `setSnapshot` with an animation name that is not registered skips `snapshot.current` and warns, instead of writing
  `paramValues` into the previous animation.

### 4.5 Drag and drop

`handle(files)` tries non-fallback handlers in reverse registration order, then built-in session handling (`.molx` and
`.molj` through `PluginCommands.State.Snapshots.OpenFile`), then fallback handlers in reverse registration order.
Session handling stays in the manager because it needs no catalog, so slim plugins keep it. The open-anything handler,
which imports `OpenFiles`, becomes the `DefaultDragAndDrop` entry with `fallback: true`, so it runs last whatever the
entry order.

## 5. What goes into an entry

A module that defines a provider exports a ready-made entry next to it. Inline records are equally valid:
`registry: [Sdf, BallAndStick, { lociLabels: [MyLabelProvider] }]`.

### 5.1 Import what you run

A provider statically imports the providers it runs, and its entry includes them. Bundling and registration then agree,
and a mistyped dependency is an import error. Registration is still required after importing: snapshots store names, and
the UI lists only registered providers.

- **Representations.** A representation entry includes the representation and the providers of its default color and
  size themes, in its scope: `BallAndStick` is
  `{ structure: { representations: [BallAndStickRepresentationProvider], themes: { color: [ElementSymbolColorThemeProvider], size: [PhysicalSizeThemeProvider] } } }`.
  These entries live in plugin-layer modules (for example `@molstar/plugin/registry/structure/ball-and-stick`), because
  graphics cannot reference `PluginRegistryEntry`. An inline record that lists a bare representation must list its
  default themes too. In development mode, the end of `init()` and each `plugin.register` warn for every registered
  representation whose default theme names are not registered in its scope; this is a check, not resolution.
- **Presets.** A preset imports the representation and theme providers it builds, passes provider objects to the
  builder, and lists their entries in its own. `BallAndStickPreset` includes `BallAndStick`. Code that runs a specific
  preset directly imports that preset's defining module, not the catalog.
- **The auto preset.** `auto` is an ordinary preset module, not a catalog. It imports the presets it chooses between
  and, through them, their representations and themes, and its entry includes theirs. Code that delegates to it (the
  extension presets and the Viewer's `ViewerAutoPreset`) imports `auto`'s defining module and accepts that this pulls in
  everything `auto` can build. Because it is not a catalog, the import-graph rule allows value imports of it.
- **Mixed props.** `createStructureRepresentationParams`, `buildRepresentation`, and the volume helpers accept a name or
  a provider in each field (`type`, `color`, `size`) and resolve each field independently, so a preset can pass
  `{ type: CartoonRepresentationProvider, color: params.theme.globalName }`. Param types are derived per field: provider
  → its params; built-in name → `BuiltInParams`; other string → `{}`. Names follow §4.3. The current "any string field
  selects the by-name path" dispatch is removed.
- **Formats.** A format entry includes the format provider, the actions that should appear for it, and, for volume,
  particle, and shape formats, the representation entries and scope themes its `visuals` step applies, passed as
  provider objects. CCP4's entry includes `Isosurface` (isosurface with its `uniform` themes). Built-in format entries
  list only actions in today's default action list (per-format lists in plan step 3), so `DefaultRegistry` registers no
  action that 5.x did not.
- **Delegation to toggleable behaviors is policy.** A preset that delegates to presets owned by a behavior that can be
  turned off (the Viewer's `ViewerAutoPreset` calling the model-archive QA and SB-NCBR presets) looks them up by id
  through the registry and skips them when absent. It does not import and include them.

### 5.2 Policy goes through config

Which preset a format or action applies after loading is policy. It goes through config, not an import; otherwise the
SDF format would import the default preset and with it cartoon and the surfaces. A fixed representation that a volume,
particle, or shape format's `visuals` step builds directly is part of that format and is imported (§5.1).

- New `PluginConfig.Structure.DefaultHierarchyPreset` (default `'preset-trajectory-default'`) replaces hard-coded
  `'default'` wherever a hierarchy preset is applied after loading.
- `PluginConfig.Structure.DefaultRepresentationPreset` (default `'preset-structure-representation-auto'`) replaces
  hard-coded `'auto'` and `PresetStructureRepresentations.auto.id`, including in the hierarchy presets. Their
  `representationPreset` param is typed with the built-in id and alias unions through `import type` and has no built-in
  default; when it is absent the preset reads the config item. Hierarchy presets then no longer value-import the
  representation-preset catalog.
- When the configured preset is not registered, loading fails with the §4.3 preset error; there is no silent fallback. A
  slim app sets `DefaultHierarchyPreset` when it does not register the default hierarchy preset, and
  `DefaultRepresentationPreset` when it does not register the auto preset.
- Other generic code that names a specific provider either imports it and its entry includes it, or reads the choice
  from config and the registries. Code that applies a provider only a toggleable behavior registers checks `has` and
  skips that part when the provider is absent.

### 5.3 Catalogs and defaults

Catalog modules assemble complete sets: the format lists, the representation and theme maps, the preset maps, the query
catalog, `PluginBehaviors`, `StateActions`, and the transformer catalog (§6.1). Only catalog modules,
default-composition modules, and apps value-import catalogs; base modules may use `import type` (§8).

`DefaultRegistry` lives in `@molstar/plugin/default-registry`. It is the concatenation of exported named entries, so
apps compose subsets by identity:

```ts
export const DefaultRegistry: readonly PluginRegistryEntry[] = [
  DefaultThemes, // every built-in theme in all three scopes, as today, plus external-structure/-volume
  DefaultStructureRepresentations,
  DefaultVolumeRepresentations,
  DefaultParticleRepresentations,
  DefaultActions, // today's full action list in today's order; before DefaultFormats (below)
  DefaultFormats, // today's order: volume, topology, coordinates, shape, particles, trajectory
  DefaultPresets,
  DefaultSelectionQueries, // today's catalog and residue queries, in today's order (§7)
  DefaultMarkdownExtensions,
  DefaultDragAndDrop,
  DefaultAnimations,
];
```

The module, these eleven export names, and their order are part of the public contract; migration examples, apps, and
the migration tool refer to them. `DefaultActions` comes before `DefaultFormats` because format entries carry actions
and actions are listed per type in registration order; this way the format entries' actions are already registered
(counted by id) and the 5.x action order is kept. `DefaultPresets` holds today's hierarchy and representation presets;
the new `BallAndStickPreset` (§12) is not added to it. The module also imports the transformer catalog (§6.1) and
`@molstar/model/script/transpilers/all` (§7) for their side effects. `DefaultPluginSpec()` in
`@molstar/plugin/default-spec` returns `{ registry: DefaultRegistry, behaviors }`; the default custom-property behaviors
stay behaviors, with unchanged runtime registration.

A theme is available only in the scopes it is listed for.

### 5.4 External color themes

`external-structure` and `external-volume` were deferred in the prototype because of two dependency problems:

1. **Package direction.** In 5.x they were part of the graphics `ColorTheme.BuiltIn` map, but they depend on the plugin:
   their `PD.ValueRef` params list structures or volumes from plugin state (`PluginStateObject`, `PluginContext`), and
   external-structure runs the plugin's `backbone` selection query. Graphics cannot import plugin, so they moved to
   `@molstar/plugin/themes/*` and fell out of `ColorTheme.createRegistry()`.
2. **Registration through `PluginContext`.** The only remaining place to add them to all six theme registries was a
   special helper in `PluginContext`. That made the context value-import an optional theme, and through it
   `state/helpers/structure-selection-query.ts`, which imports the `StateTransforms` facade, whose `transforms/model`
   imports the query module back: the evaluation-order cycle the facade's lazy getters work around (issue #1791).

The registry design removes both:

- **Owner.** The providers stay in `@molstar/plugin/themes/external-structure` and `external-volume`. Graphics catalogs
  never mention them, so there is no upward package edge.
- **Registration is data.** An `ExternalColorThemes` entry next to them lists both providers under
  `structure.themes.color`, `volume.themes.color`, and `particles.themes.color`, matching 5.x, where every scope's
  registry contained them. `DefaultThemes` includes it; other apps list it when they want the themes. `PluginContext`
  and other base modules never import them.
- **No cycle.** external-structure imports `backbone` from `state/queries/structure/structure.ts` (§7), which imports
  only core and model script code. With `applyBuiltInSelection` removed and the facade gone, neither theme reaches a
  transform module, `PluginContext`, or a catalog; `PluginContext` remains a type-only import in the `PD.ValueRef`
  getters. The import-graph check asserts this (plan step 4).
- **Snapshots.** A snapshot that uses an external theme in a plugin that does not list the entry falls back to the
  registry default with a warning (§4.3).

This replaces the temporary ledger entry and closes the checklist item.

## 6. Transformers and snapshots

### 6.1 Global registration

Transformers stay in the global `StateTransformer` registry: importing a defining module registers its transformers by
id. A transformer is available to snapshots if any imported module defines it. After the split, several built-in
transformers are imported by nothing but their own module (for example `ImportString`, `ImportJson`, and `ParseJson`),
so a transformer catalog module, `state/transforms/catalog.ts`, value-imports every built-in transform and format module
and exports nothing. `@molstar/plugin/default-registry` imports it, so the default plugin and the Viewer can restore any
5.x snapshot. Slim apps import leaf modules instead; tools that need every transformer import the catalog.

### 6.2 Snapshot validation

Today `PluginState.setSnapshot` stops the animation and applies structure-component options and the behavior tree before
a missing transformer throws in the data tree, leaving the plugin half-updated, and
`PluginStateSnapshotManager.setStateSnapshot` clears the snapshot list first. Validation runs in two phases, because
provider names can only be checked after the snapshot's behaviors have registered their providers. Both functions read
ids and names only; transition frames carry their own data trees.

- **`PluginState.validateSnapshotTransformers(snapshot)` (fatal).** It collects transformer ids from the behavior tree,
  the data tree, and every `transition.frames[i].data`, and throws one error listing the ids missing from the global
  registry. `setSnapshot` calls it first, before `animation.stop()`; `setStateSnapshot` calls it for every entry before
  `clear()`. A failing snapshot leaves the plugin unchanged.
- **`PluginState.reportUnregisteredNames(snapshot)` (warning).** `setSnapshot` calls it after the behavior tree is
  applied and before the data tree is applied. It collects the type, color, and size names from
  `StructureRepresentation3D`, `VolumeRepresentation3D`, and `ParticlesRepresentation3D` params in the data tree and
  transition frames, checks them with `has(name)` in the matching scope, and warns for each missing name, naming the
  provider and scope and saying the registry default replaces it. It does not name a specific substitute, which depends
  on the data (`getApplicableTypes`). The data-tree update then normalizes as today. The check cannot be fatal:
  behaviors enabled by the snapshot register providers when the behavior tree is applied, and 5.x snapshots naming
  extension themes the Viewer does not ship degrade rather than fail. An empty registry in a scope the snapshot uses
  fails that cell (§4.3).
- **Risk: custom properties.** A snapshot naming a custom property whose behavior is absent still fails during data-tree
  restore (`customModelProperties.get` throws), and behaviors the snapshot enables make these names impossible to
  pre-check. The restore error must name the property.

### 6.3 Module layout

Remove the `StateTransforms` facade and its lazy getters, and split the seven transform modules by functionality (target
tree in the [plan](../plans/plugin-composition.md#step-1-module-splits)):

- Generic operations live under `state/transforms/<area>/`, with subgroups where an area is large.
- Format-specific parse and conversion transformers move next to their format provider and entry, one module per format
  under `state/formats/<family>/` with a `catalog.ts` per family; for example, `formats/trajectory/sdf.ts` holds
  `TrajectoryFromSDF`, `SdfProvider`, and the `Sdf` entry. Transformers shared by several formats live in a shared
  module.
- Non-transformer helpers move out of large modules; format category constants move to small modules so importing a
  constant does not load a format family.
- Every import of a module must be justified by what its users need. Exact file names are settled during the split.

Transformer ids, names, and the `DeflateData` id typo stay unchanged. Identity checks such as
`transform.transformer === X` keep working because each transformer is still created once.

### 6.4 Actions

`transformer.toAction()` creates a new action with a new UUID on each call today, so `remove(transformer)` never finds
it and adding a transformer twice adds duplicates. In 6.0 it returns a cached action, so a transformer and its action
share one identity and id. `StateActionManager` counts registrations by action id; `remove` takes effect at zero.
Extensions that add an action in `register()` and remove it in `unregister()` no longer remove an action the spec also
lists.

## 7. Module-split rules

A module that defines a provider does not also hold a catalog or a registry class with a preload, and base modules reach
optional functionality only through ids, config, `has`, and `import type`.

- **Presets.** The preset modules split into a types-and-helpers module (`PresetProvider` interfaces, `CommonParams`,
  `reprBuilder`, `presetStaticComponent`, `updateFocusRepr`) and one module per preset (or small group) exporting the
  provider and its entry. `PresetStructureRepresentations` and `PresetTrajectoryHierarchy` become catalogs. The builders
  import types only. Layout: `state/builder/structure/representation-presets/` (`types.ts`, one module per
  representation preset, `catalog.ts`) and `state/builder/structure/hierarchy-presets/` (`types.ts`, one module per
  hierarchy preset, `crystal-symmetry.ts` shared by `unitcell` and `supercell`, `catalog.ts`).
- **Focus representation.** Presets and `StructureComponentManager` reach the focus representation behavior only through
  a small id module, `plugin.state.hasBehavior`/`updateBehavior`, and `import type` for its params. They do not name
  representations or themes for it, and they update the behavior only when it is present instead of inserting it (plan
  step 3).
- **Behaviors.** `BuiltInPluginBehaviors` (static command wiring) moves to its own module. `behavior.ts` loses its
  `export *` and becomes only the `PluginBehaviors` catalog. Base modules import `PluginBehavior` from
  `@molstar/plugin/behavior/behavior` ([architecture §4.4](architecture.md#44-no-barrel-files)).
- **Actions.** `StateActions` becomes a catalog module; library code imports actions from `state/actions/*`.
- **Selection queries.** Grouped by functionality like the transforms, using the existing `StructureSelectionCategory`
  values as the subgroups. Each group module exports its queries and an entry with them under
  `structure.selectionQueries`:

  ```text
  state/queries/structure/
    query.ts       # StructureSelectionQuery type and constructor, StructureSelectionCategory
    registry.ts    # StructureSelectionQueryRegistry (separate from query.ts so the preload does not form an import cycle)
    common.ts      # entity and residue tests shared by the groups
    basic.ts       # all, current
    type.ts        # polymer, protein, nucleic, water, ion, lipid, branched, ligand, coarse, and their hidden
                   # internal variants (ligandPlusConnected, branchedPlusConnected, ...)
    structure.ts   # trace, backbone, sidechain, sidechainWithTrace, helix, beta
    bond.ts        # disulfideBridges, nosBridges
    residue.ts     # nonStandardPolymer, ring, aromaticRing, and the amino-acid and nucleic-base residue queries
    manipulate.ts  # surroundings, surroundingLigands, surroundingAtoms, complement, covalentlyBonded, ...
    dynamic.ts     # ResidueQuery, ElementSymbolQuery, EntityDescriptionQuery and the get*Queries helpers the UI
                   # builds per structure; never registered
    catalog.ts     # the StructureSelectionQueries map
  ```

  Registration only controls what the selection UI offers. Code that runs a query (presets, static components, the
  external-structure theme, the component manager) imports the query value from its group module and does not need it
  registered. The residue queries are created once at module scope, so identity counting deduplicates them.
  `DefaultSelectionQueries` lists the individual queries in today's order rather than concatenating the group entries,
  so the UI order is unchanged. The unused `applyBuiltInSelection` is removed; with it and `structure-component.ts`
  importing group modules, no query module imports a transform, which removes the facade cycle. The component manager
  picks its default query by identity, not position (plan step 3).

- **Component manager options.** The interactions code leaves `StructureComponentManager`, but the `interactions` option
  slot stays so it round-trips through snapshot JSON; the Interactions behavior registers the handler that applies it,
  and without the behavior the option is preserved but ignored.
- **Markdown extensions.** `BuiltInMarkdownExtension` moves out of the manager module into catalog entries. The `query`
  extension is its own entry and evaluates scripts in whatever languages the app has enabled.
- **Script languages.** The core `@molstar/model/script/script` module (the `Script` type, `is`, `areEqual`, `Info`, and
  `toExpression`/`toLoci`/`toQuery`/`getStructureSelection`) handles `mol-script` directly and reaches other languages
  only through a module-global transpiler table, like `DefaultQueryRuntimeTable`. The table lives in
  `@molstar/model/script/transpile`, whose `parse(lang, str)` reads it; `transpile.ts` no longer imports
  `transpilers/all.ts`, which stops exporting `_transpiler`. Importing `@molstar/model/script/transpilers/<lang>`
  (`pymol`, `vmd`, `jmol`) at the top of the app or spec file registers that language for every plugin on the page;
  `@molstar/model/script/transpilers/all` imports all three. Code that calls `parse` directly imports the languages it
  parses the same way. Nothing goes in the spec, and languages that are not imported are not bundled.
  `@molstar/plugin/default-registry` and the Viewer import `transpilers/all`. A script in a language that is not enabled
  throws `Script language '<x>' is not available in this build` (in a snapshot, that cell fails), and the script-param
  UI offers only enabled languages.
- **Graphics registries.** `StructureRepresentationRegistry`, `VolumeRepresentationRegistry`,
  `ParticleRepresentationRegistry`, and `ThemeRegistry` keep their classes without preloads. The built-in maps move to
  catalog modules such as `@molstar/graphics/repr/structure/catalog` and `@molstar/graphics/theme/color/catalog`.

## 8. Built-in name types

The name types keep their names. Each is derived in its catalog module from the value list, and other modules import it
with `import type`. Format catalogs become provider arrays, which removes the key/name duplication:

```ts
// state/formats/trajectory/catalog.ts
export const BuiltInTrajectoryFormats = [MmcifProvider, CifCoreProvider, PdbProvider /* ... */] as const;
export type BuiltInTrajectoryFormat = (typeof BuiltInTrajectoryFormats)[number]['name'];

// state/builder/structure.ts
import type { BuiltInTrajectoryFormat } from '@molstar/plugin/state/formats/trajectory/catalog';
```

- **Format name types:** `BuiltInTrajectoryFormat`, `BuiltInCoordinatesFormat`, `BuiltInTopologyFormat`,
  `BuiltInParticlesFormat`, and the misspelled `BuildInVolumeFormat` and `BuildInShapeFormat`, which stay as deprecated
  aliases of the added `BuiltInVolumeFormat` and `BuiltInShapeFormat`. The types move from `state/formats/<family>` to
  `state/formats/<family>/catalog`.
- **Preset name types:** `BuiltInTrajectoryHierarchyPresetId`/`Alias` and
  `BuiltInStructureRepresentationPresetId`/`Alias` (§4.3).
- **Namespace types keep their paths.** `ColorTheme.BuiltIn`/`BuiltInParams`, `SizeTheme.BuiltIn`/`BuiltInParams`,
  `StructureRepresentationRegistry.BuiltIn`/`BuiltInParams`, and the volume and particle equivalents stay namespace
  types; the owning module type-imports its catalog in the same package (for example `theme/color.ts` declares
  `type BuiltIn = keyof typeof BuiltInColorThemes` from `./color/catalog.js`). This type-only cycle inside one package
  is allowed by [architecture §3.3](architecture.md#33-plugin-initialization-cycles), so type use sites stay unchanged.
  The namespace values (`ColorTheme.BuiltIn`, `SizeTheme.BuiltIn`, `StructureRepresentationRegistry.BuiltIn`, ...) are
  removed; value users import providers directly.
- **Catalog checks.** Representation and theme catalogs stay name-keyed maps because `BuiltInParams<T>` indexes by key.
  A compile-time check replaces today's constructor check:
  `namedCatalog<T extends { [K in keyof T]: { readonly name: K } }>(t: T): T`.
- **Boundary.** An import-graph check rejects value imports of catalog modules outside catalog modules,
  default-composition modules, and apps. Type imports must not point to a higher package; model cannot type-import a
  plugin catalog. Declarations reference catalog declarations, so type-checking loads them; bundles do not.
- **Classification.** Catalog names are not uniform, so the check reads a checked-in manifest of catalog and
  default-composition modules. Apps are everything under `apps/*`, `examples/*`, `cli/*`, `servers/*`, and `smoke/*`;
  any other module is a base module. `MVSRuntimeRegistry` (§11) lives in a catalog module of `packages/mvs/runtime`
  listed in the manifest, which is how MVS consumes catalogs as a library package.
- **Typed names do not prove registration.** Lookup reports a provider missing from the plugin (§4.3). If fast types are
  adopted in v7, the derived types become explicit unions checked against the catalogs.

## 9. Decoupling `PluginContext`

Base modules must not value-import catalogs or optional functionality:

- Registries and managers have no built-in preloads (§7); `DataFormatRegistry` no longer imports the format lists;
  `context.ts` imports `BuiltInPluginBehaviors` from its own module; the drag-and-drop fallback and the snapshot
  transition are no longer implicit (§4.4, §4.5).
- `parseTrajectory(blob)` moves into the mmCIF entry's module or becomes format-neutral. The Cube format's structure
  path imports the transforms it uses.
- `DownloadStructure` belongs in `DefaultActions`, not in a format's entry. It derives its sources and URL formats from
  the registered formats (plan step 3) and its preset values from config (§5.2).
- `@molstar/plugin/spec` keeps only `PluginSpec`, `PluginRegistryEntry`, and the `PluginSpec.Action`/`Behavior` helpers.
- Library helpers that construct plugins (`PluginViewModel`, `PluginUIViewModel`, and their hooks in
  `@molstar/plugin-extension`) require an explicit spec; the hooks' `spec: (defaultSpec) => spec` form is removed. The
  Viewer global exposes defaulting wrappers at the app boundary.
- Base entry points (`@molstar/plugin/context`, `@molstar/plugin/spec`, `@molstar/plugin-ui`, `@molstar/plugin-ui/spec`)
  must not value-import default specs or catalogs; the import-graph check enforces this.

## 10. UI

- `DefaultPluginUISpec` moves from `@molstar/plugin-ui/spec` to `@molstar/plugin-ui/default-spec`
  ([architecture §9.1](architecture.md#91-import-map)); `@molstar/plugin-ui/spec` keeps only the `PluginUISpec` type.
  `createPluginUI` requires a spec, and `index.ts` stops importing the default.
- The base `Plugin`/`ControlsWrapper` no longer falls back to `DefaultStructureTools`. The full tools component moves to
  `@molstar/plugin-ui/default-spec` and is set as `DefaultPluginUISpec().components.structureTools`. Without
  `structureTools`, the base layout renders a minimal list (structure source, components, measurements) that reaches no
  catalog. A spec that spreads `DefaultPluginUISpec()` but replaces `components` must spread the default components to
  keep the full tools.
- Quick styles resolve presets by id through the registry (`DefaultRepresentationPreset`, then fixed ids) and hide
  buttons whose preset is not registered, instead of importing the preset map.
- The volume controls take the representation type from the volume registry or config and import transformers from leaf
  modules; volume streaming follows §5.2.
- Follow-up: a UI registry on `PluginUISpec` for structure tools, import controls, and generic representation controls.
  These currently live on core `PluginContext` as untyped maps (`customStructureControls`, `customImportControls`,
  `genericRepresentationControls`) holding React components in a non-UI package; they move to `PluginUIContext`, and
  `customParamEditors` (already a typed `PluginUISpec` field) joins the registry. Until then they stay
  behavior-registered.

## 11. Extensions, MVS, and apps

- **Extensions** register providers imperatively in behavior `register()`/`unregister()` today and keep their behaviors.
  A behavior may replace hand-written registration with `plugin.register(entry)`, built from plugin state if needed, and
  call the undo in `unregister()`. An extension can also export a plain entry for providers that never need toggling.
  Custom properties and resolvers stay imperative (§3.3). Extensions are base modules for the import-graph check.
- **MVS.** Only the provider part of the MVS `Registrables` record (`mvs/runtime/src/behavior.ts`) becomes an entry,
  which `MolViewSpec.register()` builds and registers with an undo; custom properties and resolvers stay imperative. MVS
  documents select representations and themes by name, so MVS is a deliberate catalog consumer: the runtime exports
  `MVSRuntimeRegistry`, listing every representation and color/size theme provider MVS can name. An MVS-only app lists
  it plus the format entries and markdown extensions MVS uses.
- **Viewer.** The registry is `[...DefaultRegistry, ViewerEntry, customFormatsEntry]`, with the custom-formats entry
  last so built-in format order and `auto()` tie-breaking are unchanged. `customFormats` keeps its `[name, provider]`
  shape with `DataFormatProvider.Unnamed` providers, each registered through `withName` (§3.1). A tuple whose name
  matches a built-in format overrides it as in 5.x without putting two providers under one name; mesoscale-explorer does
  the same (plan §3.1).
- **Viewer API values stay.** `Viewer.create`, its options, and the Viewer's loading methods keep their accepted values.
  `loadTrajectory({ preset: 'all-models' })` keeps working through the built-in preset aliases (§4.3);
  `LoadTrajectoryParams.preset` is typed with the hierarchy id and alias unions.
- **Viewer global.** `molstar.lib.plugin` keeps `StateTransforms` with its 5.x keys and members as an app-level object
  literal assembled at the final bundle boundary ([architecture §4.4](architecture.md#44-no-barrel-files)), and adds
  `DefaultRegistry`. Its other library values follow the 6.0 library API (plan §3.1). This carve-out is recorded in
  [architecture §9.3](architecture.md#93-compatibility-contract); `Viewer.create` and its options are unchanged.
- **Other apps, examples, and tools** replace post-`init()` clearing with an explicit registry, add `DefaultRegistry`
  (or a narrower list) to literal specs, import default specs from `default-spec`, replace default entries instead of
  overriding removed spec fields, and use preset ids. Tools that enumerate the built-in set import the transformer
  catalog and build their context from `DefaultPluginSpec()`. Per-app details are in the
  [plan](../plans/plugin-composition.md#3-per-consumer-migration).

## 12. Acceptance

```ts
import type { PluginSpec } from '@molstar/plugin/spec';
import { PluginConfig } from '@molstar/plugin/config';
import { Sdf } from '@molstar/plugin/state/formats/trajectory/sdf';
import { DefaultHierarchyPreset } from '@molstar/plugin/state/builder/structure/hierarchy-presets/default';
import { BallAndStickPreset } from '@molstar/plugin/state/builder/structure/representation-presets/ball-and-stick';
import { createPluginUI } from '@molstar/plugin-ui';
import { renderReact18 } from '@molstar/plugin-ui/react18';

const spec: PluginSpec = {
  registry: [Sdf, DefaultHierarchyPreset, BallAndStickPreset],
  behaviors: [/* highlight, select, camera */],
  config: [[PluginConfig.Structure.DefaultRepresentationPreset, 'preset-structure-representation-ball-and-stick']],
};
const plugin = await createPluginUI({ target: document.getElementById('app')!, spec, render: renderReact18 });
const data = await plugin.builders.data.download({ url: '/ligand.sdf' }, { state: { isGhost: true } });
const trajectory = await plugin.builders.structure.parseTrajectory(data, 'sdf');
await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'preset-trajectory-default');
```

Module paths are illustrative. `BallAndStickPreset` (`preset-structure-representation-ball-and-stick`) is a new preset
written for this target. The example relies on the base layout's minimal structure tools (§10) and enables no script
language.

The example must render an SDF ligand while the import graph and bundle exclude the modules below. Verify with the
import-graph check (including the UI modules and the modules split in §7), esbuild metafile output for both a split and
a single-file build, and a rendering smoke test. The check names the excluded modules concretely; the acceptance step
pins the final list, starting from:

| Exclusion              | Modules                                                                                          |
| ---------------------- | ------------------------------------------------------------------------------------------------ |
| mmCIF parser           | `@molstar/model/formats/structure/mmcif`, `@molstar/plugin/state/formats/trajectory/mmcif`       |
| CCP4 parser            | `@molstar/io/reader/ccp4/*`, `@molstar/plugin/state/formats/volume/ccp4`                         |
| Cartoon                | `@molstar/graphics/repr/structure/representation/cartoon`                                        |
| Volume representations | `@molstar/graphics/repr/volume/{direct-volume,isosurface,slice,dot,segment}`                     |
| Other presets          | Every preset module except the default hierarchy and ball-and-stick presets; the preset catalogs |
| Query catalog          | The selection-query catalog module (`StructureSelectionQueries`)                                 |
| Script transpilers     | `@molstar/model/script/transpilers/{all,pymol,vmd,jmol}` and the parsers under them              |
| MP4 export             | `@molstar/mp4-export-extension`                                                                  |

Beyond the slim target, the default compositions must register the same providers as the step-0 baseline (apart from
listed intentional differences), existing snapshots must restore in the Viewer, and the error and warning behavior of
§4.3 must hold. The full checklist is in the [plan](../plans/plugin-composition.md#6-acceptance-checklist).

## 13. Compatibility ledger entries

Record in [breaking-v6-changes.md](../plans/breaking-v6-changes.md) when implemented:

- `new PluginContext(spec)` starts with empty registries; a literal spec that did not spread the default spec adds
  `DefaultRegistry` to keep 5.x behavior (§3.1).
- `PluginSpec.actions`, `animations`, and `customFormats` are removed; `PluginContext` throws when they are present
  (§3.1).
- `DataFormatProvider` requires `name`; `DataFormatRegistry.add(name, provider)` with a differing name registers a named
  copy and warns; `BuiltIn*Formats` become provider arrays instead of `[name, provider]` tuples (§3.1, §4.3, §8).
- Format name types, category constants, providers, and transformers move (§6.3, §8); `BuiltInVolumeFormat` and
  `BuiltInShapeFormat` are added, and `BuildIn*` remain as deprecated aliases.
- `StateTransforms` is removed, except on the Viewer global (§6.3, §11).
- `PluginBehaviors`, `StateActions`, `StructureSelectionQueries`, `PresetStructureRepresentations`, and
  `PresetTrajectoryHierarchy` become catalogs; the `BuiltIn` namespace values are removed and the types remain (§7, §8).
- `DefaultPluginSpec` and `DefaultPluginUISpec` move to `default-spec` modules; `createPluginUI` requires a spec; the
  base UI, and a spec that replaces `components` without spreading the defaults, get the minimal structure tools (§10).
- `registry.default` may be `undefined`. Unregistered representation and theme names are still substituted; builders,
  param helpers, and snapshot restore now warn, naming the provider and scope. An empty registry is an error (§4.3).
- Representation and theme registries have `has(nameOrProvider)`; removing an unknown provider is a no-op; duplicate
  registration is reference-counted, and a different provider under an existing key throws (§4.2).
- Presets resolve through the registry by `id` or the new optional `alias`; built-in presets keep their short keys as
  aliases. An unresolved string throws instead of returning `undefined`. `applyPreset` is typed by the built-in id and
  alias unions, and `PluginConfig.Structure.DefaultHierarchyPreset` is added (§4.3, §5.2).
- `toAction()` is cached (§6.4).
- `PluginAnimationManager.current` can be `undefined`; `AnimateStateSnapshotTransition` is no longer registered
  implicitly (§4.4).
- `PluginDragAndDropEntry` replaces bare handlers in entries; `addHandler` with a different function under an existing
  name throws instead of replacing it (§4.2).
- `PluginViewModel`, `PluginUIViewModel`, and their hooks require an explicit spec (§9).
- `molstar.lib.plugin` values other than `StateTransforms` follow the 6.0 API (specs with the removed fields throw,
  literal specs get empty registries), and `DefaultRegistry` is added; `Viewer.create` and its options, including
  `customFormats`, are unchanged (§11).
- PyMOL, VMD, and Jmol scripts require a top-level import of their transpiler module or `transpilers/all`; the default
  registry and the Viewer import `all` (§7).
- `StructureComponentManager.setOptions` no longer inserts the focus representation behavior, and the `interactions`
  option takes effect only with the Interactions behavior (§7).
- The temporary "external color themes unavailable by default" entry is resolved (§5.4).

## 14. Open questions

- The UI contribution registry (§10) is designed after the plugin registry lands.
- Whether the development-mode default-theme check (§5.1) should become an error once in-repo entries are clean.
