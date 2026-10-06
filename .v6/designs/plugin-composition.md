# Mol\* 6.0: plugin composition

Design draft, 2026-10-06. This replaces the composition details in
[architecture §6](architecture.md#6-plugin-composition) and the `PluginFeature` proposal in the
[summary](summary.md#plugin-composition). It describes the target; none of it is implemented yet.

## 1. Goals and constraints

- **Include only what is imported and listed.** An application's plugin contains the formats, representations, themes,
  presets, and other providers that its spec lists. Those values come from static imports, so the bundle contains them
  and nothing else.
- **Cheap plugin creation.** Given a spec, `PluginContext` registers the listed entries in order. It does not resolve
  dependencies, look up providers by string name, or load catalogs.
- **Static builds.** Composition uses only static imports and synchronous registration. Single-file bundles (classic
  IIFE or one ESM file) tree-shake the same way as split builds. No API requires selecting providers by string name;
  name-based selection is acceptable only where an app has deliberately imported a full catalog.
- **Behaviors stay as they are.** Behaviors remain the mechanism for runtime-toggleable, parameterized functionality.
  The new registry entries do not live in the behavior state tree.
- **Compatibility.** Transformer ids, snapshot JSON, provider names, `PluginSpec.Action`/`Behavior`, and the
  `BuiltInTrajectoryFormat`-style name types remain.

## 2. Current coupling

Constructing `new PluginContext(spec)` loads the full built-in catalog even with an empty spec. `DefaultPluginSpec` adds
only actions, behaviors, and animations on top. No spec field can remove a built-in provider.

| Edge                                                          | Mechanism                                                                    | What it loads                                               |
| ------------------------------------------------------------- | ---------------------------------------------------------------------------- | ----------------------------------------------------------- |
| `context.ts` → representation registries                      | Constructors add every `BuiltIn` entry                                       | 17 structure, 5 volume, 4 particle representations          |
| `context.ts` → `ColorTheme`/`SizeTheme.createRegistry()`      | `ThemeRegistry` requires a built-in map; created for three scopes            | 37 color and 6 size themes, three times                     |
| `context.ts` → `DataFormatRegistry`                           | Constructor adds all built-in formats; each format module imports the facade | 37 formats, ~130 transformers, all IO readers/model formats |
| `context.ts` → `StructureBuilder` → preset builders           | Constructors add all presets; lookup checks the static map before registry   | 5 hierarchy and 11 representation presets                   |
| `context.ts` → `StructureSelectionQueryRegistry`              | Constructor preload                                                          | 65 queries                                                  |
| `context.ts` → `behavior.ts`                                  | `BuiltInPluginBehaviors` shares a module with all dynamic behaviors          | All dynamic behaviors, including eight custom-property ones |
| `context.ts` → `StructureComponentManager`                    | Direct value imports                                                         | `InteractionsProvider`, focus representation behavior       |
| `context.ts` → `MarkdownExtensionManager`, `DragAndDrop`      | Constructor preload; fallback handler imports `OpenFiles`                    | 14 markdown extensions; all formats again                   |
| `createPluginUI` → `DefaultPluginUISpec`                      | Fallback when no spec is passed                                              | The full default spec and volume-streaming UI               |
| Facade cycle: `transforms/model` → selection queries → facade | `StateTransforms` lazy getters work around it (issue #1791)                  | —                                                           |

Other findings that affect the design:

- About 17 call sites use `registry.default`, which is the first registered entry (themes: first after sorting). An
  empty registry would make them crash.
- Generic code contains format-specific logic. `builders.structure.parseTrajectory(blob)` is hard-wired to mmCIF;
  `DownloadStructure` builds its format list from a static array; the Cube volume format builds a structure through the
  model transforms.
- A snapshot that names an unimported transformer throws during data-tree restoration, after the behavior snapshot has
  already been applied. The plugin is left half-updated.
- `transformer.toAction()` creates a new action with a new UUID on each call. `ActionManager.remove(transformer)` never
  finds it, and adding the same transformer twice adds duplicates.
- Extensions register providers imperatively in behavior `register()`/`unregister()`. The MVS runtime's `Registrables`
  record is the declarative form of the same pattern.

## 3. The spec's registry

`PluginSpec` gains a `registry` field: a list of declarative entries that `PluginContext` registers during `init()`.

```ts
interface PluginSpec {
  registry?: PluginRegistryEntry[];
  behaviors: PluginSpec.Behavior[];
  config?: [PluginConfigItem, unknown][];
  canvas3d?: PartialCanvas3DProps;
  layout?: { initial?: Partial<PluginLayoutStateProps> };
}

interface PluginRegistryEntry {
  themes?: {
    structure?: { color?: ColorTheme.Provider[]; size?: SizeTheme.Provider[] };
    volume?: { color?: ColorTheme.Provider[]; size?: SizeTheme.Provider[] };
    particles?: { color?: ColorTheme.Provider[]; size?: SizeTheme.Provider[] };
  };
  representations?: {
    structure?: StructureRepresentationProvider[];
    volume?: VolumeRepresentationProvider[];
    particles?: ParticleRepresentationProvider[];
  };
  formats?: DataFormatProvider[];
  presets?: {
    hierarchy?: TrajectoryHierarchyPresetProvider[];
    structure?: StructureRepresentationPresetProvider[];
  };
  selectionQueries?: StructureSelectionQuery[];
  lociLabels?: LociLabelProvider[];
  markdownExtensions?: MarkdownExtension[];
  dragAndDrop?: PluginDragAndDropHandler[];
  actions?: PluginSpec.Action[];
  animations?: PluginStateAnimation[];
}
```

Field types are illustrative; exact provider types come from their owning modules. `PluginRegistryEntry` lives in
`@molstar/plugin/spec` and refers to graphics/model types with `import type` only.

The spec keeps only what is not a registration: `behaviors` (stateful, parameterized, part of the behavior state tree),
`config`, `canvas3d`, and `layout`. The former spec-level `actions`, `animations`, and `customFormats` move into
entries:

| 5.x spec field               | 6.0                                           |
| ---------------------------- | --------------------------------------------- |
| `actions: [...]`             | `registry: [{ actions: [...] }]`              |
| `animations: [...]`          | `registry: [{ animations: [...] }]`           |
| `customFormats: [[name, p]]` | `registry: [{ formats: [p] }]`, with `p.name` |

These fields are removed rather than kept as deprecated aliases: the migration is mechanical, and code that spreads
`DefaultPluginSpec().actions` fails immediately instead of silently registering nothing. Record the change in the
compatibility ledger and handle literal spec objects in the migration tool.

`DataFormatProvider` gains a required `name`, like representation and theme providers, so entries list providers
directly and registries key them consistently. Today the name exists only in the `[name, provider]` tuples of the
built-in lists and `customFormats`.

### 3.1 Rules

- **Pure data.** An entry contains provider values only. It cannot contain behaviors, config, or functions run at
  registration. This keeps entries from turning back into a feature system with its own lifecycle.
- **No resolution.** Entries have no `requires`. Each entry lists everything it directly needs (§5). What is registered
  is exactly what the spec lists.
- **Order.** `init()` registers entries in spec order, and within an entry in a fixed order: themes, representations,
  formats, presets, selection queries, loci labels, markdown extensions, drag-and-drop handlers, actions, animations.
  Registration happens after managers and builders exist and before behaviors are applied, so behaviors can rely on the
  registered providers.
- **Not serialized.** Entries are not part of the state tree or snapshots. They are never toggled at runtime.
- **One place for registrations.** Everything registered from the spec comes from `spec.registry`; there are no parallel
  spec-level lists.

### 3.2 Runtime registration

`PluginContext` exposes `register(entry: PluginRegistryEntry): () => void`. `init()` uses it for the spec's entries, and
the returned function undoes exactly that registration. Behaviors can use it in `register()` to replace hand-written
registration code; the MVS runtime's `Registrables` loop becomes a call to it. Existing imperative registration keeps
working, so extensions migrate when convenient.

## 4. Registries

All plugin registries start empty: structure/volume/particle representations, the color/size theme registries for each
scope, data formats, hierarchy and representation presets, selection queries, markdown extensions, drag-and-drop
handlers, and animations.

- **Identity and reference counting.** Registering the same provider object again increments a count. Registering a
  different provider under an existing name is an error. Unregistering decrements the count and removes the provider
  only at zero. This makes overlapping entries safe (two presets that both import ball-and-stick) and stops a behavior's
  `unregister()` from removing a provider that the spec also registered.
- **Explicit defaults.** `registry.default` returns the configured default if it is registered, otherwise the first
  registered applicable provider, otherwise throws an error naming the registry. Review each of the ~17 current call
  sites; most should pass the requested name and use the default only as a fallback.
- **Theme registries.** `ThemeRegistry` no longer requires a built-in map. `ColorTheme.BuiltIn`/`SizeTheme.BuiltIn` move
  to catalog modules (§7).
- **Preset lookup.** Preset builders resolve names through their registry only. Remove the static-map-first
  `resolveProvider` and the fixed `defaultProvider`; the default preset comes from config.
- **Format lookup.** `DataFormatRegistry` keeps priority-then-registration-order resolution for `auto()`.
  `dataFormats.get(name)` for an unregistered name throws an error that names the format and says it is not registered
  in this plugin.

## 5. What goes into an entry

A module that defines a provider can export a ready-made entry next to it. Plain inline records are equally valid:

```ts
registry: [
  Sdf,
  { representations: { structure: [SpacefillRepresentationProvider] } },
],
```

### 5.1 Import what you run

A provider statically imports the providers it directly uses, and its entry includes them. This makes bundling and
registration agree:

- A representation preset imports the representation and theme providers it builds, passes provider objects (not name
  strings) to the builder, and lists them in its entry. Registering the preset's entry registers everything it needs. A
  mistyped dependency becomes an import error.
- A format's entry includes its format provider and the actions that should appear for it in the UI.
- Registration is still required after importing. Snapshots store representation/theme names and
  `StructureRepresentation3D` resolves them through the plugin's registry; the UI lists only registered types.

### 5.2 Policy goes through config

Choices about what to show after loading are not imports. A format's `visuals` step applies the configured hierarchy
preset by name, and the hierarchy preset applies the configured representation preset
(`PluginConfig.Structure.DefaultRepresentationPreset`). If that preset is not registered, loading fails with a clear
error. Otherwise the SDF format would import the default preset and with it cartoon and the surface representations.

Apply the same rule to other generic code that names a specific representation or theme today, such as the focus
representation behavior (`ball-and-stick`, `physical`, `uniform`) and volume-streaming's direct use of the isosurface
provider. Either the code imports the provider it runs and its entry includes it, or it reads the choice from config.

### 5.3 Catalogs and defaults

Catalog modules assemble complete sets: for example `DefaultRegistry` (an array of entries), the built-in format list,
and the built-in representation/theme maps. Only default-spec modules and apps that deliberately want the full set
import them as values. Runtime modules may import catalogs only with `import type` (§7).

```ts
export const DefaultPluginSpec = (): PluginSpec => ({
  registry: DefaultRegistry, // includes the default actions and animations
  behaviors: [...],
});
```

The default custom-property behaviors (secondary structure, interactions, accessible surface area, and so on) remain
behaviors in `DefaultPluginSpec`; their runtime registration of themes and representations is unchanged.

## 6. Transformers

### 6.1 Global registration

Transformers stay in the global `StateTransformer` registry: importing a defining module registers its transformers by
id, as today. They are not listed in the spec. A transformer is available to snapshots if any imported module defines
it.

Snapshot loading checks every transformer id in the behavior and data trees against the global registry before applying
any part of the snapshot. If ids are missing, it throws one error listing them and leaves the plugin unchanged.

### 6.2 Module layout

Remove the `StateTransforms` facade and its lazy getters, and split the seven transform modules by functionality, with
subgroups:

```text
state/transforms/
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
  trajectory/{mmcif,cif-core,pdb,gro,xyz,lammps,mol,sdf,mol2}.ts
  coordinates/{dcd,xtc,trr,nctraj,lammps}.ts
  topology/{psf,prmtop,top}.ts
  volume/{ccp4,dsn6,mtz,dx,cube,density-server,segmentation,structure-factors}.ts
  shape/{ply,obj,vtp}.ts
  particles/{star,tbl,ndjson,em,simularium,mmcif-assembly}.ts
```

- Generic operations stay under `transforms/<area>/`.
- Format-specific parse and conversion transformers move next to their format provider and entry. For example,
  `formats/trajectory/sdf.ts` holds `TrajectoryFromSDF`, `SdfProvider`, and the `Sdf` entry.
- Move non-transformer helpers out of large modules: `VolumeRepresentation3DHelpers`, `getTrajectory`, `getBoxMesh`, and
  the measurement-data helpers.
- Format category constants (`TrajectoryFormatCategory` and siblings) move to small modules so importing a constant does
  not load a format family.
- Exact file names are settled during the split; the rule is that every import of a module is justified by what its
  users need, and no module imports a catalog as a value.

Transformer ids, names, and the `DeflateData` id typo stay unchanged. Identity checks such as
`transform.transformer === X` keep working because each transformer is still created once.

### 6.3 Actions

`transformer.toAction()` returns a cached action, so the same transformer always maps to the same action. Registering an
action twice is ignored by identity, and `ActionManager.remove(transformer)` works.

## 7. Built-in name types

`BuiltInTrajectoryFormat`, `BuiltInVolumeFormat`, `BuiltInShapeFormat`, `BuiltInCoordinatesFormat`,
`BuiltInTopologyFormat`, and the representation/theme name types remain. Each is derived in its catalog module from the
value list, and other modules import it with `import type`:

```ts
// state/formats/trajectory/catalog.ts: imports every trajectory format
export const BuiltInTrajectoryFormats = [
  ['mmcif', MmcifProvider],
  ['sdf', SdfProvider],
  // ...
] as const;
export type BuiltInTrajectoryFormat = (typeof BuiltInTrajectoryFormats)[number][0];

// state/builder/structure.ts
import type { BuiltInTrajectoryFormat } from '@molstar/plugin/state/formats/trajectory/catalog';
```

With `verbatimModuleSyntax` the type import is erased from emitted JavaScript, so bundles never contain the catalog. The
catalog keeps name-keyed tuples (or a name-keyed map) so the names stay literal types; a check verifies that each key
matches its provider's `name`, as `StructureRepresentationRegistry.BuiltIn` does today.

- The boundary check rejects value imports of catalog modules outside default-spec modules and apps.
- Declarations of consumers reference the catalog declarations; type-checking loads them, but bundles do not.
- Type imports must not point to a higher package; for example model cannot type-import a plugin catalog.
- A typed name does not prove registration. Runtime lookup reports a format or provider missing from this plugin.
- Generic defaults such as `R extends StructureRepresentationRegistry.BuiltIn` in the representation-parameter helpers
  switch to the same `import type` form.
- If fast types are adopted in v7, these derived types become explicit unions checked against the catalogs.

## 8. Decoupling `PluginContext`

Base modules must not import catalogs or optional functionality:

- Remove the built-in preloads from the representation, theme, data-format, selection-query, preset, and
  markdown-extension registries.
- Split `behavior.ts` so `BuiltInPluginBehaviors` (static command wiring) does not share a module with the dynamic
  behaviors.
- Move the interactions-specific code out of `StructureComponentManager`, and make its use of the focus representation
  behavior and selection queries explicit.
- Move the drag-and-drop fallback handler (which imports `OpenFiles`) into a registry entry included by the default
  registry.
- Remove the implicit `AnimateStateSnapshotTransition` registration when no animations are given; MVS registers the
  animation it needs in its entry.
- Make `parseTrajectory(blob)` format-neutral or move it into the mmCIF entry; build `DownloadStructure` format options
  from the registry; give the Cube format's structure path an explicit dependency on the transforms it uses.
- Base entry points (`@molstar/plugin/context`, `@molstar/plugin/spec`, `@molstar/plugin-ui`) must not import default
  specs or catalogs as values. Enforce this with the import-graph check.

## 9. UI

- `createPluginUI` requires an explicit spec; `DefaultPluginUISpec` stays in `@molstar/plugin-ui/default-spec`.
- UI code that hard-codes catalog entries (quick styles, volume and particle controls, `DefaultStructureTools`) follows
  §5: import what it runs, or read choices from config and the registries.
- Follow-up: give `PluginUISpec` the same idea for UI contributions (custom parameter editors, structure tools, import
  controls) that currently live on `PluginContext` as untyped maps. Design it after the plugin registry lands.

## 10. Extensions, MVS, and apps

- Extensions keep their behaviors. A behavior may replace hand-written registration with `plugin.register(entry)` and
  call the returned function in `unregister()`. An extension can also export a plain entry for providers that never need
  toggling.
- The MVS runtime's `Registrables` loop becomes `plugin.register`. Its loader maps format names directly to transformers
  and hard-codes representation names; review these against §5 so an MVS-only app registers what MVS needs.
- The Viewer composes `DefaultPluginUISpec` plus its extensions as today. Its `Viewer.create` options are unchanged:
  `customFormats` keeps its `[name, provider]` shape, and the Viewer turns it into a registry entry, using each tuple's
  name for the provider. The removal of `PluginSpec.customFormats` does not affect the Viewer API.

## 11. Acceptance

The slim-plugin target from architecture §6.5 becomes:

```ts
import type { PluginSpec } from '@molstar/plugin/spec';
import { PluginConfig } from '@molstar/plugin/config';
import { Sdf } from '@molstar/plugin/state/formats/trajectory/sdf';
import { DefaultHierarchyPreset } from '@molstar/plugin/state/builder/structure/hierarchy/default';
import { BallAndStickPreset } from '@molstar/plugin/state/builder/structure/presets/ball-and-stick';
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
await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');
```

Module paths are illustrative. The example must render an SDF ligand while the import graph and bundle exclude the mmCIF
and CCP4 parsers, cartoon and volume representations, the other presets, and MP4 export. Verify with the import-graph
check, esbuild metafile output for both a split and a single-file build, and a rendering smoke test. Also verify:

- `DefaultPluginSpec` and the Viewer register the same providers, in an order that keeps today's defaults.
- Existing snapshots restore in the Viewer; a snapshot needing an unimported transformer fails before changing state.
- Loading an unregistered format or preset produces a clear error.

## 12. Compatibility ledger entries

Record in [breaking-v6-changes.md](../plans/breaking-v6-changes.md) when implemented:

- `new PluginContext(spec)` starts with empty registries; `spec.registry` supplies providers.
- `StateTransforms` is removed; transformer modules move as in §6.2.
- `createPluginUI` requires a spec.
- `registry.default` semantics and preset lookup change as in §4.
- `PluginSpec.actions`, `animations`, and `customFormats` are removed; use registry entries (§3).
- `DataFormatProvider` requires `name`.
- The external color themes return as registry entries (see the [checklist](../plans/checklist.md)).

## 13. Implementation order

Each step keeps the build and the full default Viewer working.

1. **Registries:** identity/reference-counted registration, explicit defaults, `ThemeRegistry` without a built-in map,
   `plugin.register`, `spec.registry`, cached `toAction()`. Move built-in maps into catalogs and register them through
   `DefaultRegistry`; behavior is unchanged.
2. **Transform split:** split modules by functionality, move format transformers next to their providers, remove the
   facade and its cycle, and migrate consumers.
3. **Context decoupling:** remove constructor preloads and implicit imports (§8); presets and generic code follow §5.
4. **Snapshot pre-validation and errors** for missing transformers, formats, and presets.
5. **Slim acceptance**, boundary checks, and migration of extensions/MVS to `plugin.register` where it simplifies them.

## 14. Open questions

- Should registry entries be allowed to carry config defaults, for example a preset entry suggesting itself as the
  default? Current answer: no; config stays in the spec.
- Where should small shared entries live, for example `DefaultHierarchyPreset`? Proposed: next to the provider they
  register.
- Which selection queries belong in `DefaultRegistry` versus being imported directly by the presets and UI that use
  them?
