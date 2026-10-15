# `PrismViewerSpec`: what a Prism viewer is made from

The one object a Prism viewer is created from, as a plugin is from its `PluginSpec` and a Seq\* viewer from its
`SeqViewerSpec`: plain data, explicit ordered lists, nothing resolved by name, a fresh default from a function. It
gathers what README §3.1–3.5 decided into one contract. Status: **accepted** (2026-10-15).

## 1. At a glance

```ts
// @molstar/prism-viewer/spec: types only (no DOM, no plugin import), like @molstar/plugin/spec
interface PrismViewerSpec {
  readonly plugin: PrismViewerPluginParts;              // §2.1: the engine's base, before the extensions' parts
  readonly options?: Partial<PrismViewerOptions>;        // §2.2: one PD definition: JS options, URL, Settings
  readonly extensions?: readonly (ViewerExtension | LazyViewerExtension)[];   // §3, in order
  readonly overrides?: Readonly<Record<PartId, Part | null>>;                 // §4: replace or remove any part by id
  readonly sequence?: PrismViewerSequenceSpec | false;   // §2.3: the one sequence viewer, or none
  readonly agent?: PrismViewerAgentSpec | false;         // §2.4: the tools and skills, and the assistant that runs them
  readonly ui?: PrismViewerUISpec;                       // §2.5: layout, theme, keybindings; ignored headless
  readonly persist?: PrismViewerPersist;                 // §2.6: light state in localStorage and sessionStorage; kept records
  readonly finalize?: (spec: PluginSpec) => PluginSpec;  // §2.1: the app's last word on the composed plugin spec
}

function DefaultPrismViewerSpec(): PrismViewerSpec;     // a new object each call; the viewer's own features (§3.1)
```

## 2. The fields

### 2.1 `plugin`: the engine's base

```ts
interface PrismViewerPluginParts {
  readonly registry: readonly PluginRegistryEntry[];          // DefaultRegistry for parity; MVSRuntimeRegistry
  readonly behaviors: readonly PluginSpec.Behavior[];         // highlight, select, focus, camera, snapshots, custom props, MVS
  readonly config?: readonly [PluginConfigItem, unknown][];
  readonly canvas3d?: PartialCanvas3DProps;
}
```

- **Exported as named parts** (`PrismViewerBehaviors`, `PrismViewerRegistry`), so an app narrows or swaps them by
  identity, as with `DefaultRegistry`. Nothing of `PluginUISpec`: no `layout`, `components`, `customParamEditors`.
- **Composed** into one `PluginSpec` by the viewer (`sources/compose-plugin.ts`), in a fixed order: registry: the base,
  then each extension's entries in list order, then the custom formats; behaviors: the base, then the extensions';
  config: the base, the extensions', then the options' (§2.2), later wins. Then `finalize` sees the result.
- **Fixed by the viewer**: `PluginConfig.State.HistoryCapacity` is 0 (undo is the viewer's, sequence.md §7); the
  resize debounce is not used (README §3.5).

### 2.2 `options`: one definition, three readers

```ts
const PrismViewerOptions = {                                   // PD params; a value of each is the spec's options
  renderer: PD.Group({ pixelScale, pickScale, transparency, resolutionMode, powerPreference, preferWebgl1,
                       allowMajorPerformanceCaveat, disableAntialiasing, illumination }),
  canvas: PD.Group({ background, focusBehavior: 'default' | 'secondary-zoom' | 'disabled' }),
  data: PD.Group({ pdbProvider, emdbProvider, volumeStreamingServer, volumeStreamingDisabled,
                   stateServer, saccharideCompIdMapType }),
  ui: PD.Group({ showSequence, showHistory, showRemoteStates, collapsed, buttons: { reset, screenshot, expand,
                 fullscreen, settings, selectionMode, animation, trajectory } }),
};
```

- **In code**: `packages/prism/viewer/src/options.ts`; `auto` (transparency, resolution, WebGL 1) leaves the plugin's
  own detection. Turned into config items and canvas props by `sources/compose-plugin.ts`.
- **One definition, three readers**: `createPrismViewer({ options })` in code, the URL (v5's parameter names kept,
  `pixel-scale`, `pdb-provider`, `hide-controls`, …; the full list is [app.md](app.md) §4), and the Settings tab, which
  shows the same params. A value becomes config items and canvas props (renderer, data), the camera behaviors chosen
  (`focusBehavior`), or hidden parts (`ui.buttons`: a `false` is `overrides: { 'core.toolbar.<id>': null }`).
- Extension settings are not here: an extension brings its own group (`ViewerExtension.settings`), shown in Settings
  under its title.

### 2.3 `sequence`: the one sequence viewer

```ts
interface PrismViewerSequenceSpec {
  readonly behaviors: readonly SequenceBehavior[];             // the viewer's (follow, secondary structure); extensions
                                                               // add theirs (§3: model-archive's pLDDT and PAE, validation)
  readonly style?: 'tracks' | 'text';                          // the initial one; Tracks by default
  readonly seq?: (base: SeqViewerSpec) => SeqViewerSpec;       // Seq*'s own spec (its behaviors, registry, render props)
}
```

`false` builds no sequence viewer (no band, no Sequence button, no sequence commands or tools). One viewer, keyed
`'default'` (sequence.md §6).

### 2.4 `agent`: tools, skills, the assistant

```ts
interface PrismViewerAgentSpec {
  readonly core?: readonly string[];                           // tool ids offered from the first turn (the viewer's default core)
  readonly tools?: readonly ToolRegistration<ViewerApi>[];     // the viewer's tools (agent.md §3); extensions add theirs (§3)
  readonly skills?: readonly ToolGroup<ViewerApi>[];           // the viewer's skills; extensions add theirs
  readonly assistant?: () => Promise<AssistantFactory>;        // the UI: import('@molstar/prism-assistant/viewer')
}
type AssistantFactory = (viewer: PrismViewer, agent: ViewerAgent) => AssistantUI;   // ViewerAgent: createViewerAgent's
                                                                                    // result (agent.md §1); AssistantUI:
                                                                                    // its tab, Ask bar, turn card, Inspector
```

- The tools, skills and the session factory (`createViewerAgent`) are DOM-free (the viewer's `agent/` folder), so
  headless hosts (evals, MCP) use them without the assistant. The disclosure tools are always in the core, and the
  toolkit is dynamic (README decision 25).
- The assistant is lazy and opt-in: the viewer's package does not depend on `@molstar/prism-assistant`; the app passes
  the import, which is its own chunk. Without it there is no Assistant tab, no Ask bar, no turn card.
- `false`: no tools at all (a viewer for an embedding host that wants none).

### 2.5 `ui`: the shell

```ts
interface PrismViewerUISpec {
  readonly layout?: Omit<ShellOptions, 'buttons'> & {          // docks, mode, compact width, compact rail; buttons are parts (§4)
    readonly tabs?: Readonly<Record<string, { dock?: DockId; order?: number }>>;
    readonly start?: string;                                   // the tab open at start (Home)
  };
  readonly theme?: UITheme;                                    // 'light' | 'dark' | 'blue' | 'auto'
  readonly keybindings?: readonly { command: string; key: KeyChord | null }[];
}
```

What shows is decided by parts and `overrides` (§4), not here: there is one way to hide a tab, a button or a command.
Headless, `ui` is ignored.

### 2.6 `persist`

```ts
interface PrismViewerPersist {
  readonly key: string;                 // the namespace of every stored key ('prism-viewer'); one per viewer on a page
  readonly preferences?: PrismStore;    // localStorage: layout, theme, settings (small, lasting)
  readonly restore?: PrismStore;        // sessionStorage: this tab's light restore after a reload (journal.md §6)
  readonly rows?: RowStore;             // IndexedDB: a session record the person keeps (journal.md §1)
}
interface PrismStore { get(key: string): unknown; set(key: string, value: unknown): void; remove(key: string): void }
```

Each is optional and absent by default: nothing is kept. The stores are small by design: a value over its budget is not
written, and the viewer says so once (journal.md §6). Every stored value carries a version (journal.md §7). Two live
viewers on a page with the same `key` are refused at creation, naming it (README decision 29).

## 3. Extensions

Plain types and helpers, not a framework: every extension lives in this repository (README §2.13), so there is no
stability promise beyond the repository's own.

```ts
interface ViewerExtension {
  readonly id: string;                                         // 'geo-export'; every part id starts with it (§4)
  readonly title: string;
  readonly plugin?: Partial<PrismViewerPluginParts>;         // registry, behaviors, config
  readonly api?: Readonly<Record<string, (ctx: ApiAreaContext) => ApiArea>>;  // its areas by id, each made once (api.md §6)
  readonly contributes?: (ctx) => readonly { area; name; operation }[];      // into core areas: built-in features only;
                                                               // `name` is `memberName(operation.id)` (`remote.list` → `remoteList`)
  readonly builtIn?: true;                                     // one of the viewer's own features (§3.1)
  readonly commands?: (api: ViewerApi) => readonly CommandDefinition[];      // made once its operations exist
  readonly annotations?: readonly AnnotationProvider[];        // what it can show of a structure (`annotations.apply`)
  readonly tools?: readonly ToolRegistration<ViewerApi>[];     // reserved: read from phase 4 (agent.md §3)
  readonly skills?: readonly ToolGroup<ViewerApi>[];           // reserved: phase 4; never core tools
  readonly sequence?: readonly SequenceBehavior[];             // reserved: phase 5
  readonly settings?: PD.Params;                               // reserved: its Settings group (plan 3.3)
  readonly ui?: (ctx: ViewerUIContext) => Promise<ViewerExtensionUI> | ViewerExtensionUI;  // its views; its own chunk
}

interface ViewerUIContext { readonly api: ViewerApi; readonly sources: ViewerSources; readonly managers: ViewerManagers }

interface ViewerExtensionUI {                                  // every id `<owner>.<kind>.<name>`, checked (§4)
  readonly sections?: readonly { id; slot: ViewerSlot; title; icon?; order?; keywords?; render: ViewRenderer }[];
  readonly tabs?: readonly { id; title; icon; dock?: DockId; order?; render: ViewRenderer }[];
  readonly overlays?: readonly { id; title; render: ViewRenderer }[];         // cards over the canvas
  readonly rows?: readonly { id; when: (node: SceneNode) => boolean; render: ViewRenderer }[];  // not shown yet
  readonly toolbar?: readonly { id; icon; title; command?: string; overlay?: PartId; group?; order?;
                                place?: 'rail' | 'canvas' }[];   // runs a command, or opens one of its overlays
}

type ViewerSlot =                                              // closed and documented: the viewer decides where things go
  | 'home.import' | 'export' | 'selection' | 'canvas-settings' | 'animation' | 'settings' | `inspector:${SceneNode['kind']}`;

interface LazyViewerExtension {
  readonly id: string;
  readonly title: string;
  readonly builtIn?: true;                                     // its stub may declare operations in core areas
  readonly stub: {                                             // every part it will bring, and what it answers for
    readonly parts: readonly { id: PartId; kind: PartKind; slot?: ViewerSlot; title?: string }[];  // placeholders; overrides
    readonly operations?: readonly { id: string; title: string }[];  // api.run of one loads the chunk and replays the call
    readonly commands?: readonly { id: string; title: string; category?: string; icon?: IconName }[];
    readonly skills?: readonly { id: string; title: string; summary: string }[];   // in the catalog; disclosing loads it
    readonly owns?: { transformers?: string[]; formats?: string[]; fileExtensions?: string[];
                      representations?: string[]; themes?: string[] };   // a session or file that needs it loads it
  };
  readonly load: () => Promise<ViewerExtension>;
  readonly preload?: boolean;                                  // loaded during creation, before `ready` (§3.1)
}
```

- **Headless-safe by shape**: `plugin`, `api`, `commands`, `skills`, `sequence` and `settings` import no DOM; `ui` is a
  function, so its views are a chunk loaded when the UI first needs them (README §3.4).
- **Commands first** (README decision 27): an extension's first job is operations with commands and tools; `ui` only
  where a view is needed.
- **Where it lives**: `@molstar/<name>-extension/prism-viewer`, depending on `@molstar/prism-viewer` (the types above,
  from `spec.ts` and `api.ts`) and prism-ui; a lazy extension's stub lives in the same package, in a module that does not
  import the rest (`/prism-viewer/stub`). Data-only extensions are `{ id, title, plugin }` (the `pluginExtension` and
  `lazyExtension` helpers come with phase 5's ports, if they earn their keep).
- **A loaded extension is checked against its stub**: an area or command the stub did not declare fails the load, and
  so does a declared operation it does not bring (a stub never replays itself); its views are checked as they register
  (an undeclared one, or one not named `<id>.<kind>.<name>`, is not shown, and said). **Activation is all or nothing**:
  a failure undoes what it did and puts the stub back; a later call tries again.
- **A built-in feature's lazy stub** declares its operations in core areas too (`session.remote.list`): stubs are
  contributed to the area, and replaced by the feature's `contributes` when it loads.
- **Owned files and sessions**: `load.files` loads the lazy extension owning a file's extension first; `session.open`
  of a `.molj` the ones owning its transformers (a `.molx`, zipped, not yet).

### 3.1 The viewer's own features: a fixed core and built-in extensions

Both (David, 2026-10-15): what everything else stands on is fixed code; the rest of the viewer's own features are
`ViewerExtension`s like any other, listed in `DefaultPrismViewerSpec().extensions`, so they can be eager or lazy, and an
app filters, reorders or overrides them as it does Mol\*'s extensions. The extension API is then used by the viewer
itself, which keeps it honest.

| Fixed core (always there; parts still removable by `overrides`)                                | Built-in features (in `DefaultPrismViewerSpec().extensions`)                               |
| ----------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------- |
| the shell: rail, status column, History, Settings, Help, theme, command menu                    | `measurements` (Selection › Measure, measurement nodes)                                     |
| the canvas and its chrome: camera, expand and fullscreen, hover label, mode switch, selection bar | `superposition` (Selection › Superpose); lazy                                               |
| loading: Home (the registry's actions, `load.action`), drag and drop (`load.files`), `.molx` / `.molj` | `snapshots` (Snapshots tab, the playback card)                                         |
| the scene: Structure tab, scene tree, inspectors, presets, quick styles, components             | `state-tree` (State tree tab); lazy                                                         |
| selection, focus and targets; the history; the agent's account                                  | `animation` (Animation card), `canvas-settings` (Canvas settings card)                      |
| the sequence viewer (configured by `sequence`, §2.3)                                            | `export-image`, `export-session` (eager); `remote-states` (behind its setting; lazy)       |
|                                                                                                 | `loaders` (database and URL loaders, `load.entry`, `load.url`, `load.mvs`, `load.snapshot`; URL parameters; `load.svs` later) |

- **One module per feature** (review F44): `@molstar/prism-viewer/features/<name>` exports the feature itself;
  `@molstar/prism-viewer/features/lazy` exports a stub for each (`LazyFeatures.superposition`), whose `load` imports the
  feature's module, so importing the lazy list pulls in none of the features' code. The default list takes the lazy
  form for the heavier ones; an app that wants everything up front lists the features themselves (or sets
  `preload: true` on a lazy entry: loaded during creation, before `ready`, still its own chunk).
- **The single-file build** inlines every lazy entry (its dynamic imports bundled into the one file), from the app's
  list; it is built by `apps/prism-viewer` (app.md §2).
- **The app's list**: `[...DefaultPrismViewerSpec().extensions, ...MolstarViewerExtensions]`; the built-ins come first,
  so an extension's sections land after the viewer's own in a slot unless ordered otherwise.

## 4. Parts and overrides

- **A part** is anything with an id the viewer registers: a tab, a section, a structure row, a toolbar item, a command,
  a tool, a skill, a sequence behavior, an API module. One scheme (review F8): `<owner>.<kind>.<name>`, where the owner
  is `core` for the fixed core and the feature's or extension's id otherwise (`core.toolbar.xr`,
  `core.sequence.follow`, `snapshots.tab.snapshots`, `geo-export.section.export`). Feature and extension ids are plain
  (`snapshots`, `geo-export`). Operation and command ids are their own (`<area>.<verb>`, operations.md §1), and an
  extension's areas start with its id (`geo-export.export`). The ids are listed in the viewer's `_spec`.
- **`overrides`** (`extensions/overrides.ts`, one rule for every part, the core's and the shell's included): `null`
  removes a part, a value replaces its fields (its id kept): `{ 'core.toolbar.reset': null, 'camera.reset': null,
  'geo-export.section.export': { title: 'Export geometry' } }`. Commands by command id; views by part id. An override
  of an id nobody registers throws, naming it: core parts and lazy stubs' parts are checked when the spec is,
  commands after creation; an eager extension's views are known only when they register. A lazy extension's
  placeholders follow overrides too. A whole extension is dropped by leaving it out of `extensions`.
- **The core's views keep short ids in the shell** (`home`, `reset`, `sequence-legend`), which are
  `core.<kind>.<id>` as parts; `ui.layout` and `view.open` take either form.
- **Dependents follow**: removing an API module leaves its commands and tools listed with the reason they cannot run
  (as Stories' modules); removing a command removes its keybinding and toolbar item.
- **Duplicates throw** at creation, naming both owners; there is no silent last-wins.

## 5. Defaults and the app

`DefaultPrismViewerSpec()` holds the viewer's base: `plugin` (the named parts), the default `sequence` and `agent`
(no assistant), its own built-in features as `extensions` (§3.1), no `persist`. Mol\*'s extensions are the app's,
since they depend on the viewer (README §3.1). `apps/prism-viewer`:

```ts
const spec: PrismViewerSpec = {
  ...DefaultPrismViewerSpec(),
  extensions: [...DefaultPrismViewerSpec().extensions!, ...MolstarViewerExtensions],   // built-ins, then Mol*'s
  options: optionsFromUrl(location.search),                    // v5's names kept
  agent: { ...DefaultPrismViewerSpec().agent, assistant: () => import('@molstar/prism-assistant/viewer') },
  persist: { key: 'prism-viewer', preferences: localStore(), restore: sessionStore(), rows: createIndexedDbRows('prism-viewer-records') },
};
const viewer = await createPrismViewer({ target: document.getElementById('app')!, spec });
// ?pdb=1tqn … : the loaders feature reads the URL once the viewer is ready and runs its API calls
```

## 6. Creating a viewer

```ts
createPrismViewer({ target, spec }): Promise<PrismViewer>;          // in a page; <prism-viewer> wraps it (shadow root)
createPrismViewerHeadless({ spec, plugin? }): Promise<PrismViewer>; // @molstar/prism-viewer/headless: no DOM; `plugin`
                                                                    // builds a HeadlessPluginContext when images are wanted
interface PrismViewer {
  readonly plugin: PluginContext;                                   // the engine: an escape hatch
  readonly api: ViewerApi;                                          // every operation (README §3.3)
  run<T>(label: string, work: (plugin: PluginContext, api: ViewerApi) => Promise<T>): Promise<T>;  // one step; `api`'s calls join it
  readonly managers: ViewerManagers;                                // commands, transactions, history, navigation, viewState, notifications
  readonly views: ViewerViews;                                      // the root view models (shell, scene, selection, …)
  readonly sources: ViewerSources;                                  // read-only states over the plugin (state.md §4)
  readonly spec: PrismViewerSpec; readonly options: PrismViewerOptionValues;  // what it was made from
  readonly state: { ready; busy; empty };                           // ReadonlyStates (README §3.6)
  readonly events: { hover; click; error };                         // non-replaying events; hover and click carry targets
  dispose(): Promise<void>;
}
```

Creation, in order: validate the spec (ids, overrides, slots), compose the plugin spec (§2.1), create and `init()` the
plugin, create the API, managers and views, activate the extensions in order (lazy ones as stubs), then, with a page,
mount the canvas and render the shell; `ready` is set last, and only then do URL calls and the light restore run
(journal.md §6), so nothing that reads the camera runs before the canvas exists. `dispose()` undoes it in reverse, and a
creation that fails part way undoes what it made before throwing (a page's persist key included). An extension id that
is a core area's id is refused. The spec's layout (`ui.layout.docks`, `tabs`) and `options.ui.buttons` seed the
`viewState` and `navigation` managers; `ui.keybindings` replace the shortcuts of the commands they name.

## 7. From v5's `ViewerOptions`

| v5                                                                         | `PrismViewerSpec`                                                          |
| -------------------------------------------------------------------------- | -------------------------------------------------------------------------- |
| `extensions`, `disabledExtensions` (names in `ExtensionMap`)               | `extensions` (a list; leave one out)                                       |
| `customFormats`                                                            | an extension's registry entry, `pluginExtension('custom-formats', { registry: [{ formats }] })`, listed last (a PD value cannot hold providers) |
| `layoutIsExpanded`, `layoutShowControls`, `layoutControlsDisplay`          | `ui.layout` (`mode`, docks)                                                |
| `layoutShowSequence`, `layoutShowLog`, `layoutShowRemoteState`, `layoutShowLeftPanel`, `collapseLeftPanel`, `collapseRightPanel` | `options.ui` (the log is History, the panels the shell's) |
| `rcsbAssemblySymmetryDefaultServerType`, `…DefaultServerUrl`, `…ApplyColors` | the assembly-symmetry extension's settings                                |
| `viewportShow*`                                                            | `options.ui.buttons` (hidden parts)                                        |
| `viewportFocusBehavior`, `viewportBackgroundColor`, `illumination`         | `options.canvas`, `options.renderer`                                       |
| renderer and server options (`pixelScale`, `pdbProvider`, …)               | `options.renderer`, `options.data`                                         |
| `config: [[item, value]]`                                                  | `plugin.config`, or `finalize`                                             |

## 8. Open

1. Which features are core and which built-in (§3.1) is a first cut; it moves as the views are wired.
2. Part ids are named as the code is written, in §4's scheme, and listed in the viewer's `_spec`.
