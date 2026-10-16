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
  readonly languages?: readonly QueryLanguageEntry[];   // §2.7: the query languages target text is read in
  readonly overrides?: Readonly<Record<PartId, Part | null>>;                 // §4: replace or remove any part by id
  readonly sequence?: PrismViewerSequenceSpec | false;   // §2.3: the one sequence viewer, or none
  readonly agent?: PrismViewerAgentSpec | false;         // §2.4: the tools and skills, and the assistant that runs them
  readonly ui?: PrismViewerUISpec;                       // §2.5: layout, theme, keybindings; ignored headless
  readonly persist?: PrismViewerPersist;                 // §2.6: light state in localStorage and sessionStorage; kept records
  readonly palettes?: PrismViewerPalettes;               // §2.8: the lists the data palettes take
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
  ui: PD.Group({ showSequence, showHistory, showRemoteStates, collapsed, homeWhenEmpty, splash, buttons: { reset, screenshot,
                 expand, fullscreen, settings, selectionMode, animation, trajectory } }),
};
```

`homeWhenEmpty` (off by default) opens Home at start while nothing is loaded; the app turns it on when its URL loads
nothing, a page that drives the viewer itself leaves it off. `splash` (on by default) is the empty canvas's splash:
light through the Mol\* star (`PrismSplash`, prism-ui), in with the wordmark once the viewer is idle on the empty
canvas, answering the pointer.

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
  readonly runtime: () => Promise<SequenceRuntime>;            // Seq* and the views over it, its own chunk:
                                                               // () => import('@molstar/prism-viewer/sequence/runtime')
  readonly behaviors: readonly SequenceBehavior[];             // the viewer's (follow, secondary structure); extensions
                                                               // add theirs (§3: model-archive's pLDDT and PAE, validation)
  readonly style?: 'tracks' | 'text';                          // the initial one; Tracks by default
  readonly seq?: (base: SeqViewerSpec) => SeqViewerSpec;       // Seq*'s own spec (its behaviors, registry, render props)
}
```

The runtime loads when a structure is first shown, or when an operation (`sequence.*`, `load.svs`) or a saved state
needs it; until then the core holds the plugin side only (`sources/sequence.ts`) and a saved state as given
(`views.sequence.viewers` is undefined, `load()` makes them). A page shows its first frame without Seq\* (the app's
initial load: 4.25 MB, 1.25 MB gzipped, from 5.46 and 1.68), and a spec that names no runtime has none of it.

`false` builds no sequence viewer (no band, no Sequence button, no sequence commands or tools). One viewer, keyed
`'default'` (sequence.md §6).

### 2.4 `agent`: tools, skills, the assistant

```ts
interface PrismViewerAgentSpec {
  readonly core?: AgentListChange<ViewerToolSpec>;   // the tools from the first turn: the viewer's CoreTools, with these added or removed
  readonly skills?: AgentListChange<ViewerSkill>;     // the viewer's Skills, with these added (an id of a default's replaces it) or removed; extensions add theirs
  readonly assistant?: () => Promise<{ createAssistant: AssistantFactory; createViewerAgent }>; // assistant(options) from '@molstar/prism-assistant/spec'
}
interface AgentListChange<T> { readonly add?: readonly T[]; readonly remove?: readonly string[] }
// The declarations (CoreTools, Skills, their types) are plain data in `agent/skills.ts`, no imports: a spec that changes
// them does not bring the toolkit (`@molstar/ai` tools, zod) into the page's startup bundle (review B10).
// The assistant creates the agent itself: its models come from the person's keys, which the assistant keeps (agent.md §4).
interface AssistantHost {
  readonly viewer: PrismViewer;
  readonly createAgent: (models: ModelSource) => ViewerAgent;  // createViewerAgent over this viewer (agent.md §1)
}
interface AssistantUI {
  readonly ui: ViewerExtensionUI;              // its tab (and Inspector), registered as an extension's views
  readonly commandMenu: CommandMenuAssistant;  // the command menu's Ask row
  readonly agentView: () => Child;             // its turn card under the status line
  readonly commands?: readonly CommandDefinition[];  // the Inspector's (Assistant trace, tools, limits), once loaded
  dispose(): void;
}
type AssistantFactory = (host: AssistantHost) => Promise<AssistantUI> | AssistantUI;
type ViewerToolSpec = { id; operation; description?; fields?; fixed?; formats? };     // a tool over an operation (agent.md §2)
type ViewerSkill = { id; title; summary; guidance; tools: ViewerToolSpec[]; feature? };  // `feature`: whose operations
```

- The tools, skills and the session factory (`createViewerAgent`) are DOM-free (the viewer's `agent/` folder), so
  headless hosts (evals, MCP) use them without the assistant. The disclosure tools are always in the core, and the
  toolkit is dynamic (README decision 25).
- The assistant is lazy and opt-in: the viewer's package does not depend on `@molstar/prism-assistant`; the app passes
  the import, which is its own chunk. The shell loads it with the agent's code (another chunk, with the AI SDK) on first
  use: until then its tab is a placeholder whose opening loads it, and the command menu's Ask row loads it and then asks
  (or opens the setup when no model is set up); the app's `?assistant` opens the tab at start (G8 perf). Without it
  there is no Assistant tab, no Ask row, no turn card.
- `false`: no tools at all (a viewer for an embedding host that wants none): `createViewerAgent` gets an empty toolkit
  but for its own tools (`get_scene`, `ask`, `undo_turn`).

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

### 2.7 `languages`: what target text is read in

```ts
interface QueryLanguageEntry { id: 'pymol' | 'vmd' | 'jmol' | 'molscript'; name: string; hint: string; transpile(text: string): Expression }
// @molstar/prism-viewer/targets/languages/{pymol,vmd,jmol,molscript,all}: one module each, so a spec bundles only what it names
```

The languages a target's query text may be written in (targets.md §1). Keywords, address lists and MolQL selectors are
always read; text in a language the spec does not list fails, naming the languages it does read. `DefaultPrismViewerSpec`
lists all four (`AllQueryLanguages`, PyMOL first, the default for text); absent, none: the transpilers (84 kB minified)
are not in the bundle (`examples/minimal-prism-viewer`). The languages are registered when a viewer is made and are
code a page has, not a per-viewer setting: two viewers on a page read what either listed. The operations' `language`
parameter keeps its four values.

### 2.8 `palettes`: the data palettes' lists

```ts
interface PrismViewerPalettes {
  readonly prism?: readonly (number | string)[];       // Prism's own; Many-Distinct (Prism) when not given
  readonly colorblind?: readonly (number | string)[];  // the colorblind-safe one; Okabe–Ito when not given
}
```

The colors the categorical themes take (chains, entities, polymers, operators, models; README §3.7), in the order
categories take them, each `0xRRGGBB` or `'#rrggbb'`. The option `canvas.palette` picks the one in use: `prism` (the
default), `classic` (Mol*'s own Many-Distinct, the themes' unchanged) or `colorblind`; the person changes it in Settings,
a page with `?palette=`, the API with `canvas.palette`. A spec switches Prism's colors with one list
(`palettes: { prism: ['#07a57b', '#e57434', …] }`) and the default choice with `options: { canvas: { palette } }`. The
choice is kept in sessions (`prism-viewer.palette`); a session without it (from v5 or before Prism's palette) opens on
the classic one.

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
  readonly skills?: readonly ViewerSkill[];                    // its tools, in skills over its operations; never core tools
  readonly sequence?: readonly SequenceBehavior[];             // its sequence behaviors (sequence.md §3)
  readonly settings?: PD.Params;                               // reserved: its Settings group (plan 3.3)
  readonly ui?: (ctx: ViewerUIContext) => Promise<ViewerExtensionUI> | ViewerExtensionUI;  // its views; its own chunk
}

interface ViewerUIContext { readonly api: ViewerApi; readonly sources: ViewerSources; readonly managers: ViewerManagers }

interface ViewerExtensionUI {                                  // every id `<owner>.<kind>.<name>`, checked (§4)
  readonly sections?: readonly { id; slot: ViewerSlot; title; icon?; order?; keywords?; render: ViewRenderer }[];
  readonly tabs?: readonly { id; title; icon; dock?: DockId; order?; render: ViewRenderer }[];
  readonly overlays?: readonly { id; title; render: ViewRenderer }[];         // cards over the canvas
  readonly pages?: readonly { id; title; icon?; header?: ViewRenderer; render: ViewRenderer }[];  // over the whole app
                                                               // (the assistant's Inspector), `navigation.setPage`
  readonly rows?: readonly { id; when: (node: SceneNode) => boolean; render: ViewRenderer }[];  // scene tree rows under
                                                               // the nodes `when` holds for (the node is the view's
                                                               // `subject`); `SceneRow` draws one as the tree's own
                                                               // (picked, its controls in the inspector, Turn on while off)
  readonly toolbar?: readonly { id; icon; title; command?: string; overlay?: PartId; group?; order?;
                                place?: 'rail' | 'canvas' }[];   // runs a command, or opens one of its overlays
}

type ViewerSlot =                                              // closed and documented: the viewer decides where things go
  | 'home.import' | 'export' | 'export.picture' | 'selection' | 'canvas-settings' | 'animation' | 'settings' | `inspector:${SceneNode['kind']}`;

interface LazyViewerExtension {
  readonly id: string;
  readonly title: string;
  readonly builtIn?: true;                                     // its stub may declare operations in core areas
  readonly stub: {                                             // every part it will bring, and what it answers for
    readonly parts: readonly { id: PartId; kind: PartKind; slot?: ViewerSlot; title?: string }[];  // placeholders; overrides
    readonly operations?: readonly LazyOperationStub[];  // api.run of one loads the chunk and replays the call; its flags
                                                         // (read, history, undoable, approval, addresses, modules) are
                                                         // the real one's, checked when it loads
    readonly commands?: readonly { id: string; title: string; category?: string; icon?: IconName }[];
    readonly skills?: readonly ViewerSkill[];                // in the catalog; disclosing loads it (tools over its operations)
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
- **An extension names its own areas** (David, 2026-10-16): the viewer does not know the extensions it will get, so the
  names come from them (`membrane` from `anvil`, `video` from `video-export`: named for what they add), and each is the
  extension's alone. An area another extension or the viewer has fails the extension's activation, naming who has it;
  only the viewer's own features contribute operations to core areas.
- **Mounted at runtime**: `viewer.extensions.mount(ext)` adds an extension the spec did not list to a running viewer
  (its plugin part on the running plugin, its areas, commands and views; a lazy one as its stub) and returns its
  removal. It is the first step to extensions loaded from elsewhere (an open registry, later); today they still live
  in this repository.

### 3.1 The viewer's own features: a small core and built-in extensions

Both (David, 2026-10-15; the core made smaller 2026-10-16, so a viewer bundles only what its spec lists): what
everything else stands on is fixed code; the rest of the viewer's own features are `ViewerExtension`s like any other,
listed in `DefaultPrismViewerSpec().extensions`, so they can be eager or lazy, and an app filters, reorders or overrides
them as it does Mol\*'s extensions. The extension API is then used by the viewer itself, which keeps it honest.

- **The scene's API is a feature too**: `structure` (the `structure`, `contents`, `component`, `representation` and
  `styles` areas), `volume` (`volume`, `particles`), `trajectory`, `annotations` (it `requires` `structure`). Their
  areas are typed as present (api.md §2); left out of a spec, reading one throws, naming the feature. A feature's
  `requires` lists the extensions it needs before it; the spec check refuses it without them. An extension's
  annotation providers are offered only where the `annotations` feature is.

| Core (always there; parts still removable by `overrides`)                                       | Built-in features (in `DefaultPrismViewerSpec().extensions`)                                |
| ------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------- |
| the shell: rail, status column, History, Settings, theme, command menu, the Export card's frame  | `structure`, `volume`, `trajectory`, `annotations`: the scene's API (areas and commands)    |
| the canvas and its chrome: camera, expand and fullscreen, hover label, mode switch, selection bar | `home` (Home), `structure-panel` (Structure tab, scene tree, inspectors, Preset, Styles), `selection-panel` (Selection tab), `help` (Help) |
| loading files: drag and drop (`load.files`), `.molx` / `.molj` sessions                          | `measurements` (Selection › Measure, measurement nodes)                                     |
| selection, focus, highlight and targets (keywords, address lists); the history; the agent's account | `superposition` (Selection › Superpose); lazy                                             |
| the sequence viewer (configured by `sequence`, §2.3)                                              | `snapshots` (Snapshots tab, the playback card); `state-tree` (State tree tab); lazy         |
|                                                                                                  | `animation` (Animation card), `canvas-settings` (Canvas settings card)                      |
|                                                                                                  | `export-image`, `export-session` (eager); `remote-states` (behind its setting; lazy)       |
|                                                                                                  | `loaders` (database and URL loaders, `load.entry`, `load.url`, `load.mvs`, `load.snapshot`; URL parameters; `load.svs` later) |

- **The panels keep the core's ids**: Home, Structure, Selection and Help are features, but their tabs are `home`,
  `structure`, `selection`, `help` in the shell (`core.tab.*` as parts), as `layout.start`, `view.open` and the
  features' commands name them. A built-in feature brings them through two hooks extensions do not have: `views` (its
  view models, made with the viewer, headless too: `viewer.views.scene`, `homeActions`; typed as present, a missing one
  throws naming its feature) and `shell` (its tabs, sections and commands in the page's shell, with the shell's
  services: the command menu at a group, opening files).
- **What a smaller spec saves** (the minimal example, `examples/minimal-prism-viewer`, with only `structure`): 359 kB of
  the start bundle minified (98 kB gzipped) against the core with every panel in it, 2026-10-16.

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
  (`snapshots`, `geo-export`). Operation and command ids are their own (`<area>.<verb>`, operations.md §1). An
  extension's area is named for what it does, usually its id (`geo-export.export`), sometimes a shorter word the agent
  and the command menu read better (`video.export` for `video-export`, `membrane.fit` for `anvil`); never a core area's.
  The ids are listed in the viewer's `_spec`.
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
(no assistant), every query language (§2.7), its own built-in features as `extensions` (§3.1), no `persist`. Mol\*'s extensions are the app's,
since they depend on the viewer (README §3.1). `apps/prism-viewer`:

```ts
const spec: PrismViewerSpec = {
  ...DefaultPrismViewerSpec(),
  extensions: [...DefaultPrismViewerSpec().extensions!, ...MolstarViewerExtensions],   // built-ins, then Mol*'s
  options: optionsFromUrl(location.search),                    // v5's names kept
  agent: { ...DefaultPrismViewerSpec().agent, assistant: assistant() },  // from '@molstar/prism-assistant/spec'
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
