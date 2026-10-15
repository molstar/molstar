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
                   stateServer, saccharideCompIdMapType, customFormats }),
  ui: PD.Group({ showSequence, showLog, showRemoteStates, collapsed, buttons: { reset, screenshot, expand,
                 fullscreen, settings, selectionMode, animation, trajectory } }),
};
```

- **One definition, three readers**: `createPrismViewer({ options })` in code, the URL (v5's parameter names kept,
  `pixel-scale`, `pdb-provider`, `hide-controls`, …; api.md), and the Settings tab, which shows the same params. A value
  becomes config items and canvas props (renderer, data), a camera behavior's params (`focusBehavior`), or hidden
  parts (`ui.buttons`: a `false` is `overrides: { 'viewer.toolbar.<id>': null }`).
- Extension settings are not here: an extension brings its own group (`ViewerExtension.settings`), shown in Settings
  under its title.

### 2.3 `sequence`: the one sequence viewer

```ts
interface PrismViewerSequenceSpec {
  readonly behaviors: readonly SequenceBehavior[];             // the viewer's (follow, pLDDT, secondary structure, PAE);
                                                               // extensions add theirs (§3)
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
  readonly skills?: readonly ToolGroup<ViewerApi>[];           // the viewer's skills; extensions add theirs (§3)
  readonly assistant?: () => Promise<AssistantFactory>;        // the UI and session: import('@molstar/prism-assistant/…')
}
```

- The tools and skills are DOM-free (`@molstar/ai-tools`; the viewer's `agent/` folder), so headless hosts (evals,
  MCP) use them without the assistant. The disclosure tools are always in the core (README §3.3).
- The assistant is lazy and opt-in: the viewer's package does not depend on `@molstar/prism-assistant`; the app passes
  the import, which is its own chunk. Without it there is no Assistant tab, no Ask bar, no turn card.
- `false`: no tools at all (a viewer for an embedding host that wants none).

### 2.5 `ui`: the shell

```ts
interface PrismViewerUISpec {
  readonly layout?: ShellOptions & {                           // docks, mode, compact width, compact rail
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
  readonly key: string;                 // the namespace of every stored key ('prism-viewer')
  readonly preferences?: PrismStore;    // localStorage: layout, theme, settings (small, lasting)
  readonly restore?: PrismStore;        // sessionStorage: this tab's light restore after a reload (journal.md §6)
  readonly rows?: RowStore;             // IndexedDB: a session record the person keeps (journal.md §1)
}
interface PrismStore { get(key: string): unknown; set(key: string, value: unknown): void; remove(key: string): void }
```

Each is optional and absent by default: nothing is kept. The stores are small by design: a value over its budget is not
written, and the viewer says so once (journal.md §6).

## 3. Extensions

Plain types and helpers, not a framework: every extension lives in this repository (README §2.13), so there is no
stability promise beyond the repository's own.

```ts
interface ViewerExtension {
  readonly id: string;                                         // 'geo-export'; every part id starts with it
  readonly title: string;
  readonly plugin?: Partial<PrismViewerPluginParts>;         // registry, behaviors, config
  readonly api?: (ctx: ViewerApiContext) => ApiModule;         // its operations, under its namespace (README §3.3)
  readonly commands?: (api: ViewerApi) => readonly CommandDefinition[];
  readonly skills?: readonly ToolGroup<ViewerApi>[];           // never core tools
  readonly sequence?: readonly SequenceBehavior[];
  readonly settings?: PD.Params;                               // its Settings group
  readonly ui?: () => Promise<ViewerExtensionUI> | ViewerExtensionUI;   // views in the viewer's slots; its own chunk
}

interface ViewerExtensionUI {
  readonly sections?: readonly { id: PartId; slot: ViewerSlot; title: string; render: ViewRenderer }[];
  readonly tabs?: readonly { id: PartId; title: string; icon: IconKey; dock?: DockId; render: ViewRenderer }[];
  readonly rows?: readonly { id: PartId; when: (node: SceneNode) => boolean; render: ViewRenderer }[];  // structure rows
  readonly toolbar?: readonly { id: PartId; command: string; icon: IconKey }[];  // canvas stack only
}

type ViewerSlot =                                              // closed and documented: the viewer decides where things go
  | 'home.import' | 'export' | 'selection' | 'canvas-settings' | 'animation' | 'settings' | `inspector:${SceneNode['kind']}`;

interface LazyViewerExtension {
  readonly id: string;
  readonly title: string;
  readonly stub: {                                             // what shows and what it answers for before it is loaded
    readonly commands?: readonly { id: string; title: string }[];
    readonly sections?: readonly { id: PartId; slot: ViewerSlot; title: string }[];
    readonly skills?: readonly { id: string; title: string; summary: string }[];
    readonly owns?: { transformers?: string[]; formats?: string[]; fileExtensions?: string[];
                      representations?: string[]; themes?: string[] };   // a session or file that needs it loads it
  };
  readonly load: () => Promise<ViewerExtension>;
  readonly preload?: boolean;                                  // loaded during creation, before `ready` (§3.1)
}

function pluginExtension(id: string, plugin: Partial<PrismViewerPluginParts>): ViewerExtension;  // data-only ones
function lazyExtension(e: LazyViewerExtension): LazyViewerExtension;
```

- **Headless-safe by shape**: `plugin`, `api`, `commands`, `skills`, `sequence` and `settings` import no DOM; `ui` is a
  function, so its views are a chunk loaded when the UI first needs them (README §3.4).
- **Where it lives**: `@molstar/<name>-extension/prism-viewer`, depending on `@molstar/prism-viewer/extension` (the
  types and helpers above) and prism-ui. Data-only extensions need no entry: the app writes
  `pluginExtension('rcsb-validation', { behaviors: [RCSBValidationReport] })`.

### 3.1 The viewer's own features: a fixed core and built-in extensions

Both (David, 2026-10-15): what everything else stands on is fixed code; the rest of the viewer's own features are
`ViewerExtension`s like any other, listed in `DefaultPrismViewerSpec().extensions`, so they can be eager or lazy, and an
app filters, reorders or overrides them as it does Mol\*'s extensions. The extension API is then used by the viewer
itself, which keeps it honest.

| Fixed core (always there; parts still removable by `overrides`)                                | Built-in features (`ViewerFeatures`, in the default list)                                   |
| ----------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------- |
| the shell: rail, status column, Log, Settings, Help, theme, command menu                        | `viewer.measurements` (Selection › Measure, measurement nodes)                              |
| the canvas and its chrome: camera, expand and fullscreen, hover label, mode switch, selection bar | `viewer.superposition` (Selection › Superpose); lazy                                        |
| loading: Home (the registry's actions), drag and drop, `.molx` / `.molj`                        | `viewer.snapshots` (Snapshots tab, the playback card)                                       |
| the scene: Structure tab, scene tree, inspectors, presets, quick styles, components             | `viewer.state-tree` (State tree tab); lazy                                                  |
| selection, focus and targets; the history; the agent's account                                  | `viewer.animation` (Animation card), `viewer.canvas-settings` (Canvas settings card)        |
| the sequence viewer (configured by `sequence`, §2.3)                                            | `viewer.export.image`, `viewer.export.session`; `viewer.remote-states` (behind its setting); lazy |
| | `viewer.loaders` (database and URL loaders over `@molstar/plugin-extension/loaders`; URL parameters) |

- **Eager or lazy is the list's choice, not the feature's code**: each built-in is exported both ways
  (`ViewerFeatures.superposition` and `ViewerFeatures.lazy.superposition`), and the default list takes the lazy form
  for the heavier ones. An app that wants everything up front lists the eager forms (or sets `preload: true` on a lazy
  entry: loaded during creation, before `ready`, still its own chunk).
- **The single-file build** inlines every lazy entry (its dynamic imports bundled into the one file), from the same
  list.
- **The app's list**: `[...DefaultPrismViewerSpec().extensions, ...MolstarViewerExtensions]`; the built-ins come first,
  so an extension's sections land after the viewer's own in a slot unless ordered otherwise.

## 4. Parts and overrides

- **A part** is anything with an id the viewer registers: a tab, a section, a structure row, a toolbar item, a command,
  a skill, a sequence behavior, an API module. Ids are `<owner>.<kind>.<name>`: the viewer's own are `viewer.*`
  (`viewer.tab.snapshots`, `viewer.toolbar.xr`, `viewer.sequence.follow`), an extension's start with its id
  (`geo-export.section.export`). The viewer's ids are listed in its `_spec` and kept stable.
- **`overrides`** replaces a part with another of its kind, or removes it with `null`, built-in or contributed alike:
  `{ 'viewer.toolbar.xr': null, 'geo-export.section.export': MyGeoExportSection }`. An override of an id nobody
  registers throws, naming it. A whole extension is dropped by leaving it out of `extensions`.
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
// ?pdb=1tqn … : the loaders feature reads the URL at start and runs its API calls
```

## 6. Creating a viewer

```ts
createPrismViewer({ target, spec }): Promise<PrismViewer>;          // in a page; <prism-viewer> wraps it (shadow root)
createPrismViewerHeadless({ spec, plugin? }): Promise<PrismViewer>; // @molstar/prism-viewer/headless: no DOM; `plugin`
                                                                    // builds a HeadlessPluginContext when images are wanted
interface PrismViewer {
  readonly plugin: PluginContext;                                   // the engine: an escape hatch
  readonly api: ViewerApi;                                          // every operation (README §3.3)
  run<T>(label: string, work: (plugin: PluginContext) => Promise<T>): Promise<T>;  // a script's history step (README §3.2)
  readonly managers: ViewerManagers;                                // commands, transactions, history, navigation, viewState, notifications
  readonly views: ViewerViews;                                      // the root view models (shell, scene, selection, …)
  dispose(): Promise<void>;
}
```

Creation, in order: validate the spec (ids, overrides, slots), compose the plugin spec (§2.1), create and `init()` the
plugin, create the API, managers and views, activate the extensions in order (lazy ones as stubs), then, with a page,
mount the canvas and render the shell. `dispose()` undoes it in reverse.

## 7. From v5's `ViewerOptions`

| v5                                                                         | `PrismViewerSpec`                                                          |
| -------------------------------------------------------------------------- | -------------------------------------------------------------------------- |
| `extensions`, `disabledExtensions` (names in `ExtensionMap`)               | `extensions` (a list; leave one out)                                       |
| `customFormats`                                                            | `options.data.customFormats` (registered last, as in v5)                   |
| `layoutIsExpanded`, `layoutShowControls`, `layoutControlsDisplay`          | `ui.layout` (`mode`, docks)                                                |
| `layoutShowSequence`, `layoutShowLog`, `layoutShowRemoteState`, collapse   | `options.ui`                                                               |
| `viewportShow*`                                                            | `options.ui.buttons` (hidden parts)                                        |
| `viewportFocusBehavior`, `viewportBackgroundColor`, `illumination`         | `options.canvas`, `options.renderer`                                       |
| renderer and server options (`pixelScale`, `pdbProvider`, …)               | `options.renderer`, `options.data`                                         |
| `config: [[item, value]]`                                                  | `plugin.config`, or `finalize`                                             |

## 8. Open

1. Which features are core and which built-in (§3.1) is a first cut; it moves as the views are wired.
2. Part ids are named as the code is written (`<owner>.<kind>.<name>`, §4) and listed in the viewer's `_spec`.
