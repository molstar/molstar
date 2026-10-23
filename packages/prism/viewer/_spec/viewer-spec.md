# `PrismViewerSpec`

What a viewer is made from (`src/spec.ts`), the extension API, parts and overrides, creation. Plain data, explicit
ordered lists, nothing resolved by name, a fresh default from a function (`DefaultPrismViewerSpec()`). The code is the
truth for how.

## 1. Fields

`plugin` (§2.1), `options` (§2.2), `extensions` (§3), `languages` (§2.7), `overrides` (§4), `sequence` (§2.3), `agent`
(§2.4), `ui` (§2.5), `persist` (§2.6), `palettes` (§2.8), `finalize` (§2.1). Types only in `spec.ts`: no DOM, no value
import of the plugin.

## 2. The fields

### 2.1 `plugin`: the engine's base

- **Named parts** (`prismViewerPluginParts()`, `sources/base.ts`): registry and behaviors the viewer's features need
  (highlight, select, focus, camera, snapshots, the default custom properties, MolViewSpec). Nothing of `PluginUISpec`
  (no `layout`, `components`, `customParamEditors`); behaviors that only drive the classic UI are left out.
- **Composed in a fixed order** (`sources/compose-plugin.ts`): registry: base, extensions' entries in list order, custom
  formats last (v5's format order and `auto()`); behaviors: base, then extensions'; config: base, extensions', options',
  later wins. `finalize(spec)` has the last word. Headless the composition is plain data and runs in Node.
- **Conflicts throw at `init()`**, naming every one; behaviors apply with `ref: transformer.id`, so the same behavior
  twice collides. An extension's plugin part goes with the extension.
- **Fixed**: `PluginConfig.State.HistoryCapacity` 0 (undo is the viewer's, [sequence.md](sequence.md) §7); no resize
  debounce (README decision 28).
- **Lazy extensions' plugin parts** (`sources/extension-plugin.ts`): `plugin.register(entry)`, `config` set on the
  running plugin (later wins), behaviors applied as `initBehaviors()` does (custom properties first); formats register
  after those of creation. A registration conflict fails the load and leaves the stub.

### 2.2 `options`: one definition, three readers

- One PD definition (`src/options.ts`: `renderer`, `canvas`, `data`, `ui` groups) read by
  `createPrismViewer({ options })`, the URL ([app.md](app.md) §4, v5's names) and the Settings tab. A value becomes
  config items and canvas props, the camera behaviors (`focusBehavior`), or hidden parts (`ui.buttons.x: false` is
  `overrides: { 'core.toolbar.x': null }`). `auto` values leave the plugin's own detection.
- `ui.undoSteps` (10, 1–100): how many of the latest steps Undo takes back; only they keep a checkpoint, older steps
  stay in History, too old to undo ([sequence.md](sequence.md) §7.4). Read when the viewer is made.
- `ui.homeWhenEmpty` (off): Home at start while nothing is loaded; the app turns it on when its URL loads nothing.
  `ui.splash` (`random`; `quantum`, `light`, `off`; the app's `?splash=`): the empty canvas's splash. `ui.hoverLabel`
  (`corner`, the default as in v5: at the canvas' bottom right, clear of the rail, a floating panel and an open card,
  above the toolbar; `cursor`: beside the pointer): where the canvas' hover label shows. A long label (Markdown, a long
  entity name) is cut at `--msp-hover-label-height` (16rem) and fades out there; it takes no pointer, so it cannot
  scroll. `ui.xr` (`auto`; `always`, `never`; classic `ShowXR`): the toolbar's XR button, only where WebXR is,
  everywhere (where it is not, clicking says so in a notification), or not at all; never where the `canvas.xr` operation
  is not. The XR params are the Canvas card's XR section.
- Extension settings are the extension's own group (`ViewerExtension.settings`, reserved), not here.

### 2.3 `sequence`: the one sequence viewer

- `runtime` is a dynamic import (Seq\* and its views, a chunk), loaded when a structure is first shown or an operation
  or a saved state needs it; until then the core holds the plugin side only (`sources/sequence.ts`). A spec without a
  runtime bundles none of Seq\*.
- `false`: no band, no Sequence button, no sequence commands or tools. One viewer, keyed `'default'`
  ([sequence.md](sequence.md) §6).

### 2.4 `agent`: tools, skills, the assistant

- `core` and `skills` are `{ add?, remove? }` changes over the viewer's defaults (an added id of a default replaces it):
  data, composes by id, needs no import of the defaults.
- The declarations are plain data in `agent/skills.ts` (no imports), so changing them keeps the toolkit (`@molstar/ai`,
  zod) out of the startup bundle. Tools, skills and `createViewerAgent` are DOM-free, for headless hosts.
- `assistant` is a lazy loader (`assistant()` from `@molstar/prism-assistant/spec`): the viewer never depends on the
  assistant package. It returns `AssistantFactory` (`(host) => AssistantUI`: its tab and Inspector as extension views,
  the command menu's Ask row, the turn card); the assistant creates the agent itself (`host.createAgent(models)`), since
  the models come from the person's keys ([agent.md](agent.md) §4). Until first use its tab is a placeholder whose
  opening loads it.
- `false`: no tools; `createViewerAgent` gets only its own.

### 2.5 `ui`: the shell

`layout` (docks, mode, compact width, `tabs` placement, `start` tab, `compactRail`), `theme`, `keybindings` (the
keymap's app layer: `{ command, key }`, `key: null` removes one; ui's commands.md §6). The person's remaps (Help) are
kept with the preferences (`persist.preferences`, key `keybindings`; `commands/keymap-store.ts`); Forget on them in
Storage also resets the keymap's person layer (`viewerStorageSources` `resetKeybindings`). The Select commands made from
the loaded structures (`select.entity.*`, `.chain.*`, `.residue.*`, `.element.*`) exist only once a structure is loaded:
an override or a keybinding of a well-formed id of those families is known at creation (`isStructureQueryId`). What
shows is decided by parts and `overrides` (§4), one way to hide anything. Ignored headless.

- **The embedded layout** (`layout.compactRail: true`; the app's `?hide-rail`): the rail is hidden behind a button in
  the canvas's corner and the canvas goes edge to edge; the button slides the rail in over the canvas, a panel opens
  beside it. For a viewer inside a page (Storybook: Viewer / With extensions, Embedded; the Seq\* introduction).
  `layout.mode: 'embedded'` is that layout as a mode, with the rest of what an embedded viewer needs
  ([shell.md](../../ui/_spec/shell.md) §3): no pinned or docked right panel, compact where narrow, and no panel at start
  (not even Home with `homeWhenEmpty`).
- **The kept layout**: with `persist.preferences`, the layout the person leaves is restored at start (§2.6); the spec's
  `ui.layout` is what Reset layout (`view.reset-layout`) goes back to. It holds the docks, their tabs, sizes, float and
  pin and the hidden tabs; the mode, the compact rail and the buttons stay the spec's (`restoreLayout`).
- **The canvas toolbar** (`ui.layout.toolbar`: `'auto' | 'expanded' | 'collapsed'`, default `auto`) is the host's
  choice, not the person's: no control, no command, not kept. `collapsed` is one button whose menu lists the tools (ui
  [shell.md](../../ui/_spec/shell.md) §2), `expanded` the full bar, `auto` the button in the compact mode (below
  `compactWidth`, an embedded layout too) and the bar otherwise, following the mode as the viewer resizes. The viewer
  gives the menu the mode choices (`modeMenuRows`: Interact and Select, the picking level as a submenu) through the
  shell's `canvasLeadMenu`.

### 2.6 `persist`

`key` (namespace of every stored key), `preferences` (localStorage: the theme and the layout), `layout` (keep the layout
in `preferences`: on, except for an embedded layout, `ui.layout.mode: 'embedded'` or `compactRail`, which is the host
page's), `restore` (sessionStorage, the light restore, [journal.md](journal.md) §6), `rows` (IndexedDB, kept records,
[journal.md](journal.md) §1). Each absent by default: nothing kept. Values are versioned and budgeted
([journal.md](journal.md) §7). Two live viewers on a page with the same `key` are refused at creation (`page.ts`, README
decision 29).

### 2.7 `languages`: what target text is read in

`QueryLanguageEntry[]` (`targets/languages/{pymol,vmd,jmol,molscript,all}`, one module each, so a spec bundles only what
it names; [targets.md](targets.md) §1). Keywords, address lists and MolQL are always read; text in an unlisted language
fails, naming those read. The default lists all four, PyMOL first (the default for text); none keeps the transpilers out
of the bundle (`examples/minimal-prism-viewer`). **Per viewer** (README decision 29): each viewer has its own set
(`QueryLanguages`, `targets/languages.ts`), made when it is made from the page-wide defaults (`registerQueryLanguages`:
none unless a host sets some) and then its spec's; one viewer's languages are not offered by another, in its `language`
params, targets or command-menu `target` field. The target functions (`compileTarget`, `parseTarget`, `compileQuery`)
take the viewer's set as `languages` and read the page-wide one without.

### 2.8 `palettes`: the data palettes' lists

`{ prism?, colorblind? }`: the colors the categorical themes take (chains, entities, polymers, operators, models), in
order. `options.representation.defaultPalette` picks `prism` (default), `classic` (Mol\*'s own) or `colorblind`
(Okabe–Ito); also Settings, `?palette=`, the `representation.defaultPalette` operation. Kept in sessions (custom state
`prism-viewer.palette`); a session without it opens on `classic` (`sources/palette.ts`).

## 3. Extensions

A `ViewerExtension` (`spec.ts`) is a plain object: `id`, `title`, `plugin` (§2.1), `api` (its areas by id, each made
once, [api.md](api.md) §6), `contributes` (operations into core areas: built-in features only), `builtIn`, `requires`,
`commands(api)`, `annotations`, `skills`, `sequence` (behaviors, [sequence.md](sequence.md) §3), `settings`, `ui(ctx)`
(its views, a chunk). Types and helpers, not a framework: every extension lives in this repository, no stability promise
beyond it (README decision 13).

- **Headless-safe by shape**: everything but `ui` imports no DOM; `ui` is a function, so views load when first needed.
- **Commands first** (README decision 27): operations with commands and tools; `ui` only where a view is needed.
- **Views go in named slots** (`ViewerSlot`, closed: Home's import, Export, Selection, canvas settings, animation,
  Settings, inspector sections), plus `tabs`, `overlays` (cards), `pages` (over the app), scene tree `rows`, `toolbar`
  items, and `selectionActions`: rows of the selection bar's menu (a title, a line, a command or `run`, by `order`; the
  clipping extension's Clip…), part ids `<ext>.selection-action.<name>`, and `nodeActions`: rows of the scene tree's
  node menus (the inspector shows the same menu) on the node kinds given and passing `when`, after the node's own and
  before Remove, by `order`; `run` gets the node (the clipping feature's Clip this…), part ids
  `<ext>.node-action.<name>`. No cards or rail buttons of its own. A tab for what is loaded (the Kinemage panel) passes
  `visible`: on the rail only while it holds, so an app with the extension shows nothing for it until it has content
  (not the person's hidden tabs). A panel with a tree of its own draws it with `SceneRow` (`children`, `depth`;
  `onActivate` for a row that is an action, a kinemage's view), so it reads as the Structure panel. A row's eye that
  turns something on (a row with controls; the viewer's Volume streaming and Unit cell too) picks it: its controls are
  in the inspector at once. A scene tree row is one row of the tree's height (`SceneRow`; the tree is virtualized,
  view-models.md §3): its controls go to the inspector, never under it. A row whose `when` holds takes a slot, so its
  view renders a row whenever `when` holds (one that renders nothing leaves a gap); `when` says what the view would
  (assembly symmetry: only a PDB entry's assembly). A `when` that follows more than the node and its structure passes
  `whenChanged`, told when it may answer otherwise (a PAE fetched is a state cell, no change of the scene): the tree
  asks it again for the rows it shows. Streaming starts in the Auto view (around the whole structure, the selection
  while there is one). A structure with several frames says how many on its row ("1,231 atoms · 10 models"; `frames` for
  a simulation's coordinates) and has a row of its own first under it (`views/structure/trajectory-row.tsx`): the frame
  shown, a scrubber, ‹ ›, play; picked, the same at full width and how it plays (mode, seconds per pass) in the
  inspector. The row plays that structure (`trajectory.play` with `structure`; without one, the operation plays every
  structure's frames, and every row says so); starting one stops what played. The canvas' trajectory bar follows the
  structure that plays (paused, the one that played last): while that structure's row is on screen (the Structure tab
  shown, the structure open, not under a card), a compact chip while frames play (name, frame, pause); otherwise the
  full bar (scrubber, ‹ ›, play), paused or not, so the controls are always at hand (another tab, the panel closed, the
  canvas expanded). For a structure in a simulation box the full bar ends with the box menu (molecules whole, in the
  box, centered on the polymers, no jumps, as in the file; `trajectory.periodic`) and a periodic images switch
  (`structure.symmetry`, the cells -1 to 1).
- **Where it lives**: `@molstar/<name>-extension/prism-viewer`, depending on `@molstar/prism-viewer` and prism-ui; a
  lazy stub in a module that imports nothing else. Data-only extensions are `{ id, title, plugin }`
  (`apps/prism-viewer/src/data-extensions.ts`).
- **Lazy** (`LazyViewerExtension`): a `stub` declares every part (ids, kinds, slots), operations with their flags,
  commands, skills, and what it `owns` (transformers, formats, file extensions, representations, themes, and a `sniff`
  for files whose extension does not tell: mmCIF holding particles); `load()` brings the real one. First use loads the
  chunk (a status task) and replays the call. `preload: true` loads during creation.
- **Checked against its stub**: an undeclared area, command or view, a declared operation not brought, or other flags
  fail the load. **Activation is all or nothing**: a failure undoes it and puts the stub back; a later call retries
  (`extensions/host.ts`).
- **Owned files and sessions**: `load.files` loads the extension owning a file's extension, or claiming its content by
  `sniff` (read only for the extensions a stub sniffs), first; `session.open` of a `.molj` those owning its
  transformers.
- **An extension names its own areas**: a name another extension or the viewer has fails activation, naming who has it;
  only built-in features contribute to core areas.
- **Mounted at runtime**: `viewer.extensions.mount(ext)` adds one to a running viewer and returns its removal.

### 3.1 The viewer's own features: a small core and built-in extensions

- **Core** (fixed code; its parts still removable by overrides): the shell (rail, status, History, Settings, theme,
  command menu, the Export card's frame), the canvas and its chrome, drag and drop and `.molx`/`.molj` sessions,
  selection, focus, highlight, targets, the history, the agent's account, the sequence viewer.
- **The selection bar** (`views/canvas/selection-bar.tsx`, `selectionModel`): Unfocus first while something is focused
  (`focus.clear`, the selection kept; beside the canvas toolbar when nothing is selected); what is selected is the
  button of a menu (hover opens it, as the canvas stack's menus; the stack's corner mark says so), and Clear. The menu:
  Focus, Show in structure tree, Make a component, the extensions' rows, More in the Selection panel. Show in structure
  tree turns on the tree's selection filter: its search's non-text part (`SceneModel.selectionFilter`, a chip before the
  text), the rows that hold the selection (structures; components and their representations it overlaps; Contents'
  entities and chains), asked per row and cached until the selection changes, followed only while on, off when the
  selection is cleared; the first component that holds it picked.
- **Built-in features** (`DefaultPrismViewerSpec().extensions`, `src/features/<name>`): `structure`, `volume`,
  `trajectory`, `annotations` (the scene's API; `requires` lists what must come before), `home`, `structure-panel`,
  `particles` (lazy: Mol*'s particles are opt-in, `ParticlesRegistry` is its plugin registry, so a viewer holds no
  particle code until a particle file, an mmCIF holding particles or a session with particles needs it),
  `selection-panel`, `help`, `loaders`, `canvas-settings`, `measurements`, `clipping` ([clipping.md](clipping.md)),
  `superposition` (lazy), `state-tree` (lazy), `remote-states` (lazy, behind its setting), `snapshots`, `animation`,
  `export-image`, `export-session`. Opt-in (not in the default spec; the Prism app adds it): `structure-move`
  ([structure-move.md](structure-move.md)).
- Built-in features' areas and operations in core areas are typed as present (`BuiltInApis`,
  `BuiltInOperationFeatures`); left out of a spec, reading one throws, naming the feature ([api.md](api.md) §2).
- **Panels keep the core's ids** (`home`, `structure`, `selection`, `help`; `core.tab.*` as parts). Built-in features
  have two hooks extensions lack: `views` (view models made with the viewer, headless too: `viewer.views.scene`) and
  `shell` (tabs, sections, commands in the page's shell).
- **One module per feature**: `@molstar/prism-viewer/features/<name>`; `features/lazy` exports stubs whose `load`
  imports the feature, so the lazy list pulls in no feature code.
- **The app's list**: `[...DefaultPrismViewerSpec().extensions, ...MolstarViewerExtensions]` (built-ins first); Mol\*'s
  extensions are the app's, since they depend on the viewer (`apps/prism-viewer/src/extensions.ts`). The single file
  inlines every lazy entry ([app.md](app.md) §2).

### 3.2 Particles

- **Opt-in, lazy** (`features/particles`, `LazyParticles`): Mol*'s `ParticlesRegistry` is the feature's plugin registry,
  `DefaultRegistry` holds none of it. The stub `owns` the particle transformers, formats and extensions (`.star`,
  `.tbl`, `.ndjson`, `.em`, `.simularium`) and sniffs `.cif`/`.bcif` for CellPack and PetWorld assemblies
  (`mmcif-sniff`); a `.molx` session is read for its transformers too. Why: most viewers never read particles.
- **Operations** (`particles.*`, operations.md): `list`, `representation`, `color`, `visibility`, `remove`, `focus`,
  `unit-cell` (as `structure.unit-cell`), `targets` (what stands for the particles: a file, a PDB id or a URL, one
  target per call), `source` (how the list is read: the format transform's params, the inspector's Source form),
  `frame`, `play`, `stop` (a trajectory's frames: viewing, no step).
- **Read at the end of the decorators**: operations and the scene read a list as its representations see it
  (`StateTree.getDecoratorRoot`: `ParticleListWithTargets` adds the targets); the hierarchy's ref holds the list below
  them.
- **Scene** (`sources/scene.ts`, by transformer id: no import of the feature): a list's row says its count and targets;
  its menu adds a representation, colors all, sets targets, focuses, shows the unit cell; a representation's row says
  its color theme and picks one; the unit cell is a row; the inspector lists source, orientations, targets, entities,
  compartments, fibers, attribute ranges, the cell, and the Source form. A trajectory's list has the frame row and the
  canvas' trajectory bar follows it (`frameInfos`, `kind: 'particles'`).
- **Default size**: particles without radii are drawn at half the typical nearest-neighbor distance
  (`getParticleSpacingRadius`, sampled, cached); an explicit size or radius wins.
- Measured (2026-10-21, headless): 100k picks load in 0.55 s, their default size in 84 ms once; a Simularium frame step
  is 1.4 ms with 20 structures loaded.

## 4. Parts and overrides

- **A part** is anything with an id the viewer registers (tab, section, row, toolbar item, command, tool, skill,
  sequence behavior, API area). Id scheme: `<owner>.<kind>.<name>`, owner `core` or the feature's or extension's plain
  id (`core.toolbar.xr`, `snapshots.tab.snapshots`, `geo-export.section.export`). Operation and command ids are
  `<area>.<verb>` ([api.md](api.md) §1).
- **`overrides`** (`extensions/overrides.ts`, one rule for every part): `null` removes it, a value replaces its fields
  (id kept). Commands by command id, views by part id. An override of an id nobody registers throws: core parts and lazy
  stubs' parts when the spec is checked, commands after creation, an eager extension's views when they register.
- The core's views keep short ids in the shell (`home`, `reset`) = `core.<kind>.<id>` as parts; `ui.layout` and
  `view.open` take either.
- **Dependents follow**: a removed area leaves its commands and tools listed with the reason they cannot run; a removed
  command takes its keybinding and toolbar item.
- **Duplicates throw** at creation, naming both owners.

## 5. Defaults and the app

`DefaultPrismViewerSpec()` (`src/default-spec.ts`): the plugin parts, the built-in features, every query language, the
sequence runtime, `agent: {}` (no assistant), no `persist`. The app (`apps/prism-viewer`) adds Mol\*'s extensions, URL
options, the assistant and `persist`; its loaders read the URL once the viewer is ready.

## 6. Creating a viewer

- `createPrismViewer({ target, spec })` (page, shadow root; `<prism-viewer>` wraps it) and
  `createPrismViewerHeadless({ spec, plugin? })` (no DOM; `plugin` makes a `HeadlessPluginContext` for images). Both
  give a `PrismViewer` (`src/viewer.ts`): `plugin` (escape hatch), `api`, `run(label, work)` (one step), `managers`,
  `views`, `sources` (read-only), `camera` (the view to follow and set: api.md §12), `spec`, `options`, `state`
  (`ready`, `busy`, `working`, `empty`: `busy` is the plugin's debounced flag, `working` anything in progress now, an
  operation from its start, a lazy extension's load, the plugin updating; the empty canvas's splash waits on `working`,
  as the plugin's flag comes 750 ms late and misses what is not plugin work), `events` (`hover`, `click`, `error`),
  `extensions`, `dispose()`. Several viewers on a page: each its own `persist.key` (§2.6, README decision 29).
- **Order**: validate the spec (`spec-check.ts`: ids, overrides, slots), compose the plugin spec, `init()` the plugin,
  make the API, managers and views, activate the extensions in order (lazy as stubs), then with a page mount the canvas
  and render the shell; `ready` last, and only then URL calls and the light restore ([journal.md](journal.md) §6).
- `dispose()` undoes it in reverse; a creation that fails part way undoes what it made (the persist key included) before
  throwing. `ui.layout` and `options.ui.buttons` seed `viewState` and `navigation`.
