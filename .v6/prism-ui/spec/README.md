# The Prism viewer: the specification

How the Prism viewer is built: Mol\*'s own UI on Prism, at feature parity with the v5 Viewer, driven by commands that a
person, a script or an agent runs alike. This spec is what the build follows; the designs in [../design](../design) are
the reasons and the mocks behind it, and where they disagree, this spec holds. Once the code exists, a small `_spec/`
next to it (as `packages/ai/_spec`) describes what was built, and this folder keeps the decisions and the order of work.

Status: **accepted, with the review's fixes** (2026-10-15). The decisions are in §2; the designs in §3 and the files
beside this one were accepted by David on 2026-10-15 and revised after the review ([review.md](review.md)); §5 lists
what is written and what is to come.

## 1. What is built

```
apps/prism-viewer          private: index.html, the entry, URL options, the extension list; built (ESM, chunks), deployed on the web
  └ @molstar/prism-viewer  public: the viewer: its state, commands, views, the extension API, the JS API, a headless entry
      ├ @molstar/prism-ui  public: the library: runtime, components, shell, commands, ink, the assistant's ports
      ├ @molstar/plugin    the engine (state tree, canvas, formats, MVS), behind the viewer's sources/ and api/
      └ @molstar/seq-*     the sequence view
  └ extensions             each ships a prism-viewer entry against the viewer's extension API (§3.1)
@molstar/prism-assistant   opt-in: the agent's UI, on prism-ui's ports
apps/prism-stories         later: the stories app (an improved Stories Untold), its own composition on prism-ui
```

## 2. Decisions

1. **prism-ui is a library, not a framework** (David, 2026-10-14). It holds what any Mol\* UI on Prism needs: the DOM
   runtime, components, the shell model and view, commands and the command menu, the canvas ink, status, theme, the
   params form, the ports the assistant attaches to. It has no opinion on how an app is composed: no modules, no
   contributions registry, no modes. The framework in [../design/contexts.md](../design/contexts.md) §3–4 and
   [../design/shell.md](../design/shell.md) §5 is not built.
2. **Opinionated apps on it.** The Prism viewer and, later, Prism stories each own their features in their own code
   (as the sequence view does, [../../seq/prism-sequence.md](../../seq/prism-sequence.md) §2).
3. **No modes** (David, 2026-10-14). The viewer and stories are separate apps. No `PrismMode`, no splash, no "Switch
   mode…"; [../design/modes.md](../design/modes.md) is superseded. A link from one app to the other is a navigation.
4. **The public package builds the viewer; the app is thin.** The reverse of v5, where `apps/viewer` is both the app and
   the library. `@molstar/prism-viewer` creates the viewer from a composition (its options and extensions) and is
   published; `apps/prism-viewer` is private: the page, URL options, the list of extensions it ships, and an ESM build
   with code splitting, deployed on the web. Packages under `packages/` are public at release; apps are not.
5. **Extensions bring their own UI, the viewer can override it** (David, 2026-10-14). An extension package ships a
   `prism-viewer` entry built against the viewer's extension API; the viewer's composition can replace or hide any part
   of an extension by id. Heavy extensions are lazy: a stub until first use, then a chunk (§3.1).
6. **Parity with v5 in function, not in form.** Everything the v5 Viewer does has a place (the inventory in
   [../analysis/current-ui.md](../analysis/current-ui.md) §5 and [../analysis/extensions.md](../analysis/extensions.md)),
   each panel designed for its task rather than ported. The `load*` methods and the URL parameters carry over; the rest
   of the JS API is designed anew (§3.2).
7. **The agent can drive the whole viewer, and the UI renders the result** (David, 2026-10-14): an agent alone can
   operate the viewer, so evals run it headless and the app can show what it did (§3.3).
8. **Commands and tools are two surfaces over one API** (David, 2026-10-15), as in Stories Untold: commands are what
   people use (forms, previews, keys), tools what the agent uses (flat inputs, progressively disclosed in skills), and
   both call the viewer's API, which is the only way state changes (§3.3). This replaces "tools are commands" in
   [../design/agent.md](../design/agent.md) §1–2.
9. **Targets are MolQL underneath** (David, 2026-10-15): permissive parsing on the way in (keywords, address lists,
   PyMOL and other query texts, schema items), MolQL JSON for everything stored and evaluated; Stories' `any-molql` is
   the start ([targets.md](targets.md)).
10. **Sessions and settings**: `.molx` and `.molj` open and save, with the viewer's own state (the sequence's) in an
    owner-keyed field of the plugin's snapshot ([sequence.md](sequence.md) §7); remote states are behind a setting; the plugin's
    settings that matter to people (v5 has many) are in the Settings tab.
11. **Deployment**: `apps/prism-viewer` is published on GitHub Pages (ESM with chunks). Beside it, a **single-file
    viewer**: one ESM file, CSS included (Prism's CSS is compiled to JS modules), that a page imports and runs with no
    build step; its lazy parts are inlined.
12. **Embedding**: `<prism-viewer>` with a shadow root, several on a page, and the embeddable parts (a viewport alone,
    for mvs-stories).
13. **Extensions are code in this repository**, not a public framework: all of them live in the Mol\* repo today, so
    the extension API is the viewer's plain types and helpers, with no stability promise or versioning beyond the
    repository's own; outside code uses the same types at its own risk.
14. **The agent's policies start from Stories Untold's**: what is preauthorized, what asks (network, uploads),
    untrusted text from files and documents, where keys are kept. Stored provider keys are always scrambled, in every scope: a
    per-browser key in IndexedDB when allowed, else an app key in the source ([agent.md](agent.md) §4).
15. **Mobile, touch and accessibility are required**: the viewer works on phones (compact mode), by keyboard and with a
    screen reader; colorblind-safe and high-contrast palettes, for the UI and for the data colors.
16. **Not ported now**: volumes and segmentations (its implementation is to be redesigned, not ported). The PAE plot is
    refreshed: Seq\* draws PAE as a dot plot already; a custom rendering may come as an extra.
17. **Order**: the plan is [plan.md](plan.md): the real layers built beside the prototype in the same packages, views
    switched to real models feature by feature, and the app becoming the viewer people get at the end.
18. **Hiding part of a structure, two ways** (David, 2026-10-15): at the entity and chain level (a contents row's eye:
    what is hidden is left out of the structure's representations) and as transparency overlays (a layer over any
    target, from opaque to hidden, keeping the geometry). Both are API operations over targets
    ([parity.md](parity.md) §8).
19. **Tests, to start**: Storybook (with a real plugin and fixtures for every data type) and headless integration tests
    over the API. Evals, smoke tests of the built app and the single-file build, and the runtime bench as a guard come
    as the code settles; much will change first.
20. **Theming, to start: tokens.** A host restyles the viewer through the `--msp-*` tokens (they cross the shadow root);
    named parts and component tokens wait until asked for.
21. **A session record** (David, 2026-10-15): the agent runs on an `AgentSession`, and the viewer keeps what happened
    (the person's runs, the agent's turns, undos) in a fiber, shown back in History; in memory by default, kept in
    IndexedDB only when the person asks, and then only as a record (no state, no undo after a reload; at most 20, each
    forgettable), and in a saved `.molx` when the person opts in; the viewer composes the AI packages' blocks, the
    library keeps no transcript ([journal.md](journal.md)).
22. **A reasonable rxjs surface everywhere** (David, 2026-10-15): every public object of the viewer exposes its state
    as `ReadonlyState`s and its happenings as non-replaying `Observable` events (§3.6), never a mutable field with a void
    `changed` event.
23. **Light state in the browser's small stores** (David, 2026-10-15): preferences in localStorage, and a per-tab light
    restore in sessionStorage (the scene without data, when it fits a small budget and its data comes from URLs), so a
    reload brings a light scene back; neither holds data ([journal.md](journal.md) §6).
24. **Selective about the plugin's managers** (David, 2026-10-15): the viewer uses the plugin's state trees, builders,
    registries and the managers that do engine work; it rewrites in its own code the ones that carry UI-era policy or
    side effects (the snapshot manager's list, toasts, layout), and calls `plugin.state` directly where a manager would
    add side effects (the history's checkpoints). The list is [state.md](state.md) §5.
25. **A dynamic toolkit** (David, 2026-10-15): tools and skills are added and removed while a session runs (an
    extension loaded, a feature removed, a registry changed); a removed tool's definitions are an eviction hint at the
    next compaction. A mechanism in `@molstar/ai-tools` and `ai-session`; the viewer decides what is registered (§4,
    [agent.md](agent.md) §3).
26. **Approvals that stay out of the way** (David, 2026-10-15): known data services are allowed; another address asks
    once with "Allow once", "Allow this domain" (for the session or in this browser) and "Deny"; uploads and clearing
    ask each time. Domain rules are shown and removable ([agent.md](agent.md) §4).
27. **Extensions register commands as the norm** (David, 2026-10-15): an extension brings operations with their commands
    and tools first, and UI only where a view is needed. Model-archive quality (pLDDT, QMEAN, PAE) and ANVIL are
    designed anew this way, not ported; v5 conveniences such as the "Automatic (w/ Annotation)" preset become commands
    (apply the annotations a structure has) rather than copies ([parity.md](parity.md) §6).
28. **No changes to the GL layer for now** (David, 2026-10-15): the WebGPU backend is being built (Alex), so Canvas3D's
    passes, render targets and canvas sizing are not touched; resizes improve at the shell and plugin level first, and the
    rest waits for that work (§3.5).
29. **Several viewers on a page** (David, 2026-10-15): keys go to the viewer that has focus (or the pointer), never to
    the window; stored state is namespaced per viewer; page-wide settings (debug and timing modes) are said to be
    page-wide; how many WebGL or WebGPU contexts a page can hold is the host's concern. One canvas split into
    viewports, explicit `Canvas3DContext` management and several state trees need a new `PluginContext`: Prism's step
    2 (§6).
30. **Operations that create things have an identity** (David, 2026-10-15): loading what is already loaded, adding a
    component that exists, says so and does nothing, unless asked for a copy (`copy: true`); where the intent is
    unclear, the agent asks the person: an `ask` tool with options and a fill-in ([operations.md](operations.md) §1,
    [agent.md](agent.md) §3).
31. **Accessibility and touch by best practice** (§3.7), with colorblind-safe and high-contrast palettes for the UI and
    the data.
32. **Stored data is versioned, collected and shown** (David, 2026-10-15): every stored value carries a version and is
    dropped when it does not match; kept records, preferences and keys are collected by age and size; a Storage view
    (beside the assistant's Trace, Tools and Limits, and without the assistant in Settings) shows what is kept where,
    with Forget ([journal.md](journal.md) §7).

## 3. Designs (accepted 2026-10-15)

The spec object these sections add up to, `PrismViewerSpec`, is laid out in [viewer-spec.md](viewer-spec.md).

### 3.1 The extension API

- **Defined by the viewer**, in `@molstar/prism-viewer/extension`; prism-ui takes no position. An extension's entry
  (`@molstar/geo-export-extension/prism-viewer`) depends on it and on prism-ui.
- **An extension is a plain object** (viewer-spec.md §3): `plugin` is the engine side (§3.4); `api` its operations,
  `commands`, `tools` and `skills` the surfaces over them (§3.3), the norm (decision 27); `ui` places views in the
  viewer's **named slots**, a closed list the viewer owns (`ViewerSlot`: Home's import actions, Export sections,
  Selection sections, canvas settings, animation, Settings groups, inspector sections), plus structure rows, canvas
  toolbar items, and tabs for large standalone views (each a tab on the rail). Cards and rail buttons of its own it
  cannot add.
- **Every part has an id** (`<owner>.<kind>.<name>`: `geo-export.section.export`). The composition overrides by id:
  `overrides: { 'geo-export.section.export': { title: 'Geometry' }, 'kinemage.tab.views': null }`. A duplicate id throws unless it is
  an override.
- **Lazy**: `lazyExtension({ id, stub, load: () => import('…/prism-viewer') })`. The stub is what shows before the load
  (command titles, a section's title, a toolbar item); first use loads the chunk as a task in the status column, swaps
  the real extension in and replays what triggered it.
- **Mol\*'s extension list lives in the app**, not in `@molstar/prism-viewer`: extensions depend on the viewer, so the
  viewer cannot depend on them (the workspace check rejects cycles). The app's list is the opinionated "full Mol\*
  viewer", also behind the single-file build's `FullPrismViewerSpec()` (app.md §2); the viewer's own features are
  built-in extensions in `DefaultPrismViewerSpec()` (viewer-spec.md §3.1). `@molstar/prism-viewer`'s dependency on
  `model-archive-extension` goes: model-archive quality becomes an extension with commands (decision 27).

### 3.2 The JS API and URL parameters

- **Loaders are plugin extension functions, imported as needed** (David, 2026-10-15): v5's `Viewer.load*` methods are
  already plain functions over the plugin (`@molstar/plugin-extension/loaders`: `loadPdb`, `loadAlphaFoldDb`,
  `loadEmdb`, `loadStructureFromUrl`, `loadMVSFromUrl`, `loadFiles`, …), so the viewer has no `load` facade that would
  pull them all in. A script imports what it uses and runs it in a history step:

  ```ts
  import { loadPdb } from '@molstar/plugin-extension/loaders';
  const viewer = await createPrismViewer({ target, spec });
  await viewer.run('Load 1TQN', (plugin) => loadPdb(plugin, '1tqn'));   // one step: origin `script`, undoable
  viewer.api;                                   // every operation, as commands and tools call them
  viewer.views;                                 // the root view models: read-only states (loaded, selection, focus)
  viewer.plugin;                                // the engine, an escape hatch
  await viewer.dispose();
  ```

- **The loaders feature** (`loaders`, a built-in extension, viewer-spec.md §3.1) wraps the same functions as API
  operations (`load.entry`, `load.url`, `load.mvs`, `load.snapshot`, …) with their commands, tools and URL parameters;
  an app that leaves it out ships none of them. A change made on the plugin outside the API (a script calling a loader
  without `run`) still shows, but it is no history step, and the steps before it can no longer be undone
  ([sequence.md](sequence.md) §7.3).

- **URL parameters**: the loaders feature reads the load parameters and runs them as API calls once the viewer is ready;
  the options read the rest at creation; the one list is [app.md](app.md) §4. v5's names keep working (`pdb`, `afdb`, `emdb`, `pdb-ihm`,
  `model-archive`, `url` + `url-format`, `structure-url`, `mvs-url`, `mvs-data`, `snapshot-url`, `snapshot-id`, and the
  renderer options `pixel-scale`, `pick-scale`, `transparency`, `illumination`, `power-preference`, …). Layout
  parameters (`hide-controls`, `collapse-left-panel`) map to the UI state. The viewer's options are one PD definition,
  read by the JS options, the URL and the Settings tab alike.

### 3.3 One API under commands and tools; state rendered by the UI

Stories Untold's layering ([../analysis/stories-untold.md](../analysis/stories-untold.md); in that repository
`designs/architecture/operations.md` and `agent.md`), carried over to the viewer:

```
views ──► view models ─┐
commands (people) ─────┤
tools in skills (agent)┼──► the API ──────────► state: the engines (Mol*, Seq*) · the viewer's · the UI's
JS API, URL ───────────┘   (origin, run log)        │
views ◄── view models ◄── sources ◄─────────────────┘
          (optimistic edits: managers.transactions)
```

- **The API is the only way state changes.** Each operation is typed, in the person's and the model's terms, addressed
  by targets (below), and **idempotent by desired state** ("show A:248-256 as sticks" twice changes nothing the second
  time, so a retried or resumed step is safe). It returns one result shape for every caller:
  `{ ok: true, did, changed, warnings?, next? } | { ok: false, code, message, hint? }`. It lives in the viewer's
  `api/` folder, which with `sources/` is the only code that imports `@molstar/plugin` (contexts.md §7).
- **Commands** (people): forms with typed fields, previews, keybindings, the command menu, over the API.
- **Tools** (the agent): flat inputs with few fields, normalized for small models, one result rendering; composed with
  `@molstar/ai-tools`' toolkit (`createToolkit`, [packages/ai/\_spec/tools.md](../../../packages/ai/_spec/tools.md)),
  as Stories does:
  - **The core**, offered from the first turn, always includes the two disclosure tools, `disclose_skills` and
    `disclose_tools`, and `list_skills` (a read, by tag). The disclosure tools' descriptions are the generated catalog,
    so the model knows every skill by its summary without loading it.
  - **Skills** (guidance plus tools; selections, volumes, export, analysis, code, …) load with `disclose_skills`, several
    in one call; a skill's extras, which its guidance names, and single tools load by name with `disclose_tools`, or when
    called. Every tool is declared from the start; where the model allows, a load is appended where it happens, so the
    cached prompt stays.
  - **Extensions add skills, not core tools**, so the core's size does not grow with them (Stories' core: 24 tools,
    about 6.8k tokens); a lazy extension's skill is listed by its stub's summary, and disclosing it loads the chunk.
  - **View tools** (focus, select, frame, open a panel) never change the structures' content; their operations are not
    history steps (`history: false`).
- **Coverage is a check, not a hope**: every API operation has a command, a tool, or both, or says why not (as Stories'
  `COMMAND_TOOLS`); a test enumerates them.
- **Targets**: one address language shared by the API, commands and tools, MolQL underneath
  ([targets.md](targets.md)).
- **Managers and views** (David's split): `PrismViewer` has `managers` (`commands`, `transactions`, `history`,
  `navigation`, `viewState`, `notifications`; generic, DOM-free, in prism-ui) and `views` (the root view models: `shell`, `agent`,
  and the viewer's `scene`, `selection`, `snapshots`, …), next to `plugin`, as `PluginContext` has its managers.
  `transactions` are **optimistic UI transactions** (any-ui's coordinator, view-models.md §5): an edit shows at once
  and is absorbed when the source reports it, or rolled back. `navigation` and `viewState` are the UI's state,
  serializable.
- **Undo is the viewer's, in one history manager** (David, 2026-10-15): every undoable step, the person's API call or
  the agent's whole turn, is one entry with a checkpoint taken before it (`plugin.state.getSnapshot` with chosen flags,
  never the snapshot manager, and the sequence's state); Undo restores the latest. Grouping, runs in flight, failures
  and changes outside the API: [sequence.md](sequence.md) §7.3.
  Mol\*'s state history is off (`HistoryCapacity` 0) and so is Seq\*'s (`history: false`); `canUndo` is not used
  ([sequence.md](sequence.md) §7).
- **Three layers of state**: the engines' (Mol\*'s state tree, camera, canvas, selection, focus, selection mode and
  granularity; Seq\*'s document and views), the viewer's own (what the sequence follows,
  the agent's turn) and the UI's (the open panel and section, the docks). What a view shows comes from read-only states fed by the sources, never from a call's return value.
  Ephemeral UI (hover, scroll, an unsaved form) is not state.
- **Headless**: `@molstar/prism-viewer/headless` creates the same viewer without a DOM, on `HeadlessPluginContext` when
  images are wanted: the same API, tools and state account (agent.md §3). As in Stories, a tool declares the API
  modules it uses (structures, view, canvas, files, …); where the host lacks one, the tool is still listed and says
  why it cannot run. prism-ui's `reactive/`, `commands/` and `forms/` stay free of the DOM, so no separate
  `@molstar/prism` core package is needed.
- **Evals**: a task runs an agent over the headless viewer; the eval checks the state account (and images where it
  matters); the run log (API calls with their arguments and origins) replays in the app to show what the agent did.
- **Different from Stories**: the viewer has no closed operation language of its own; its "document" is the plugin's
  state tree, which any transformer can change, plus the sequence's SeqViewSpec. The API's operations compile to
  state-tree updates and sequence patches, and undo restores snapshots, not inverse operations. Stories gets a
  schema-based document for its whole state; the viewer's agent is narrower.
- **The agent's first scope**: change the background and canvas, apply presets and styles, build MolQL and SeqQL
  selections, fetch structures (by id, by URL), focus and select, open panels. More skills follow the parity list.

### 3.4 The plugin spec

`createPrismViewer` builds a plain v6 `PluginSpec` ([../../designs/plugin-composition.md](../../designs/plugin-composition.md))
from four layers, in a fixed order, with nothing resolved by name; the code is in the viewer's `sources/`.

1. **The viewer's base**: the registry entries and behaviors its own features need (highlight, select, focus, camera,
   snapshot controls, the focus representation, the default custom properties; MolViewSpec with `MVSRuntimeRegistry`;
   `DefaultRegistry` for parity), exported as named parts so an app can swap or narrow them. Behaviors that only drive
   the classic UI are left out; nothing of `PluginUISpec` (no `layout`, `components`, `customParamEditors`).
2. **The extensions' `plugin` parts**, `{ registry?, behaviors?, config? }`, in the list's order. Data-only extensions
   (rcsb, pdbe, dnatco, anvil, sb-ncbr, wwpdb, g3d, backgrounds) need no `/prism-viewer` entry: the app wraps them,
   `pluginExtension('rcsb-validation', { behaviors: [RCSBValidationReport] })`; their presets and themes reach the UI
   through the generic commands and tools.
3. **The options**, as config items and canvas props.
4. **The app's last word**: `finalize: (spec) => spec`.

- **Order**: registry base → extensions → custom formats last (format order and `auto()` as in v5); config: later wins
  (`PluginConfigManager` sets in order), so options win over extension defaults.
- **Conflicts throw at `init()`**, atomically, naming every one; behaviors apply with `ref: transformer.id`, so the same
  behavior twice collides. An extension's plugin part goes with the extension: leaving the extension out of the list drops it.
- **Slim viewers**: with a narrower base, features whose providers are absent hide, and their commands and tools say
  why (plugin-composition.md §5.2).
- **Lazy extensions**: on load, `plugin.register(entry)` (its undo kept), its `config` items set on the running plugin
  (`plugin.config.set`, later wins as at creation), and the behaviors applied to the behavior tree as `initBehaviors()`
  does, custom properties first. Its formats register after the ones loaded at creation (a later `auto()` tie-break);
  a registration conflict fails the load and says why, leaving the stub. A stub declares every part it will bring
  (ids, kinds, slots, operations, tools) and what it owns (transformer ids, format names and file extensions,
  representation and theme names), so overrides are checked at creation and a session or file that needs it loads it
  first (viewer-spec.md §3, review F6, F45).
- **Headless**: the composition is plain data and runs in Node. An extension's `plugin`, `api`, `commands` and
  `skills` stay DOM-free; its views come through `ui: () => import(…)`, a chunk of their own.

### 3.5 The canvas: smooth resizes, few draws

Today the plugin resizes after a 50 ms debounce and a trailing 100 ms throttle (`context.ts`, `input.resize`), so the
canvas is stretched for up to ~150 ms while a window or a panel moves. The plugin already asks Canvas3D to resize in its
next frame (`requestResize`); the cost is sizing the canvas element (`resizeCanvas` / `setCanvasSize` in
`graphics/src/canvas3d/util.ts`), which clears it, so the frame after a size change is blank until the next draw
(review F20). Canvas3D and the GL layer are not touched for now (decision 28). The plan:

1. **Most layout changes are not resizes** (now, shell and plugin level). The canvas fills the viewer and the panels
   float over it (shell.md §2). A pinned panel, the sequence band or a card does not shrink the canvas: the camera's
   framing keeps the scene in the visible area instead, by moving the camera's target and position as the band's
   framing does (`paddedFocus`, `viewer/src/sources/framing.ts`; not `camera.viewOffset`, which is the jitter passes'
   channel, review F3). Only the viewer's own size (the window, fullscreen, the host's element) resizes the canvas.
   The shell's pinned panel and rail stop shrinking the canvas (shell.md §2 changes with this).
2. **Fewer waits** (now, plugin level): Prism's plugin replaces the debounce and throttle with one resize per animation
   frame while sizes change, and a last one when they settle; measured before and after, kept only if it does not add
   blank frames. As built (`sources/canvas-host.ts`): a ResizeObserver on the viewer's element sizes the canvas and
   draws it in the same callback (after layout, before paint), and an unchanged size is not set again; checked in the
   browser: at each width the canvas matched its host and its center was drawn, not cleared.
3. **Later, with the WebGPU work**: sizing the canvas inside the render frame (so the cleared canvas is drawn in the same
   frame), cheap frames while a resize is going on (no illumination or multi-sample accumulation), and render targets
   in size steps. Designed with Alex once the backend lands.

### 3.6 State and events: one rxjs surface

Where something can be watched, it is watched the same way (decision 22), following prism-ui's `reactive/` and what
Stories learned (analysis/stories-untold.md §2):

- **State is a `ReadonlyState<T>`**: an `Observable<T>` with a synchronous `value`, replaying its current value to a new
  subscriber, emitting only on change (equality-suppressed, shallow for records), never erroring (an error is a value).
  Writers stay private; readers get the readonly view.
- **Happenings are events**: non-replaying `Observable`s (`Subject` inside), for what happened rather than what is: a
  run finished, a turn ended, an extension loaded. An outcome is never kept as state (a replayed outcome shows in the
  wrong place).
- **Streams are made once and kept** on their owner (a `ReadonlyState` per field, an `events` object), so subscribing
  twice costs nothing and the identity is stable for views.
- **No other shapes**: no callbacks or listener lists in public APIs, no mutable field plus a void `changed`; where the
  plugin has one, the viewer's sources adapt it (state.md §4).

The surface, by owner:

| Owner                    | State                                                                                       | Events                                              |
| ------------------------ | ------------------------------------------------------------------------------------------- | --------------------------------------------------- |
| `PrismViewer`            | `ready`, `busy`, `empty` (nothing loaded)                                                   | `error`                                             |
| `api`                    | `running` (the runs in progress), `definitions` (the operations, as extensions come and go) | `started`, `finished` (each with its `ApiResult`)   |
| `managers.history`       | `steps` (undoable, newest first), `canUndo`, `undoBlocked` (why not), `running`, `checkpointing` (on or off for a large scene) | `undone`, `checkpointingOff`              |
| `managers.transactions`  | the projections over each view model's base                                                 | `rejected` (a projection rolled back)               |
| `managers.navigation`    | `tab`, `card`, `section`, `page`                                                            | —                                                   |
| `managers.viewState`     | `layout` (docks, sizes, pinned, compact)                                                    | —                                                   |
| `managers.notifications` | `current`, `tasks`, `chips`                                                                 | `posted`                                            |
| `views.*`                | each view model's fields (scene nodes as records, selection, focus, snapshots, …)           | per view model where it has happenings              |
| the sequence viewer      | `style`, `colors`, `structure`, `rail`, `edits`, `legend`                                                     | —                                                   |
| extensions               | `loaded` (ids), per lazy entry `status` (stub, loading, loaded, failed)                     | `activated`, `failed`                               |
| the session record       | `sessions` (recent), a window's `entries`                                                   | `appended`                                          |
| the agent                | `turn`, `view`, `settings`                                                                  | `turnEnded`                                         |

Seq\* already follows the same split (`SeqState<T>`, `SeqEvent<T>`, `Live<T>` for pointer-rate values), so the viewer
passes its states through where their shape fits. The AI packages do not use rxjs (they import nothing of it); the
assistant adapts their hooks and observers into this surface.

### 3.7 Accessibility and touch

Best practice, checked per view (decisions 15, 31):

- **WCAG 2.2 AA** for the UI: contrast from the tokens (and a high-contrast theme), visible focus, targets of at least
  24 px (44 px in compact mode, the canvas's own buttons 38 px beside the molecule), no meaning by color alone; colorblind-safe palettes for data colors (structure themes,
  the sequence, legends) as a choice in Settings.
- **Keyboard**: every view usable without a pointer (the command menu is a combobox, trees are ARIA trees with arrow
  keys, cards trap and restore focus); the canvas is reached through the scene tree and targets, so picking,
  measuring and focusing have a path that needs no sight and no pointer (select by target, measure by targets, focus
  by target).
- **Screen readers**: one polite live region for results (`ApiResult.did`), progress and the agent's finished replies
  (not every streamed token); the canvas has a text alternative (the scene's outline and the selection).
- **Touch**: what hover gives elsewhere is a long press (the hover label, the sequence's tooltip and cards); ordered
  picks (measurements, superposition) list their points and take them one tap at a time, or by target; gestures follow
  the trackball's and are listed in Help.
- **Reduced motion**: `prefers-reduced-motion` turns off the shell's transitions and the camera's animated framing
  (pinning, focus) in favor of cuts.

## 4. Changes outside Prism

What the spec asks of other packages, collected from the sections above; each is additive unless it says otherwise.

| Package                      | Change                                                                                             | From             |
| ---------------------------- | -------------------------------------------------------------------------------------------------- | ---------------- |
| `@molstar/plugin`            | `PluginState.registerCustomState(key, { get, set })`: owner-keyed state written into every snapshot and applied on every restore (the snapshot manager and sessions included) | sequence.md §7 |
| `@molstar/plugin`            | a behavior holding the latest registered set (registry changes as a state)                         | state.md §3      |
| `@molstar/plugin`, graphics  | every path that changes canvas props emits `settingsUpdated` (a bug where one does not)            | state.md §3      |
| `@molstar/plugin`            | Prism's plugin resizes once a frame instead of after a debounce (no GL change); the rest after the WebGPU work | §3.5 |
| `@molstar/core` (assets)     | a URL-keyed, ref-counted, read-through byte cache in `AssetManager`                                | sequence.md §2   |
| `@molstar/ai-tools`, `ai-session` | a dynamic toolkit: add and remove tools and skills between turns; removed definitions as eviction hints for compaction (a mechanism; the viewer decides) | decision 25 |
| `@molstar/ai-fibers`         | a listing of fibers by owner (kept records are listed from it), or the viewer keeps its own index | journal.md §1 |
| `@molstar/prism-ui`          | the style registry adopts sheets into shadow roots as well as documents; keyboard handling scoped to a viewer's root | app.md §3, decision 29 |
| `@molstar/plugin-extension`  | the loaders split from the React parts (the package declares `@molstar/plugin-ui`), so the viewer depends on no classic UI (review F43) | §3.2 |
| Seq\* (`@molstar/seq-*`)     | none for data (its `Loaders` port suffices); `history: false` on the viewer's document port        | sequence.md §2, §7 |
| Stories Untold → this repo   | `any-molql` ported into `packages/prism/viewer/src/targets/`; its agent policies as the start      | targets.md §4, §2.14 |

## 5. Files

| File                                 | Content                                                                                    | State   |
| ------------------------------------ | ------------------------------------------------------------------------------------------ | ------- |
| [viewer-spec.md](viewer-spec.md)     | `PrismViewerSpec`: fields, extensions and the viewer's own features, parts and overrides, creation | written |
| [targets.md](targets.md)             | the address language: keywords, addresses, query texts, schema items; MolQL underneath     | written |
| [nodes.md](nodes.md)                 | node selectors: cells of the state tree by CSS-like text, over Mol*'s `StateSelection`    | written |
| [sequence.md](sequence.md)           | the sequence view and the 3D view: shared data, sequence behaviors, the link, MVS and SVS, history | written |
| [state.md](state.md)                 | the plugin's state, area by area: its source, view model, story, what is missing           | written |
| [parity.md](parity.md)               | v5 → Prism, feature by feature: where it lives, its API area, its state                    | written |
| [operations.md](operations.md)       | the API by area: operations, targets, results, undo; the commands over them                | written |
| [agent.md](agent.md)                 | tools and skills, the core, the account, policies, headless                                | written |
| [library.md](library.md)             | prism-ui's public surface and what it leaves to apps                                       | written |
| [journal.md](journal.md)             | the session record: a fiber of runs and turns, History, the agent's session                | written |
| [review.md](review.md)               | the multi-agent review of 2026-10-15 and how each finding was folded in                    | written |
| [api.md](api.md)                     | the viewer's API: one object per area (the managers and builders done again), operations as their methods, the types around them, grouping, extensions, the code skill | accepted |
| [app.md](app.md)                     | `apps/prism-viewer` and the single-file viewer: builds, chunks, deployment, embedding, URL | written |
| [plan.md](plan.md)                   | the order of work: parallel by layer, switched per feature; tracks and seven phases        | accepted |

Folded in: the viewer's composition and the extension API are [viewer-spec.md](viewer-spec.md); the command set is in
[operations.md](operations.md); the JS API is §3.2 and [app.md](app.md).

## 6. Later: step 2

Out of this spec's scope, decided to come after it (David, 2026-10-15):

- **A new `PluginContext`**: explicit `Canvas3DContext` management, one canvas split into viewports (Canvas3D supports
  viewport splitting; v5 did it with several plugins sharing a context), several state trees in one viewer, and the
  plugin's managers that the viewer rewrites (decision 24) moved into it. Too many changes alongside this spec, and it
  meets the WebGPU work (decision 28).
- **The canvas's own resize path** (§3.5, step 3), with the WebGPU backend.
