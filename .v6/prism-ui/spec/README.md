# The Prism viewer: the specification

How the Prism viewer is built: Mol\*'s own UI on Prism, at feature parity with the v5 Viewer, driven by commands that a
person, a script or an agent runs alike. This spec is what the build follows; the designs in [../design](../design) are
the reasons and the mocks behind it, and where they disagree, this spec holds. Once the code exists, a small `_spec/`
next to it (as `packages/ai/_spec`) describes what was built, and this folder keeps the decisions and the order of work.

Status: **draft** (2026-10-15). The decisions are in §2; the designs in §3 and the files beside this one were accepted
by David on 2026-10-15; §5 lists what is written and what is to come.

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
17. **Order**: the plan (`plan.md`) is written after this spec is approved.
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

## 3. Designs (accepted 2026-10-15)

The spec object these sections add up to, `PrismViewerSpec`, is laid out in [viewer-spec.md](viewer-spec.md).

### 3.1 The extension API

- **Defined by the viewer**, in `@molstar/prism-viewer/extension`; prism-ui takes no position. An extension's entry
  (`@molstar/geo-export-extension/prism-viewer`) depends on it and on prism-ui.
- **An extension is a plain object**: `{ id, plugin?, api?, commands?, skills?, ui?, settings? }`. `plugin` is the
  engine side (§3.4); `api` its operations, `commands` and `skills` the two surfaces over them (§3.3); `ui` places views in the viewer's **named slots**, a closed,
  documented list the viewer owns: Export sections, structure rows, inspector sections, Home's import actions, tabs,
  settings groups, canvas toolbar items. The viewer decides where things may go; an extension cannot add a rail button
  by itself.
- **Every part has an id** (`geo-export.export-section`). The composition overrides by id:
  `overrides: { 'geo-export.export-section': MyGeoExport, 'kinemage.views': null }`. A duplicate id throws unless it is
  an override.
- **Lazy**: `lazyExtension({ id, stub, load: () => import('…/prism-viewer') })`. The stub is what shows before the load
  (command titles, a section's title, a toolbar item); first use loads the chunk as a task in the status column, swaps
  the real extension in and replays what triggered it.
- **The default list lives in the app**, not in `@molstar/prism-viewer`: extensions depend on the viewer, so the viewer
  cannot depend on them (the workspace check rejects cycles). The app's list is the opinionated "full Mol\* viewer";
  `@molstar/prism-viewer`'s own dependency on `model-archive-extension` moves into that extension's entry.

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

- **The loaders feature** (`viewer.loaders`, a built-in extension, viewer-spec.md §3.1) wraps the same functions as API
  operations (`load.entry`, `load.url`, `load.mvs`, `load.snapshot`, …) with their commands, tools and URL parameters;
  an app that leaves it out ships none of them. A change made on the plugin outside the API (a script calling a loader
  without `run`) still shows, and it ends any held agent checkpoint, but it is no history step.

- **URL parameters** become API calls at start, handled by the loaders feature. v5's names keep working (`pdb`, `afdb`, `emdb`, `pdb-ihm`,
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
  - **View tools** change only the view (focus, select, open a panel), never the engine's state.
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
  the agent's whole turn, is one entry with a checkpoint taken before it (the plugin's snapshot through its snapshot
  manager, and the sequence's choices); Undo restores the latest.
  Mol\*'s state history is off (`HistoryCapacity` 0) and so is Seq\*'s (`history: false`); `canUndo` is not used
  ([sequence.md](sequence.md) §7).
- **Three layers of state**: the engines' (Mol\*'s state tree, camera, canvas, selection, focus, selection mode and
  granularity: what a plugin snapshot holds; Seq\*'s document and views), the viewer's own (what the sequence follows,
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
  behavior twice collides. The plugin part is overridden per extension only (replaced or dropped).
- **Slim viewers**: with a narrower base, features whose providers are absent hide, and their commands and tools say
  why (plugin-composition.md §5.2).
- **Lazy extensions**: on load, `plugin.register(entry)` (its undo kept) and the behaviors applied to the behavior tree
  as `initBehaviors()` does, custom properties first. A stub declares what it owns (transformer ids, format names and
  file extensions, representation and theme names), so opening a session or a file that needs it loads it first.
- **Headless**: the composition is plain data and runs in Node. An extension's `plugin`, `api`, `commands` and
  `skills` stay DOM-free; its views come through `ui: () => import(…)`, a chunk of their own.

### 3.5 The canvas: smooth resizes, few draws

Today the plugin resizes after a 50 ms debounce and a trailing 100 ms throttle (`context.ts`, `input.resize`), so the
canvas is stretched for up to ~150 ms while a window or a panel moves. The plan, in order:

1. **Most layout changes are not resizes.** The canvas fills the viewer and the panels float over it (shell.md §2). A
   pinned panel, a docked sequence band or a card does not shrink the canvas: the camera's projection is shifted
   instead (`camera.viewOffset`, as the band's framing does with `paddedFocus`), so the scene stays centered in the
   visible area, and pinning is a camera change animated in the render loop, with no reallocation. Only the viewer's
   own size (the window, fullscreen, the host's element) resizes the canvas.
2. **A resize happens in the frame.** Every resize event calls `canvas3d.requestResize()`, which resizes inside the next
   render (`render`, `resizeRequested`): at most one resize and one draw a frame, never a draw of its own. The plugin's
   debounce and throttle go (for Prism's plugin; the classic UI keeps them).
3. **A frame during a resize is cheap.** The cost is `passes.updateSize()`, which reallocates the draw, pick,
   multi-sample and illumination targets. While a resize is going on, frames skip the costly passes (illumination,
   multi-sample accumulation) and one full-quality draw follows when the size has settled (no change for a few frames).
4. **If that is not enough**: targets are allocated in size steps (64 px, say) and drawn into a sub-viewport of them
   (Canvas3D already draws a viewport inside its drawing buffer), so most frames of a drag reallocate nothing.

Steps 2–4 are Canvas3D and plugin changes, measured in the browser's performance panel before and after (a window drag
over a large structure with illumination on): frames per second, longest task, GPU memory churn.

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
| `api`                    | `running` (the runs in progress)                                                            | `started`, `finished` (each with its `ApiResult`)   |
| `managers.history`       | `steps` (undoable, newest first), `canUndo`, `checkpointing` (on or off for a large scene)  | `undone`                                            |
| `managers.transactions`  | the projections over each view model's base                                                 | `rejected` (a projection rolled back)               |
| `managers.navigation`    | `tab`, `card`, `section`, `page`                                                            | —                                                   |
| `managers.viewState`     | `layout` (docks, sizes, pinned, compact)                                                    | —                                                   |
| `managers.notifications` | `current`, `tasks`, `chips`                                                                 | `posted`                                            |
| `views.*`                | each view model's fields (scene nodes as records, selection, focus, snapshots, …)           | per view model where it has happenings              |
| the sequence viewer      | `style`, `theme`, `structure`, `legend`                                                     | —                                                   |
| extensions               | `loaded` (ids), per lazy entry `status` (stub, loading, loaded, failed)                     | `activated`, `failed`                               |
| the session record       | `sessions` (recent), a window's `entries`                                                   | `appended`                                          |
| the agent                | `turn`, `view`, `settings`                                                                  | `turnEnded`                                         |

Seq\* already follows the same split (`SeqState<T>`, `SeqEvent<T>`, `Live<T>` for pointer-rate values), so the viewer
passes its states through where their shape fits. The AI packages do not use rxjs (they import nothing of it); the
assistant adapts their hooks and observers into this surface.

## 4. Changes outside Prism

What the spec asks of other packages, collected from the sections above; each is additive unless it says otherwise.

| Package                      | Change                                                                                             | From             |
| ---------------------------- | -------------------------------------------------------------------------------------------------- | ---------------- |
| `@molstar/plugin`            | `PluginState.Snapshot.custom?: Record<string, unknown>`, kept by the snapshot manager and sessions | sequence.md §7   |
| `@molstar/plugin`            | a behavior holding the latest registered set (registry changes as a state)                         | state.md §3      |
| `@molstar/plugin`, graphics  | every path that changes canvas props emits `settingsUpdated` (a bug where one does not)            | state.md §3      |
| `@molstar/plugin`, graphics  | resizes through `requestResize`; cheap frames while resizing; targets in size steps if needed      | §3.5             |
| `@molstar/core` (assets)     | a URL-keyed, ref-counted, read-through byte cache in `AssetManager`                                | sequence.md §2   |
| AI packages                  | none: the viewer depends on `ai-session`, `ai-fibers`, `ai-components`, `ai-tools` and composes them (journal.md) | journal.md |
| `@molstar/plugin-extension`  | none: the loaders are used as they are                                                             | §3.2             |
| Seq\* (`@molstar/seq-*`)     | none for data (its `Loaders` port suffices); `history: false` on the viewer's document port        | sequence.md §2, §7 |
| Stories Untold → this repo   | `any-molql` ported into `packages/prism/viewer/src/targets/`; its agent policies as the start      | targets.md §4, §2.14 |

## 5. Files

| File                                 | Content                                                                                    | State   |
| ------------------------------------ | ------------------------------------------------------------------------------------------ | ------- |
| [viewer-spec.md](viewer-spec.md)     | `PrismViewerSpec`: fields, extensions and the viewer's own features, parts and overrides, creation | written |
| [targets.md](targets.md)             | the address language: keywords, addresses, query texts, schema items; MolQL underneath     | written |
| [sequence.md](sequence.md)           | the sequence view and the 3D view: shared data, sequence behaviors, the link, MVS and SVS, history | written |
| [state.md](state.md)                 | the plugin's state, area by area: its source, view model, story, what is missing           | written |
| [parity.md](parity.md)               | v5 → Prism, feature by feature: where it lives, its API area, its state                    | written |
| [operations.md](operations.md)       | the API by area: operations, targets, results, undo; the commands over them                | written |
| [agent.md](agent.md)                 | tools and skills, the core, the account, policies, headless                                | written |
| [library.md](library.md)             | prism-ui's public surface and what it leaves to apps                                       | written |
| [journal.md](journal.md)             | the session record: a fiber of runs and turns, History, the agent's session                | written |
| [app.md](app.md)                     | `apps/prism-viewer` and the single-file viewer: builds, chunks, deployment, embedding, URL | written |
| `plan.md`                            | the order of work, each step buildable and reviewable                                      | after approval |

Folded in: the viewer's composition and the extension API are [viewer-spec.md](viewer-spec.md); the command set is in
[operations.md](operations.md); the JS API is §3.2 and [app.md](app.md).
