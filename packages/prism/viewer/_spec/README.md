# The Prism viewer: the specification

`@molstar/prism-viewer`: Mol\*'s viewer on Prism, its new UI, over a real plugin. This spec is for working on or with
the code without reading all of it: what the package holds, how a call flows through it, the contracts a host or an
extension relies on, where the code is, and the decisions below. Code comments cite these files' sections (`api.md §2`,
`journal.md §6`, `README decision 30`). The code is the truth for how. The base it builds on:
[prism-ui's \_spec](../../ui/_spec/README.md); the assistant: [its \_spec](../../assistant/_spec/README.md); what is
missing: [../../TODO.md](../../TODO.md).

| File                             | What it covers                                                                                                |
| -------------------------------- | ------------------------------------------------------------------------------------------------------------- |
| [viewer-spec.md](viewer-spec.md) | `PrismViewerSpec`: fields, extensions, built-in features, parts and overrides, creation                       |
| [api.md](api.md)                 | the API: operations, areas, run path, origins, groups, extensions' areas, commands over operations, coverage  |
| [operations.md](operations.md)   | every operation of the default spec: id, title, flags, params (generated from the code; a test keeps it true) |
| [targets.md](targets.md)         | the address language: keywords, addresses, query texts; MolQL underneath                                      |
| [nodes.md](nodes.md)             | node selectors over the state tree                                                                            |
| [marks.md](marks.md)             | marks on the scene: highlights and what views draw over it                                                    |
| [state.md](state.md)             | the plugin's state area by area: sources, view models                                                         |
| [journal.md](journal.md)         | the session record, History, what is kept in the browser                                                      |
| [sequence.md](sequence.md)       | the sequence viewer (Seq\*) beside the 3D view, the link, the history's rules                                 |
| [agent.md](agent.md)             | the agent: tools, skills, the account, policies, headless and evals                                           |
| [parity.md](parity.md)           | v5 → Prism, feature by feature                                                                                |
| [app.md](app.md)                 | `apps/prism-viewer`, the single file, embedding, URL parameters                                               |

## Entries

| Module                                        | What it is                                                                                         |
| --------------------------------------------- | -------------------------------------------------------------------------------------------------- |
| `page.ts`: `createPrismViewer`                | the viewer in a page: a shadow root on the target, the shell, the canvas; one per `persist.key`    |
| `headless.ts`: `createPrismViewerHeadless`    | the same viewer with no DOM (Node, tests, evals): the API, sources, history, record, agent         |
| `default-spec.ts`: `DefaultPrismViewerSpec()` | the base: the plugin's named parts, the built-in features (eager and lazy), no assistant           |
| `spec.ts`                                     | the types: `PrismViewerSpec`, `ViewerExtension`, `LazyViewerExtension`, parts, slots, the agent's  |
| `api.ts`, `api/`                              | the API: areas of operations, the runner, results, origins                                         |
| `viewer.ts`                                   | `createViewerCore`: the plugin, the API host, the history, the record, the sources, the extensions |

The app (`apps/prism-viewer`) adds Mol\*'s extensions (`FullPrismViewerSpec()`), options from the URL, `persist` and the
assistant (`@molstar/prism-assistant/viewer`, its own chunk), and builds the single file `prism-viewer.js`.

## Folders

| Folder                 | Holds                                                                                                                                                                        | May touch the plugin |
| ---------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------- |
| `api/`                 | the core areas (`load`, `structure`, `component`, `representation`, `select`, `focus`, `camera`, `canvas`, `history`, …), the runner (`run.ts`), targets as params           | yes                  |
| `sources/`             | what the views and the agent read: the scene, the selection, the canvas, the pointer, the sequence's structure, checkpoints, the light restore, the palette, the canvas host | yes                  |
| `features/`            | built-in features as extensions: `operations.ts` and `source.ts` (may touch the plugin), `index.ts`, `ui.tsx` (its own chunk)                                                | operations, sources  |
| `extensions/`          | the extension host (activation all or nothing, lazy stubs and their check, sequence behaviors, skills) and the one override rule                                             | no                   |
| `commands/`            | the command menu's commands over operations (`operationCommand`)                                                                                                             | no                   |
| `model/`               | view models (DOM-free): Home, the scene tree, the inspector, presets, quick styles                                                                                           | no                   |
| `views/`, `ui/`        | the views (Prism DOM JSX) and the shell's composition (`ui/shell.tsx`: tabs, toolbar, cards, the assistant's mount)                                                          | no                   |
| `sequence/`            | the sequence viewer (Seq\*) over the structure the 3D view shows; its behaviors (`behaviors.ts`)                                                                             | no                   |
| `targets/`             | target strings (`"A:10-20"`, `"ligand"`, MolQL, PyMOL text) to loci and back                                                                                                 | no                   |
| `agent/`               | the agent: tools over operations, skills, the account, policies, `createViewerAgent`, `runEval`                                                                              | no                   |
| `journal/`             | the session's record (runs, turns, undos), keeping it in IndexedDB, History's model                                                                                          | no                   |
| `persist/`, `storage/` | the light stores (versioned, budgeted), the Storage view's model and sources                                                                                                 | no                   |
| `testing/`             | mocks for the stories (Storybook); not the viewer's code                                                                                                                     | —                    |

`src/_test/boundary.test.ts` checks the plugin rule; `src/_test/coverage.test.ts` that every operation has a command, a
tool, or a reason it is neither.

## How a call flows

1. **Someone calls an operation**: a command (the person), a tool (the agent, with its turn), the URL (`link`), a
   script: `api.run(id, values)` or `api.<area>.<op>(values)`, each through `api.as(origin)`.
2. **The runner** (`api/run.ts`) checks the operation's modules and `when`, merges the values over its params, and runs
   it as a **history step** (a checkpoint before; one step per call, or one per agent turn as a group) unless it is a
   read or `history: false`. A result that changed nothing (`unchanged`) leaves no step.
3. **The operation** changes the plugin (a state update, a manager call) and returns `ApiResult.ok(did, changed)` or
   `fail(code, message, hint)`. Failures are results, not throws.
4. **The run is finished**: `api.events.finished` carries it with its origin; the **record** writes it
   (`journal/record.ts`); the agent's notes take the person's runs; the shell saves a file an agent's run made.
5. **The sources** see the plugin change (one subscription each, coalesced per frame) and the views follow.

Undo restores a checkpoint (the data tree, canvas props, behaviors, custom state such as the sequence's and the
palette's; never the camera or the selection); Redo puts back the state Undo left until a new step.

## Contracts

- **Operations** (`api/area.ts`, `OperationDefinition`): `id` `<area>.<verb>`, `params` (PD), flags `read`,
  `history: false` (view state), `undoable: false` (cannot be restored: an export, an upload), `approval` (`network` |
  `upload` | `destructive`), `addresses(values)` (what a network call reaches; absent: the agent asks), `modules`,
  `when`, `surface` (UI-only or agent-only, with a reason).
- **Extensions** (`spec.ts`): `{ id, title, plugin?, api?, commands?, annotations?, sequence?, skills?, ui? }`; a lazy
  one is a stub (its parts, operations with their flags, commands, skills) and `load()`. Part ids are
  `<owner>.<kind>.<name>`; overrides replace or remove any part or command by id.
- **Sequence behaviors** (`sequence/behaviors.ts`): `attach({ shown, patch, offer })`: patches onto the document `shown`
  names, menus and legend blocks offered.
- **The agent's policy** (`agent/policies.ts`): reads and view changes go; changes while the turn can be undone go;
  `destructive`, `upload` and `undoable: false` ask; network by `addresses` (known services, the request's hosts, the
  domains allowed).
- **Stored values** (`persist/stores.ts`): an envelope `{ v, d }` with the app's version, a budget per store; the
  record's codec is `prism-viewer.record.v1`.

## Tests

Headless with `src/_test/harness.ts` (`headlessViewer({ files: [cbs()] })`), no network (stubbed `fetch`), Vitest in
`_test/` folders; DOM tests opt in with jsdom. Live agent runs are `*.live.ts` (`pnpm test:live`), never unit tests. The
packed single file is checked by `smoke/browser/prism`.

## Decisions

Numbered as code cites them (`README decision 30`); numbers of dropped decisions are not reused.

- **1. prism-ui is a library, not a framework**: runtime, components, shell, commands, ink, status, theme, params form,
  the assistant's ports; no modules, no contributions registry, no modes. Apps compose it their own way.
- **2. Opinionated apps on it**: the viewer (and later stories) owns its features in its own code.
- **3. No modes**: the viewer and stories are separate apps; a link between them is a navigation.
- **4. The public package builds the viewer; the app is thin**: `@molstar/prism-viewer` makes the viewer from a spec;
  `apps/prism-viewer` (private) is the page, URL options, its extension list and the build.
- **5. Extensions bring their own UI; the viewer overrides it by id**; heavy ones are lazy
  ([viewer-spec.md](viewer-spec.md) §3–4).
- **6. Parity with v5 in function, not form**: every v5 feature has a place, each panel designed for its task
  ([parity.md](parity.md)).
- **7. The agent can drive the whole viewer and the UI renders the result**: evals run it headless.
- **8. Commands and tools are two surfaces over one API**, the only way state changes ([api.md](api.md) §11).
- **9. Targets are MolQL underneath**: permissive text in, MolQL JSON stored and evaluated ([targets.md](targets.md)).
- **10. Sessions**: `.molx`/`.molj` open and save; the viewer's own state rides in owner-keyed custom state of the
  plugin's snapshot ([sequence.md](sequence.md) §7); remote states behind a setting.
- **11. Deployment**: the app (ESM with chunks) on GitHub Pages, and a single-file ESM viewer, CSS included, lazy parts
  inlined ([app.md](app.md)).
- **12. Embedding**: `<prism-viewer>` with a shadow root (host CSS stops there), several on a page (`src/element.ts`).
- **13. Extensions are code in this repository**: the extension API is the viewer's plain types, no stability promise
  beyond the repository's.
- **15. Mobile, touch and accessibility are required** ([accessibility.md](../../ui/_spec/accessibility.md)).
- **16. Not ported**: Volumes & Segmentations (to be redesigned).
- **18. Hiding part of a structure, two ways**: entity and chain level (`contents.visibility`, left out of the
  components) and transparency layers over any target (`representation.transparency`, geometry kept) ([api.md](api.md)
  §2).
- **19. Tests**: stories with a real plugin and fixtures for every data type; headless integration tests over the API; a
  smoke of the built single file.
- **20. Theming is tokens**: a host restyles through `--msp-*` (they cross the shadow root).
- **21. A session record**: runs, turns and undos kept in a thread, shown in History; in memory by default, kept only
  when the person asks, never with state ([journal.md](journal.md)).
- **22. One rxjs surface** (`src/viewer.ts`, `api/`, `sources/`): state is a `ReadonlyState<T>` (synchronous `value`,
  replays current, emits on change only, never errors; writers private); happenings are non-replaying `Observable`
  events (a run finished, a turn ended), never kept as state; streams made once and kept on their owner; no callbacks,
  listener lists or mutable field plus `changed` in public APIs (sources adapt the plugin's). Seq\*'s `SeqState` passes
  through; the AI packages use no rxjs, the assistant adapts them.
- **23. Light state in the browser's small stores**: preferences in localStorage, a per-tab light restore in
  sessionStorage; neither holds data ([journal.md](journal.md) §6).
- **24. Selective about the plugin's managers**: use the state trees, builders, registries and engine managers; rewrite
  those with UI-era policy (snapshot list, toasts, layout); call `plugin.state` directly where a manager adds side
  effects ([state.md](state.md) §5).
- **25. A dynamic toolkit**: tools and skills added and removed while a session runs (`@molstar/ai` `tools/`,
  `session/`); the viewer decides what is registered ([agent.md](agent.md) §3).
- **26. Approvals stay out of the way**: known data services allowed; another address asks (once, this domain, deny);
  uploads and clearing ask each time ([agent.md](agent.md) §4).
- **27. Extensions register commands first**: operations with commands and tools, UI only where a view is needed; v5
  conveniences become commands ([parity.md](parity.md) §6).
- **28. No changes to the GL layer** (WebGPU is being built): resizes improve at the shell and plugin level. Panels
  float over the canvas and do not resize it (framing keeps the scene in the visible area, `sources/framing.ts`); only
  the viewer's own size does. A ResizeObserver on the viewer's element sizes and draws the canvas in the same callback,
  and an unchanged size is not set again (`sources/canvas-host.ts`).
- **29. Several viewers on a page**: keys go to the viewer that has focus (`ownsKey`, never the window); stored state is
  namespaced by `persist.key`, two live viewers with one key are refused (`src/page.ts`).
- **30. Operations that create things have an identity**: what exists is said and not made again unless `copy: true`
  ([api.md](api.md) §1).
- **31. Colorblind-safe and high-contrast palettes** for the UI and the data
  ([accessibility.md](../../ui/_spec/accessibility.md), [viewer-spec.md](viewer-spec.md) §2.8).
- **32. Stored data is versioned, collected and shown**: a version on every stored value, collection by age and size, a
  Storage view with Forget ([journal.md](journal.md) §7).
