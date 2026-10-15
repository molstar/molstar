# The Prism viewer: the specification

`@molstar/prism-viewer`: Mol\*'s viewer on Prism, its new UI, over a real plugin. This spec is for working on or with
the code without reading all of it: what the package holds, how a call flows through it, the contracts a host or an
extension relies on, and where the code is. Why the design is so is in
[.v6/prism-ui/spec](../../../../.v6/prism-ui/spec/README.md) (the decisions, the order of work, the reviews); code
comments cite its sections (`operations.md §2`, `journal.md §6`). The code is the truth for how.

| File                           | What it covers                                                                                                |
| ------------------------------ | ------------------------------------------------------------------------------------------------------------- |
| [operations.md](operations.md) | every operation of the default spec: id, title, flags, params (generated from the code; a test keeps it true) |

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
