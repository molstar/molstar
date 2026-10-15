# Plan: building the Prism viewer

The order of work for the spec in this folder. Each step is buildable and reviewable on its own (a PR, or a commit on
the branch for small ones), keeps `pnpm test`, `npx tsc -b` and Storybook green, and leaves the mocks working until
the real thing replaces them. Status: **accepted** (David, 2026-10-15); in progress on `v6-prism-impl-proto`.

## 1. The approach: parallel by layer, switched per feature

David proposed a parallel re-implementation of the prototype in `packages/prism/*`, switched over when ready. Agreed,
with the line drawn by layer rather than by package, because most of `packages/prism` is not prototype:

| Code (lines, 2026-10-15)                                                     | What it is                                    | In the plan                                                  |
| ---------------------------------------------------------------------------- | --------------------------------------------- | ------------------------------------------------------------ |
| prism-ui `dom/`, `reactive/`, `commands/`, `forms/`, `components/` (~8k)     | the library, real and tested                  | kept; evolved in place (managers, styles, keys)              |
| prism-ui `shell/`, `views/` (~4.3k)                                          | the shell and its views, over real models     | kept; `ShellModel` split into the managers (step 1.2)        |
| prism-viewer `views/` (~2k), `sequence/` and `sources/` (~1.9k)              | views over the `model/` contracts; the first real slice | kept; views move from mock models to real ones, feature by feature |
| prism-viewer and prism-ui `testing/` (~5k), `model/` contracts (~0.7k)       | the mocks and the contracts they implement    | kept as Storybook's data while the real models arrive; retired per feature |
| prism-assistant (~4.9k)                                                      | the assistant's views over mock models        | kept; wired to `createViewerAgent` in phase 4                |
| what does not exist: `api/`, the composition, features, `agent/`, the app    | —                                             | built new, beside the mocks                                  |

- **Parallel**: the real layers (`api/`, `sources/`, `features/`, `agent/`, `createPrismViewer`, `apps/prism-viewer`)
  are written next to the mocks in the same packages; nothing in the mocks is edited to make the real thing work.
- **Switched per feature**: a view moves from its mock model to the real one when that feature's area is built
  (its Storybook story gets a real-plugin variant, the mock variant stays until the feature is done), so the app grows
  one working feature at a time and there is no big-bang switch inside the packages.
- **The one switch at the end** is the app: `apps/prism-viewer` becomes the viewer people get (phase 7), and the
  classic `apps/viewer` stays beside it for a major version.
- **Why not a separate package or folder tree**: the views are the valuable part of the prototype and are already
  written against contracts; a second copy would drift from the first within weeks, and Storybook would have to show
  both.

## 2. Tracks that run beside the phases

These touch other packages, are small and additive, and unblock later phases; each is its own PR, delegated where
mechanical:

| Track                                                                                       | Needed by   |
| ------------------------------------------------------------------------------------------- | ----------- |
| A1. `PluginState.registerCustomState` (sequence.md §7.1), with tests on snapshots and sessions | phase 2     |
| A2. The registry-state behavior (state.md §3)                                               | phase 3     |
| A3. Canvas props paths that do not emit `settingsUpdated`: find and fix (state.md §3)       | phase 1     |
| A4. Loaders split from `@molstar/plugin-extension`'s React parts (README §4)                | phase 3     |
| A5. A URL-keyed, ref-counted, read-through byte cache in `AssetManager` (sequence.md §2)    | phase 2     |
| A6. A dynamic toolkit in `@molstar/ai-tools` and `ai-session`: add and remove tools and skills between and during turns, eviction hints (README decision 25) | phase 4 |
| A7. A fiber listing by owner in `@molstar/ai-fibers`, or the viewer's own index (journal.md §1) | phase 4     |
| A8. Port `any-molql` from Stories Untold into `packages/prism/viewer/src/targets/` (targets.md §4) | phase 2 |

No track touches Canvas3D or the GL layer (README decision 28).

## 3. Phases

### Phase 1: the skeleton, end to end

The smallest real viewer: created from a spec, a plugin behind it, one operation through the run path, shown in the
app and run headless.

1. **prism-viewer's spec types and composition**: `PrismViewerSpec`, `DefaultPrismViewerSpec()`, `composePluginSpec`
   (viewer-spec.md §2.1, README §3.4), `createPrismViewer` and `createPrismViewerHeadless` (viewer-spec.md §6), the
   creation order and `ready`.
2. **prism-ui managers**: `managers/` with `history` (steps, groups, the `HistoryPort`), `navigation` and `viewState`
   (split from `ShellModel`, which keeps registering parts), `notifications` (today's `StatusModel`); the style registry
   adopting into shadow roots; keyboard listeners on the viewer's root (library.md §3).
3. **The API's frame**: `ApiArea`, `OperationFn`, `api.run` (when, params, history step, transactions, the record
   hook), `api.as(origin)`, `api.group` (api.md §3–5); the sources helpers (state.md §4); the history's checkpoint port
   over `plugin.state` (sequence.md §7.2).
4. **First areas**: `canvas.background`, `camera.reset`, `load.files`, `structure.get`; the shell over the real managers
   with the canvas mounted; Home's drop target.
5. **`apps/prism-viewer`**: the page, the entry, an ESM build with splitting (`scripts/esbuild/app.mjs` gains the
   mode: `molstar.format: "esm"`), previewed as a static build on the developer's machine (David, 2026-10-15: no
   Pages deployment for now; `pnpm --filter @molstar/prism-viewer-app build`, served from `apps/prism-viewer/build`).
6. **Tests**: headless integration tests (create, load a fixture, an operation, undo, with no canvas); a real-plugin
   Storybook story of the shell.

Accepted when: the app opens, a dropped file loads and shows, the background changes and undoes, and the same runs
headless in a test.

### Phase 2: the scene, selection and the sequence

1. **Targets** (track A8 first): `TargetResolver`, the address and keyword forms, PyMOL text, schema items, "what it
   selected" (targets.md).
2. **Scene areas**: `structure` (presets, model, assembly), `contents` (with visibility, README decision 18),
   `component`, `representation` (with the selection layers), `styles`, `annotations`; the Structure tab's scene tree
   and inspectors over real view models (records per node, `VirtualList`).
3. **Selection and focus**: `select`, `focus`, `highlight`, the selection mode and granularity; the selection bar, the
   mode switch, the hover label; `events.hover` and `click`.
4. **The sequence on the API**: the link's writes through operations (origin `link`), `load.url` for its misses, sequence
   behaviors with stable refs, the edits overlay, the sequence's state in snapshots (track A1) (sequence.md §2–5, §7.1).
5. **Stories with data**: a real-plugin variant per data type (structures, models and trajectories, assemblies, maps
   and streamed density, particles, shapes) (README decision 19), particles in the scene tree.

Accepted when: everything in parity.md §2 and §3 marked core is built over a real plugin and covered by a headless test.

### Phase 3: features as built-in extensions

1. **The extension API**: `ViewerExtension`, stubs that declare every part and operation, lazy loading and replay,
   overrides by part id, one module per feature (viewer-spec.md §3–4).
2. **The features**, one PR each, the lighter first: `canvas-settings`, `snapshots` (with playback), `export-image`,
   `export-session`, `animation` and `trajectory`, `measurements`, `state-tree`, `loaders` (the URL table, app.md §4),
   `remote-states` (its client out of the classic view), `superposition` (its algorithms out of the classic view).
3. **Settings and History**: the Settings tab over the options and the plugin's settings; History over the record's
   memory rows (journal.md §3).

Accepted when: parity.md §1–§5 are built or dropped, every feature lazy or eager by the list alone, and the coverage
test passes (every operation has a command or a tool, or says why).

### Phase 4: the agent

1. **`createViewerAgent`** (agent.md §1), `paramsShape` moved into the viewer, the core tools and the first skills
   (selections, styles, canvas), `ask`; on the dynamic toolkit (track A6).
2. **Policies**: approvals per address with domain rules, the turn's undoability, quoted data (agent.md §4).
3. **The record**: the writer, entries, keeping a record in IndexedDB, the index (track A7), `.molx` opt-in
   (journal.md).
4. **The assistant on it**: prism-assistant's views over the real agent; key storage (agent.md §4); the Storage view
   (journal.md §7).
5. **Evals, a first one**: a task on the headless viewer, checked on the state account.

Accepted when: an agent can do the first scope (README §3.3) in the app and headless, a turn undoes as one step, and a
kept record survives a reload.

### Phase 5: Mol\*'s extensions

1. **Commands first** (README decision 27): model-archive quality (pLDDT, QMEAN; the pLDDT track and the PAE dot plot
   as sequence behaviors; the viewer's dependency removed), `annotations.apply` for what a structure carries.
2. **Ports**: assembly-symmetry, geo-export, model-export, video-export (lazy), zenodo (lazy), kinemage, ANVIL; the
   data-only ones through `pluginExtension`.
3. **The app's list and `FullPrismViewerSpec()`**; the single-file build (app.md §2).

Accepted when: parity.md §6 is built or dropped, and the single file runs on a page with no build step.

### Phase 6: parity pass and polish

Accessibility and touch (README §3.7), compact mode, colorblind-safe and high-contrast palettes, the light restore
(journal.md §6), versions and collection of stored data (journal.md §7), several viewers on a page (README decision 29),
plugin-level resize improvements (README §3.5, steps 1–2), the docs (`apps/docs`), performance budgets.

### Phase 7: the switch

`apps/prism-viewer` becomes the viewer people get, when: parity.md has no row left missing; the tests and the first
evals pass; the performance budgets hold; the classic `apps/viewer` stays published beside it for a major version.
Then the mocks that are left are retired, and the `_spec` folders next to the code are written (README, top).

## 4. Lanes and review gates

Two lanes by kind of work (David, 2026-10-15), and gates where contracts settle, so the cheaper lane fans out only on
patterns the expensive one has fixed.

### 4.1 Lanes

| Lane                         | Takes                                                                                                  |
| ---------------------------- | ------------------------------------------------------------------------------------------------------ |
| **Opus** (design and integration) | contracts and cross-cutting code: the composition, the API's frame and run path, the history and checkpoints, the managers split, the extension API and lazy loading, the sequence link, the agent's session and policies; the first of each pattern (the first area, the first feature, the first extension port); each panel's UX pass; integrating the Sonnet lane's PRs; the gate reviews |
| **Sonnet** (well-specified work) | work with a pattern to follow and a spec section to cite: the side tracks A1–A5, A7, A8; areas after the first (phase 2.2–2.3); features after the first (phase 3.2); extension ports after the first (phase 5.2); real-data stories and fixtures; headless tests per area; moves, renames and docs |

- **A Sonnet task's brief** names the spec section, the pattern to copy (a file that already does it), the files it may
  touch, and the tests that prove it; its verification is scoped to those tests.
- **Parallel Sonnet tasks run in worktrees on disjoint files**; two tasks never edit the same file. Opus merges and
  resolves.
- **Escalation**: a Sonnet task that finds the spec silent or wrong stops and reports instead of deciding; Opus (or
  David) decides and the spec is updated first.

### 4.2 Per phase

| Phase | Opus                                                                                 | Sonnet (after Opus's first)                                               |
| ----- | ------------------------------------------------------------------------------------ | ------------------------------------------------------------------------- |
| 1     | 1.1 composition, 1.2 managers split, 1.3 the API's frame, history, sources helpers   | 1.4 the first areas once 1.3 lands; 1.5 the app's build; 1.6 tests        |
| 2     | 2.1 the target resolver over the port, the first scene area, 2.4 the sequence on the API | the other scene areas, selection and focus areas, 2.5 stories per data type |
| 3     | 3.1 the extension API, the first feature (`canvas-settings`), 3.3 History             | the other features, one task each                                         |
| 4     | 4.1 `createViewerAgent` and the core tools, 4.2 policies                             | the skills after the first, 4.3 the record, the Storage view              |
| 5     | 5.1 model-archive (a new design), the first port                                     | the other ports, the data-only list, the single-file build                |
| 6     | accessibility and touch design, the resize work                                      | palettes, the light restore, stored-data versions, docs                   |
| 7     | the switch                                                                           | —                                                                         |

### 4.3 Review gates

A gate is a review of what landed since the last one, before the next fan-out; it blocks the next phase until its
blockers are fixed. Run as a small workflow: three or four reviewers by lens (correctness against the spec, contracts
and boundaries, tests and coverage, UX where views changed), then a dedup, then **batched** verification (two or three
verifiers over batches of findings, not one per finding), then a short report; findings are fixed in the lane that owns
the code.

| Gate | When                                         | What it freezes                                                                |
| ---- | -------------------------------------------- | ------------------------------------------------------------------------------ |
| G1   | after phase 1.3, before the first Sonnet areas | `ApiArea`, `OperationFn`, `api.run`, the `HistoryPort`, the managers, the sources helpers |
| G2   | end of phase 1                               | the composition and creation order, the app's build, the headless test harness |
| G3   | end of phase 2                               | targets, the scene areas, the sequence link                                    |
| G4   | after phase 3.1, before the features fan out | the extension API, stubs and lazy loading, part ids                            |
| G5   | end of phase 3                               | parity of the core and features; the coverage test                            |
| G6   | end of phase 4                               | the agent's session, tools, policies and record                                |
| G7   | end of phase 5                               | the extension ports, the single-file build                                     |
| G8   | before the switch (phase 7)                  | parity, performance budgets, accessibility                                     |

Between gates, every Sonnet PR gets one Opus read before it merges (not a workflow), and every PR updates the spec
where it moves it, and parity.md's state column.
