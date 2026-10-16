# Prism: to do

What is missing, and what was considered and deferred. One entry each, with a pointer; an entry leaves when it is done
(its design or plan says so then). Gathered 2026-10-16 from the branch review
([2026-10-16-branch-review.md](../../.v6/prism-ui/reviews/2026-10-16-branch-review.md)), the specs
([.v6/prism-ui/spec](../../.v6/prism-ui/spec/README.md)) and the code. The AI package (`@molstar/ai`) has its own:
[../ai/TODO.md](../ai/TODO.md); Seq\*'s is [../sequences/todo.md](../sequences/todo.md).

## Decisions open

- **One schema vocabulary; reconsider zod** (David). The workspace describes and checks data three ways: Mol\*'s
  `ParamDefinition` (operations, options, Settings: `viewer/src/options.ts`, the API areas), io-ts in MolViewSpec
  (`mvs/builder/src/tree/generic/field-schema.ts`, the builder and the runtime) and zod in `@molstar/ai` (the agent's
  tools: field shapes from each operation's `ParamDefinition`, `viewer/src/agent/params.ts`, turned into zod, then into
  the JSON Schema the model sees and the validator). The question is the workspace's, not the agent's: one vocabulary
  that gives types, validation with good messages and JSON Schema (for the model, for MolViewSpec's published schema,
  for SeqViewSpec), or each kept where it is. zod is never in the page's startup bundle (the tool declarations are plain
  data, `viewer/src/agent/skills.ts`; the assistant is a chunk, `assistant/src/spec.ts`), but it is ~192 kB minified in
  the assistant's chunk (English only, as `scripts/esbuild/app.mjs` builds it; ~453 kB with every locale, 95 kB
  gzipped), and the AI SDK's providers import it, so it leaves the page only with providers of our own. Candidates: JSON
  Schema straight from `ParamDefinition` with a small validator over that closed vocabulary; zod everywhere (v4 emits
  JSON Schema; MolViewSpec off io-ts); a smaller validator. The harness side, slice 12b: [../ai/TODO.md](../ai/TODO.md).
- **The splash: the light or the light as a wave** (2026-10-16). The empty canvas runs the light (`prism-field.ts`, with
  its air, `fluid.ts`); the light as a wave (`quantum-field.ts`: Schrödinger by split-step FFT, Bohmian specks) is in
  Storybook only (Components / Splash: Quantum, Side by side). Pick one, or offer both
  (`ui.splash: 'light' | 'quantum'`).
- **Duplicates** (review B11): `StatusModel` and `NotificationsManager`; `ShellModel.layout` beside `viewState.save()`;
  `persist.key` that namespaces nothing; three ways to run Preset and Styles in `StructureTab`.
- **The query languages' registry is page-wide** (review B12): two viewers on a page offer each other's languages.
  Recommended: page-wide until pages with several viewers matter (plan.md phase 6).
- **Size** (review C14–C18; C13 was a no-op: esbuild minifies identifiers under `minify: true`): MolViewSpec in the
  startup bundle (~197 kB, a lazy MVS feature with stub format providers); the Seq\* chunk (one default for now,
  everything but streaming; its floor of ~721 kB is aligners, SeqQL, streams and palettes:
  [../sequences/todo.md](../sequences/todo.md)); eager extensions in the app (~190 kB, lazy stubs owning their formats);
  panels shown on demand as static imports (command menu, History, Settings with the params form, volume inspectors);
  `sideEffects` declared almost nowhere.
- **Build budgets** (2026-10-16, app.md §1): the `index.js` budget no longer bounds the start (the start code is in
  shared chunks since `lib`'s lazy import): a budget over the start closure (65 files, 4.4 MB, 1.29 MB gzipped); lower
  the 8 MB `prism-viewer.js` budget (6.2 MB today).
- **`molstar.lib`** (`apps/prism-viewer/src/lib.ts`): also from the single file `prism-viewer.js` (+32 kB gzipped on
  every page that imports it)? `StateTransforms` for full v5 parity (architecture.md §9.3 keeps it for the classic
  Viewer)?
- **The night's choices, keep or change** (review D): built-in areas typed as present and throwing when left out; panels
  keeping their tab ids; History and Settings in the core; the runner re-checking `exists` only for loads; the
  kept-records index growing (journal.md §1).

## Missing: the viewer

- A fast path in the runner for pointer-rate view operations; measure first (review, "Not done").
- `project()` and `selectionItems()` in `viewer/src/targets/evaluate.ts` read atomic-only properties (unused; B9).
- The drop overlay and Home read `sources.formats` once, when the shell mounts (review, "Not done").
- The PAE row's size before the matrix is read assumes a full ModelCIF PAE (review, "Not done").
- Node selectors: `:has(...)`, positions (`:first`, `:nth(n)`), kinds for shapes and particle lists; `node` on the
  structure-scoped operations ([nodes.md §5](../../.v6/prism-ui/spec/nodes.md)).
- Themes and presets that fetch as they apply carry no addresses (only known services today); the pLDDT track reads
  B-factors where the 3D theme reads ModelCIF scores (plan.md, G7 "Left").
- The record: a run's targets and `changed`, the turn's model, reply deltas, who undid a step, one entry per Undo turn,
  the index's size and folding, the record inside a `.molx` ([journal.md §1–2](../../.v6/prism-ui/spec/journal.md)).
- The sequence's caches and its edits overlay ([sequence.md](../../.v6/prism-ui/spec/sequence.md) §2, §7).
- No redo in the history (`ui/src/managers/history.ts`).
- The app: hashed entry names, the single file published to Pages beside the app (`pnpm deploy:remote` publishes the app
  at `prism/` and Storybook at `ui-storybook/`), a CI smoke of it ([app.md §1–2](../../.v6/prism-ui/spec/app.md)).
- Phase 6 and the switch: accessibility and touch, compact mode, several viewers on a page, performance budgets, and
  parity.md's missing rows ([plan.md](../../.v6/prism-ui/spec/plan.md) phases 6–7).
- Stories, its own app on prism-ui (there are no modes; `ModeChooser` is gone), is not built.

## Missing: the agent and the assistant

- Skills: SeqQL selections and the sequence skill, `viewer-code` (TypeScript against the API, api.md §7), opening panels
  as tools ([agent.md §3](../../.v6/prism-ui/spec/agent.md)).
- `ask` parking the turn with the question in the turn card; today it ends the turn
  ([agent.md §3](../../.v6/prism-ui/spec/agent.md)).
- Policies: steps made unsafe by a change outside the API when undo is off; the domain rules listed in Settings; data
  quoted in tool results ("Loaded 1CBS.") as it is in the state block ([agent.md §4](../../.v6/prism-ui/spec/agent.md)).
- What the agent sees: colors, snapshots and the canvas in the state block; a state version for repeated reads
  ([agent.md §5](../../.v6/prism-ui/spec/agent.md)).
- Without a page: the eval's run log replayed in the app; an MCP server over the headless viewer
  ([agent.md §6](../../.v6/prism-ui/spec/agent.md)).
- The setup: models asked from each provider's `/models`, an own endpoint's models listed
  (`assistant/src/viewer/settings-model.ts`).
- The command menu's assist: every request is a turn; reading words into one command (`assistant/src/viewer/ui.tsx`).

## Considered and deferred

- `agent.skills` as a function of the defaults: `{ add, remove }` was chosen (data, composes by id, needs no import of
  the defaults; [viewer-spec.md §2.4](../../.v6/prism-ui/spec/viewer-spec.md)).
- Keys in `assistant(options)`: a page never passes a key; keys are the person's, typed into the Assistant tab
  ([agent.md §4](../../.v6/prism-ui/spec/agent.md)). Dev keys from `.env.local` serve only the live test runs
  (`pnpm test:live`); a host with its own models passes `provided` (a server-side adapter).
- Field edges in the default themes: kept as they are (David, 2026-10-16; plan.md, G8).
