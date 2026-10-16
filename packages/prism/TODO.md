# Prism: to do

What is missing, and what was considered and deferred. One entry each, with a pointer; an entry leaves when it is done
(its design or plan says so then). Gathered 2026-10-16 from the branch review
([2026-10-16-branch-review.md](../../.v6/prism-ui/reviews/2026-10-16-branch-review.md)), the specs
([.v6/prism-ui/spec](../../.v6/prism-ui/spec/README.md)) and the code. The AI packages have their own:
[../ai/TODO.md](../ai/TODO.md).

## Decisions open

- **Reconsider zod** (David). The agent's tools are field shapes (`paramsShape`, `viewer/src/agent/params.ts`, from each
  operation's `ParamDefinition`) that `@molstar/ai-tools` turns into zod, then into the JSON Schema the model sees and
  the validator. zod is never in the page's startup bundle (the declarations are plain data,
  `viewer/src/agent/skills.ts`; the assistant is a chunk, `assistant/src/spec.ts`), but it is ~192 kB minified in the
  assistant's chunk (English only, as `scripts/esbuild/app.mjs` builds it; ~453 kB with every locale, 95 kB gzipped).
  Alternatives: JSON Schema straight from `ParamDefinition` (operations already have it) or from the field shapes, with
  a small validator over that closed vocabulary; a smaller validator. The AI SDK's providers import zod too, so it
  leaves the page only with adapters of our own. The harness side: [../ai/TODO.md](../ai/TODO.md).
- **What is public** (review B9): the package root as the embedder's entry, a curated export list instead of `./*`,
  `viewer.sources` read-only, the unused ported targets code removed (`viewer/src/targets/index.ts`).
- **Duplicates** (review B11): `StatusModel` and `NotificationsManager`; `ShellModel.layout` beside `viewState.save()`;
  `persist.key` that namespaces nothing; three ways to run Preset and Styles in `StructureTab`.
- **The query languages' registry is page-wide** (review B12): two viewers on a page offer each other's languages.
  Recommended: page-wide until pages with several viewers matter (plan.md phase 6).
- **Size** (review C13–C18): identifier minification off (~280 kB gzipped on the app); MolViewSpec in the startup bundle
  (~197 kB, a lazy MVS feature with stub format providers); the Seq\* chunk built from Seq\*'s whole default spec (a
  Prism composition); eager extensions in the app (~190 kB, lazy stubs owning their formats); panels shown on demand as
  static imports (command menu, History, Settings with the params form, volume inspectors); `sideEffects` declared
  almost nowhere.
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
- The app: hashed entry names, publishing the single file to Pages, a CI smoke of it
  ([app.md §1–2](../../.v6/prism-ui/spec/app.md)).
- Phase 6 and the switch: accessibility and touch, compact mode, several viewers on a page, performance budgets, and
  parity.md's missing rows ([plan.md](../../.v6/prism-ui/spec/plan.md) phases 6–7).
- The story mode (`stories`) is not built.

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
