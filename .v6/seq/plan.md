# Plan: the sequence viewer (Seq\*) into Mol\* v6, then Prism over it

Supersedes the order of work in [../prism-ui/handoff.md](../prism-ui/handoff.md) (its "How the user works" and
"Practical notes" still hold). Background:
[../prism-ui/analysis/harness-and-seqstar.md](../prism-ui/analysis/harness-and-seqstar.md).

## 1. The order

1. **Port Seq\*** from `stories-untold/lib/seqstar` (§2): the layout stays, the review and the fixes come after. The
   move is one way: from then on Seq\* changes here.
2. **Review, fix, merge** this branch. With Seq\* it is 300k+ lines, so it merges before the Prism viewer is written.
3. **Prism over Seq\*** (§3), on a branch of its own.
4. **Then the agent in Prism**: the handoff's former step 2 (command metadata, the selection param, the core viewer
   commands, the checkpointer, `AgentModel` over `AgentSession`, History over fibers), decided with the user.

## 2. The port

Seq\* today (2026-10-13): `core` (71k lines, no dependencies), `viewer` (113k: canvas, DOM and React UIs, rxjs),
`nightingale` (5k, the Nightingale elements), `molstar-poc` (5k, Mol\* v5), `scripts` (the spec and export checks),
`storybook` (13k: all 31 stories and the docs pages, React), `_spec` (the spec, part generated), `todo.md`. Tests run on
Vitest: a node tier, a Chromium tier (Playwright), heap guards and benchmarks.

| Step                   | What                                                                                                                                                                                                                                                                                                  |
| ---------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Layout                 | `packages/seq/{core,viewer,nightingale,molstar,scripts,storybook}`, with `_spec/` and `todo.md` beside them; Seq\*'s `AGENTS.md` adapted (the repository's rules, not Stories' shared-tree ones)                                                                                                      |
| The repository's rules | relative imports to `.js` (Seq\* writes `.ts`), exports in the `molstar-src` / `types` / `import` shape with a `lib` build, the repository's Biome formatting, the inventory, tsconfig references, `noRestrictedImports` per package (from Seq\*'s `layers.grit` and `core/_test/layers.test.ts`)     |
| msgpack                | Seq\*'s copy (`core/src/model/msgpack`) is Mol\*'s plus MessagePack extension types (`STORIES-UNTOLD(M2)`, encode and decode): those go into `@molstar/io`'s msgpack (an addition, in the changelog), and `seq-core` imports `@molstar/io`, its one dependency (the spec's "no dependencies" says so) |
| `molstar-poc`          | becomes `@molstar/seq-molstar`, its 43 v5 imports (`molstar/lib/...`) moved to the v6 packages (`model`, `io`, `query-language`, `plugin`, `mvs/runtime`)                                                                                                                                             |
| Tests                  | stay on Vitest, under their own script (`pnpm test:seq`, as `jest:ai`); the browser tier and the heap guards as scripts of their own, not in `pnpm test`                                                                                                                                              |
| Storybook              | `packages/seq/storybook` registered in `storybook/packages.ts` under a `Seq*` group; its React stories mounted by a decorator of that entry (the Storybook is `web-components-vite`); its paths into other packages' `_test` become exports (`./testing`)                                             |
| Dependencies           | to the catalog, pinned exactly: the Nightingale elements, Vitest and its browser provider, Storybook's React pieces                                                                                                                                                                                   |
| Docs                   | the root `AGENTS.md` table, `.v6/plans/checklist.md`, the changelog                                                                                                                                                                                                                                   |

Mechanical steps go to Sonnet subagents at low or medium effort, each checked against the tests its change touches.

## 3. Prism over Seq\*

- **The context**: Seq\*'s view context on Prism's context, at least in the viewer mode (the story mode may differ), per
  [contexts.md](../prism-ui/design/contexts.md).
- **The UI**: written in prism-ui's own framework, not Seq\*'s React UI. Where it lives (prism-ui, or the viewer mode's
  package) is decided then, under prism-ui's opt-in boundary.
- **The mocks**: the sequence strip replaced by the real viewer. Start with a plain text view or a single track with a
  toggle. No viewer frame: its theme follows the canvas background.

## 4. Follow-ups

An item leaves this list in the commit that does it. Each says what triggers it. Seq\*'s own deferred features stay in
`packages/seq/todo.md`.

| Item                                                                                                                           | Trigger                                      |
| ------------------------------------------------------------------------------------------------------------------------------ | -------------------------------------------- |
| One test runner: Seq\*'s Vitest and the AI packages' ESM Jest folded together (the checklist's "Jest as ESM")                  | the Jest-to-ESM work                         |
| Seq\*'s story support vendors Prism's tokens (`support/seqstar-prism`): use prism-ui's own                                     | the port's review                            |
| Seq\*'s React UI (`viewer/react`): kept for React hosts, or dropped once Prism has its own                                     | Prism's sequence UI built                    |
| The subscription guard (`prism/ui/src/_test/subscriptions.test.ts`) covers prism-ui only; the other Prism packages want it too | a leak found there, or a shared check script |
| The model catalog's prices are secondhand (ported from Stories, with sources and dates)                                        | credits shown to people                      |
| The Inspector's mock Trace link always lands on the latest turn                                                                | a session feeding the Inspector              |
| Whether browsers encode AAC and Opus at Gemini's 24 kHz (a runtime check falls back to WAV)                                    | narration in the story mode                  |
