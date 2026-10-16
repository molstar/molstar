# The migration, in slices

Each slice leaves `main`-quality state: it builds, its tests pass, and the docs it moves are updated in the same commit
(`AGENTS.md`, "Keep the docs true"). Order: cut first (less to move), then collapse, then the renames and the providers
rewrite inside the one package. The owner's decisions (README, open questions) are all taken (2026-10-16); the last slice is a review gate.

## Checks

| Name      | Commands                                                                                                  |
| --------- | --------------------------------------------------------------------------------------------------------- |
| **ai**    | `pnpm test:unit packages/ai`, `npx tsc -b packages/ai/*` (after the collapse: `npx tsc -b packages/ai`)    |
| **prism** | `pnpm test:unit packages/prism`, `npx tsc -b packages/prism/viewer packages/prism/assistant`              |
| **repo**  | `npx tsc -b`, `pnpm lint`, `pnpm check:workspace`, `npx biome format` / `npx prettier --check` on the touched files only (other sessions edit this checkout: never `pnpm format` on everything) |
| **build** | `pnpm build:lib`, `pnpm build:apps` (the assistant's lazy provider chunks), the minimal viewer's bundle size (`examples/minimal-prism-viewer/README.md:107` records today's AI share) |
| **stories** | the assistant's stories still load (`packages/prism/assistant/src/_stories/assistant.stories.tsx` imports `scriptedModel`) |

Never `pnpm test:live` (real models, money) unless the owner asks.

## Slices

| #   | Slice                                                                                       | Touches                                                                                                     | Checks                       | After   |
| --- | ------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------- | ---------------------------- | ------- |
| 1   | **Cut the conformance suite to the essentials** **done** (`bce8b8a8c`) (conformance.md §2): delete dropped cases, P and S files, `meter.ts`, the IndexedDB run; keep `runHarnessConformance` for now on the 59, on today's API | `packages/ai/conformance/src/cases/**`, `components/src/harness/_test/conformance-indexeddb.test.ts`, `_spec/testing.md` | ai, repo | —       |
| 2   | **The basic in-memory harness** **done** (`ce4401c0d`) (harness.md; decided 2026-10-16): `createHarness` replaces the host, the runner port, the storage port and the supervisor; the harness's IndexedDB (persistence.md §2) goes with them; the session over it (§5); synchronous hashing with `@molstar/core`'s MurmurHash3 over bytes and the narrow Biome exception for it (harness.md §7) with the clock's digest tracking gone; blobs stay references; the 59 cases and the `createSession` fixture ported onto it; `tools/_test/turn.test.ts` and `trace/_test/helpers.ts` follow the fixture | `packages/ai/{harness,components,session,conformance}`, `tools/src/_test/`, `trace/src/_test/`; prism-viewer `agent/session.ts:9,29–30,121–122`; prism-assistant `viewer/inspector-sources.ts:97`; `_spec/{harness,session,testing,README}.md`, `.v6/ai/{README,harness}.md` | ai, prism, repo, stories | 1, 3 |
| 3   | **Delete dead code** **done** (`6362aaf7b`) (dead-code.md §1): `runOnce`, session prices, `callTool`, retention inspection, generators, small exports (narration and the fiber deltas stay; `./zod` is 12a) | `packages/ai/{session,harness,tools,conformance}`, `_spec/{session,harness}.md` | ai, repo | — |
| 4   | **Kept records out of `@molstar/ai`** **done** (`2824c7a57`) (decided: option A, persistence.md §3): `createIndexedDbRows` moves into prism-viewer; `@molstar/ai-components` is then empty and goes | `components/src/fibers/indexeddb.ts` → `packages/prism/viewer/src/journal/`; `record.ts:774`, `ui/history.tsx:92`, two tests; `.v6/prism-ui/spec/journal.md`; the components package's plumbing | ai, prism, repo | 2 |
| 5   | **Fiber cleanup pages out** **done** (`f2cd5fc51`) (threads.md §3; the deltas stay, decided): `nextFiberCleanupPage`, `cleanupPageSize` and the limits validated but unused | `fibers/src/{kernel,limits}.ts`, their tests, `conformance/src/fibers/generators.ts`, `_spec/fibers.md` | ai, prism | — |
| 6   | **Collapse into `@molstar/ai`** **done** (`0b5463ea9`; deviations in packaging.md "As done") (packaging.md): `git mv` into `src/<folder>/`, one `package.json` with the exports map, relative imports (196 lines), the conformance cases into `harness/_test/`, the 57 external import sites, plumbing (`pnpm-workspace.yaml`, `inventory.json`, root and Prism `tsconfig.json`s, Prism `package.json`s, `biome.json` folder rules, `opt-in-boundary.test.ts`), `pnpm install`; resolve the name collisions (§3) | everything under `packages/ai` (narration joins as `src/narration/`, `mediabunny` a dependency), 35 files in `packages/prism`, the plumbing, `apps/docs/.../assistant.md` | ai, prism, repo, build, stories | 1–5 |
| 7   | **Fibers → threads** **done** (`20ed4f3cb`) (threads.md §2)                                                         | `packages/ai/src/threads/`, `testing/`, prism-viewer `journal/`, `storage/sources.ts`, `ui/history.tsx`, `viewer.ts`; `_spec/fibers.md` → `threads.md`; `.v6/prism-ui/spec/{journal,README,plan}.md` | ai, prism, repo | 6 |
| 8   | **Providers** **done** (`d75fdb019`; deviations in providers.md "As done") (providers.md): new `src/providers/` behind the existing port, the tests cut (§4), then the renames (`Provider`, `createAnthropicProvider`, …) and the assistant's `settings-model.ts` | `packages/ai/src/providers/`, `harness/` (port name), `session/` (`adapterFor`), prism-assistant `viewer/settings-model.ts` and tests, prism-viewer `spec.ts`, `agent/_test/live-models.ts`, `_spec/models.md`, docs page | ai, prism, repo, build | 6 (native tool additions kept, decided) |
| 9   | **Duplicates and legacy paths** **done** (`c9da6ff06`; deviations in duplicates.md "As done") (duplicates.md, persistence.md §5): markers, estimates, one `Clock`, one limits source, pre-label fallbacks | `packages/ai/src/{trace,session,tools,threads}` | ai, prism | 6 |
| 10  | **The spec, rewritten for the simplified code** **done** (`3de7599e5`): `packages/ai/_spec/*` and `.v6/ai/*` (README packages table → one package and its folders; decisions 7–9; `harness.md` §1, §14, §18; `plan.md` "Where it comes from" gets a line on this simplification), `packages/ai/AGENTS.md`, `packages/ai/TODO.md` (its SQLite, data-migration and suite-waits entries no longer apply as written); this folder marked done | docs only                                                                                                   | `npx prettier --check` on them | 6–9 |
| 11  | **Leave the serial Vitest project** **done** (`870bafbde`; measured below) (conformance.md §4): `serial` and its comment removed from `vitest.config.ts`, the AI tests join `unit`; no file needed isolation. Before: the serial project ran after everything else, a 5.3 s tail (53 files, 515 tests; the whole run 71 and 78 s on a quiet machine). After: 9 whole-repository runs and 13 `packages/ai` runs under load (load average 30 to 70, one full run alongside), all green; the AI files are about 4.4 s of 505 s of file time, so no tail. Whole-run times were dominated by machine load (95 to 170 s; an old/new A/B on the same machine: old 121, 98, 144 s, new 89, 141, 164 s), so the saving is the 5 s tail, not visible in them | `vitest.config.ts`                                                                                          | several full `pnpm test:unit` runs | 1, 6 |
| 12a | **zod isolated in `packages/ai`** **done** (`70f84c488`) (tool-schemas.md §0, the owner's goal): the `./zod` entry and `zod.ts` gone, one Biome override restricting zod everywhere but `packages/ai/**`, the folder rules inside | `biome.json`, `packages/ai/src/tools/`, `_spec/tools.md`, `.v6/ai/tools-and-session.md` §1 | ai, prism, repo (`pnpm lint` over the whole repository) | 6 (or before, on today's packages) |
| 12b | **Tool schemas without zod** (tool-schemas.md §3; optional, **deferred** 2026-10-16): field shape → JSON Schema and the check in our code, an open enum for the disclosure tools; first a one-off test that today's schemas for every viewer tool equal the new generator's | `packages/ai/src/tools/` | ai, prism | 12a |
| 13  | **Review gate** (below): the whole branch reviewed and measured against this analysis before it lands | review notes and `results.md` in this folder; fixes go back to the slice they belong to | all of them, below | every other slice (12b only if it was done) |

Slices 1–5 each also update the docs they touch; slice 10 is the pass that makes the spec read as one design rather
than a list of edits.

## What can run in parallel

- **1, 3 and 5** touch disjoint files (cases / dead modules / the fiber kernel's cleanup pages) and can run at once,
  each in its own worktree; all three land before 6 starts.
- **2** is the one design-heavy slice: it waits for 1 (the cases it ports) and 3 (3 deletes `callTool`, the retention
  inspection, `runOnce` and prices from the same packages; or 3's harness and session parts are done inside 2). It
  touches two Prism files: coordinate with the sessions working there.
- **4** touches `packages/prism/viewer`, where other sessions work (the prism API cleanup): coordinate it with them, or
  do it inside 6.
- **6 is a single big-bang slice**: nothing else may touch `packages/ai` or the AI import sites in `packages/prism`
  while it runs, and it lands in one commit (a half-collapsed workspace does not build). It is mechanical: a scripted
  `git mv` plus import rewrites, then the hand-made parts (exports map, collisions, Biome rules).
- **7 and 8** both follow 6 and touch different folders (`threads/` and prism-viewer's journal; `providers/` and
  prism-assistant's settings): they can run in parallel. **9** touches `trace/`, `session/`, `tools/`: parallel with 7
  and 8 except where 9 changes the session's model limits (with 8).
- **12a** and **12b** follow 6 and touch only `src/tools/` and `biome.json`: parallel with 7–9.
- **10** last of the work; **11** whenever 1 and 6 are in; **13**, the review gate, after everything.

## Who implements which slice

Most slices are mechanical once the design files here say what moves; Sonnet can do those, each from a brief that
names the files, the lists (case names, symbols, rename maps) and the checks, with an Opus review of the diff before it
lands. The design-heavy ones stay with Opus.

| #   | Slice                         | Model                  | Why                                                                                                                                       |
| --- | ----------------------------- | ---------------------- | ----------------------------------------------------------------------------------------------------------------------------------------- |
| 1   | cut the conformance suite     | **Sonnet**             | deletion by a list (conformance.md §2 names every kept and dropped case); the checks say whether it worked                                 |
| 2   | the basic in-memory harness   | **Opus**               | a rewrite of the runner's loop and the session over it, with judgment on every case ported (harness.md §3–5); the one slice that designs  |
| 3   | dead code                     | **Sonnet**             | deletions listed with paths (dead-code.md §1)                                                                                            |
| 4   | kept records into prism-viewer | **Sonnet**             | one file moved, four import sites, two tests; coordination with the Prism sessions is the hard part, not the code                       |
| 5   | fiber cleanup pages out       | **Sonnet**             | two functions and their limits, located (threads.md §3)                                                                                  |
| 6   | collapse into `@molstar/ai`   | **Sonnet**, Opus for the exports map and folder rules | scripted `git mv` and import rewrites (packaging.md §2–4 list every site); the exports map, the name collisions and the Biome folder rules are decisions already written down but worth a second look |
| 7   | fibers → threads              | **Sonnet**             | a rename map with a stated scope (threads.md §2), including what must not be renamed (the molecular fibers)                              |
| 8   | providers                     | **Opus**               | a new layer: the `Message` → `ModelMessage` converter, stream mapping, reasoning replay, cache options per provider (providers.md §2–3)   |
| 9   | duplicates and legacy paths   | **Sonnet**             | small, each item located (duplicates.md, persistence.md §5)                                                                              |
| 10  | the spec rewritten            | **Opus**               | the design written as one coherent document in the repository's voice, not a list of edits                                               |
| 11  | leave the serial project      | **Sonnet** (or Haiku)  | a config change and repeated full runs                                                                                                   |
| 12a | zod isolated in `packages/ai` | **Sonnet** (or Haiku)  | an entry and a file removed, one lint override                                                                                           |
| 12b | tool schemas without zod      | **Sonnet**             | the behavior is pinned by `input.test.ts` and the one-off equivalence test (tool-schemas.md §3); the code is ~180 lines over a closed vocabulary |
| 13  | the review gate               | **Opus**, then the owner | judgment over the whole diff, and the owner's sign-off                                                                                  |

So 10 of 14 slices can go to Sonnet (the collapse with an Opus look at its exports and rules), and Opus keeps 2, 8, 10
and the review gate, 13.

## The review gate (slice 13)

Nothing of this simplification lands on `master` before it passes; a finding goes back to the slice it belongs to, and
the gate runs again on what changed.

1. **The whole repository green**, on a clean checkout of the branch: `pnpm install --frozen-lockfile`, `npx tsc -b`,
   `pnpm lint`, `pnpm format:check`, `pnpm check:workspace`, `pnpm test:unit`, `pnpm build`, the browser smoke
   (`node smoke/run.mjs browser`), and the CI workflow (`.github/workflows/node.yml`) green on the pushed branch.
2. **A code review of the branch diff** against its base: `/code-review` at high effort in a fresh Opus session (the
   owner may also run `/code-review ultra`, the billed multi-agent review, which only the owner can start). Every
   finding fixed or written down with a reason.
3. **The invariants this analysis set, checked by grep and recorded**:
   - no `@molstar/ai-` package name left in code or configuration (history in changelogs excepted);
   - one package: `packages/ai/package.json` with explicit exports, no `./*`, `"sideEffects": false`;
   - no `zod` import outside `packages/ai` (the Biome override fails the lint if there is one), none in the public types;
   - in `packages/ai`, no Mol\* import except `@molstar/core/data/util/hash-functions`, no `crypto.subtle`, no
     `node:*`, no rxjs;
   - no `fiber` identifier left in `packages/ai` or the viewer's journal (the molecular fibers elsewhere untouched);
   - no IndexedDB in `packages/ai`; the kept records' rows in prism-viewer;
   - every `harness.md §n` / `tools-and-session.md §n` citation in code points at a section that exists.
4. **The numbers against this analysis**, in `results.md` here: lines removed per category (README "The categories,
   sized"), the AI tests' time locally and in CI (today ~22.5 s and ~52 s, conformance.md §1), the conformance case
   count (59), the bundles (the minimal viewer's AI share, `examples/minimal-prism-viewer/README.md:107`; where zod and
   each provider land in the Prism viewer app's chunks). A large miss is explained, not just reported.
5. **The spec read against the code**: `packages/ai/_spec/*` and `.v6/ai/*` describe the simplified package, the
   decisions in the README are reflected, and `packages/ai/AGENTS.md` and `TODO.md` are current.
6. **The app**: the Prism viewer app starts, the Assistant runs a scripted turn (the stories, `scriptedModel`), an
   approval parks and resumes, Stop works, the Inspector shows the turn's trace, and a kept record survives a reload.
   One live turn on a real model only if the owner asks (it costs money).
7. **The owner's sign-off.**

## Later (not in these slices)

- **Persistence back**: no longer "add a backend": with the storage port, leases and the supervisor gone (slice 2), a
  durable or remote harness is a redesign. It can start from git (`packages/ai/{harness,components,conformance}` at
  `e91cf9bc0`): the host protocol, the durable core, and the dropped cases (conformance.md §2 "Dropped").
