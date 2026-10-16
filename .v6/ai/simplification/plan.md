# The migration, in slices

Each slice leaves `main`-quality state: it builds, its tests pass, and the docs it moves are updated in the same commit
(`AGENTS.md`, "Keep the docs true"). Order: cut first (less to move), then collapse, then the renames and the providers
rewrite inside the one package. The owner's decisions (README, open questions) gate the slices marked **(decision)**.

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
| 1   | **Cut the conformance suite to the essentials** (conformance.md §2): delete dropped cases, P and S files, `meter.ts`, the IndexedDB run; keep `runHarnessConformance` for now on the 59, on today's API | `packages/ai/conformance/src/cases/**`, `components/src/harness/_test/conformance-indexeddb.test.ts`, `_spec/testing.md` | ai, repo | —       |
| 2   | **The basic in-memory harness** (harness.md; decided 2026-10-16): `createHarness` replaces the host, the runner port, the storage port and the supervisor; the harness's IndexedDB (persistence.md §2) goes with them; the session over it (§5); synchronous hashing with the repository's MurmurHash3 over bytes (harness.md §7) with the clock's digest tracking gone; blobs stay references; the 59 cases and the `createSession` fixture ported onto it; `tools/_test/turn.test.ts` and `trace/_test/helpers.ts` follow the fixture | `packages/ai/{harness,components,session,conformance}`, `tools/src/_test/`, `trace/src/_test/`; prism-viewer `agent/session.ts:9,29–30,121–122`; prism-assistant `viewer/inspector-sources.ts:97`; `_spec/{harness,session,testing,README}.md`, `.v6/ai/{README,harness}.md` | ai, prism, repo, stories | 1, 3 |
| 3   | **Delete dead code** (dead-code.md §1): narration **(decision)**, `runOnce` **(decision)**, session prices, `callTool`, retention inspection, generators, small exports, `./zod` entry **(decision)** | `packages/ai/{narration,session,harness,tools,conformance}`; for narration also `inventory.json`, `tsconfig.json`, `biome.json` (`mediabunny` stays in the catalog), `prism/ui/src/_test/opt-in-boundary.test.ts`, `_spec/narration.md` | ai, repo | — |
| 4   | **Kept records out of `@molstar/ai`** **(decision, option A)** (persistence.md §3): `createIndexedDbRows` moves into prism-viewer; `@molstar/ai-components` is then empty and goes | `components/src/fibers/indexeddb.ts` → `packages/prism/viewer/src/journal/`; `record.ts:774`, `ui/history.tsx:92`, two tests; `.v6/prism-ui/spec/journal.md`; the components package's plumbing | ai, prism, repo | 2 |
| 5   | **Fiber deltas and cleanup out** **(decision)** (threads.md §3)                              | `fibers/src/{kernel,host,domain,limits}.ts`, their tests, `conformance/src/fibers/`, `_spec/fibers.md`      | ai, prism                     | —       |
| 6   | **Collapse into `@molstar/ai`** (packaging.md): `git mv` into `src/<folder>/`, one `package.json` with the exports map, relative imports (196 lines), the conformance cases into `harness/_test/`, the 57 external import sites, plumbing (`pnpm-workspace.yaml`, `inventory.json`, root and Prism `tsconfig.json`s, Prism `package.json`s, `biome.json` folder rules, `opt-in-boundary.test.ts`), `pnpm install`; resolve the name collisions (§3) | everything under `packages/ai`, 35 files in `packages/prism`, the plumbing, `apps/docs/.../assistant.md` | ai, prism, repo, build, stories | 1–5 |
| 7   | **Fibers → threads** (threads.md §2)                                                         | `packages/ai/src/threads/`, `testing/`, prism-viewer `journal/`, `storage/sources.ts`, `ui/history.tsx`, `viewer.ts`; `_spec/fibers.md` → `threads.md`; `.v6/prism-ui/spec/{journal,README,plan}.md` | ai, prism, repo | 6 |
| 8   | **Providers** (providers.md): new `src/providers/` behind the existing port, the tests cut (§4), then the renames (`Provider`, `createAnthropicProvider`, …) and the assistant's `settings-model.ts` | `packages/ai/src/providers/`, `harness/` (port name), `session/` (`adapterFor`), prism-assistant `viewer/settings-model.ts` and tests, prism-viewer `spec.ts`, `agent/_test/live-models.ts`, `_spec/models.md`, docs page | ai, prism, repo, build | 6; the tool-addition decision |
| 9   | **Duplicates and legacy paths** (duplicates.md, persistence.md §5): markers, estimates, one `Clock`, one limits source, pre-label fallbacks | `packages/ai/src/{trace,session,tools,threads}` | ai, prism | 6 |
| 10  | **The spec, rewritten for the simplified code**: `packages/ai/_spec/*` and `.v6/ai/*` (README packages table → one package and its folders; decisions 7–9; `harness.md` §1, §14, §18; `plan.md` "Where it comes from" gets a line on this simplification), `packages/ai/AGENTS.md`, `packages/ai/TODO.md` (its SQLite, data-migration and suite-waits entries no longer apply as written); this folder marked done | docs only                                                                                                   | `npx prettier --check` on them | 6–9 |
| 11  | **Leave the serial Vitest project** if the AI tests pass in parallel under load (conformance.md §4) | `vitest.config.ts`                                                                                          | several full `pnpm test:unit` runs | 1, 6 |

Slices 1–5 each also update the docs they touch; slice 10 is the pass that makes the spec read as one design rather
than a list of edits.

## What can run in parallel

- **1, 3 and 5** touch disjoint files (cases / dead modules / the fiber kernel) and can run at once, each in its own
  worktree. 3's narration removal edits `inventory.json`, `tsconfig.json` and `biome.json`, which 6 rewrites: land 3
  before 6 starts.
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
- **10** last; **11** whenever 1 and 6 are in.

## Later (not in these slices)

- **Persistence back**: no longer "add a backend": with the storage port, leases and the supervisor gone (slice 2), a
  durable or remote harness is a redesign. It can start from git (`packages/ai/{harness,components,conformance}` at
  `e91cf9bc0`): the host protocol, the durable core, and the dropped cases (conformance.md §2 "Dropped").
