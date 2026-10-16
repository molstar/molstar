# Results: the review gate (slice 13)

The simplification as built, measured against this analysis, 2026-10-16, on `v6-prism-impl-proto` at `6b08f12ba`. The
gate is [plan.md](plan.md) § "The review gate"; each item below is one of its seven.

## 1. The whole repository green

Twice, each on a clean checkout of the branch tip (a fresh worktree): first at `3de7599e5` (before the review fixes),
then at `6b08f12ba` (after).

| Step                                | `3de7599e5`             | `6b08f12ba`                                       |
| ----------------------------------- | ----------------------- | ------------------------------------------------- |
| `pnpm install --frozen-lockfile`    | ok                      | ok                                                |
| `npx tsc -b`                        | ok                      | ok                                                |
| `pnpm lint`                         | ok                      | ok                                                |
| `pnpm format:check`                 | ok                      | **fails** on `packages/prism/ui/src/components/icon/wordmark.tsx:38` |
| `pnpm test:unit`                    | ok: 8,825 passed        | 8,865 passed, **1 unhandled error** (exit 1)      |
| `pnpm build`                        | ok                      | ok                                                |
| `pnpm check:workspace`              | ok                      | ok                                                |
| `node smoke/run.mjs browser`        | ok                      | ok                                                |

Both failures at `6b08f12ba` come from other sessions' Prism work landed in between (the Home heading as glass, the
splash), not from this simplification: no commit of it touches those files. The unhandled error is
`ReferenceError: OffscreenCanvas is not defined`, attributed to
`packages/prism/viewer/src/features/home/_test/load-shows-structure.test.tsx`. The lead session was told.

Not run: CI on the pushed branch (it needs a push, the owner's call). A Prism structure test
(`packages/prism/viewer/src/sources/_test/representation-quick.test.ts`) timed out once under heavy load and passed
3 of 3 runs on its own: flaky under load, not AI.

## 2. Code review of the branch

An Opus review of the 13 commits (`1fd60f0af`…`3de7599e5`): 1 high, 4 medium, 16 low, plus 26 low in the spec read
(item 5). Verified with probe tests before reporting. Every finding was fixed, each with a test that fails before and
passes after (where there is behavior to test), in 12 commits `6218b3799`…`6b08f12ba` (plan.md slice 13 maps each to
its findings):

- **High**: the person's own OpenAI-compatible endpoints could not run a turn (no model limits; predates the work).
  The owner decided: optional context-window and max-output fields in the own-endpoint setup, passed as `modelLimits`,
  defaults 32,768 and 4,096 (`6b08f12ba`). The form was checked in Storybook.
- **Medium**: the session's output clamp never ran with a real provider (a slice 9 regression: any model with an output
  maximum under 16,384 failed every turn); the harness kept its first provider for the whole conversation; a resumed
  turn lost its parked step's plan (tool attempts reported as model calls); a park inside a nested run re-ran finished
  sibling calls with `prior: 'none'` (a write applied twice).
- **Low**: reset with no completed turn; parallel calls overshooting `toolCalls`; the `endBy` watcher outliving its
  step; `./testing/threads` not packable; `./trace` pulling tools and zod through `session/`; two wall clocks; a
  class-based provider losing its methods; repeated "shown" events after a resume; noisy test names; comments citing
  documents or sections that do not exist.

The review also confirmed what holds: atomic steps (every write after the last await; stop beats a settlement), no race
or unhandled rejection in stop / park / resume / reset, nested failures stay local, projection pure and digests stable,
the providers' replay paths, the OpenAI fetch rewrite, error mapping with the key redacted, the MurmurHash3 fix.

## 3. The invariants

Checked by `git grep` at `6b08f12ba`; all hold:

- no `@molstar/ai-` package name in code, configuration or the lockfile (only in history and this folder);
- one package, `packages/ai/package.json`: `@molstar/ai`, `"sideEffects": false`, explicit exports (`.`, `./harness`,
  `./providers`, `./providers/{anthropic,openai,google}`, `./models`, `./trace`, `./testing`, `./testing/threads`,
  `./narration`, `./narration/mediabunny`), no `./*`;
- no `zod` import outside `packages/ai` (a repository-wide Biome ban; within, only `providers/` and `tools/`), none in
  the public types;
- in `packages/ai/src`, one Mol\* import (`@molstar/core/data/util/hash-functions`, in `messages/hash.ts`), no
  `crypto.subtle`, no `node:*`, no rxjs;
- no `fiber` identifier in `packages/ai` or the viewer's journal (the molecular fibers elsewhere untouched);
- no IndexedDB in `packages/ai/src`; the kept records' rows are `packages/prism/viewer/src/journal/indexeddb-rows.ts`;
- every `harness.md §n` / `tools-and-session.md §n` citation names an existing section (198 checked in slice 10, the
  review's B14 comments fixed in `809d662a1`).

## 4. The numbers

**Lines** (`.ts`/`.tsx` under `packages/ai`, `wc -l`; source excludes `_test/`):

| | Before (`e91cf9bc0`) | After slice 10 (`3de7599e5`) | After the fixes (`6b08f12ba`) |
| --- | ---: | ---: | ---: |
| Source | 25,921 | 13,887 | 14,747 |
| Tests | 11,250 | 12,032 | 11,990 |
| **Total** | **37,171** | **25,919** (−30%) | **26,737** (−28%) |

The conformance cases count as source before (they were the conformance package's `src`) and as tests after
(`harness/_test/`). Against the README's estimate (~15,500 of ~38,400 lines, ~40%) the cut is smaller, by category:

| Category | Estimated | Measured |
| --- | --- | --- |
| harness + components (durability, storage, supervisor) | ~2,750 source | 5,530 → 2,888: −2,642 |
| conformance (219 → 59 cases, meter, suite) | ~5,650 | 8,914 → 3,799 (`harness/_test` 2,482 with the fix round's cases; `testing/` 1,317 with the threads' suite): −5,115 |
| providers | ~1,400 source, ~2,250 tests | 2,630 → 2,256 source (−374), 3,488 → 2,865 tests (−623) |
| dead code | ~700 / ~350 | ~500 in slice 3 (`runOnce`, prices, `callTool`, …), 323 with `zod.ts` and its test in 12a |
| package plumbing | 10 → 1 `package.json` and `tsconfig.json`, 9 Biome blocks → 9 folder blocks | as estimated |

The providers are the gap: the rewrite kept the strict attachment checks, the stream and error mapping and each
provider's replay code as they were, rather than leaning on the AI SDK for them (providers.md "As done"). The fix round
added ~850 source lines back (the own-endpoint form, reservation at admission, results kept across a nested park) and
18 harness cases.

**Packages**: 10 → 1. **Conformance**: 219 cases on two backends → 59 cases as plain harness tests (77 tests in
`harness/_test/` after the fixes).

**Test time**: the AI tests were ~22.5 s of file time locally and a ~52 s serial tail in CI (run 38695888995). Now 53–55
files, 515–540 tests, ~3–4.5 s of file time on a quiet machine (10.4 s while a full repository run loaded the
machine), in the parallel `unit` project: no serial tail (slice 11; 9 full runs and 13 AI runs under a load average of
30–70, all green). CI not measured (no push).

**Bundles** (`pnpm build`, the Prism viewer app): the startup files (`index.js` 6 KB, `prism-viewer.js` 6.3 MB) contain
neither zod nor the AI SDK; zod is in one 195 KB shared lazy chunk; each provider is its own lazy chunk (Anthropic
138 KB, Google 150 KB, OpenAI 222 KB). Not measured before the work, so no comparison.

## 5. The spec against the code

The review's B1–B13 (retry counts, what runs next, the hash's input, `segments` → `prefix`, a missing
`CompactionDomain` member, what observations carry, how labels and markers are read, the line-count claim, missing
slice hashes, the tool groups' fields, direct runs' argument check, a stale packaging note, stale catalog comments) are
fixed in `38d6e6d55`; B14 (own endpoints need limits) with the high finding. Links and section citations resolve.

## 6. The app

In Storybook (`Viewer/Assistant`: the real viewer and the real `AgentSession` over a scripted model, after slice 10, and
the own-endpoint form after the fixes):

- a turn completes: two tool calls run, the surface is shown, the card offers Undo turn; no console errors;
- an approval parks the turn ("Waiting for you", Deny / Allow once / Allow example.org); Allow once resumes it and the
  tool runs (the load then fails on the story's unreachable URL, as scripted to go on: the story's model has one step);
- Stop ends a running turn as stopped, with its usage;
- the Inspector shows the turn's trace;
- the own-endpoint form shows the optional context window and max output with their defaults.

Not done: a kept record surviving a reload (covered by the viewer's journal tests, not checked in the browser), and a
live turn on a real model (the owner did not ask; it costs money).

Seen on the way, not AI: a load from an unreachable URL surfaces in the viewer as
`Cannot read properties of undefined (reading 'cell')` instead of a plain network failure.

## 7. The owner's sign-off

Pending.
