# `packages/ai` simplification: the analysis

Read-only analysis of the agent harness (`packages/ai/*`) for a drastic simplification, 2026-10-16, on
`v6-prism-impl-proto` at `e91cf9bc0`. No code was changed. The target is the owner's:

- `packages/ai` becomes one package, `@molstar/ai`;
- "fibers" become "threads";
- model adapters become **providers**, a thin layer over the AI SDK;
- the spec (`.v6/ai/*`, `packages/ai/_spec/*`) describes the simplified code;
- IndexedDB persistence goes; only the in-memory store stays;
- the conformance suite keeps only the essentials, so CI gets faster.

Line counts are `wc -l` over `.ts`/`.tsx` (no `cloc` on this machine). Source excludes `_test/`. Test timings are from
one local run (`vitest run packages/ai --project serial`, Apple silicon) and one CI run (GitHub Actions run
`38695888995`, "Native tests", 2026-10-15). Where a claim is **inferred** rather than read, it says so.

| File                               | What it covers                                                                               |
| ---------------------------------- | -------------------------------------------------------------------------------------------- |
| [packaging.md](packaging.md)       | ten packages into one: the folders, the exports map, name collisions, every caller to update |
| [threads.md](threads.md)           | fibers → threads: the rename map, what of the fiber machinery an in-memory thread still needs |
| [providers.md](providers.md)       | model adapters → providers over the AI SDK: what the layer does today, what stays, what goes |
| [harness.md](harness.md)           | the basic in-memory harness (decided): what stays, what goes, its API, the session over it   |
| [persistence.md](persistence.md)   | IndexedDB removed: what goes, the viewer's kept records, the durability machinery behind it  |
| [conformance.md](conformance.md)   | the suite: what it costs in CI, the essential cases, what goes                                |
| [dead-code.md](dead-code.md)       | exports nobody imports, speculative features, legacy fallbacks                               |
| [duplicates.md](duplicates.md)     | the same concept in two or three places                                                      |
| [plan.md](plan.md)                 | the migration in slices that keep the build green, checks per slice, what runs in parallel   |

## Today's packages

| Package                      | Purpose                                                                      | Src lines | Test lines | Workspace deps                       | Other deps                                         | Imported from outside `packages/ai` by                                                   |
| ---------------------------- | ---------------------------------------------------------------------------- | --------: | ---------: | ------------------------------------ | -------------------------------------------------- | ---------------------------------------------------------------------------------------- |
| `@molstar/ai-agent`          | `Message` and parts, `ModelAdapter` port, `hashMessage`, `digest`, `Json`     |       377 |        144 | none                                 | —                                                  | prism-viewer (`Json`), prism-assistant (`ModelAdapter`, `ModelFrame`, `JsonObject`)        |
| `@molstar/ai-harness`        | contracts, host protocol, runner, drives, tools step, nested runs            |     4,450 |         83 | agent                                | —                                                  | prism-viewer (`HarnessStorage` type only, `agent/session.ts:9`)                           |
| `@molstar/ai-components`     | memory and IndexedDB storage (durable core), in-process supervisor, IndexedDB rows for fibers | 1,080 | 791 | agent, harness, fibers | dev: fake-indexeddb | prism-viewer: only `createIndexedDbRows` (`journal/record.ts:774`, `ui/history.tsx:92`, two tests) |
| `@molstar/ai-model-adapters` | Anthropic, OpenAI (Responses, chat), Google on the AI SDK; the model catalog |     2,624 |      3,417 | agent                                | `ai`, `@ai-sdk/{anthropic,google,openai,openai-compatible}`, zod | prism-assistant (`settings-model.ts`: three lazy deep imports and `listedModels`), prism-viewer live tests, docs |
| `@molstar/ai-tools`          | tools, field shapes, results, bindings, repeats, nesting, the toolkit (skills) |   1,715 |      2,147 | agent, harness                       | zod                                                | prism-viewer (`agent/{toolkit,tools,params,session}.ts`)                                  |
| `@molstar/ai-session`        | `AgentSession`, its drive, budgets, composition, compaction, `runOnce`, prices |   2,745 |      2,263 | agent, harness, components, tools    | —                                                  | prism-viewer (`AgentSession`, types), prism-assistant (types, `compactionFor`)            |
| `@molstar/ai-trace`          | `readTrace` and the Inspector's views; `./format`                            |     1,338 |      1,046 | agent, harness                       | —                                                  | prism-assistant (`inspector/*`, `viewer/inspector-sources.ts`)                            |
| `@molstar/ai-fibers`         | owner-scoped append-only transcripts: kernel, host, row-store port           |     2,027 |        723 | none                                 | —                                                  | prism-viewer (`journal/record.ts`, `journal/kept-records.ts`, `storage/sources.ts`, `viewer.ts`) |
| `@molstar/ai-conformance`    | the harness and fiber suites, the storage meter, test doubles (`./testing`)  | 8,914 (all test infrastructure) | 0 | agent, harness, fibers   | **vitest, fast-check as dependencies**             | prism tests and stories (`scriptedModel`, `ModelStep` from `./testing`)                   |
| `@molstar/ai-narration`      | text to speech (Gemini, OpenAI), WAV, Mediabunny AAC/Opus                    |       651 |        636 | none                                 | mediabunny                                         | **nobody** (only the opt-in boundary test lists its name)                                 |
| **Total**                    |                                                                              | **25,921** | **11,250** |                                     |                                                    | 57 import sites in 35 files (`packaging.md` §4)                                           |

Plus the spec: `packages/ai/_spec/` 870 lines, `.v6/ai/` 639 lines. Package plumbing: 10 `package.json` (460 lines), 10
`tsconfig.json` (173), ten Biome blocks (`biome.json` lines 90–600, 511 lines), ten inventory entries, ten root
`tsconfig.json` references.

## The target: one package, `@molstar/ai`

```
packages/ai/                  @molstar/ai (pnpm-workspace.yaml: 'packages/ai' instead of 'packages/ai/*')
  package.json  tsconfig.json  AGENTS.md  _spec/
  src/
    index.ts                  the host's surface: session, tools, threads, messages
    messages/                 was ai-agent: Message, parts, hashing, Json (the port moves to providers/)
    harness/                  was ai-harness, rewritten as a basic in-memory loop (harness.md); ai-components goes
      _test/                  the essential conformance cases, as plain tests (conformance.md)
    providers/                was ai-model-adapters, rewritten (providers.md); models.ts is the catalog
    tools/                    was ai-tools
    session/                  was ai-session (compaction/ inside)
    trace/                    was ai-trace (format.ts inside)
    threads/                  was ai-fibers, renamed and trimmed (threads.md)
    testing/                  was ai-conformance/testing + the createSession fixture
```

Exports map (packaging.md §2 has the reasoning):

| Entry                    | Holds                                                                            | Who needs it                                       |
| ------------------------ | -------------------------------------------------------------------------------- | -------------------------------------------------- |
| `.`                      | session, tools, threads, messages, the drive and harness contracts a host names  | prism-viewer, prism-assistant                      |
| `./harness`              | `createHarness`, `defaultDrive`, `drive`, the reader                             | tests, helpers that drive a harness directly       |
| `./providers`            | `createAnthropicProvider`, `createOpenAIProvider`, `createGoogleProvider`        | prism-assistant (lazy), live tests                 |
| `./providers/anthropic`, `./providers/openai`, `./providers/google` | one provider each, so a lazy import loads one SDK | prism-assistant `settings-model.ts:207–218` |
| `./models`               | the model catalog (data, no dependencies)                                        | prism-assistant `settings-model.ts:10`             |
| `./trace`                | `readTrace`, the views, the formatters                                           | prism-assistant inspector                          |
| `./testing`              | `scriptedModel`, the simulated clock and ids, message builders, `createSession`  | tests and stories, here and in Prism               |

No `./*` wildcard: every deep import that exists today is one of the entries above.

## The categories, sized

| Category                   | What                                                                                     | Lines removed (src / test)       | Lines moved or rewritten           | File            |
| -------------------------- | ---------------------------------------------------------------------------------------- | -------------------------------- | ---------------------------------- | --------------- |
| Rename only                | fibers → threads; `Fiber*` → `Thread*`, the `fiber` table and `fiber:` scope prefix      | 0                                | ~1,750 occurrences in ~40 files    | threads.md      |
| Collapse                   | ten packages into folders of one; one `package.json`, tsconfig, Biome block, inventory entry | ~1,000 of plumbing (package.json, tsconfig, Biome) | all 37,171 lines move; 57 external import sites, 196 internal import lines (126 files) rewritten | packaging.md |
| Replace                    | model adapters → providers: retained neutral payload, its decode and native projection go | ~1,700 / ~2,500 (inferred estimate) | ~900 new provider lines (inferred) | providers.md    |
| Remove: IndexedDB          | the durable core, IndexedDB storage and rows, their tests and conformance run            | 557 / 210                        | the viewer's kept records (decision) | persistence.md |
| Remove: durability (decided) | host protocol, storage port, leases, launches, supervisor, schema version, crash re-runs, cache, segment heads, steering | ~3,050 / ~1,300 beyond the conformance cut (inferred) | the session's takeover path | persistence.md §4 |
| Simplify: conformance      | 219 harness cases → 59; no second backend; no performance suite or storage meter        | ~5,950 of 8,914 (to ~2,950)      | the cases left become `harness/_test/` | conformance.md |
| Dead or speculative        | narration, `runOnce`, session prices, `callTool`, fiber deltas, `./zod`, pre-label fallbacks | ~1,600 / ~1,250             | —                                  | dead-code.md    |
| Duplicates                 | markers, token estimates, prices, canonical JSON, two storage ports, three model-limit sources | ~150 (and two concepts fewer) | —                                | duplicates.md   |

Totals, if every recommendation is taken (the durability machinery included, as decided), counting each line once:
about **18,000 of 38,400 lines go** (source, tests and package plumbing; ~47%), the package count goes from ten to one,
and the AI tests' serial tail in CI goes from about 52 s to an estimated 15 s (conformance.md §4), or to almost nothing
if they can leave the serial Vitest project.

## Top findings

1. **Nobody persists the harness.** `createIndexedDbStorage` has no caller outside `packages/ai`; the viewer never passes
   `storage` (`packages/prism/viewer/src/agent/session.ts:121`, and no caller sets it). Removing the harness's IndexedDB
   breaks nothing outside.
2. **The fibers' IndexedDB rows are used**: the viewer's kept records (History's Keep, Home's list) store a fiber in
   `createIndexedDbRows('prism-viewer-records')` (`journal/record.ts:774`, `ui/history.tsx:92`). Removing all IndexedDB
   removes that feature unless the rows move into prism-viewer (130 lines). This is the one removal with a visible cost.
3. **The conformance suite is most of the AI tests' time.** Locally 18 of 22.5 s; in CI 34 of ~52 s (IndexedDB 25.9 s,
   memory 8.5 s), all in the `serial` Vitest project that runs after everything else with `fileParallelism: false`
   (`vitest.config.ts:17`). Dropping the IndexedDB run alone saves ~26 s of CI; the performance (P) and liveness (S)
   cases are most of the rest.
4. **Narration has no consumer** (651 + 636 lines, the `mediabunny` dependency), and neither do `runOnce`, the
   session's `prices.ts`, the catalog's `dollarsOf`, `callTool`, the fiber deltas and the `./zod` entry.
5. **The adapter layer converts twice**: `Message` → a neutral retained payload (`common/payload.ts`) → AI SDK
   `ModelMessage` (`common/stream.ts`) → HTTP, plus a fetch-rewriting transport for OpenAI tool additions
   (`openai/disclosure-transport.ts`). A provider can project straight to the SDK's own message shape.
6. **One package removes copies the package rules forced**: trace copies the session's markers (`trace/src/read.ts:50`
   vs `session/src/compose.ts:19`), and pricing and token estimates exist twice or three times (duplicates.md).
7. **A single root entry has name collisions** to settle first: `dollarsOf`, `ModelPrice`, `estimatedTokens`,
   `REQUEST_MARKER`, `FOCUS_MARKER`, `ToolSpec`, `ToolOutcome` (packaging.md §3).
8. **The test doubles ship vitest as a dependency** (`conformance/package.json`): the suite imports vitest. Once the
   cases are `_test/` files, `./testing` needs neither vitest nor fast-check.
9. **The session uses a fraction of the harness's protocol** (harness.md §1): one item per turn, no steering, no
   interrupts, no withdraw, one runner; the host protocol, leases, launches and the supervisor serve durability and
   multiple runners only. The basic in-memory loop keeps the journal, the event log, the drive and the tools step,
   and the trace and the Inspector read it unchanged.

## Open questions for the owner

1. **Kept records** (finding 2): move the IndexedDB rows into prism-viewer (recommended: the viewer owns its one
   persistent store; `@molstar/ai` stays memory-only), drop Keep for now, or keep IndexedDB rows in `@molstar/ai`?
2. ~~Durability machinery~~: **decided 2026-10-16, removed**: "just basic in memory harness for now for driving an
   app" ([harness.md](harness.md)).
3. **Narration**: delete (recommended: no consumer; the story mode can bring it back from git), a `./narration` entry
   of `@molstar/ai` (then `mediabunny` becomes its dependency), or its own package?
4. **Tool additions mid-conversation** (providers.md §3): keep the native `tool_addition` / `additional_tools` path for
   the two catalog models with the quirk (Claude Sonnet 5.5, GPT-6 Luna), or accept that disclosing a tool changes the
   request's tool list (a cache break at the tools segment)?
5. **Fiber deltas** (threads.md §3): streamed byte-offset deltas, delta pages and gap detection have no production
   user (the viewer settles whole payloads). Remove (recommended)?
6. **Names**: `ModelAdapter` → `Provider`, `model-catalog` → `models`, `@molstar/ai/threads` or threads on the root
   entry? (packaging.md §2 proposes names.)
7. **`./zod`**: another session is adding a "reconsider zod" note to `packages/ai/TODO.md`; drop the entry here
   (unused outside) and decide zod with that note?
8. **Folder rules**: the Biome per-package import rules become per-folder rules on relative paths (inferred to be
   expressible with `noRestrictedImports` patterns; not verified), or a check in `scripts/workspace/`?
9. **`runOnce`** (the command menu's assistant, not built): delete now, or keep for that feature?
10. ~~Hashing~~: **decided 2026-10-16, synchronous** SHA-256 in plain TypeScript, same output as today
    ([harness.md](harness.md) §7).
11. ~~Blobs~~: **decided 2026-10-16, keep references**: request and result bodies stay blob refs, the trace unchanged
    ([harness.md](harness.md) §4).
