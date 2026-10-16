# Duplicated concepts

Most of these exist because the package rules forbade the natural import (`biome.json`: "The trace reads the harness
only", "Fibers import nothing of the workspace"). In one package they can be one definition.

## 1. The request markers

`REQUEST_MARKER = 'Request: '` and `FOCUS_MARKER = 'Focus: '` are defined in `session/src/compose.ts:19,21` and again
in `trace/src/read.ts:50,52`, because the trace may not import the session. With labels on every entry (`.v6/ai/README.md`
decision 3) the markers are read only for entries written before labels (persistence.md §5), so with no persistence
the trace's copies can go entirely; else the trace imports `session/compose.ts`.

## 2. Token estimates

Three "four characters a token" estimates: `model-adapters/src/common/estimate.ts:20` (`estimateTextTokens`, with a
profile), `trace/src/read.ts:152` (`estimatedTokens`), and the compaction's sizes (`session/src/compaction/plan.ts:460`,
`estimatedTokens` of an edit, which sums per-entry sizes). One text estimate in `messages/` (the providers' admission
estimate stays the providers', as it can be calibrated); the compaction one renamed (`editTokens`, packaging.md §3).
~15 lines.

## 3. Prices

`ModelPrice` and `dollarsOf` twice: the catalog's (`model-adapters/src/model-catalog.ts:35,241`, per entry) and the
session's (`session/src/prices.ts`, per table, plus `creditsOf`). Neither is called (dead-code.md §1). Keep the
catalog's data and, if any, its `dollarsOf`; delete the session's. 48 + 34 lines.

## 4. A model's limits, in three places

The catalog (`limits` per entry), each adapter's `modelLimits` override for a custom endpoint
(`_spec/models.md` § The adapters), and the session's `ModelChoice.contextWindow` / `maxOutput` and
`AgentSessionOptions.contextWindow` (`session/src/types.ts:117,119`, `_spec/session.md:33,50`). The session then
reconciles them (`session/src/session.ts:276–295`, `capabilitiesOf`, which also forces `appendsTools: false` at :295).
Proposal: one source, the provider's `capabilities()` (catalog, else `modelLimits`); the session reads it and drops its
own window options. Size small (inferred ~40 lines), the gain is one answer to "how big is this model's window".

## 5. Two storage ports, two clocks, two JSON encodings

| Concept         | Harness                                                   | Fibers                                                 |
| --------------- | --------------------------------------------------------- | ------------------------------------------------------ |
| storage port    | `HarnessStorage` (`harness/src/contracts/storage.ts`)     | `RowStore` (`fibers/src/storage.ts`)                   |
| memory backend  | `createMemoryStorage` (272 lines)                         | `createMemoryRows` (37)                                |
| IndexedDB       | `createIndexedDbStorage` on `createDurableStorage`        | `createIndexedDbRows`                                  |
| suite           | `runHarnessConformance` (219 cases)                       | `runFiberConformance` (12)                             |
| clock           | `Clock` (`harness/src/contracts/runner.ts:37`)            | `FiberClock` (`fibers/src/domain.ts:124`, `{ now() }`) |
| canonical JSON  | `canonicalJson` (`agent/src/hash.ts:107`, RFC 8785)        | `canonicalData` (`fibers/src/data.ts:58`)              |
| branded ids     | `HarnessId`, `EntryId`, `CallId` …                        | `FiberId`, `FiberEntryId` …                            |

Threads stay independent of the harness by design (`_spec/fibers.md`: "Independent of the harness"), which is worth
keeping; but in one package they can share `messages/`' `Json`, `canonicalJson` and one `Clock` type. With IndexedDB
gone the two ports, backends and suites stop being parallel infrastructure (persistence.md, threads.md §3).

## 6. "Segments", two meanings

The harness stores heads as `Segment`s (content-addressed slices of a head in size classes,
`harness/src/head/index.ts`); adapters report a request's `segments` (short digests of tools, system and each message,
`model-adapters/src/common/adapter.ts:43`) for the trace's prefix breaks. Same word, unrelated things, both on the
model-call path (`harness/src/runner/model.ts:85`). Rename the request's to `prefix` or `parts` in the providers
rewrite.

## 7. Packages that only exist as a layer

`@molstar/ai-agent` (377 lines) holds types and hashing so that harness and adapters need not import each other; and
`@molstar/ai-components` exists so that "only one package touches the platform" (`.v6/ai/README.md` decision 8), which
with IndexedDB gone means it holds a memory store and the in-process supervisor, both platform-free. As folders,
`messages/` and `harness/` take them (packaging.md §2).

## As done (slice 9, 2026-10-16)

- **§1 markers**: the trace still needs `REQUEST_MARKER` (to tell a request entry's state from its words, by the
  label's prompt), so it imports it from `session/compose.ts` (the one file its Biome block allows); its copies and
  `FOCUS_MARKER` are gone. Every reading of an entry by its text or marker is gone (persistence.md §5): the toolkit's
  `groupsIn`, the session's `catchUp`, compaction's `isRequest`, `requestTextOf`, `isCarryOver` and `latestGuidance`,
  the trace's `defaultRequestOf` and `defaultWindowPart`. Without a label an entry is no request, no guidance, no
  skill.
- **§2 estimates**: `estimatedTokens(text)` (four characters a token) is in `messages/estimate.ts`, used by the trace,
  the harness's context and output sizes, and the scripted model. The compaction's per-edit estimate is `editTokens`;
  its `sizeOf` (three characters a token, cautious on purpose) and the providers' calibrated `estimateTextTokens` stay
  theirs.
- **§4 limits**: nothing outside `packages/ai` set `ModelChoice.contextWindow` / `maxOutput` or
  `AgentSessionOptions.contextWindow` (Prism's settings model builds its providers without `modelLimits`, so an own
  endpoint's models, known only by `modelLimits`, fail the turn with the provider's error; it did already, as no model
  choice gave a window). They are removed; the session reads the provider's `capabilities()` and `contextWindow()` is
  200k until the first call has asked. An own endpoint's limits, when the assistant adds them, go to the provider's
  `modelLimits`.
- **§5 clock and JSON**: `Clock` (with `sleep`) is defined in `messages/` and re-exported by the harness's
  `contracts/runner.ts`; threads read `Pick<Clock, 'now'>` (`ThreadClock` is gone) and use `canonicalJson` of
  `messages/` (`canonicalData` and the unused `strictCanonicalData` are gone). The threads' Biome block allows
  `messages/`; the strict-JSON checks (`assertStrictJsonData`) stay.
- **§6 segments**: the providers' request digests are `ProjectedRequest.prefix` (the model call's evidence and the
  trace follow); no other `segments` remains in `packages/ai`.
