# Threads: fibers renamed, and what an in-memory thread needs

## 1. What a fiber is today

`@molstar/ai-fibers` (2,027 source lines, 723 test lines in `fibers/src/_test/`, plus the 705-line suite at
`conformance/src/fibers/` and 187 lines of runners and row tests in `components/src/fibers/_test/`): an owner-scoped,
append-only list of entries with idempotent commands, streamed entries and bounded window reads, independent of the
harness (`packages/ai/_spec/fibers.md`).

| File (`fibers/src/`) | Lines | What                                                                                      |
| -------------------- | ----: | ----------------------------------------------------------------------------------------- |
| `kernel.ts`          |   663 | pure decisions per command (`prepareFiber*`), window and delta projection, cleanup pages   |
| `host.ts`            |   605 | `createFiberHost`: load rows, run the kernel, persist, retry `RowConflict`, notify, watch  |
| `domain.ts`          |   212 | records, commands, results, `FiberPort`                                                   |
| `limits.ts`          |   123 | `FiberLimits` (15 fields, no defaults), hard bounds, relation checks                       |
| `state.ts`           |   116 | copy-on-write drafts and the diff to row writes                                           |
| `index.ts`           |    95 | 79 exported symbols                                                                       |
| `data.ts`            |    73 | strict JSON checks, `cloneData`, canonical JSON                                           |
| `storage.ts`         |    65 | `RowStore` port, `RowConflict`, `FIBER_TABLES`                                            |
| `ids.ts`             |    38 | branded ids                                                                               |
| `memory.ts`          |    37 | `createMemoryRows`                                                                        |

**The one consumer** is prism-viewer's session record (`packages/prism/viewer/src/journal/record.ts`, 874 lines). It
calls `create`, `append`, `begin`, `heartbeat`, `settle`, `abort`, `window` and `forget` (`git grep "host\.\|target\."`
in `record.ts`), never `appendDelta`, `deltas`, `remove`, `subscribe` or `watch`. Its rows are `createMemoryRows()` by
default (`record.ts:239`) and `createIndexedDbRows('prism-viewer-records')` for a kept record (`record.ts:774`,
`ui/history.tsx:92`).

## 2. The rename map

Scope: `packages/ai/fibers`, `packages/ai/conformance/src/fibers`, `packages/ai/components/src/fibers`, and in
prism-viewer `journal/`, `storage/sources.ts`, `ui/history.tsx`, `viewer.ts`. About 1,750 case-insensitive occurrences
of "fiber" in 41 files (`git grep -io fiber` over those paths). **Not** in scope: the molecular "fibers" of
`packages/{graphics,model,plugin}`, `apps/mesoscale-explorer` and the viewer's particle representation
(`packages/prism/viewer/src/api/particles.ts:89`, `commands/volume.ts:57`); a blind replace would break them.

| Today                                                     | Proposed                                   |
| --------------------------------------------------------- | ------------------------------------------ |
| `@molstar/ai-fibers`                                      | `@molstar/ai` (`src/threads/`)             |
| `createFiberHost`, `FiberHost`, `FiberHostOptions`        | `createThreadHost`, `ThreadHost`, `ThreadHostOptions` |
| `FiberId`, `fiberId()`                                    | `ThreadId`, `threadId()`                   |
| `FiberEntry`, `FiberEntryId`, `fiberEntryId()`            | `ThreadEntry`, `ThreadEntryId`, `threadEntryId()` (not `EntryId`: the harness has one) |
| `FiberOwnerId`, `fiberOwnerId()`, `FiberCommandId`, `fiberCommandId()` | `ThreadOwnerId`, `threadOwnerId()`, `ThreadCommandId`, `threadCommandId()` |
| `FiberLimits`, `createFiberLimits`, `FIBER_LIMIT_HARD_BOUNDS`, `FIBER_CONFORMANCE_LIMITS` | `ThreadLimits`, `createThreadLimits`, `THREAD_LIMITS` (§3), — |
| `FIBER_PAYLOAD_VERSION`, `FiberPayloadEnvelope`           | `THREAD_PAYLOAD_VERSION`, `ThreadPayload`  |
| `FiberRows`, `FiberTable`, `FIBER_TABLES`                 | `ThreadRows`, `ThreadTable`, — (unused)    |
| `Fiber*Command`, `Fiber*Result`, `FiberWindow`, `FiberFailure`, `FiberPort` | `Thread*`                |
| table `'fiber'`, scope prefix `'fiber:'` (`host.ts:82`)   | `'thread'`, `'thread:'`                    |
| `_spec/fibers.md`                                         | `_spec/threads.md`                         |

**Stored data**: renaming the table and the scope prefix makes records kept under the old names unreadable. v6 is
unreleased and `.v6/ai/README.md` decision 9 says there is no data migration, so this is acceptable; if the kept
records move into the viewer (persistence.md §3) the viewer may keep its old names instead.

## 3. Trimming, now that threads live in memory

What the fiber machinery exists for, and whether an in-memory thread still needs it:

| Mechanism                                                            | Where                                                     | Needed in memory?                        | Size                     |
| -------------------------------------------------------------------- | --------------------------------------------------------- | ---------------------------------------- | ------------------------ |
| **Streamed deltas**: `appendDelta`, byte offsets, delta pages, gap detection | `kernel.ts:278–353, 585–650`, `host.ts:597+`, domain types, `deltas` table | no production user (§1)            | ~250 src, ~250 test (inferred from the case list) |
| **Cleanup pages** (`nextFiberCleanupPage`, `cleanupPageSize`)        | `kernel.ts:651–663`, limits; "validated but not used" (`_spec/fibers.md` Gotchas) | no                            | ~40                      |
| **Optimistic concurrency** across hosts: generations, `RowConflict`, retry ×8, `generation()` re-checks | `host.ts:131–197`, `storage.ts` | only with a shared store (another tab) | ~80              |
| **Watch** (another handle's writes)                                  | `host.ts:103–125, 511–520`, `RowStore.watch`             | only with IndexedDB (BroadcastChannel)   | ~40                      |
| **Receipts** (idempotent commands across reloads, cross-fiber reuse) | kernel `replayOrConflict`, `receipt` table                | a retry inside one process is still useful; across reloads no | keep, simplify |
| **Row diffing** (copy-on-write drafts → row writes)                  | `state.ts`                                                | only to write rows; memory could hold the state | 116 if the row port goes |
| **Inactivity aborts and heartbeat**                                  | kernel, host                                              | yes: the viewer's turn entries heartbeat (`record.ts:329`) | keep        |
| **Limits** (15 fields, hard bounds, relations)                       | `limits.ts`                                               | partly: the viewer uses hard bounds with `maxActiveStreams: 8` (`record.ts:208`) | ~60 of 123 |

Two shapes, by open question 1 (kept records). **Decided** (2026-10-16): A, with the kept records' IndexedDB rows in
prism-viewer; the streamed deltas **stay** (a later streaming record may need them), only the unused cleanup pages go:

- **A. Keep the row port** (`RowStore`) because the viewer's kept records write IndexedDB rows: rename, drop the
  cleanup pages, keep the deltas, the concurrency and the watch for the IndexedDB rows wherever they live. About **40
  source lines** go (the decided shape).
- **B. No row port**: if nothing persists threads, the host holds `ThreadState` in a `Map`; `state.ts`, `storage.ts`,
  `memory.ts`, the conflict loop and the watch go. About **700 source lines** go (inferred), and the viewer's
  `FiberRows` uses (`journal/kept-records.ts:9`, `storage/sources.ts:9`, `viewer.ts:11`, the spec's `persist.rows`) are
  removed with Keep.

Both keep the API the viewer calls (§1) unchanged apart from names.

## 4. Tests

- `fibers/src/_test/kernel.test.ts` (494), `host.test.ts` (163), `review-cursor-1.test.ts` (66): stay (the delta cases too).
- `conformance/src/fibers/` (705: `conformance.ts` 539, `fixtures.ts` 95, `generators.ts` 47) runs on memory and
  IndexedDB rows (`components/src/fibers/_test/conformance-{memory,indexeddb}.test.ts`). With one backend it becomes
  `threads/_test/thread.test.ts`; its three unused exports (`FIBER_INVALID_ENVELOPES`, `fiberCleanupCursorPages`,
  `fiberEnvelopeIsValid`, `generators.ts:23–47`) go.
- `components/src/fibers/_test/rows.test.ts` (129) compares memory and IndexedDB rows: goes with the IndexedDB rows, or
  moves with them.

## 5. What breaks outside

prism-viewer: `journal/record.ts` (imports at :9–25, 101 occurrences), `journal/kept-records.ts`, `storage/sources.ts`,
`viewer.ts`, `ui/history.tsx:92`, tests `journal/_test/{record,kept-records}.test.ts`, `storage/_test/sources.test.ts`.
Docs: `.v6/prism-ui/spec/journal.md` (9 mentions), `spec/README.md:377`, `spec/plan.md:45` (A7, "a fiber listing by
owner"), `handoff.md`, `design/agent.md`, `.v6/ai/analysis/viewer-mode.md` §6, `.v6/plans/breaking-v6-changes.md`,
`CHANGELOG.md`, `examples/minimal-prism-viewer/README.md:107`.
