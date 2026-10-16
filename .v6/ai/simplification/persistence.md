# Persistence: IndexedDB out, memory only

## 1. What exists

Two storage ports, each with a memory and an IndexedDB backend, all in `@molstar/ai-components` except the fibers'
memory rows:

| Backend                                   | File                                           | Lines | Used in production by                                      |
| ----------------------------------------- | ---------------------------------------------- | ----: | ---------------------------------------------------------- |
| harness memory storage                    | `components/src/harness/storage/memory.ts`     |   272 | `AgentSession` by default (`session/src/session.ts:395`)   |
| harness durable core (`createDurableStorage`, `DurableBackend`) | `components/src/harness/storage/durable.ts` | 221 | only `indexeddb.ts`                                     |
| harness IndexedDB (`createIndexedDbStorage`, `createIndexedDbBackend`, BroadcastChannel) | `components/src/harness/storage/indexeddb.ts` | 206 | **nobody** |
| fiber memory rows (`createMemoryRows`)    | `fibers/src/memory.ts`                         |    37 | prism-viewer's session record (`journal/record.ts:239`)    |
| fiber IndexedDB rows (`createIndexedDbRows`) | `components/src/fibers/indexeddb.ts`        |   130 | prism-viewer's **kept records** (`journal/record.ts:774`, `ui/history.tsx:92`) |

The viewer's agent accepts a `storage` (`packages/prism/viewer/src/agent/session.ts:29,121`) but nothing passes one
(`git grep "createViewerAgent("`: `ui/shell.tsx:518`, `agent/eval.ts:94`, the docs), so every harness in the
repository lives in memory today.

**Migrations**: there are none. `.v6/ai/README.md` decision 9 ("No data migration") and the schema version
(`SCHEMA_VERSION = 1`, `harness/src/contracts/records.ts`; a claim from another version blocks the turn,
`stale-schema`) are the whole story. So "remove its migrations" removes nothing; the schema-version check is option B
(§4).

## 2. What goes (the harness's IndexedDB)

- `components/src/harness/storage/durable.ts` (221), `indexeddb.ts` (206); their exports
  (`components/src/index.ts:9–16`: `createDurableStorage`, `DurableBackend`, `RecordKind`, `WriteSet`,
  `createIndexedDbBackend`, `createIndexedDbStorage`).
- `memory.ts` keeps the event-log helpers `durable.ts` imports from it (`addEvent`, `createLog`, `deepEqual`,
  `firstAfter`, `lastSeq`, `durable.ts:26`); they become private.
- `components/src/harness/_test/conformance-indexeddb.test.ts` (43): the whole second run of the suite, 220 cases,
  **13.9 s locally, 25.9 s in CI** (conformance.md §1).
- `fake-indexeddb` from `components/package.json:44`.
- Session option docs: `storage: createIndexedDbStorage('my-app-agent')` (`packages/ai/_spec/session.md:28`), "a reload
  with IndexedDB storage and the same `harnessId` finds the conversation" (`_spec/README.md` § How it binds, step 7),
  `_spec/harness.md` § Storage, `.v6/ai/harness.md` §14 and §18 ("IndexedDB in the page"), `.v6/ai/README.md` packages
  table and decision 8.
- The session test that reopens a harness on shared storage (`session/src/_test/drive.test.ts:186–199`, `h-legacy`)
  shares one `createMemoryStorage()` between two sessions (`:178`); it stays.

Nothing outside `packages/ai` breaks.

## 3. The decision: the viewer's kept records

The fibers' IndexedDB rows back a real feature: History's **Keep** copies the session record into
`'prism-viewer-records'`, and Home lists kept records from a small index in the same rows
(`.v6/prism-ui/spec/journal.md` §1, "Listing kept records"; `journal/kept-records.ts`; `storage/sources.ts`; the spec's
`persist.rows` / the record's `keepIn`). Tests: `journal/_test/record.test.ts:354–365`,
`storage/_test/sources.test.ts:111`.

| Option                                                        | `@molstar/ai`                         | prism-viewer                                                         |
| ------------------------------------------------------------- | ------------------------------------- | -------------------------------------------------------------------- |
| **A. Move the rows into prism-viewer** (recommended)          | memory only; keeps the `ThreadRows` port (threads.md §3 A) | gains `journal/indexeddb-rows.ts` (130 lines) and `rows.test.ts`' IndexedDB half; imports change at `record.ts:774`, `history.tsx:92`, two tests |
| B. Drop Keep for now                                          | memory only, no row port (threads.md §3 B) | Keep, the kept-records index, Home's list and `persist.rows` go; journal.md §1 rewritten |
| C. Keep the IndexedDB rows in `@molstar/ai`                   | one IndexedDB file stays (`threads/indexeddb.ts`) | unchanged apart from import paths                                    |

A meets "only the in-memory store stays" for the library and keeps the feature: the viewer owns its one persistent
store. It touches files other sessions are editing (`packages/prism/**`), so it is its own slice (plan.md).

## 4. Option B: the durability machinery behind the storage port

IndexedDB was the reason for much of the harness's protocol: a turn that survives a reload, a runner in another tab
taking over, two tabs writing one harness. With memory only, a harness dies with its page. What exists only for that
(read in the code; sizes partly inferred):

| Mechanism                                   | Where                                                                                       | Lines (approx.)  |
| ------------------------------------------- | ------------------------------------------------------------------------------------------- | ---------------- |
| Supervisor (liveness, takeover after a lease or launch window) | `components/src/harness/supervisor/in-process.ts` (233), `harness/src/contracts/supervisor.ts` (33), `conformance/src/testing/supervisor.ts` (48) | 314 |
| Leases and launch keys: `claim`, `renew`, `release`, `check`, `launching`/`live` runner state | `harness/src/host/runner.ts` (`hold`, `currentLease`), `host/client.ts` (`launchResult`, `withLaunch`), `runner/index.ts` (21 lines name them), `contracts/{ports,records}.ts` | ~400 (inferred) |
| Schema version and `blocked: 'stale-schema'` | `contracts/records.ts`, `host/{index,runner,writer}.ts`, `runner/index.ts` (21 mentions)   | ~60              |
| Runner cache with a spill under `cacheBytes` (to avoid store reads) | `harness/src/runner/state.ts` `WorkingSet` (`:91–262`), `contracts/inspect.ts`           | ~200             |
| Heads as content-addressed segments in size classes (constant-cost appends to a store) | `harness/src/head/index.ts` (88), `Segment`/`HeadRecord` records, the `segments` store | ~200 with callers |
| Tombstones (a removed id is never reused)   | `contracts/storage.ts`, `host/index.ts`, `memory.ts`                                        | ~40              |
| The session's takeover path (`interrupted()`, resume after the old lease, `leaseMs`) | `session/src/session.ts:114–115,396–397,705–771`                                  | ~100             |

Taken together ~1,200–1,800 source lines and, beyond what conformance.md already cuts, ~600 test lines
(`components/src/harness/_test/supervisor.test.ts` 173, the remaining T/S cases on leases). The loop would become: a
session submits, a runner in the same process runs the turn to its end or its park, and storage is a plain in-memory
journal and event log.

**Recommendation: keep it for this round.** The owner's words are IndexedDB removed "for now"; with the storage port,
leases and the supervisor intact, persistence comes back as one backend plus the conformance cases on it. Stripping
them is a redesign of the harness (its §4, §12, §14, §15 in `.v6/ai/harness.md`) and would have to be rebuilt for a
durable or remote host (§15 "Remote deployments"). Open question 2.

## 5. Legacy fallbacks that go with persistence

Code that reads entries "written before labels" exists because stored conversations could predate labels (2026-10-14,
`.v6/ai/README.md` decision 3). With nothing persisted there are no such entries:
`session/src/compaction/domain.ts:50`, `compaction/plan.ts:129`, `session/src/session.ts:361`,
`tools/src/toolkit.ts:330`, `trace/src/read.ts:56`, and the test `tools/src/_test/toolkit.test.ts:374`. Small (inferred
~40–80 lines), but each is a second code path.
