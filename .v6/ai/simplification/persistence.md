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
`stale-schema`) are the whole story. So "remove its migrations" removes nothing; the schema-version check goes with the
durability machinery (§4).

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

## 4. The durability machinery behind the storage port: removed too (decided)

IndexedDB was the reason for much of the harness's protocol: a turn that survives a reload, a runner in another tab
taking over, two tabs writing one harness. With memory only, a harness dies with its page. The owner decided
(2026-10-16) that this goes as well: "just basic in memory harness for now for driving an app". What goes, what stays
and the new shape are in [harness.md](harness.md): the host protocol, the storage port and its memory backend, leases
and launches, the supervisor, schema versions, crash re-runs, the runner's cache and spill, segment heads, steering and
interrupts. About 2,750 of the harness's and components' ~4,950 source lines go, and ~300 of the session's.

Consequence for later: persistence, a remote runner or a second tab are no longer "add a backend"; they are a redesign
that can start from git (`packages/ai/{harness,components}` at `e91cf9bc0`) and the dropped conformance cases.

## 5. Legacy fallbacks that go with persistence

Code that reads entries "written before labels" exists because stored conversations could predate labels (2026-10-14,
`.v6/ai/README.md` decision 3). With nothing persisted there are no such entries:
`session/src/compaction/domain.ts:50`, `compaction/plan.ts:129`, `session/src/session.ts:361`,
`tools/src/toolkit.ts:330`, `trace/src/read.ts:56`, and the test `tools/src/_test/toolkit.test.ts:374`. Small (inferred
~40–80 lines), but each is a second code path.
