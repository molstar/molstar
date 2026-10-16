# Threads (`src/threads/` of `@molstar/ai`)

A thread is an owner-scoped, append-only list of entries with idempotent commands. An entry is appended settled, or
streams byte-offset deltas and then settles or aborts. Reads are bounded windows by position and delta pages with gap
detection. It is what the harness lacks for a person-facing history: the journal is what the model sees (compaction
drops turns from it) and turns have no records with a prompt and a reply. Independent of the harness: it imports only
`messages/` (`canonicalJson`, the `Clock` type) of the other folders and has its own storage port. Background:
`.v6/ai/analysis/viewer-mode.md` §6.

**Status**: built and conformance-tested. A building block an application composes: the Prism viewer's record is a
thread (`packages/prism/viewer/src/journal/record.ts`): every change the viewer's API runs and every agent turn (begun
at the prompt, settled at its end) is an entry; History reads them in windows, and a kept record goes to the viewer's
IndexedDB rows.

## Files (`src/threads/`)

| File         | What                                                                                                 |
| ------------ | ---------------------------------------------------------------------------------------------------- |
| `domain.ts`  | the records, commands, results and `ThreadPort` (the ten operations)                                 |
| `kernel.ts`  | pure decisions: `prepareThread*` per command against the state given, `projectThreadWindow`/`Deltas` |
| `host.ts`    | `createThreadHost`: loads state from rows, runs the kernel, persists, retries conflicts, notifies    |
| `state.ts`   | `ThreadState` (one thread as loaded), copy-on-write drafts, the diff to row writes                   |
| `storage.ts` | the row-store port `RowStore`, `RowConflict`; `memory.ts`: `createMemoryRows`                        |
| `limits.ts`  | `ThreadLimits`, `createThreadLimits`, `THREAD_LIMIT_HARD_BOUNDS`, `THREAD_CONFORMANCE_LIMITS`        |
| `ids.ts`     | branded ids (`ThreadId`, `ThreadEntryId`, `ThreadOwnerId`, `ThreadCommandId`) and their constructors |
| `data.ts`    | strict JSON checks and `cloneData`; fingerprints and byte sizes use `messages/`' `canonicalJson`     |

## Using it

```ts
const threads = createThreadHost({ rows: createMemoryRows() }); // or a RowStore of the app's (Prism: IndexedDB)
// ids are branded: threadId('…'), threadEntryId, threadOwnerId, threadCommandId (non-empty strings)
await threads.create({ commandId, threadId, ownerId, limits: createThreadLimits({ ...THREAD_LIMIT_HARD_BOUNDS, maxActiveStreams: 4 }) });
const envelope = (value: unknown) => ({ codec: 'my.turn.v1', version: THREAD_PAYLOAD_VERSION, value }); // strict JSON
await threads.append({ commandId, threadId, entryId, descriptor: envelope({ kind: 'turn' }), payload: envelope(turn) });
const page = await threads.window({ threadId, start, limit: 20, direction: 'backward' });
```

- **Commands** (`ThreadPort`): `create`, `append` (settled), `begin` (streaming), `appendDelta` (`start` must equal the
  stream's `nextOffset`), `heartbeat`, `settle` (the caller gives the whole final payload; deltas are dropped), `abort`,
  `remove`; reads `window`, `deltas`. Commands return `{ status: 'accepted', value, commandId }` or a
  `ThreadFailure { status: 'conflict' | 'invalid' | 'not-found' | 'closed', reason }`; `heartbeat` returns a bare
  `ThreadDeltaResult | ThreadFailure`; the reads return the page or a `ThreadFailure` (`not-found`: `thread-missing`,
  `entry-missing`; `invalid`: a limit out of range, an entry over the window's bytes).
- **Idempotent** by `commandId`: the same command again returns its stored result; the same id with other input is
  `conflict` (`command-reused`). `heartbeat` has no receipt.
- **Payloads** are envelopes `{ codec, version: 1, value }` of strict JSON. A delta's offset advances by the canonical
  JSON byte size of its whole envelope, not of streamed text.
- **Abort reasons**: `caller-interrupted`, `inactivity`, `limit-exceeded`, `host-closed`. A delta over a stream's limits
  is accepted with `kind: 'aborted'` (`limit-exceeded`) and the entry becomes terminal.
- **Inactivity** is checked at the start of every operation, reads included (there is no timer): a stream idle for
  `inactivityTimeoutMs` is aborted, and that is persisted.
- **Writes**: an accepted command has been written to the rows (`rows.apply` before it returns). Within a host, commands
  are serial; across hosts or tabs, optimistic concurrency on a per-scope generation (`RowConflict`, retried up to 8
  times). No leases.
- **Observing**: `subscribe(listener)` for changes through this host, `watch(listener)` for other handles' (through
  `RowStore.watch`; memory rows have none, the viewer's IndexedDB rows post on a `BroadcastChannel`). `forget(threadId)`
  drops the rows.

## Storage

`RowStore<T>`: `load(scope)`, `apply(scope, expected, writes)` (atomic; generation `expected + 1` or `RowConflict`),
`drop(scope)`, `generation(scope)`, optional `watch`. A thread is the scope `thread:<id>` with tables `thread`, `entry`,
`stream`, `deltas`, `receipt`, each row `{ json }`. The one backend here is `createMemoryRows`; a store that outlives
the page is the application's. The viewer's IndexedDB rows are such a store (`createIndexedDbRows(name)`,
`packages/prism/viewer/src/journal/indexeddb-rows.ts`), and its tests run the threads' conformance
(`runThreadConformance` of `@molstar/ai/testing/threads`, the suite in `src/testing/threads/`) on them.

## Limits

`ThreadLimits` has no defaults: all 14 fields are given (start from `THREAD_LIMIT_HARD_BOUNDS`). It has entry, batch,
window and delta-page sizes, active streams (count, bytes, batches), heartbeat, inactivity and aborted-preview
retention. `THREAD_LIMIT_HARD_BOUNDS` are the ceilings (e.g. entry 1 MB, 32 active streams, window 256 entries,
inactivity 1 h); relations are checked (heartbeat < inactivity, batch ≤ entry, window ≥ 2 × entry …).

## Gotchas

- `create` refuses limits beyond `THREAD_LIMIT_HARD_BOUNDS` (`invalid`, `limits`); build them with
  `createThreadLimits(limits)` to find out early.
- Reads can write (an inactivity abort).
- A removed thread stays fenced (`removed`); creating the same id again conflicts. `forget` drops it.
- Receipts are kept until `remove`; `deltaFlush*` and `heartbeatIntervalMs` are validated but the host does not use them
  (batching and heartbeat pacing are the caller's).
- Cross-thread `commandId` reuse is caught only for threads loaded in this host: derive command ids from the thread.
