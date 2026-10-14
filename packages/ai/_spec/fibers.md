# Fibers (`@molstar/ai-fibers`)

A fiber is an owner-scoped, append-only list of entries with idempotent commands. An entry is appended settled, or
streams byte-offset deltas and then settles or aborts. Reads are bounded windows by position and delta pages with gap
detection. It is what the harness lacks for a person-facing history: the journal is what the model sees (compaction
drops turns from it) and turns have no records with a prompt and a reply. Independent of the harness: it imports nothing
of the workspace and has its own storage port. Background: `.v6/ai/analysis/viewer-mode.md` §6.

**Status**: built and conformance-tested; no production code uses it yet. A building block an application composes: for
example an app writing an entry per turn from the session's hooks, and its History reading them in windows.

## Files (`fibers/src/`)

| File         | What                                                                                               |
| ------------ | -------------------------------------------------------------------------------------------------- |
| `domain.ts`  | the records, commands, results and `FiberPort` (the ten operations)                                |
| `kernel.ts`  | pure decisions: `prepareFiber*` per command against the state given, `projectFiberWindow`/`Deltas` |
| `host.ts`    | `createFiberHost`: loads state from rows, runs the kernel, persists, retries conflicts, notifies   |
| `state.ts`   | `FiberState` (one fiber as loaded), copy-on-write drafts, the diff to row writes                   |
| `storage.ts` | the row-store port `RowStore`, `RowConflict`, `FIBER_TABLES`; `memory.ts`: `createMemoryRows`      |
| `limits.ts`  | `FiberLimits`, `createFiberLimits`, `FIBER_LIMIT_HARD_BOUNDS`, `FIBER_CONFORMANCE_LIMITS`          |
| `ids.ts`     | branded ids (`FiberId`, `FiberEntryId`, `FiberOwnerId`, `FiberCommandId`) and their constructors   |
| `data.ts`    | strict JSON checks, `cloneData`, canonical JSON (fingerprints and byte sizes)                      |

## Using it

```ts
const fibers = createFiberHost({ rows: createIndexedDbRows('my-app-fibers') }); // or createMemoryRows()
// ids are branded: fiberId('…'), fiberEntryId, fiberOwnerId, fiberCommandId (non-empty strings)
await fibers.create({ commandId, fiberId, ownerId, limits: createFiberLimits({ ...FIBER_LIMIT_HARD_BOUNDS, maxActiveStreams: 4 }) });
const envelope = (value: unknown) => ({ codec: 'my.turn.v1', version: FIBER_PAYLOAD_VERSION, value }); // strict JSON
await fibers.append({ commandId, fiberId, entryId, descriptor: envelope({ kind: 'turn' }), payload: envelope(turn) });
const page = await fibers.window({ fiberId, start, limit: 20, direction: 'backward' });
```

- **Commands** (`FiberPort`): `create`, `append` (settled), `begin` (streaming), `appendDelta` (`start` must equal the
  stream's `nextOffset`), `heartbeat`, `settle` (the caller gives the whole final payload; deltas are dropped), `abort`,
  `remove`; reads `window`, `deltas`. Commands return `{ status: 'accepted', value, commandId }` or a
  `FiberFailure { status: 'conflict' | 'invalid' | 'not-found' | 'closed', reason }`; `heartbeat` returns a bare
  `FiberDeltaResult | FiberFailure`; the reads return the page or a `FiberFailure` (`not-found`: `fiber-missing`,
  `entry-missing`; `invalid`: a limit out of range, an entry over the window's bytes).
- **Idempotent** by `commandId`: the same command again returns its stored result; the same id with other input is
  `conflict` (`command-reused`). `heartbeat` has no receipt.
- **Payloads** are envelopes `{ codec, version: 1, value }` of strict JSON. A delta's offset advances by the canonical
  JSON byte size of its whole envelope, not of streamed text.
- **Abort reasons**: `caller-interrupted`, `inactivity`, `limit-exceeded`, `host-closed`. A delta over a stream's limits
  is accepted with `kind: 'aborted'` (`limit-exceeded`) and the entry becomes terminal.
- **Inactivity** is checked at the start of every operation, reads included (there is no timer): a stream idle for
  `inactivityTimeoutMs` is aborted, and that is persisted.
- **Durability**: an accepted command has been written (`rows.apply` before it returns). Within a host, commands are
  serial; across hosts or tabs, optimistic concurrency on a per-scope generation (`RowConflict`, retried up to 8 times).
  No leases.
- **Observing**: `subscribe(listener)` for changes through this host, `watch(listener)` for other handles' (through
  `RowStore.watch`; the IndexedDB rows post on a `BroadcastChannel`, memory rows have none). `forget(fiberId)` drops the
  rows.

## Storage

`RowStore<T>`: `load(scope)`, `apply(scope, expected, writes)` (atomic; generation `expected + 1` or `RowConflict`),
`drop(scope)`, `generation(scope)`, optional `watch`. A fiber is the scope `fiber:<id>` with tables `fiber`, `entry`,
`stream`, `deltas`, `receipt`, each row `{ json }`. Backends: `createMemoryRows` (here), `createIndexedDbRows(name)`
(`@molstar/ai-components`). SQLite later.

## Limits

`FiberLimits` has no defaults: all 15 fields are given (start from `FIBER_LIMIT_HARD_BOUNDS`). It has entry, batch,
window and delta-page sizes, active streams (count, bytes, batches), heartbeat, inactivity and aborted-preview
retention, cleanup page size. `FIBER_LIMIT_HARD_BOUNDS` are the ceilings (e.g. entry 1 MB, 32 active streams, window 256
entries, inactivity 1 h); relations are checked (heartbeat < inactivity, batch ≤ entry, window ≥ 2 × entry …).

## Gotchas

- `create` refuses limits beyond `FIBER_LIMIT_HARD_BOUNDS` (`invalid`, `limits`); build them with
  `createFiberLimits(limits)` to find out early.
- Reads can write (an inactivity abort).
- A removed fiber stays fenced (`removed`); creating the same id again conflicts. `forget` drops it.
- Receipts are kept until `remove`; `cleanupPageSize`, `deltaFlush*` and `heartbeatIntervalMs` are validated but the
  host does not use them (batching and heartbeat pacing are the caller's).
- Cross-fiber `commandId` reuse is caught only for fibers loaded in this host: derive command ids from the fiber.
