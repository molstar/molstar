# The harness (`@molstar/ai-harness`)

Runs agent turns durably: the model is called, the tools it asks for run, and the loop goes on until the turn ends,
surviving reloads, crashes and runners that move. Mechanism only; every choice is the drive's. Design and rationale:
[.v6/ai/harness.md](../../../.v6/ai/harness.md) (the `§` numbers below are its sections).

## Files

| Path (`harness/src/`)     | What                                                                                                               |
| ------------------------- | ------------------------------------------------------------------------------------------------------------------ |
| `contracts/records.ts`    | every stored record: `Execution`, `ActiveTurn`, `Head`/`HeadRecord`/`Segment`, events, `Usage`, `Limits`           |
| `contracts/ports.ts`      | `HarnessHost`, `HarnessReader`, `HarnessClient`, `RunnerPort`, `Move`, `HeadEdit`, `HarnessError`                  |
| `contracts/drive.ts`      | `Decisions` (each kind's payload and answers), `Drive`, `DriveContext`, `DriveHandlers`                            |
| `contracts/runner.ts`     | `RunnerOptions`, `Runner`, `ToolBindings`/`ToolHandler`/`ToolFrame`, `RunContext`, `RunSpec`, `Clock`, `IdSource`  |
| `contracts/storage.ts`    | `HarnessStorage`, `StorageTx`, `EventLog`, `RecordStore`                                                           |
| `contracts/observe.ts`    | `HarnessObservable`, `Observation`                                                                                 |
| `contracts/supervisor.ts` | `Supervisor`, `LaunchReport`, `SupervisorPolicy`                                                                   |
| `host/`                   | `createHarnessHost`: the protocol operations over storage (`client.ts`, `runner.ts`, `writer.ts`)                  |
| `runner/`                 | `createRunner` (`index.ts`: claim, serve, commit), `model.ts` (model step), `run.ts` (shared step ops), `state.ts` |
| `tools/`                  | `plan.ts` (bind calls to definitions), `step.ts` (the tools step), `journal.ts` (a tool's journal)                 |
| `nested/`                 | `nest`: nested runs inside a call                                                                                  |
| `drive/`                  | `defaultDrive`, `drive(handlers, fallback)`, `legal`                                                               |
| `head/`                   | `buildMove`, `emptyHead`: the segment scheme                                                                       |

## The data

- **Entry**: a `Message` (`ai-agent`) addressed by `hashMessage` (`EntryId`). Immutable, deduplicated, never deleted
  while the harness exists. Its optional `label` (what it is: a request, guidance) is part of it; the harness stores it
  and never reads it, and never counts it in a size (`WorkingSet.tokensOf`).
- **Head**: an immutable list of entry ids (`Head { id, entries }`), stored as a `HeadRecord` of content-addressed
  `Segment`s in binary size classes (`head/index.ts`): an append is near constant, a head reads in ~log n segments.
- **Pointer**: `scope → head`. `MAIN` (`'main'`) is the turn's conversation and lives on the Execution; tools and
  helpers name others. A pointer moves only in a commit, by a `Move { pointer, from, head, segments }` (compare-and-set
  on `from`), built with `buildMove`.
- **Projectable**: every stored head has each tool call answered in the next user message; an assistant message with
  calls enters a head only with its results.
- **Event log**: the only evidence, numbered by `seq`. `StepEvent` (how a step ended: `calls` with attempts and usage,
  `moves`, `consumed` items, `compaction`, `next` step or `end`, `summary` on the ending one) and `ControlEvent`
  (`submitted`, `withdrawn`, `stop-requested`, `stopped`, `resumed`, `reset`, `claimed`, `blocked`, `released`,
  `launched`). Turns have no record of their own: `client.turns()` reads them from the log. Large bodies (a call's
  request, its result, a tools step's assistant message) are **blobs**, `ref = digest(body)`, put before the commit.
- **Execution**: the one mutable record per harness: `seq`, `main` (head id), `inbox: Item[]`, `runner: RunnerState`
  (`none` | `launching(key, until)` | `live(runnerId, leaseUntil, …)`) and `turn?: ActiveTurn` (configuration, budget,
  `pending` step, `park`, `stopRequested`, `blocked`, usage).
- **What runs next** is `execution.turn.pending` (from the last step event), never derived from a head.

## Handles and operations

`createHarnessHost(storage, { clock?, leaseMs = 30_000, launchMs = 60_000 })` → `HarnessHost`: `create(id, { key })`,
`remove(id)` (tombstones the id, never reused, and drops its records), `client(id)`, `runnerPort(id)`. Apps use the
**client**, runners the **runner port**; an app can never commit and a runner never submit, stop or reset.

| Client (`HarnessClient`)                          | What                                                                                                     |
| ------------------------------------------------- | -------------------------------------------------------------------------------------------------------- |
| `submit({ key, messages, options?, interrupt? })` | an inbox item; idempotent by key; `options` is the app's (the session puts the turn configuration there) |
| `withdraw({ key })`                               | removes a waiting item                                                                                   |
| `stop({ key, turnId \| item, reason? })`          | durable Stop; naming an item withdraws it if it still waits, else stops the turn it began                |
| `resume({ key, turnId })`                         | re-runs a parked step, asking the drive again                                                            |
| `reset({ key, pointers?, reason? })`              | ends any turn, moves pointers (default `main` to its head at the last completed turn); deletes nothing   |
| `check()`                                         | returns a launch key if work waits and no live runner holds it (expired lease)                           |
| reads (`HarnessReader`)                           | `execution`, `heads`, `pointers`, `entries`, `events`, `steps`, `turns`, `blobs`, `watch`                |

Runner port adds `claim({ launch, runnerId, schemaVersion })`, `putBlobs`, `begin`, `commit`, `renew`, `release`,
`headRecords`, `segments`, `lastCompaction(pointer)`. Every operation is one storage transaction; keyed ones replay by
key (the same key with different input is `invalid`). A commit is refused unless the lease is the runner's, the step is
the pending one and each move starts from its pointer's current head.

Operations that create work (`submit`, `stop`, `resume`, `reset`, `check`, `release`) return `{ seq, launch? }`: a
**launch key** when no runner is on the harness. Whoever gets it hands it to a `Supervisor.ensure(harnessId, launch)`,
which starts one runner per key; `runner.run(launch)` claims and serves.

**Errors**: one `HarnessError { code, retry, detail }`; codes `invalid`, `head-changed`, `inbox-changed`, `stale` (lease
lost), `stale-schema`, `pointer-busy`, `unknown-outcome` (retry with the same key), `integrity`. A model or tool
attempt's failure is not an operation error: it is an
`AttemptFailure { code, message, transient, effect, retryAfterMs? }` in the step's evidence, and it opens a `failure`
decision.

## The runner

`createRunner(options: RunnerOptions)`: `port`, `runnerId`, `models: ModelBindings` (`resolve(target) → adapter`),
`tools: ToolBindings` (`resolve(definition) → handler`), `drive`, `configure(items) → TurnConfiguration`, optional
`price(configuration, usage)`, `clock`, `ids`, `limits { signalEveryMs, handOffAfterMs, depth = 4, cacheBytes }`.
`run(launch)` returns a `RunnerExit` (`released`, `handed-off`, `lost`, `stale-launch`, `stale-schema`, `aborted`);
`observe()` streams `Observation`s; `abort()` is local, not a durable Stop.

The loop (`runner/index.ts` `serve`):

1. **Claim** with the launch key (lease `leaseMs`). Every operation renews it; `renew` goes out only after about a third
   of the lease in silence. Only a `stale` refusal or expiry loses the lease, never a failed read.
2. **No active turn**: ask `next-turn` (payload: context size, last compaction; the inbox is in the context). `begin`
   takes items, an optional compaction edit, and `configure(items)` gives the turn's model, tool set and limits; all in
   one commit. `release` when the drive answers so.
3. **Model step** (`runner/model.ts`): project the `main` head through the adapter with the available tools, admit
   (estimate + reserved output within the window, else `capacity`; budget, else `budget`), stream. With tool calls it
   settles as evidence only and names a tools step carrying the **plan** (each call bound to the definition it was
   planned against). Without calls the reply is appended to `main` and the turn completes at the boundary.
4. **Tools step** (`tools/step.ts`): `repair` unbound calls, `permission` for each call not preauthorized, `schedule`
   into stages, run handlers (retries per `failure`), then append in **one move**: the assistant message, one user
   message with every result, then contributions.
5. **Boundary** between rounds: `boundary` (proceed, steer waiting items in, compact, or end when completing).
6. Each step ends in one commit with an ending: `settled`, `parked` (permission answered `park`; the runner releases),
   `interrupted` (an interrupt item), `abandoned` (Stop, `fail`, time). Within a step: stop > fail > park > interrupt. A
   crash loses the step and it re-runs (`prior: 'possible'` on its calls).

A turn runs in memory; the store is touched at commits. Everything but the Execution and pointers is immutable, so the
runner caches freely (`WorkingSet`).

## Decisions (`contracts/drive.ts`)

A `Decision { kind, path, payload }` (`path`: call ids down to the nested run that asks) goes to the innermost drive,
then outward. A `Drive.decide` returns an answer or `'pass'`; `drive(handlers, fallback)` composes. `'stop'` and
`{ fail }` are legal for every kind; a drive that throws (`drive-error`), passes all the way out or answers illegally
(`legal()`, `invalid-answer`) fails the run: the turn at depth 0, only the nested run inside a call.

| Kind           | Payload                                          | Answers                                                                                                | `defaultDrive`                                                                         |
| -------------- | ------------------------------------------------ | ------------------------------------------------------------------------------------------------------ | -------------------------------------------------------------------------------------- |
| `permission`   | `call`                                           | `'allow'`, `{ deny: text }`, `{ park: { note } }`                                                      | allow                                                                                  |
| `schedule`     | `calls`                                          | `{ stages: CallId[][], maxConcurrency? }` (every call once)                                            | one stage per call                                                                     |
| `budget`       | `dimension`, `budget`, `need?`                   | `{ extend: Limits }`                                                                                   | fail                                                                                   |
| `failure`      | `call?`, `attempt?`, `failure`, `unstarted?`     | `{ retry }`, `{ replace: config }` (model), `{ report: { skip? } }` (tool), `{ compact }` (projection) | transient: retry ×3 with back-off (not before `retryAfterMs`); tool: report; else fail |
| `model-result` | `finish` (`output-limit`, `filtered`), `message` | `'accept'`, `{ replace: config }`                                                                      | accept                                                                                 |
| `repair`       | `call`, `reason`, `unstarted`                    | `{ report: { parts?, contributions?, skip? } }`                                                        | report                                                                                 |
| `boundary`     | `size`, `completing`, `lastCompaction?`          | `'proceed'`, `{ steer: keys }`, `{ compact: HeadEdit }`, `'end'` (completing)                          | steer waiting items when completing                                                    |
| `interrupt`    | `item`, `stepKind`                               | `'abandon'`, `'skip'` (tools step only)                                                                | skip tools / abandon model                                                             |
| `next-turn`    | `size`, `lastCompaction?`                        | `{ begin: { keys, compact? } }`, `'release'`                                                           | begin with every waiting item                                                          |

`failure` without `call` and with `attempt` is a model attempt; without `attempt` it is a projection failure
(`capacity`, `unrepresentable`, or `provider`: no adapter, or `project` threw). `skip { text, calls? }` names only
`unstarted` calls (default: all). The `DriveContext` gives `harnessId`, `turnId`, `depth`, `budget`, `head`,
`configuration`, `inbox` (depth 0), `signal`, `journal.read(scope)` and `run(spec)` (a nested run inside the open step;
how a drive compacts with a model).

## Tools, as the harness sees them

- A turn's `ToolSet { enabled, preauthorized }`. A tool is **available** when enabled and loaded in the head: a
  `tool-definition` part with `deferred: false` (`deferred: true` only declares it). Each call is planned against the
  last loaded definition before it; an unavailable tool, a missing handler or non-object arguments make it **unbound**
  (`repair`, reasons `unknown-tool`, `not-loaded`, `not-enabled`, `no-handler`, `invalid-arguments`). Schema validation
  is the handler's.
- **Handler**: `ToolBindings.resolve(definition) → ToolHandler`, `invoke(RunContext)` yields `progress` and `usage`
  frames, then exactly one `result { status, parts, contributions?, effect?, transient? }`. A frame after the result or
  no result is a `stream-protocol` failure.
- **`RunContext`**: `harnessId`, `turnId`, `callId`, `path`, **`key`** (names the intended effect: the call id at top
  level, saved before any effect; a digest of parent key, tool, arguments and count in a nested run), **`prior`**
  (`'possible'` on a retry or re-run: a non-idempotent tool looks for its earlier effect), `attempt`, `definition`,
  `arguments`, `signal`, `journal`, `run(spec)`.
- **`journal`**: `read(scope)` any pointer (sees this step's staged moves); `append`/`edit` on a pointer the tool holds,
  never `main` (`pointer-busy`); staged until the step commits.
- **Contributions**: messages a result adds after the results message (a skill's guidance, tool definitions); one
  identical to an entry already in the head is not appended again.

## Nested runs (`nested/index.ts`)

`context.run(spec: RunSpec)` from a handler or a drive: `input`, `tools` (with their definitions), `model`, `limits`,
optional `drive`, and where the context starts (`pointer`: continued and moved; `from`: a fork of a head; neither:
empty). Same engine, its own budget (not charged to the parent's limits, but in the turn's usage), its calls recorded
under the parent call (`parentCallId`), nothing durable until the parent step commits: a re-run parent re-runs it. Not
steerable; depth limit 4 (a deeper run fails with `reason: 'depth'`); it may not write `main`. Ends `completed` with its
last `output` or `failed`; a stop or a park ends the whole step.

## Budgets and projection

- `Limits` over `modelCalls`, `toolCalls`, `uncachedInput`, `cacheRead`, `cacheWrite`, `output`, `reasoning`, `cost`
  (with the runner's `price`), `endBy` (a time). Consumption per turn and per nested run; each ended attempt is charged
  once, estimated where usage is missing. A short dimension opens `budget`.
- Projection is pure (`adapter.project({ messages, tools, configuration, capabilities })`); appending keeps the
  projected prefix byte-stable so provider caches hold (see [models.md](models.md)).

## Observation

`runner.observe()`: `model.delta`, `tool.progress`, `attempt.started`, `attempt.ended`, `decision.opened`,
`decision.answered`, and `event` (each committed event). Best effort, ids not payloads. `client.watch()` is
level-triggered: it emits the current `seq`, then each change, so a missed emission loses nothing; re-read on each.

## Storage (`contracts/storage.ts`)

`HarnessStorage`: `transact(harnessId, run(tx))` (serializable per harness; all writes or none), `putBlobs`/`getBlobs`
(outside transactions), `watch(harnessId)`, `drop`. `StorageTx` has record stores (`execution`, `messages`, `pointers`,
`heads`, `segments`), the `events` log with its indexes (`after`, `byStep`, `byTurn`, `endedTurns`, `byKey`,
`lastCompaction`) and `store` (tombstones). The host holds all protocol logic; a backend only stores and indexes.
Foreign JSON (arguments, results) is stored as an opaque string.

Backends in `@molstar/ai-components`: `createMemoryStorage()` (standalone: a serial queue over cloned records) and
`createIndexedDbStorage(name)` (one database; other tabs hear writes on a `BroadcastChannel`). IndexedDB sits on
`createDurableStorage(backend: DurableBackend)`: a transaction reads through its staged writes, then applies them
atomically guarded by the harness's write version, and re-runs on a conflict. Nothing may be awaited inside an IndexedDB
transaction, so hashing happens before. A new backend implements `DurableBackend` (or `HarnessStorage` directly) and
runs the conformance suite ([testing.md](testing.md)). SQLite is planned, not built.

`createInProcessSupervisor({ client, startRunner(launch), clock?, marginMs = 1000, …SupervisorPolicy })`
(`components/src/harness/supervisor/in-process.ts`): starts a runner per launch key, watches the Execution, and calls
`check()` when a launch window or lease runs out, so a turn whose owner died is taken up unasked.

## Schema version

`SCHEMA_VERSION = 1`. A claim from another version blocks the turn (`blocked: 'stale-schema'`); stop or reset it. There
is no migration.

## Gotchas

- `hashMessage` and `digest` are async (WebCrypto). Tests on a simulated clock must let hashing settle
  ([testing.md](testing.md)).
- `main` is the harness's: tools and nested runs cannot write it; reset and compaction move it.
- Don't read "what runs next" from heads, and don't put bookkeeping in `Message`s: a `label` says what an entry is,
  never what happened to it (that is the event log's).
- A `HeadEdit` is `{ entries: (EntryId | Message)[] }`: kept hashes and new messages, in order; it must stay
  projectable.
- Stop wins over a settlement in flight; a parked turn holds no runner until `resume`.
