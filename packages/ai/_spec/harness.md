# The harness (`src/harness/` of `@molstar/ai`)

Runs agent turns in memory: the model is called, the tools it asks for run, and the loop goes on until the turn ends or
parks. Mechanism only; every choice is the drive's. One loop in the page that calls it, one turn at a time, nothing kept
beyond that page's memory. Design and rationale: [.v6/ai/harness.md](../../../.v6/ai/harness.md) (the `§` numbers code
comments cite are its sections).

## Files

| Path (`src/harness/`)  | What                                                                                                                                                            |
| ---------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `contracts/records.ts` | the records: `Execution`, `ActiveTurn`, `Head`, `Pointer`, events, `Usage`, `Limits`                                                                            |
| `contracts/harness.ts` | `Harness`, `HarnessReader`, `HarnessOptions`, `TurnHandle`, `TurnOutcome`, `HeadEdit`, `HarnessError`                                                           |
| `contracts/drive.ts`   | `Decisions` (each kind's payload and answers), `Drive`, `DriveContext`, `DriveHandlers`                                                                         |
| `contracts/runner.ts`  | `Clock` (re-exported from `messages/`), `IdSource`, `ToolBindings`/`ToolHandler`/`ToolFrame`, `ToolJournal`, `RunContext`, `RunSpec`                            |
| `contracts/observe.ts` | `HarnessObservable`, `Observation`                                                                                                                              |
| `loop/`                | `createHarness` (`index.ts`: turns, steps, writes, reads), `model.ts` (model step), `run.ts` (shared step ops), `state.ts` (`Store`, `Collector`, exits, usage) |
| `tools/`               | `plan.ts` (bind calls to definitions), `step.ts` (the tools step), `journal.ts` (a tool's journal)                                                              |
| `nested/`              | `nest`: nested runs inside a call                                                                                                                               |
| `drive/`               | `defaultDrive`, `drive(handlers, fallback)`, `legal`                                                                                                            |
| `observe/`             | `Subject`, the hot observable behind `observe()` and `watch()`                                                                                                  |

## The data

- **Entry**: a `Message` (`messages/`) addressed by `hashMessage` (`EntryId`, 32 hex digits). Immutable, deduplicated,
  never deleted while the harness exists. Its optional `label` (what it is: a request, guidance) is part of it; the
  harness stores it and never reads it, and never counts it in a size (`Store.tokensOf`).
- **Head**: an immutable list of entry ids (`Head { id, entries }`), `id` the digest of the list.
- **Pointer**: `scope → head`. `MAIN` (`'main'`) is the turn's conversation and lives on the Execution; tools and nested
  runs hold others. A pointer moves only at the end of a step.
- **Projectable**: every head has each tool call answered in the next user message; an assistant message with calls
  enters a head only with its results.
- **Event log**: the only evidence, numbered by `seq` (from 1). `StepEvent` (how a step ended: `calls` with attempts and
  usage, `moves`, `compaction`, `next` step or `end`, `summary` on the ending one, `begun` on a begin step) and
  `ControlEvent` (`stopped`, `resumed`, `reset`). Turns have no record of their own: `turns()` reads them from the log.
  Large bodies (a call's request, its result, a tools step's assistant message) are **blobs**, `ref = digest(body)`,
  written with the step.
- **Execution**: the reader's snapshot: `seq`, `main` (head id), `turn?: ActiveTurn` (configuration, budget, `pending`
  step, `park`, usage).
- **What runs next** is `execution.turn.pending` (from the last step event), never derived from a head.
- **State** is one object in memory (`loop/state.ts` `Store`): entries, heads, pointers, the event log, blobs, the
  active turn. A step collects its writes (`Collector`) and applies them at its end in one synchronous assignment: a
  stop, a failure or a reset mid-step leaves the journal as it was.

## The harness

`createHarness({ models, tools, drive, configure?, price?, clock?, ids?, limits: { depth = 4 }?, id? })` → `Harness`:
`models: ProviderBindings` (`resolve(target) → provider`, asked once per configuration in each `turn` and each `resume`,
so a provider the host replaced between them is the one used), `tools: ToolBindings` (`resolve(definition) → handler`),
the `drive`, `configure(messages)` (a turn's configuration when `turn` gives none), `price(configuration, usage)`.

| Operation                            | What                                                                                                                                                                               |
| ------------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `turn({ messages, configuration? })` | begins a turn; refused (`invalid`) while one is active. Returns `{ turnId, done }`; `done` resolves with the `TurnOutcome`: `{ turnId, end: TurnEnd }` or `{ turnId, parked }`     |
| `stop(reason?)`                      | the running step is aborted and recorded `abandoned`, end `stopped`; a parked turn ends at once (`stopped`)                                                                        |
| `resume()`                           | re-runs a parked turn's step, asking the drive again (`resumed`); calls that succeeded before the park keep their results (below); refused (`invalid`) unless one is parked        |
| `reset(heads?, reason?)`             | ends any turn as `reset` (a running step's writes are dropped) and moves pointers: the given heads, else `main` to its head at the end of the last completed turn; deletes nothing |
| `observe()`                          | live observations and each event as it is written                                                                                                                                  |
| reads (`HarnessReader`)              | `id`, `execution`, `heads`, `pointers`, `entries`, `events`, `steps`, `turns`, `blobs`, `watch`                                                                                    |

The reads are async (a reader may be remote: the trace and the Inspector take a `HarnessReader`). `emptyHead()` is the
empty head (the same id in every harness; `forget` resets `main` to it).

**Errors**: one `HarnessError { code, detail }`; codes `invalid` (a turn while one is active or without a configuration,
a resume with none parked, an unknown head, entry or blob, an unprojectable edit) and `pointer-busy` (a tool writing
`main` or another writer's pointer). A model or tool attempt's failure is not an operation error: it is an
`AttemptFailure { code, message, transient, effect, retryAfterMs? }` in the step's evidence, and it opens a `failure`
decision.

## The loop (`loop/index.ts`)

1. **Begin step**: ask `turn-start` (payload: context size, last compaction); a `compact` edit applies first, then the
   turn's messages are appended to `main`. One step event (`stepKind: 'begin'`, `begun.configuration`).
2. **Model step** (`loop/model.ts`): project the `main` head through the provider with the available tools, admit
   (estimate + reserved output within the window, else `capacity`; budget, else `budget`), stream. With tool calls it
   settles as evidence only and names a tools step carrying the **plan** (each call bound to the definition it was
   planned against). Without calls the reply is appended to `main` and the turn completes at the boundary.
3. **Tools step** (`tools/step.ts`): `repair` unbound calls, `permission` for each call not preauthorized, `schedule`
   into stages, run handlers (retries per `failure`), then append in **one move**: the assistant message, one user
   message with every result, then contributions.
4. **Boundary** between rounds: `boundary` (proceed, compact, or end when completing).
5. Each step ends with one event and an ending: `settled`, `parked` (permission answered `park`; `done` resolves
   `parked`), `abandoned` (Stop, `fail`, time; an unexpected error fails the turn too). Within a step: stop > fail >
   park.

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
| `boundary`     | `size`, `completing`, `lastCompaction?`          | `'proceed'`, `{ compact: HeadEdit }`, `'end'` (completing)                                             | proceed                                                                                |
| `turn-start`   | `size`, `lastCompaction?`                        | `'proceed'`, `{ compact: HeadEdit }`                                                                   | proceed                                                                                |

`failure` without `call` and with `attempt` is a model attempt; without `attempt` it is a projection failure
(`capacity`, `unrepresentable`, or `provider`: no provider, or `project` threw). `skip { text, calls? }` names only
`unstarted` calls (default: all). The `DriveContext` gives `harnessId`, `turnId`, `depth`, `budget`, `head`,
`configuration`, `signal`, `journal.read(scope)` and `run(spec)` (a nested run inside the open step; how a drive
compacts with a model).

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
  level; a digest of parent key, tool, arguments and count in a nested run), **`prior`** (`'possible'` on a retry, or
  when the call had started before its step parked: a non-idempotent tool looks for its earlier effect), `attempt`,
  `definition`, `arguments`, `signal`, `journal`, `run(spec)`.
- **`journal`**: `read(scope)` any pointer (sees this step's staged moves); `append`/`edit` on a pointer the tool holds,
  never `main` (`pointer-busy`); staged until the step ends.
- **Contributions**: messages a result adds after the results message (a skill's guidance, tool definitions); one
  identical to an entry already in the head is not appended again.

## Nested runs (`nested/index.ts`)

`context.run(spec: RunSpec)` from a handler or a drive: `input`, `tools` (with their definitions), `model`, `limits`,
optional `drive`, and where the context starts (`pointer`: continued and moved; `from`: a fork of a head; neither:
empty). Same engine, its own budget (not charged to the parent's limits, but in the turn's usage), its calls recorded
under the parent call (`parentCallId`), nothing written until the parent step ends. Depth limit 4 (a deeper run fails
with `reason: 'depth'`); it may not write `main`. Ends `completed` with its last `output` or `failed`; a stop or a park
ends the whole step. On resume the parked step runs again: a top-level call that succeeded before the park and holds no
pointer keeps its result and is neither asked about nor run again; every other call that had started (the parked call, a
sibling the park cut short, a writer whose writes went with the step, their nested calls) runs again with
`prior: 'possible'` (harness.md §5).

## Budgets and projection

- `Limits` over `modelCalls`, `toolCalls`, `uncachedInput`, `cacheRead`, `cacheWrite`, `output`, `reasoning`, `cost`
  (with the harness's `price`), `endBy` (a time). Consumption per turn and per nested run; each ended attempt is charged
  once, estimated where usage is missing. A short dimension opens `budget`. Admission reserves what it admits until the
  attempt is charged, so calls running in parallel cannot overshoot a limit (a decision's `need` includes what the
  others hold). `endBy` passing opens one `budget` decision, shared by the step's watcher and its admissions; the step
  settles only once it is answered, so its `extend` is the step's.
- Projection is pure (`provider.project({ messages, tools, configuration, capabilities })`); appending keeps the
  projected prefix byte-stable so provider caches hold (see [models.md](models.md)).

## Observation

`harness.observe()`: `model.delta`, `tool.progress`, `attempt.started`, `attempt.ended`, `decision.opened`,
`decision.answered`, and `event` (each event as it is written). Best effort, ids not payloads. `watch()` is
level-triggered: it emits the current `seq`, then each change, so a missed emission loses nothing; re-read on each.

## Gotchas

- `hashMessage` and `digest` are synchronous: `@molstar/core`'s MurmurHash3_x86_128 (seed 42) over UTF-8 bytes, an id
  within one page's memory, not a cryptographic hash. A crafted message (text from a file in a tool result) could
  collide with an earlier entry and be taken for it; low for a single person's viewer session.
- `main` is the harness's: tools and nested runs cannot write it; reset and compaction move it.
- Don't read "what runs next" from heads, and don't put bookkeeping in `Message`s: a `label` says what an entry is,
  never what happened to it (that is the event log's).
- A `HeadEdit` is `{ entries: (EntryId | Message)[] }`: kept hashes and new messages, in order; it must stay
  projectable.
- Stop wins over a settlement in flight; a parked turn holds nothing until `resume` (or `stop`, `reset`).
- Nothing outlives the page: a reload ends the conversation.
