# The harness (`src/harness/` of `@molstar/ai`)

Runs agent turns in memory: the model is called, the tools it asks for run, and the loop goes on until the turn ends or
parks. Mechanism only; every choice is the drive's. Generic: it knows nothing of Mol\* or Prism, which hosts it (§18).
One loop in the page that calls it, one turn at a time, nothing kept beyond that page's memory.

Code comments cite the numbered sections (`harness.md §6`); the numbering is fixed, so §10, §14 and §15 stay as one-line
sections.

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

## 1. Parts and homes

| Part                      | Owns                                                                  | Where                                        |
| ------------------------- | --------------------------------------------------------------------- | -------------------------------------------- |
| Harness (`createHarness`) | the state in memory; runs steps, binds and runs tools, asks the drive | `harness/`                                   |
| Drive                     | every choice                                                          | the app; `defaultDrive` in `harness/`        |
| Provider                  | projection and streaming for one provider's API, over the AI SDK      | `providers/`, behind the port of `messages/` |
| Tool handler              | one tool's effects                                                    | the app (`tools/` makes handlers)            |

One home per concern, kept nowhere else: what a model sees is the **journal** (§2); evidence (requests, attempts,
results, usage) is the **event log** (§3); what runs next is the active turn's pending step (§4), never read from a
head; a call's authority is the definition it was planned against (§8); every choice is the drive's (§7).

## 2. The journal

- **Entry**: a `Message` (`messages/`) addressed by `hashMessage` (`EntryId`, 32 hex digits). Immutable, deduplicated,
  never deleted while the harness exists.
- **Label** (`{ kind, key?, data? }`, optional): what the entry is (a request, guidance), for the code that reads the
  journal (the session, compaction, the trace). Part of the entry (hashed: the same text under another label is another
  entry); never sent to a provider, never counted in a size (`Store.tokensOf`), never read by the harness. It says only
  what is true of the content; what happened when is the event log's.
- **Head**: an immutable list of entry ids (`Head { id, entries }`), `id` the digest of the list.
- **Pointer**: `scope → head`. `MAIN` (`'main'`) is the turn's conversation and lives on the Execution; tools and nested
  runs hold others. A pointer moves only at the end of a step, by an append or an edit (`HeadEdit`: kept hashes plus new
  messages; compaction).
- **Projectable**: every head has each tool call answered in the next user message; an assistant message with calls
  enters a head only with its results.
- A **contribution** identical to an entry already in the head is not appended again (a skill's guidance is kept once).

## 3. The event log

The only evidence, numbered by `seq` (from 1). `StepEvent`: how a step ended (`calls` with attempts and usage, nested
calls by `parentCallId`; `moves`; `compaction`; `next` step or `end`; `summary` on the ending one; `begun` on a begin
step). `ControlEvent`: `stopped`, `resumed` (a parked turn), `reset`. Turns have no record of their own: `turns()` reads
them from the log. A model call's evidence is what reproduces it (head, configuration, tools, provider revision, payload
digest), not its payload.

## 4. State and the Execution

- **State**: one object in memory (`loop/state.ts` `Store`): entries, heads, pointers, the event log, blobs, the active
  turn. A **blob** (a call's request or result, a tools step's assistant message) is kept by reference,
  `ref = digest(body)`, written with its step.
- **Execution**: the reader's snapshot: `seq`, `main` (head id), `turn?: ActiveTurn` (configuration, budget, `pending`
  step, `park`, usage).
- **What runs next** is `execution.turn.pending` (from the last step event), never derived from a head.
- One loop, in the page that calls it; no runner, lease or supervisor apart from it.

## 5. Turns

- At most one active turn (running or parked). `turn({ messages, configuration? })` begins one (refused, `invalid`,
  while one is active) and returns `{ turnId, done }`; `done` resolves with the `TurnOutcome`:
  `{ turnId, end: TurnEnd }` or `{ turnId, parked }`. A turn ends `completed`, `failed`, `stopped` or `reset`.
- **Stop** (`stop(reason?)`): the running step is aborted and recorded `abandoned` with its attempts, end `stopped`; a
  parked turn ends at once. Stop wins over a settlement in flight.
- **Reset** (`reset(heads?, reason?)`): ends any turn as `reset` (a running step's writes are dropped) and moves
  pointers: the given heads, else `main` to its head at the end of the last completed turn (the empty head before any).
  Deletes nothing; never runs on its own.
- **Park**: a `permission` answered `park` parks the step (`done` resolves `parked`). `resume()` re-runs the step and
  asks the drive again (`resumed`); refused unless one is parked. A top-level call that succeeded before the park and
  holds no pointer keeps its result and is neither asked about nor run again (a park inside a nested run ends the whole
  step); every other call that had started (the parked call, a sibling the park cut short, a writer whose writes went
  with the step, their nested calls) runs again with `prior: 'possible'`.
- A prompt that arrives during a turn is the host's to hold until the turn ends: each prompt is its own turn (§10).

## 6. Steps

1. **Begin step**: ask `turn-start` (payload: context size, last compaction); a `compact` edit applies first, then the
   turn's messages are appended to `main`. One step event (`stepKind: 'begin'`, `begun.configuration`).
2. **Model step** (`loop/model.ts`): project the `main` head through the provider with the available tools, admit (§11),
   stream. With tool calls it settles as evidence only and names a tools step carrying the **plan** (each call bound to
   the definition it was planned against, §8). Without calls the reply is appended to `main` and the turn completes at
   the boundary.
3. **Tools step** (`tools/step.ts`): `repair` unbound calls, `permission` for each call not preauthorized, `schedule`
   into stages, run handlers (retries per `failure`), then append in **one move**: the assistant message, one user
   message with every result, then contributions.
4. **Boundary** between rounds: `boundary` (proceed, compact, or end when completing).

A model step plus its tools step is a **round**. A step collects its writes (`Collector`: entries, heads, blobs, its
event) and applies them at its end in one synchronous assignment: a stop, a failure or a reset mid-step leaves every
pointer where it was. **Endings**: `settled`, `parked`, `abandoned` (Stop, `fail`, time; an unexpected error fails the
turn too); within a step, stop > fail > park. **Retries** wait in process: a model retry sends the same request, a tool
retry the same key with `prior: 'possible'`.

## 7. Decisions and drives

A `Decision { kind, path, payload }` (`path`: call ids down to the nested run that asks) goes to the innermost drive,
then outward. A `Drive.decide` returns an answer or `'pass'`; `drive(handlers, fallback)` composes (an app writes only
its overrides over `defaultDrive()`). `'stop'` and `{ fail }` are legal for every kind; a drive that throws
(`drive-error`), passes all the way out or answers illegally (`legal()`, `invalid-answer`) fails the run: the turn at
depth 0, only the nested run inside a call.

| Kind           | Payload                                          | Answers                                                                                                | `defaultDrive`                                                                                             |
| -------------- | ------------------------------------------------ | ------------------------------------------------------------------------------------------------------ | ---------------------------------------------------------------------------------------------------------- |
| `permission`   | `call`                                           | `'allow'`, `{ deny: text }`, `{ park: { note } }`                                                      | allow                                                                                                      |
| `schedule`     | `calls`                                          | `{ stages: CallId[][], maxConcurrency? }` (every call once)                                            | one stage per call                                                                                         |
| `budget`       | `dimension`, `budget`, `need?`                   | `{ extend: Limits }`                                                                                   | fail                                                                                                       |
| `failure`      | `call?`, `attempt?`, `failure`, `unstarted?`     | `{ retry }`, `{ replace: config }` (model), `{ report: { skip? } }` (tool), `{ compact }` (projection) | transient: up to 3 attempts (2 retries) with back-off (not before `retryAfterMs`); tool: report; else fail |
| `model-result` | `finish` (`output-limit`, `filtered`), `message` | `'accept'`, `{ replace: config }`                                                                      | accept                                                                                                     |
| `repair`       | `call`, `reason`, `unstarted`                    | `{ report: { parts?, contributions?, skip? } }`                                                        | report                                                                                                     |
| `boundary`     | `size`, `completing`, `lastCompaction?`          | `'proceed'`, `{ compact: HeadEdit }`, `'end'` (completing)                                             | proceed                                                                                                    |
| `turn-start`   | `size`, `lastCompaction?`                        | `'proceed'`, `{ compact: HeadEdit }`                                                                   | proceed                                                                                                    |

`failure` without `call` and with `attempt` is a model attempt; without `attempt` it is a projection failure
(`capacity`, `unrepresentable`, or `provider`: no provider, or `project` threw). `skip { text, calls? }` names only
`unstarted` calls (default: all). The `DriveContext` gives `harnessId`, `turnId`, `depth`, `budget`, `head`,
`configuration`, `signal`, `journal.read(scope)` and `run(spec)` (a nested run inside the open step; how a drive
compacts with a model).

## 8. Tools

- A turn's `ToolSet { enabled, preauthorized }`. A tool is **available** when enabled and loaded in the head: a
  `tool-definition` part with `deferred: false` (`deferred: true` only declares it). Each call is **planned** against
  the last loaded definition before it; an unavailable tool, a missing handler or non-object arguments make it
  **unbound** (`repair`, reasons `unknown-tool`, `not-loaded`, `not-enabled`, `no-handler`, `invalid-arguments`). Schema
  validation is the handler's, so it can normalize what a small model sends.
- **Handler**: `ToolBindings.resolve(definition) → ToolHandler`, `invoke(RunContext)` yields `progress` and `usage`
  frames, then exactly one `result { status, parts, contributions?, effect?, transient? }`. A frame after the result or
  no result is a `stream-protocol` failure.
- **`RunContext`**: `harnessId`, `turnId`, `callId`, `path`, **`key`** (names the intended effect: the call id at top
  level; a digest of parent key, tool, arguments and count in a nested run), **`prior`** (`'possible'` on a retry, or
  when the call had started before its step parked: a non-idempotent tool looks for its earlier effect), `attempt`,
  `definition`, `arguments`, `signal`, `journal`, `run(spec)`.
- **`journal`**: `read(scope)` any pointer (sees this step's staged moves); `append`/`edit` on a pointer the tool holds,
  never `main` (`pointer-busy`); staged until the step ends.
- **Contributions**: messages a result adds after the results message (a skill's guidance, tool definitions), deduped
  (§2).

## 9. Nested runs

`context.run(spec: RunSpec)` from a handler or a drive (`nested/index.ts`): `input`, `tools` (with their definitions),
`model`, `limits`, optional `drive`, and where the context starts (`pointer`: continued and moved; `from`: a fork of a
head; neither: empty). Same engine, its own budget (not charged to the parent's limits, but in the turn's usage), its
calls recorded under the parent call (`parentCallId`), nothing written until the parent step ends; a retried parent call
runs it again. Depth limit 4 (a deeper run fails with `reason: 'depth'`); it may not write `main`. Ends `completed` with
its last `output` or `failed`; a stop or a park ends the whole step (resume: §5).

## 10. Inbox and steering

None: no steering, interrupts or withdrawal; each prompt is its own turn (§5).

## 11. Projection, caching and budgets

- Projection is pure (`provider.project({ messages, tools, configuration, capabilities })`); appending keeps the
  projected prefix byte-stable so provider caches hold. Providers report the request's `prefix`, so a trace can say
  where a cache broke ([models.md](models.md)).
- **Admission**: a request goes only if its estimate plus the reserved output fits the window (else `capacity`) and the
  budget allows it (else `budget`).
- **Budgets**: `Limits` over `modelCalls`, `toolCalls`, `uncachedInput`, `cacheRead`, `cacheWrite`, `output`,
  `reasoning`, `cost` (with the harness's `price`), `endBy` (a time); consumption per turn and per nested run. Each
  ended attempt is charged once, estimated where usage is missing. A short dimension opens `budget`. Admission reserves
  what it admits until the attempt is charged, so parallel calls cannot overshoot a limit (a decision's `need` includes
  what the others hold). `endBy` passing opens one `budget` decision, shared by the step's watcher and its admissions;
  the step settles only once it is answered, so its `extend` is the step's.

## 12. Operations

`createHarness({ models, tools, drive, configure?, price?, clock?, ids?, limits: { depth = 4 }?, id? })` → `Harness`:
`models: ProviderBindings` (`resolve(target) → provider`, asked once per configuration in each `turn` and each `resume`,
so a provider the host replaced between them is the one used), `tools: ToolBindings`, the `drive`, `configure(messages)`
(a turn's configuration when `turn` gives none), `price(configuration, usage)`.

Operations: `turn`, `stop`, `resume`, `reset` (§5), `observe` (§13), and the reads (`HarnessReader`: `id`, `execution`,
`heads`, `pointers`, `entries`, `events`, `steps`, `turns`, `blobs`, `watch`). The reads are async (a reader may be
remote: the trace and the Inspector take a `HarnessReader`). `emptyHead()` is the empty head (the same id in every
harness; `forget` resets `main` to it). Each operation checks its precondition and refuses with `invalid` (§16): a turn
while one is active or with no configuration, a resume with none parked, a reset to an unknown head, a read of an
unknown head, entry or blob, an edit that leaves a head unprojectable. `stop` without an active turn does nothing.

## 13. Observation

`harness.observe()`: `model.delta`, `tool.progress`, `attempt.started`, `attempt.ended`, `decision.opened`,
`decision.answered`, and `event` (each event as it is written). Best effort; ids and small values (a delta's text, a
progress value), not records. `watch()` is level-triggered: it emits the current `seq`, then each change, so a missed
emission loses nothing; re-read on each. Views (the trace, cost) are projections that store ids.

## 14. Storage

None: the state is §4's one object, and a reload ends the conversation.

## 15. Remote deployments

None: a reader may be remote (the reads are async), the loop may not.

## 16. Errors

One operation error, `HarnessError { code, detail }`: `invalid` (§12) and `pointer-busy` (a tool writing `main` or
another writer's pointer). A model or tool attempt's failure is not an operation error: it is an
`AttemptFailure { code, message, transient, effect, retryAfterMs? }` (`effect`: `none`, `applied`, `unknown`) in the
step's evidence, and it opens a `failure` decision.

## 17. Out of scope

Persistence (a conversation that outlives the page), more than one active turn, agent profiles, UI stores, compaction
strategies (the harness has the edit; the session has the strategy), garbage collection of entries. Open ones:
[../TODO.md](../TODO.md).

## 18. How Prism hosts it

`@molstar/prism-viewer`'s `src/agent/` (`.v6/prism-ui/spec/agent.md`):

- **Tools are the viewer's operations**: a PD param set becomes a tool's field shapes (`paramsShape`); a tool runs the
  operation as a person's command would, with the origin `agent`.
- **The drive** is the session's with the viewer's rules: approval through the Assistant (park, then resume with the
  person's answer), reads together and writes one at a time, compaction at the boundary.
- **One undo step per turn**: the session's `onTurn` opens one history group, which the turn's changes join.
- **No storage** for the harness: the viewer's record of the conversation is a thread ([threads.md](threads.md)). Model
  keys stay in the host.

## Decisions

- **Policy in the drive, mechanism in the harness.** Flow, tool execution, retries and accounting are guaranteed; every
  choice is a decision with a replaceable default, so an app writes only its overrides.
- **One home per concern** (§1): a change that moves what a model sees out of the journal, evidence out of the event
  log, or a choice out of the drive is the wrong change.
- **Content-addressed, immutable entries, no metadata bag**: bookkeeping never changes a `Message`. Code finds an entry
  by its label, never by its text or markers; a found entry's text may be read for display or carry-over.
- **One loop, one turn at a time, writes at the step's end**: a stop, failure or reset mid-step cannot leave a half
  write. No steering, no interrupts: each prompt waits for its own turn.
- **Nested runs live in their parent's step**: nothing of them is written before it ends; a retried parent runs them
  again (durable nested runs: deferred, [../TODO.md](../TODO.md)).
- **Level-triggered watch**: a reader re-reads on each `seq` it hears, so it misses nothing.
- **In memory only**: no storage port, durability, leases or supervisor; a reload ends the conversation. What an app
  keeps is its own (Prism's record is a thread). Persistence had no host using it and was a large share of the code and
  most of the test time.
- **Synchronous hashing**: entries, heads and blobs are addressed with `@molstar/core`'s MurmurHash3-128 over UTF-8
  bytes: an id within one page's memory, not a cryptographic hash, so nothing waits on WebCrypto.

## Gotchas

- `hashMessage` and `digest`: MurmurHash3_x86_128 (seed 42), `digest` over the UTF-8 bytes of canonical JSON,
  `hashMessage` over a tagged structural encoding of the message (`encodeMessage`, never JSON text). A crafted message
  (text from a file in a tool result) could collide with an earlier entry and be taken for it; low for a single person's
  viewer session.
- `main` is the harness's: tools and nested runs cannot write it; reset and compaction move it.
- Don't read "what runs next" from heads, and don't put bookkeeping in `Message`s: a `label` says what an entry is,
  never what happened to it (that is the event log's).
- A `HeadEdit` is `{ entries: (EntryId | Message)[] }`: kept hashes and new messages, in order; it must stay
  projectable.
- A parked turn holds nothing until `resume` (or `stop`, `reset`).
