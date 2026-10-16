# The harness

The harness runs agent turns in memory: a model is called, the tools it asks for run, and the loop goes on until the
turn ends or parks. It is generic: it knows nothing of Mol\* or Prism, which hosts it (§18). It is one loop in the page
that calls it, one turn at a time, and keeps nothing beyond that page's memory: what an app needs to drive itself, and
no more (decision 9, [README](README.md)).

Code comments cite this document's sections (`harness.md §6`); the numbering is fixed, so §10, §14 and §15 stay as
one-line sections. The working spec, with the files and signatures, is
[packages/ai/\_spec/harness.md](../../packages/ai/_spec/harness.md).

## 1. Parts and homes

| Part                      | Owns                                                                     | Where                                  |
| ------------------------- | ------------------------------------------------------------------------ | -------------------------------------- |
| Harness (`createHarness`) | the state in memory; runs steps, binds and runs tools, asks the drive    | `harness/` of `@molstar/ai`            |
| Drive                     | every choice                                                             | the app; `defaultDrive` in `harness/`  |
| Provider                  | projection and streaming for one provider's API, over the AI SDK         | `providers/`, behind the port of `messages/` |
| Tool handler              | one tool's effects                                                       | the app (`tools/` makes handlers)      |

`messages/` holds the shapes both sides share: `Message` and its parts, tool definitions, `ModelConfiguration`, the
`Provider` port, `hashMessage` and `digest`.

Each concern has one home and is kept nowhere else: what a model sees is the **journal** (§2); evidence (requests,
attempts, results, usage) is the **event log** (§3); what runs next is read from the last step event, never from a head;
a call's authority is the definition it was planned against (§8); every choice is the drive's (§7).

## 2. The journal

- An **entry** is a `Message`, addressed by `hashMessage`: `@molstar/core`'s MurmurHash3_x86_128 over a structural byte
  encoding (never JSON text), synchronous. Identical messages are one entry; entries are immutable and never deleted
  while the harness exists.
- A **head** is an immutable list of entry hashes, addressed by their digest. A **pointer** names a current head
  (`main` is the turn's; tools and nested runs hold others). A pointer moves only at the end of a step, by an
  **append** or an **edit** (kept hashes plus new messages: compaction).
- Every head is **projectable**: each tool call has its result in the next user message. An assistant message with calls
  enters a head only together with its results.
- A **contribution** identical to an entry already in the head is not appended again (a skill's guidance is kept once).
- A message's optional **label** (`{ kind, key?, data? }`) says what the entry is, for the code that reads the journal
  (the session, compaction, the trace). It is part of the entry (hashed: the same text under another label is another
  entry), never sent to a provider, never counted in a size, and never read by the harness. It says only what is true
  of the content; what happened when (a turn, a step, a time) is the event log's.

## 3. The event log

The only evidence store, numbered by `seq`. A **step event** records how a step ended: its calls and their attempts
(times, outcomes, usage; nested calls by `parentCallId`), its moves, its compactions, and its `next` step or the turn's
`end`. **Control events** record what ends or goes on outside a step: a parked turn `stopped` or `resumed`, a `reset`.
Turns have no records of their own: the event that ends one carries its summary. Large bodies are blobs (§4); a model
call's evidence is what reproduces it (head, configuration, tools, provider revision, payload digest), not its payload.

## 4. State and the Execution

- The **state** is one object in memory: entries, heads, pointers, the event log, blobs and the active turn. A
  **blob** (a call's request or result, a tools step's assistant message) is kept by reference, `digest(body)`, and
  written with its step.
- The **Execution** is the reader's snapshot: `seq`, `main`, and the active turn (its pending step, park, budget,
  usage).
- One loop runs in the page that calls it, one turn at a time; there is no runner, lease or supervisor apart from it.

## 5. Turns

- At most one active turn (running or parked). `turn({ messages, configuration })` begins one: its **begin step** asks
  the drive `turn-start` (where a drive compacts first), then appends the messages to `main`. A turn ends
  `completed`, `failed`, `stopped` or `reset`.
- **Stop** aborts the running step, which is recorded `abandoned` with its attempts; a parked turn ends at once. Stop
  wins over a settlement in flight.
- **Reset** ends any turn (a running step's writes are dropped) and moves pointers back (by default `main` to its head
  at the end of the last completed turn, the empty head before any); it deletes nothing and never runs on its own.
- **Park**: a `permission` answered `park` parks the step; `resume` re-runs the step and asks the drive again. A
  top-level call that succeeded before the park (a park inside a nested run ends the whole step) and holds no pointer
  keeps its result and does not run again; every other call that had started runs again with `prior: "possible"`.
- A prompt that arrives during a turn is the host's to hold until the turn ends: each prompt is its own turn (§10).

## 6. Steps

- A **round** is a model step, plus a tools step for its calls. The model step projects, admits (window and budget) and
  streams; with calls it settles as evidence only and names the tools step, without calls it appends to `main`. The
  tools step repairs unbound calls, asks `permission`, schedules, runs, then appends the assistant message, one user
  message with every result, and the contributions, in one move.
- A step builds its writes (entries, heads, blobs, its event) as it goes and applies them at its end, in one
  synchronous assignment: a stop, a failure or a reset mid-step leaves every pointer where it was.
- **Endings**: `settled`, `parked`, `abandoned` (Stop, `fail`, time). Within a step, stop wins over fail and park.
- **Retries** wait in process: a model retry sends the same request, a tool retry the same key with `prior: "possible"`.

## 7. Decisions and drives

A decision is `{ kind, path, payload }` (`path`: the call ids down to the nested run that asks). The kinds:
`permission`, `schedule`, `budget`, `failure`, `model-result`, `repair`, `boundary`, `turn-start`. It goes to the
innermost drive, then outward; the app's drive answers or passes to its fallback (`drive(handlers, defaultDrive())`).
`stop` and `fail` are legal everywhere; a drive that throws fails the turn. `defaultDrive` answers every kind (allow,
one stage per call, retry transient failures up to 3 times, accept, proceed…), so an app writes only its overrides.

## 8. Tools

- A turn's tool set names the `enabled` and `preauthorized` tools; a tool is available when enabled and loaded in the
  head (a `tool-definition` part). Each call is **planned** against the last loaded definition before it; an unavailable
  tool, a missing handler or non-object arguments make it **unbound** (`repair`). Schema validation is the handler's, so
  it can normalize what a small model sends.
- **The handler** (`ToolBindings.resolve(definition) → invoke(context)`) yields `progress` and `usage` frames, then
  exactly one `result` (`succeeded` or `failed`, parts, optional contributions, `effect`, `transient`).
- **Keys**: a call's `key` names its intended effect (the call id at top level; a digest in a nested run), and
  `prior: "possible"` says an earlier attempt may have run (a retry, or a run before the step parked), so a tool that
  is not idempotent can look for its effect.
- **A tool's own context**: `journal.read` any pointer; `append` and `edit` on one it holds (never `main`), staged until
  the step ends.

## 9. Nested runs

A handler or a drive calls `run(spec)` to run an agent inside its call: its own input, tools, model, budget and optional
drive, starting from a pointer, a head or empty. The same engine, its own budget (never charged to the parent's limits),
its calls recorded under the parent call, nothing written until the parent step ends; a retried parent call runs it
again. Depth 4.

## 10. Inbox and steering

None: no steering, interrupts or withdrawal; each prompt is its own turn (§5).

## 11. Projection, caching and budgets

- Projection is pure (head, configuration, available tools), and appending keeps the projected prefix byte-stable, so
  provider caches hold. Providers report the request's `prefix` so a trace can say where a cache broke.
- A request is admitted only if its estimate plus the reserved output fits the window (else `capacity`).
- **Budgets** are limits plus consumption per turn and per nested run, over calls, tokens by kind, `cost` (with an app
  price) and `endBy`. Every ended attempt is charged once; estimates where usage is missing. Admission reserves what
  it admits until the attempt is charged, so parallel calls cannot overshoot a limit; `endBy` passing is asked about
  once per step, and the step settles only once that is answered.

## 12. Operations

The harness's operations are `turn`, `stop`, `resume` and `reset` (§5), `observe` (§13), and the reads (`execution`,
`heads`, `pointers`, `entries`, `events`, `steps`, `turns`, `blobs`, `watch`). Each checks its precondition and refuses
with `invalid` (§16): a turn while one is active or with no configuration, a resume with none parked, a reset to an
unknown head, a read of an unknown head, entry or blob, an edit that leaves a head unprojectable. `stop` without an
active turn does nothing.

## 13. Observation

`harness.observe()` emits live observations (`model.delta`, `tool.progress`, attempts, decisions) and, once committed,
events, as each is written. Live ones are best effort and carry ids, not payloads. `watch()` is level-triggered (it says
only that `seq` moved), so a missed emission loses nothing. Views (the trace, cost) are projections that store ids.

## 14. Storage

None: the state is the one object of §4, and a reload ends the conversation.

## 15. Remote deployments

None: a reader may be remote (the reads are async), the loop may not.

## 16. Errors

One operation error, `HarnessError { code, detail }`: `invalid` (§12) and `pointer-busy` (a tool writing `main` or
another's pointer). An attempt's failure is not an operation error: it has a code, a `transient` hint, an `effect`
(`none`, `applied`, `unknown`) and, when the provider said, `retryAfterMs`, and opens `failure`.

## 17. Out of scope

Persistence (a conversation that outlives the page, several runners: [plan.md](plan.md) "Later"), more than one active
turn, agent profiles, UI stores, compaction strategies (the harness has the edit; the session has the strategy),
garbage collection of entries.

## 18. How Prism hosts it

The Prism viewer's agent is built ([../prism-ui/spec/agent.md](../prism-ui/spec/agent.md); the design:
[../prism-ui/design/contexts.md](../prism-ui/design/contexts.md) §4–5,
[../prism-ui/design/agent.md](../prism-ui/design/agent.md)):

- **Tools are the viewer's operations**: a PD param set becomes a tool's field shapes (`paramsShape`), and a tool runs
  the operation as a person's command would, with the origin `agent`.
- **The drive** is the session's, with the viewer's rules: approval through the Assistant (park, then resume with the
  person's answer), reads scheduled together and writes one at a time, compaction at the boundary.
- **One undo step per turn**: the session's `onTurn` opens one history group, which the turn's changes join.
- **Storage**: none for the harness; the viewer keeps its own record of the conversation as a thread
  ([../../packages/ai/\_spec/threads.md](../../packages/ai/_spec/threads.md)). Model keys stay in the host.
