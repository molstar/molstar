# The harness

The harness runs agent turns in memory: a model is called, the tools it asks for run, and the loop goes on until the
turn ends. It is generic: it knows nothing of Mol\* or Prism, which host it (§18). Code comments cite this document's
sections (`harness.md §6`); the numbering is fixed. Until 2026-10-16 it ran turns durably (storage, leases, a
supervisor, a host protocol); that went, "just basic in memory harness for now for driving an app"
([simplification/harness.md](simplification/harness.md)), and the sections below say what stays.

## 1. Parts and homes

| Part                        | Owns                                                                      | Folder of `@molstar/ai`      |
| --------------------------- | ------------------------------------------------------------------------- | ---------------------------- |
| Harness (`createHarness`)   | the state in memory; runs steps, binds and runs tools, asks the drive     | `harness/`                   |
| Drive                       | every choice                                                              | the app; defaults in harness |
| Provider                    | projection and streaming for one provider, over the AI SDK                | `providers/`                 |
| Tool handler                | one tool's effects                                                        | the app                      |

`messages/` (of `@molstar/ai`) holds the model side's shapes: `Message` and its parts, tool definitions, `ModelConfiguration`, the
`Provider` port, `hashMessage`, `digest`.

Each concern has one home and is kept nowhere else: what a model sees is the **journal** (§2); evidence (requests,
attempts, results, usage) is the **event log** (§3); what runs next is read from the last step event, never from a head;
a call's authority is the definition it was planned against (§8); every choice is the drive's (§7).

## 2. The journal

- An **entry** is a `Message`, addressed by `hashMessage` (`@molstar/core`'s MurmurHash3_x86_128 over a structural byte
  encoding, never JSON text; synchronous). Identical messages are one entry; entries are immutable and never deleted.
- A **head** is an immutable list of entry hashes, addressed by their digest. A **pointer** names a current head
  (`main` is the turn's; tools and nested runs hold others). A pointer moves only at the end of a step, by an
  **append** or an **edit** (kept hashes plus new messages: compaction); a step's writes apply together, so a stop or
  a failure mid-step leaves every pointer where it was.
- Every head is **projectable**: each tool call has its result in the next user message. An assistant message with calls
  enters a head only together with its results.
- A **contribution** identical to an entry already in the head is not appended again (a skill's guidance is kept once).
- A message's optional **label** (`{ kind, key?, data? }`) says what the entry is, for the code that reads the journal
  (a session, compaction, a trace). It is part of the entry (hashed: the same text under another label is another
  entry), never sent to a provider, never counted in a size, and never read by the harness. It says only what is true
  of the content; what happened when (a turn, a step, a time) is the event log's.

## 3. The event log

The only evidence store, numbered by `seq`. A **step event** records how a step ended: its calls and their attempts
(times, outcomes, usage; nested calls by `parentCallId`), its moves, the items it consumed, its compactions, and its
`next` step or the turn's `end`. **Control events** record what ends or goes on outside a step: a parked turn
`stopped` or `resumed`, a `reset`. Turns have no records: the event that ends one carries its summary. Large bodies are blobs; a model call's
evidence is what reproduces it (head, configuration, tools, provider revision, payload digest), not its payload.

## 4. The Execution, runners and the supervisor

- The **Execution** is the reader's snapshot: `seq`, `main`, and the active turn (its pending step, park, budget, usage).
- Runners, leases, launches, the supervisor and schema versions went (2026-10-16): the loop runs in the page that
  calls it, one turn at a time.

## 5. Turns

- At most one active turn (running or parked). `turn({ messages, configuration })` begins one: its **begin step** asks
  the drive `turn-start` (where a drive compacts first), then appends the messages to `main`. A turn ends
  `completed`, `failed`, `stopped` or `reset`.
- **Stop** aborts the running step, which is recorded `abandoned` with its attempts; a parked turn ends at once. Stop
  wins over a settlement in flight.
- **Reset** ends any turn (a running step's writes are dropped) and moves pointers back (by default `main` to its head
  at the end of the last completed turn); it deletes nothing and never runs on its own.
- **Park**: a `permission` answered `park` parks the step; `resume` re-runs the step and asks the drive again.

## 6. Steps

- A **round** is a model step, plus a tools step for its calls. The model step projects, admits (window and budget) and
  streams; with calls it settles as evidence only and names the tools step, without calls it appends to `main`. The
  tools step repairs unbound calls, asks `permission`, schedules, runs, then appends the assistant message, one user
  message with every result, and the contributions, in one move.
- A step builds its writes (entries, heads, blobs, its event) as it goes and applies them at its end, in one
  synchronous assignment.
- **Endings**: `settled`, `parked`, `abandoned` (Stop, `fail`, time). Within a step, stop wins over fail and park.
- **Retries** wait in process: a model retry sends the same request, a tool retry the same key with `prior: "possible"`.

## 7. Decisions and drives

A decision is `{ kind, path, payload }` (`path`: the call ids down to the nested run that asks). The kinds:
`permission`, `schedule`, `budget`, `failure`, `model-result`, `repair`, `boundary`, `turn-start` (until 2026-10-16
also `interrupt` and `next-turn`, and `boundary` could steer waiting items in). It goes
to the innermost drive, then outward; the app's drive answers or passes to its fallback
(`drive(handlers, defaultDrive())`). `stop` and `fail` are legal everywhere; a drive that throws fails the turn.
`defaultDrive` answers every kind (allow, serial, retry transient failures up to 3 times, accept, proceed…), so an app
writes only its overrides.

## 8. Tools

- A turn's tool set names the `enabled` and `preauthorized` tools; a tool is available when enabled and loaded in the
  head (a `tool-definition` part). Each call is **planned** against the last loaded definition before it; an unavailable
  tool, a missing handler or non-object arguments make it **unbound** (`repair`). Schema validation is the handler's, so
  it can normalize what a small model sends.
- **The handler** (`ToolBindings.resolve(definition) → invoke(context)`) yields `progress` and `usage` frames, then
  exactly one `result` (`succeeded` or `failed`, parts, optional contributions, `transient`).
- **Keys**: a call's `key` names its intended effect (the call id at top level; a digest in a nested run), and
  `prior: "possible"` says an earlier attempt may have run, so a tool that is not idempotent can look for its effect.
- **A tool's own context**: `journal.read` any pointer; `append` and `edit` on one it holds (never `main`), staged until
  the step ends.

## 9. Nested runs

A handler or a drive calls `run(spec)` to run an agent inside its call: its own input, tools, model, budget and optional
drive, starting from a pointer, a head or empty. The same engine, its own budget (never charged to the parent), its
calls recorded under the parent call, nothing written until the parent step ends. Depth 4.

## 10. Inbox and steering

Gone (2026-10-16): no inbox, steering, interrupts or `withdraw`; each prompt is its own turn (§5).

## 11. Projection, caching and budgets

- Projection is pure (head, configuration, available tools), and appending keeps the projected prefix byte-stable, so
  provider caches hold. Providers report the request's `prefix` so a trace can say where a cache broke.
- A request is admitted only if its estimate plus the reserved output fits the window (else `capacity`).
- **Budgets** are limits plus consumption per turn and per nested run, over calls, tokens by kind, `cost` (with an app
  price) and `endBy`. Every ended attempt is charged once; estimates where usage is missing.

## 12. Protocol operations

Gone (2026-10-16): the host protocol (client and runner handles, keyed idempotent operations, commit batches); the
harness's operations are `turn`, `stop`, `resume` and `reset` (§5).

## 13. Observation

`harness.observe()` emits live observations (`model.delta`, `tool.progress`, attempts, decisions) and, once committed,
events, as each is written. Live ones are best effort and carry ids, not payloads. `watch()` is level-triggered (it says only that `seq`
moved), so a missed emission loses nothing. Views (the trace, cost) are projections that store ids.

## 14. Storage and cost

Gone (2026-10-16): the storage port, its memory and IndexedDB backends, segment heads and the cost cases; the state is one
object in memory (simplification/harness.md §4).

## 15. Remote deployments

Gone (2026-10-16), with the host protocol (§12).

## 16. Errors

One operation error, `HarnessError { code, detail }` (`invalid`: a turn while one is active, a resume with none
parked, an unknown head; `pointer-busy`). An attempt's failure is not an operation error: it has a code, a `transient`
hint, an `effect` (`none`, `applied`, `unknown`) and, when the provider said, `retryAfterMs`, and opens `failure`.

## 17. Out of scope

Durability (persistence, reloads, several runners; simplification/plan.md "Later"), more than one active turn, agent
profiles, UI stores, compaction strategies (helpers only), garbage collection of entries.

## 18. How Prism hosts it

Prism's side is designed in [../prism-ui/design/contexts.md](../prism-ui/design/contexts.md) §4–5 and
[../prism-ui/design/agent.md](../prism-ui/design/agent.md), and is not built yet:

- **Tools are commands**: one registry; a command's params become the tool's JSON Schema, its metadata says whether it
  reads, needs approval, and which capabilities it uses. A handler runs `execute(id, args, { source, signal })` and maps
  the result's effect to the tool result.
- **The drive** is the app's: approval through the Assistant (park, then resume with the person's answer), reads
  scheduled together and writes one at a time, compaction at the boundary.
- **One undo step per turn** through the mode's checkpointer; every change carries its origin (`agent`).
- **Storage**: none; the harness lives in the page's memory, and a reload ends the conversation. Model keys stay in the
  host.
