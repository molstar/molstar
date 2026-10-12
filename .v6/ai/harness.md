# The harness

The harness runs agent turns durably: a model is called, the tools it asks for run, and the loop goes on until the turn
ends, surviving reloads, crashes and runners that move. It is generic: it knows nothing of Mol\* or Prism, which host it
(§18). Code comments cite this document's sections (`harness.md §6`); the numbering is fixed.

## 1. Parts and homes

| Part                                | Owns                                                       | Package                      |
| ----------------------------------- | ---------------------------------------------------------- | ---------------------------- |
| Host (`createHarnessHost(storage)`) | durable state and the protocol operations, each atomic     | `@molstar/ai-harness`        |
| Runner (`createRunner`)             | runs steps, binds and runs tools, asks the drive, commits  | `@molstar/ai-harness`        |
| Drive                               | every choice                                               | the app; defaults in harness |
| Supervisor                          | liveness: while there is work, exactly one runner is on it | one per runtime              |
| Model adapter                       | projection and streaming for one provider                  | `@molstar/ai-model-adapters` |
| Tool handler                        | one tool's effects                                         | the app                      |

`@molstar/ai-agent` holds the model side's shapes: `Message` and its parts, tool definitions, `ModelConfiguration`, the
`ModelAdapter` port, `hashMessage`, `digest`.

Each concern has one home and is kept nowhere else: what a model sees is the **journal** (§2); evidence (requests,
attempts, results, usage) is the **event log** (§3); what runs next is read from the last step event, never from a head;
a call's authority is the definition it was planned against (§8); liveness is the runner state and the supervisor (§4);
every choice is the drive's (§7).

## 2. The journal

- An **entry** is a `Message`, addressed by `hashMessage` (SHA-256 over a structural byte encoding, never JSON text;
  WebCrypto, so async). Identical messages are one entry; entries are immutable and never deleted.
- A **head** is an immutable list of entry hashes. A **pointer** names a current head (`main` is the turn's; apps and
  tools name others). A pointer moves only in a commit, by an **append** or an **edit** (kept hashes plus new messages:
  compaction), compare-and-set on its old head.
- Every head is **projectable**: each tool call has its result in the next user message. An assistant message with calls
  enters a head only together with its results.
- A **contribution** identical to an entry already in the head is not appended again (a skill's guidance is kept once).

## 3. The event log

The only evidence store, numbered by `seq`. A **step event** records how a step ended: its calls and their attempts
(times, outcomes, usage; nested calls by `parentCallId`), its moves, the items it consumed, its compactions, and its
`next` step or the turn's `end`. **Control events** record the other operations (`submitted`, `stopped`, `reset`,
`claimed`…). Turns have no records: the event that ends one carries its summary. Large bodies are blobs; a model call's
evidence is what reproduces it (head, configuration, tools, adapter revision, payload digest), not its payload.

## 4. The Execution, runners and the supervisor

- One mutable **Execution** per harness: `seq`, `main`, the inbox, the runner state and the active turn.
- **Runner state** is `none`, `launching(key)` or `live(lease)`. An operation that creates work returns a launch key
  when no runner is on it; `claim` with that key makes the runner live. Every operation renews the lease, so a busy
  runner sends no heartbeat. A failed read is not a lost lease; only a `stale` refusal or expiry is. A turn written by
  another schema version is blocked (stop or reset it; there is no migration).
- **The supervisor** starts one runner per launch key and checks when a launch window or lease runs out, so a turn whose
  owner died is taken up without anyone asking. The in-process one is in `ai-components`.

## 5. Turns

- At most one active turn. A runner with no turn asks the drive `next-turn`; `begin` creates it in one commit. A turn
  ends `completed`, `failed`, `stopped` or `reset`.
- **Stop** is durable and any client may call it, naming a turn or an item (cancelling a prompt cannot race a turn that
  is beginning). With a live runner it is a request: the runner aborts and records the step `abandoned`. Stop wins over
  a settlement in flight.
- **Reset** ends any turn and moves pointers back (by default `main` to its head at the last completed turn); it deletes
  nothing and never runs on its own.
- **Park**: a `permission` answered `park` parks the step and releases the runner; `resume` re-runs the step and asks
  the drive again.

## 6. Steps

- A **round** is a model step, plus a tools step for its calls. The model step projects, admits (window and budget) and
  streams; with calls it settles as evidence only and names the tools step, without calls it appends to `main`. The
  tools step repairs unbound calls, asks `permission`, schedules, runs, then appends the assistant message, one user
  message with every result, and the contributions, in one move.
- A turn runs in memory; the store is touched at commits only, and everything stored but the Execution and pointers is
  immutable, so a runner caches it freely.
- **Endings**: `settled`, `parked`, `interrupted`, `abandoned` (Stop, `fail`, time); a crash loses the step and it
  re-runs. Within a step, stop wins over fail, park and interrupt.
- **Retries** wait in process: a model retry sends the same request, a tool retry the same key with `prior: "possible"`.

## 7. Decisions and drives

A decision is `{ kind, path, payload }` (`path`: the call ids down to the nested run that asks). The kinds:
`permission`, `schedule`, `budget`, `failure`, `model-result`, `repair`, `boundary`, `interrupt`, `next-turn`. It goes
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
- **Keys**: a call's `key` names its intended effect (the call id at top level, saved before any effect; a digest in a
  nested run), and `prior: "possible"` says it may have run before, so a tool that is not idempotent can look for its
  earlier effect.
- **A tool's own context**: `journal.read` any pointer; `append` and `edit` on one it holds (never `main`), staged until
  the step commits.

## 9. Nested runs

A handler or a drive calls `run(spec)` to run an agent inside its call: its own input, tools, model, budget and optional
drive, starting from a pointer, a head or empty. The same engine, its own budget (never charged to the parent), its
calls recorded under the parent call, nothing durable until the parent step commits: a re-run parent re-runs it whole.
Not steerable; depth 4.

## 10. Inbox and steering

`submit` adds an item (key, messages, turn options, optional interrupt), idempotent by key, at any time from any client.
Items are consumed once: by a begin, a steer at a round boundary, or an interrupt. An interrupt aborts the model step,
or skips or abandons the tools step's calls.

## 11. Projection, caching and budgets

- Projection is pure (head, configuration, available tools), and appending keeps the projected prefix byte-stable, so
  provider caches hold. Adapters report the request's segments so a trace can say where a cache broke.
- A request is admitted only if its estimate plus the reserved output fits the window (else `capacity`).
- **Budgets** are limits plus consumption per turn and per nested run, over calls, tokens by kind, `cost` (with an app
  price) and `endBy`. Every ended attempt is charged once; estimates where usage is missing.

## 12. Protocol operations

Client operations: `submit`, `withdraw`, `stop`, `resume`, `reset`, `check`. Runner operations: `claim`, `putBlobs`,
`begin`, `commit`, `renew`, `release`. Every operation is atomic; keyed ones are idempotent by key. A commit is refused
unless the lease is the runner's, the step is the pending one and each move starts from its pointer's current head. The
host trusts the runner beyond these checks; a host behind a trust boundary validates at that boundary.

## 13. Observation

`runner.observe()` emits live observations (`model.delta`, `tool.progress`, attempts, decisions) and, once committed,
events. Live ones are best effort and carry ids, not payloads. `watch()` is level-triggered (it says only that `seq`
moved), so a missed emission loses nothing. Views (the trace, cost) are projections that store ids.

## 14. Storage and cost

- A backend implements `HarnessStorage`: one serializable transaction per harness over the Execution, events, messages,
  pointers, heads and segments; blobs outside it; a watch. Protocol logic is the host's; a backend only stores and
  indexes. Foreign JSON is stored as an opaque string.
- Heads are stored as segments in size classes, so appending is near-constant and a head reads in about log _n_ reads.
- Backends: memory (tests, helpers) and IndexedDB (the page). A SQLite backend for servers and headless hosts is
  planned, not built. One durable core under the backends runs a transaction against its staged writes and retries on a
  conflict.
- **Cost is measured, not timed** (the conformance suite's performance cases): a commit costs a constant plus what it
  adds; cold work is linear in the selected head only.

## 15. Remote deployments

Apps and runners use separate handles, `HarnessClient` and `RunnerPort`, each under its own credential when remote: an
app can never commit, a runner never submit, stop or reset. Operations are serializable request and response; an
`unknown-outcome` is retried with the same key.

## 16. Errors

One operation error, `HarnessError { code, retry, detail }` (`invalid`, `head-changed`, `inbox-changed`, `stale`,
`stale-schema`, `pointer-busy`, `unknown-outcome`, `integrity`). An attempt's failure is not an operation error: it has
a code, a `transient` hint and an `effect` (`none`, `applied`, `unknown`), and opens `failure`.

## 17. Out of scope

Durable nested runs, more than one active turn, agent profiles, UI stores, compaction strategies (helpers only), garbage
collection of entries, migrating turns across schema versions.

## 18. How Prism hosts it

Prism's side is designed in [../prism-ui/design/contexts.md](../prism-ui/design/contexts.md) §4–5 and
[../prism-ui/design/agent.md](../prism-ui/design/agent.md), and is not built yet:

- **Tools are commands**: one registry; a command's params become the tool's JSON Schema, its metadata says whether it
  reads, needs approval, and which capabilities it uses. A handler runs `execute(id, args, { source, signal })` and maps
  the result's effect to the tool result.
- **The drive** is the app's: approval through the Assistant (park, then resume with the person's answer), reads
  scheduled together and writes one at a time, compaction at the boundary.
- **One undo step per turn** through the mode's checkpointer; every change carries its origin (`agent`).
- **Storage**: IndexedDB in the page; a SQLite backend for a headless or server host is planned, not built. Model keys
  stay in the host.
