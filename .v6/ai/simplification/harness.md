# The harness: a basic in-memory loop for driving an app

**Decided** (David, 2026-10-16): the durability machinery goes too, "just basic in memory harness for now for driving
an app". This replaces the analysis's first recommendation, to keep it for this round (persistence.md §4). What follows is what that means in the
code: what the harness keeps because the session, the tools and the trace depend on it, what goes, the shape of the
new API, and the callers.

## 1. What the harness is for, from its one driver

`AgentSession` (`session/src/session.ts`) is the only code that drives a harness, and prism-viewer drives only the
session (`packages/prism/viewer/src/agent/session.ts:116`). What it uses (`git grep "client\.\w\+"` in
`session/src`, `trace/src`, `packages/prism`):

| Uses                                                                                       | Does not use                                                         |
| ------------------------------------------------------------------------------------------ | -------------------------------------------------------------------- |
| `submit` (one item per turn), `stop`, `resume` (a parked step), `reset` (`forget()`), `check` (takeover after a reload) | `withdraw`, interrupt items, steering (`boundary` answered `steer`): "each prompt is its own turn" (`session/src/drive.ts:356`) |
| reads: `execution` (`turn`, `main`, `inbox`, `runner`), `turns`, `steps`, `heads`, `entries`, `watch` | `pointers`, `remove`, more than one runner, remote handles           |
| the runner's `observe()` for live hooks                                                    | `runner.abort()` as a local stop                                     |

The trace (`trace/src/read.ts:437–478`) and the Inspector (`prism/assistant/src/viewer/inspector-sources.ts:97`) read
through `HarnessReader`: `execution`, `events`, `turns`, `entries`, `heads`, `blobs`, `watch`.

## 2. What stays

These are the harness's value for an app and stay as they are (`.v6/ai/harness.md` sections):

- **The journal** (§2): entries are `Message`s addressed by hash, deduplicated (a skill's guidance once, E2
  skill-twice), labels hashed but never sent or counted; `main` plus pointers that tools and nested runs hold; heads stay
  projectable (every call answered); a step's writes apply at its end, so a stop leaves `main` as it was before the round.
- **The event log** (§3): step events with calls, attempts, usage and failures; the trace is built from it.
- **Turns and steps** (§5, §6): one active turn; model step, tools step, rounds; endings `settled`, `parked`,
  `abandoned`; stop > fail > park; park and resume for approvals; reset.
- **Decisions and the drive** (§7): `permission`, `schedule`, `budget`, `failure`, `model-result`, `repair`,
  `boundary` (without `steer`), and a turn-start decision (today `next-turn`) where the session compacts before a turn.
  `defaultDrive`, `drive(handlers, fallback)`, `legal()`.
- **Tools** (§8): planning against the loaded definition, repair of unbound calls, the frame protocol, contributions,
  a tool's own journal; `key` and `prior` on a retry.
- **Nested runs** (§9): compaction and nested tools use them (`session/src/compaction/run.ts`, `tools/src/nest.ts`).
- **Projection, admission, budgets** (§11), **observation** (§13).

## 3. What goes

| Mechanism                                                                                       | Files                                                                                                                       | Lines (today)     |
| ----------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------- | ----------------- |
| The host protocol: `HarnessHost`, `HarnessClient` / `RunnerPort` split, `begin`/`commit` batches, keyed idempotent operations (`key`, `inputDigest`, replay), `unknown-outcome` | `harness/src/host/{index,client,runner,writer}.ts`, `contracts/ports.ts` (most)                         | 919 + ~150        |
| The storage port and its backend: `HarnessStorage`, `StorageTx`, `RecordStore`, the event-log indexes, tombstones | `contracts/storage.ts`, `components/src/harness/storage/memory.ts`                                                 | 74 + 272          |
| Leases and launches: `claim`, `renew`, `release`, `check`, `LaunchKey`, `RunnerState`, `signalEveryMs`, heartbeat, `handOffAfterMs` | `runner/index.ts` (`read`, `send`, `commitOp`, `examine`, `refresh`, `watch`, `heartbeat`, `release`, `serve`), `contracts/records.ts` | ~400 of 762 |
| The supervisor                                                                                  | `components/src/harness/supervisor/in-process.ts`, `contracts/supervisor.ts`                                               | 233 + 33          |
| Schema versions: `SCHEMA_VERSION`, `blocked: 'stale-schema'`                                     | `contracts/records.ts`, `host/*`, `runner/index.ts`                                                                         | ~60               |
| Crash re-runs of a step (`prior: 'possible'` after a crash; nested keys stable across re-runs)  | `runner/index.ts` `step(turn, rerun)`, `harness/src/nested/index.ts` key digests                                           | ~60               |
| The runner's cache and spill (`WorkingSet`, `cacheBytes`, retention inspection)                 | `runner/state.ts:91–262`, `contracts/inspect.ts`                                                                           | ~190              |
| Heads as segments in size classes (`HeadRecord`, `Segment`, `buildMove`, `Move.from` compare-and-set) | `head/index.ts`, records, callers in `runner/state.ts`, `host/runner.ts`                                              | 88 + ~100         |
| The inbox's steering and interrupts: interrupt items, `withdraw`, the `interrupt` decision, `steer` answers, the `interrupted` ending | `contracts/drive.ts`, `runner/index.ts` (`steered`, `interrupting`), `tools/step.ts` (skip on interrupt) | ~150 (inferred) |
| `callTool`                                                                                      | `tools/call.ts`                                                                                                             | 56                |
| Error codes `stale`, `stale-schema`, `unknown-outcome`, `integrity`, `inbox-changed`, `head-changed` | `contracts/ports.ts`                                                                                                   | —                 |

Of the harness's 4,450 source lines plus the components' 505 (memory storage and supervisor), an estimated **~2,200
remain** (inferred, from the file sizes above and §4's sketch): about 2,750 fewer.

## 4. The shape, proposed

```ts
// @molstar/ai/harness
const harness = createHarness({
  models,                // ModelBindings: resolve(target) → provider
  tools,                 // ToolBindings
  drive,                 // Drive
  clock?, ids?, price?, limits?: { depth? },
});

const turn = harness.turn({ messages, configuration }); // refused while a turn is active
turn.done;               // Promise<TurnEnd | { parked }>: ends completed, failed, stopped, or parks
harness.stop(reason?);   // the active turn: aborts its step, ends it stopped
harness.resume();        // a parked turn: re-runs the parked step, asking the drive again
harness.reset(heads?);   // ends any turn; main back to the last completed turn's head (or the given heads)
harness.observe();       // live observations and each event, as today's runner.observe()

// the reader, as today (the trace and the Inspector need no change):
harness.execution(); harness.events({ after, limit }); harness.steps(ids); harness.turns({ before, limit });
harness.entries(ids); harness.heads(ids); harness.blobs(refs); harness.watch();
```

- **State** is one object in memory: entries (`Map<EntryId, Message>`), heads (`Map<HeadId, EntryId[]>`), pointers,
  the event log (an array; `seq` its index), blobs (`Map`), the active turn. No transactions: a step builds its new
  heads and events and applies them at its end in one synchronous assignment, so a stop or a failure mid-step leaves
  nothing (§2's "a step's writes apply at its end").
- **The loop** is the runner's `step` / `round` / `exitStep` logic (`runner/index.ts:345–575`, `runner/model.ts`,
  `runner/run.ts`, `tools/*`), called directly by `turn()`; no claim, no lease, no serve loop.
- **`Execution`** stays as the reader's snapshot (`main`, `turn` with `pending`, `park`, `budget`, `usage`), without
  `runner`, `inbox`, `schemaVersion`.
- **Events**: `StepEvent` as today minus `commitKey` / `inputDigest`; control events `turn-started`, `stopped`,
  `resumed`, `reset` (`submitted`, `withdrawn`, `claimed`, `launched`, `released`, `blocked` go).
- **Blobs** keep their refs (decided, David, 2026-10-16): request and result bodies stay in a blob map addressed by
  `digest(body)`, so `readTrace` and the Inspector stay unchanged. Only `putBlobs` before a commit goes: the step puts
  its blobs with the rest of its writes.

## 5. The session over it

`AgentSession` (945 lines) loses the setup of host, client, port and supervisor (`session.ts:394–449`), the watch-based
wait for a turn's result (`:790–835`, which polls `execution()` and `turns()` until the submitted key began and
ended), the takeover path (`interrupted()`, `resume()` without decisions after a lease, `:713–771`) and the options
`storage`, `harnessId`, `leaseMs` (`:110–115`). `turn()` awaits `harness.turn(...).done`; `resume(hooks, decisions)`
calls `harness.resume()`. Estimated ~300 lines fewer (inferred).

Callers outside:

- `packages/prism/viewer/src/agent/session.ts:9,29–30,121–122` (the `storage` and `harnessId` options and the
  `HarnessStorage` type): removed.
- `packages/prism/assistant/src/viewer/inspector-sources.ts:97` (`client.watch()`): the session exposes the harness's
  reader instead of a `HarnessClient` (`session.client()`, `_spec/session.md` § Turns); same methods.
- `packages/prism/assistant/src/viewer/agent-model.ts`: types only (`TurnResult`, `ToolEvent`, …), unchanged.
- Docs: no promise to take back. The user docs already say a reload ends the record and every turn
  (`apps/docs/src/content/docs/v6/prism-viewer/assistant.md:133`); `.v6/prism-ui/design/agent.md` and `spec/agent.md`
  mention reloads only for key scopes. The AI spec does promise it (`_spec/README.md` step 7, `_spec/session.md`
  `interrupted()`), and slice 10 rewrites that.

## 6. Tests

- `components/src/harness/_test/supervisor.test.ts` (173) and `testing/supervisor.ts` (48): go.
- `components/src/harness/_test/memory.test.ts` (375, the storage port's behavior): goes; a short state test (~100)
  covers dedup, heads and the event log reads the trace uses.
- `harness/src/_test/head.test.ts` (83, the segment scheme, fast-check): goes; `fast-check` then leaves the AI package
  if S1 is gone too (conformance.md).
- Conformance: 11 of the 70 cases the first pass kept test mechanisms that go (keys, commit protocol, compare-and-set moves,
  inbox and interrupts); 59 stay (conformance.md §2).
- Session tests that reopen a harness on shared storage (`session/src/_test/drive.test.ts:176–200`,
  `resume.test.ts:36` "a new runtime reads the harness and finishes the interrupted turn") are dropped.

## 7. Hashing: synchronous (decided)

**Decided** (David, 2026-10-16): hashing becomes synchronous. Today `hashMessage` and `digest` are WebCrypto SHA-256
(`agent/src/hash.ts:96–121`), so async (`_spec/harness.md` Gotchas): every entry id, payload digest and blob ref is an
`await`, the simulated clock patches `crypto.subtle.digest` to count digests in flight (`testing/clock.ts:21–54, 170–174`)
and gives them real time before calling a deadlock, and tests must `settle()`. Those WebCrypto races are also why the
suite retried cases and runs in the serial Vitest project (`.v6/ai/plan.md` "The suite waits on conditions",
`vitest.config.ts:14–16`).

- **The function**: SHA-256 in plain TypeScript in `messages/hash.ts` (~70 lines, inferred: a standard
  implementation), with no dependency. The AI code may not import Mol\*'s own hashes (`packages/core/src/data/util/hash-functions.ts`
  has FNV variants), and an entry id is content addressing that deduplicates, so a cryptographic-strength hash keeps
  collisions out of the question. Same algorithm, same output: ids and digests do not change, the RFC 8785 vectors
  (conformance C1 digest-jcs) still hold, and a test checks it against WebCrypto.
- **Also**: it works outside a secure context (`crypto.subtle` needs one: a viewer served over plain http from a
  non-localhost host has none).
- **What changes**: `hashMessage`, `digest` and `encodeMessage` callers drop `await` (about 50 call sites in
  `packages/ai`, `git grep "hashMessage(\|digest("`); the step code that hashes before applying gets simpler; the
  attachments' byte digest (`model-adapters/src/common/attachments.ts:51`) uses the same function;
  `testing/clock.ts` loses the digest tracking (~40 lines) and `settle()` its digest wait; the P11 hash-once case counts
  calls through an injected counter instead of spying on `crypto.subtle.digest`
  (`conformance/src/cases/performance.ts:480–482`).
- **Cost**: a JS SHA-256 is slower than WebCrypto on large inputs (inferred); entries are messages of kilobytes and
  payload digests are per model call, so it should not show. To check in the slice: hash a long conversation's
  payload and time it.
