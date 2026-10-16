# Conformance: the essentials

## 1. What it costs today

`@molstar/ai-conformance` is 8,914 lines, all test infrastructure:

| Part                                         | Lines | What                                                                       |
| -------------------------------------------- | ----: | -------------------------------------------------------------------------- |
| `src/cases/*.ts` (13 files + `helpers.ts`)   | 7,213 | 219 harness cases (217 + 2 `it.todo`), ids J T H E X N B C S O P V R       |
| `src/fibers/`                                |   705 | the fibers' suite (12 cases), fixtures, generators                         |
| `src/testing/`                               |   602 | the doubles: `scriptedModel`, `SimulatedClock`, ids, supervisor, tools, messages |
| `src/session.ts`                             |   186 | `createSession`, the harness fixture the cases and other packages' tests use |
| `src/meter.ts`                               |   158 | `meterStorage`: counts reads and writes at the storage port (the P cases)  |
| `src/suite.ts`, `index.ts`                   |    50 | `runHarnessConformance`: every case against one backend                    |

It runs twice, from `components/src/harness/_test/`: on memory (`conformance.test.ts`) and on fake IndexedDB
(`conformance-indexeddb.test.ts`). Every AI test file sits in the `serial` Vitest project, which runs after all other
projects with `fileParallelism: false` because the cases "count and wait for rounds of real work (timers, WebCrypto,
fake IndexedDB)" (`vitest.config.ts:14–17, 79–81`).

| Run                                                   | Local (Apple silicon) | CI (run 38695888995, "Native tests") |
| ----------------------------------------------------- | --------------------: | -----------------------------------: |
| `conformance-indexeddb.test.ts` (220)                 |               13.9 s  |                               25.9 s |
| `conformance.test.ts` (219)                           |                4.1 s  |                                8.5 s |
| `components/src/harness/_test/supervisor.test.ts` (6) |                1.5 s  |                                1.4 s |
| all other AI test files (48)                          |                3.0 s  |       ~16 s of the tail (with imports) |
| **the AI tests**                                      |  **22.5 s** (30 s wall) |  **~52 s** of a 200 s Vitest run (07:03:38 → 07:04:30) |

CI's whole "Native tests" job took 4 min; the AI serial tail is about a quarter of its Vitest time.

Inside the memory run, by group (count / time): P 33 / 2,098 ms, S 7 / 895 ms, J 11 / 432 ms (J4's two long-history
cases 331 ms), V 37 / 183 ms, E 25 / 131 ms, X 20 / 77 ms, B 8 / 68 ms, T 34 / 60 ms, R 12 / 51 ms, N 12 / 36 ms,
C 10 / 26 ms, H 8 / 21 ms, O 2 / 6 ms. On IndexedDB, P alone is 7.3 s and J 4.4 s.

## 2. What is essential, once there is one backend

The suite's design is "the contract every storage backend runs" (`_spec/testing.md`). With memory the only backend that
reason goes: the cases become the harness's own tests, and the question is which behaviors a host relies on. The rule
used here:

- **keep** one case per clause of the contract a host in one page relies on: each decision kind's legal answers and
  default, each step ending and their precedence, the journal's invariants (compare-and-set moves, `main` read-only,
  projectable heads, dedup of contributions, labels never counted), idempotency by key, stop / park / resume / reset,
  nested runs (depth, pointers, usage), budgets (charge once, extend, `endBy`), projection (digest stable, prefix, window
  admission, unrepresentable), and the regressions that are not about storage or racing runners;
- **drop** what only a durable or shared store, or a second runner, makes reachable: backend faults
  (`unknown-outcome`, late faults, integrity of rewrites), lease and launch races, takeover after expiry, schema
  versions, tombstones, storage cost (the P cases and `meterStorage`), liveness under a supervisor (S), and crash
  re-runs. With memory only a lease never expires with work pending unless a runner dies inside the page (persistence.md
  §4); those paths are unreachable in production today. If durability returns, these cases come back from git with the
  backend (plan.md, "Later").

### Kept (70 of 219)

| Group | Kept | Cases                                                                                                                                                       |
| ----- | ---: | ----------------------------------------------------------------------------------------------------------------------------------------------------------- |
| J     | 4/11 | J2 stale-move, J3 stop-during-tools, J3 interrupt-model-step, J4 hash-stable                                                                               |
| T     | 10/34 | T1 begin-while-active, T1 wrong-step, T2 same-key-same-result, T2 same-key-other-input, T3 begin-order, T4 stop-live-runner, T4 stop-beats-settlement, T4 stop-item-waiting, T5 reset-default, T6 catch-up |
| H     | 3/8  | H1 contribution-like-answer, H1 contribution-with-calls, H1 redeclared                                                                                     |
| E     | 10/25 | E1 invalid-answer, E2 round-order, E2 no-calls-round, E2 skill-twice, E2 label-not-counted, E4 retry-same-digest, E4 retry-after, E4 replace-new-call, E5 completing-steer, E6 park-resume |
| X     | 9/20 | X1 denied, X1 unbound-report, X1 failure-report, X1 invalid-arguments, X1 skip-rest, X1 no-result, X2 park-no-effects, X3 key-across-retry, X5 sibling-failure |
| N     | 5/12 | N1 nested-records, N2 depth-exceeded, N2 pointer-busy (main), N2 tool-journal, N3 nested-usage (the running budget excludes it)                            |
| B     | 4/8  | B1 missing-usage, B1 cache-categories, B2 token-shortfall, B2 end-by-mid-step                                                                              |
| C     | 5/10 | C1 same-input-same-digest, C1 append-prefix (after a model-and-tool round), C2 over-window, C2 unrepresentable, C2 unknown-tool-call                       |
| S     | 0/7  | —                                                                                                                                                           |
| O     | 1/2  | O1 throwing-subscriber                                                                                                                                      |
| P     | 1/33 | P11 hash-once ("no round hashes more than the first"): hashing is real CPU in memory too                                                                   |
| V     | 12/37 | V1 not-loaded, V2 drive-throws, V3 reset-clears-park, V4 nested-fail, V4 nested-stop, V4 nested-pass, V5 resume-not-parked, V6 tool-calls-limit, V7 schedule-invalid, V7 default-model-retry, V7 default-tool-retry, V8 filtered |
| R     | 6/12 | R2 interrupt-already-stored, R3 drive-run-in-model-result, R4 stop-refuses-exit, R5 extend-in-interrupted-step, R6 interrupted-stored-marks, R12 no-adapter |

The case names are the suite's (`conformance/src/cases/*.ts`); the list is a proposal to argue with, case by case.
Borderline drops worth a second look: T4 stop-item-begun and stop-parked-with-items, E6 nested-park, X1 skip-calls,
N2 fork, V7 max-concurrency, R8–R11 (interrupt and withdraw races inside one runner, which memory can still hit).

### Dropped (149), by reason

| Reason                                     | Cases                                                                                                                         |
| ------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------- |
| storage cost and the meter                 | P (32 of 33), J2 many-pointers, J4 long-history, J4 edit-then-append                                                          |
| liveness, leases, launches, supervisor     | S (all 7), T1 expired-lease, T1 zombie, T4 stop-then-expiry, stop-while-launching, stop-right-after-claim, stop-before-release, only-abandon-after-request, T6 renew-silent, two-events, E4 retry-backoff, R7 stale-with-expired-lease |
| backend faults and integrity               | J1 rewrite-refused, T2 unknown-outcome, T2 late-fault, R1 cold-read-fault                                                     |
| schema version, tombstones                 | T1 schema-differs, T2 remove-tombstone, V3 reset-ends-blocked                                                                 |
| crash re-runs                              | J3 crash-after-model-step, T3 steer-across-crash, E1 crash-before-decision, E3 crash-nothing, X3 key-across-crash, N1 nested-crash, B1 usage-once, B2 extend-persisted (lost to a crash) |
| a second case of a clause already kept     | the rest (e.g. the second C1 append-prefix, the second N3, the second V7 fail-beats-park, E1's commit-shape checks, which `.v6/ai/README.md` decision 4 says the host does not police) |
| `it.todo`                                  | H1 import-rule (Biome checks it), the S3 policy-hooks todo                                                                    |

## 3. Where the cases go, and sizes

- The 70 cases become `packages/ai/src/harness/_test/*.test.ts` (one file per group or per theme), calling
  `createSession` from `testing/` directly. No `ConformanceFixture`, no `runHarnessConformance`, no `suite.ts`, and
  `vitest`/`fast-check` leave the shipped code (`conformance/package.json` lists them as dependencies today).
- `meter.ts` (158) goes with P; `createSession`'s metering (`session.ts:30`) goes with it.
- `testing/` stays the public `./testing` entry: `scriptedModel` (163), `SimulatedClock` (197), ids, messages, tools;
  `testing/supervisor.ts` (48, `simulatedSupervisor`) goes with S unless option A's supervisor test keeps it.
- `components/src/harness/_test/supervisor.test.ts` (173, 1.5 s): keep two of six (one runner per launch key; a lost
  runner relaunched) while the supervisor stays (persistence.md §4).
- The fibers' suite becomes the threads' tests (threads.md §4).

| Part                         | Today | After (estimate) |
| ---------------------------- | ----: | ---------------: |
| harness cases + helpers      | 7,213 | ~2,300 (70 cases at the ~30 lines a case the non-P files average) |
| fibers suite                 |   705 | ~350 (no delta cases, no generators) |
| meter, suite, index          |   208 | 0                |
| fixture (`session.ts`)       |   186 | ~150             |
| testing doubles              |   602 | ~550             |
| **total**                    | 8,914 | **~3,350**       |

## 4. CI after

- The IndexedDB run goes: −25.9 s in CI.
- P, S and J4 were 3.3 of the memory run's 4.1 s; 70 cases of the rest take an estimated 0.4 s locally, ~1 s in CI.
- Estimated AI tail in CI: ~52 s → **~15 s** (the 46 other files' share stays, minus the files deleted with narration,
  `compatibility.test.ts`, the IndexedDB rows). Inferred, from the per-file timings above.
- The `serial` project may then go: its reason was the P, S and IndexedDB cases losing races under load
  (`vitest.config.ts:14–16`; `.v6/ai/plan.md` "The suite waits on conditions"). If the remaining files pass in the
  parallel `unit` project under load (to verify: several full runs), the AI tests stop adding a serial tail at all.
  Inferred, not measured.
