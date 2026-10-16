# Testing (`@molstar/ai-conformance`)

Tests live in `_test/` folders next to the code and run in the repository's one Vitest run:

```sh
pnpm test:unit packages/ai             # all AI packages
pnpm test:unit packages/ai/session     # one package
```

No test calls a real model or a network: models are scripted, provider HTTP is replayed from fixtures. Anything that
calls a real model costs money; ask before running it.

## The conformance suite

`runHarnessConformance({ name, storage: () => HarnessStorage })` runs the 59 essential cases against a backend. It is
the harness's contract: a change to a contract (a record, an operation, a decision, an ending) comes with a case here.
The one runner is `components/src/harness/_test/conformance.test.ts` (memory); the IndexedDB run is gone, and with it
the cases only a durable or shared store makes reachable: backend faults, leases and launch races, supervisor liveness,
schema versions, tombstones, crash re-runs and storage cost (what was cut, and why:
`.v6/ai/simplification/conformance.md`).

| File (`conformance/src/cases/`) | Ids | Cases | What                                                                             |
| ------------------------------- | --- | ----: | -------------------------------------------------------------------------------- |
| `journal.ts`                    | J   |     2 | stop during tools leaves `main` as it was; the message hash is stable            |
| `turns.ts`                      | T   |     5 | begin while active, Stop (a live runner, before a settlement), reset, catch-up   |
| `homes.ts`                      | H   |     3 | negative cases: each concern stays in its home                                   |
| `steps.ts`                      | E   |     9 | step endings, rounds, retries, park and resume                                   |
| `tools.ts`                      | X   |     9 | denied, unbound, failed, skipped calls, keys, parked permission, sibling failure |
| `nested.ts`                     | N   |     5 | nested runs: records, depth, pointers, usage                                     |
| `budgets.ts`                    | B   |     4 | charging, token shortfall, `endBy`                                               |
| `projection.ts`                 | C   |     5 | digest, prefix, window admission, unrepresentable, unknown tool                  |
| `observation.ts`                | O   |     1 | observers cannot affect a turn                                                   |
| `performance.ts`                | P   |     1 | P11: no round hashes more than the first                                         |
| `coverage.ts`                   | V   |    12 | unbound calls, drive errors, reset, nested decisions, limits, defaults, filtered |
| `regressions.ts`                | R   |     3 | bugs found, each pinned                                                          |

Case names start with their id and cite the design section (`'R4 stop-refuses-exit: … (harness.md §5, §6)'`).

`@molstar/ai-conformance/fibers`: `runFiberConformance({ name, create, durable? })`, the fibers' suite, run on memory
and IndexedDB rows in `components/src/fibers/_test/`.

## Test doubles (`@molstar/ai-conformance/testing`)

| Export                                                                    | What                                                                                                                                                                        |
| ------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `scriptedModel(steps, options?)`                                          | a `ModelAdapter` playing `ModelStep`s (`text`, `calls`, `reasoning`, `finish`, `usage`, `fail`, `waitMs`, `deltas`); `.bindings`, `.projections`, `.invocations`, `.push()` |
| `scriptedTools(specs)`                                                    | `ToolBindings` from functions or specs; records calls                                                                                                                       |
| `SimulatedClock`, `settle()`                                              | virtual time (`advance`, `runUntil`); `settle` lets WebCrypto hashing finish                                                                                                |
| `simulatedIds()`, `simulatedSupervisor()`                                 | deterministic ids; a supervisor that starts once per launch key, optionally losing launches (`failLaunch`)                                                                  |
| `user`, `system`, `assistant`, `toolDefinition`, `declareTools`, `textOf` | message builders                                                                                                                                                            |

`createSession(fixture, setup)` (from the main entry) wires host, client, port, scripted model and tools on a simulated
clock for harness-level tests: `send(prompt)`, `run(launch)`, `events()`, `main()`, `turns()`.

## Testing at each level

- **Harness change**: a conformance case with `createSession`, scripted model and tools, and assertions on events,
  `main()` and turns.
- **Tools**: `tools/src/_test/`: `tool.test.ts` checks tools directly (`directRun`, handlers), `turn.test.ts` runs them
  in harness turns (`createSession`).
- **Session and drive**: `session/src/_test/helpers.ts` has a toy host (a counter with a version), its toolkit and an
  `AgentSession` over `scriptedModel`; tests drive `turn`, `resume`, `stop` and check `TurnResult` and hooks.
- **Adapters**: SSE fixtures and a recording `fetch` (`model-adapters/src/_test/support.ts`); `compatibility.test.ts`
  pins request bodies by hash, so a change to wire bytes must be deliberate: update `retained-wire.json` and raise the
  adapter's `revision`.
- **Trace**: `trace/src/_test/helpers.ts` runs turns on `createSession` (with segmenting model bindings) and reads them
  back with `readTrace`.

## Pitfalls

- Hashing is WebCrypto (off the microtask queue). On a simulated clock, wait with `settle()` or `runUntil`, not a fixed
  number of ticks; the clock counts in-flight digests as pending work.
- Some cases still wait a fixed number of event-loop rounds and retry under load (`.v6/ai/plan.md`, "The suite waits on
  conditions"); prefer waiting on the condition in new cases.
- Tests import packages by name (`@molstar/ai-conformance/testing`); Vitest resolves them to sources through the
  `molstar-src` condition.
