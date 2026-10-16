# Testing (`@molstar/ai`)

Tests live in `_test/` folders next to the code and run in the repository's one Vitest run:

```sh
pnpm test:unit packages/ai                 # the whole package
pnpm test:unit packages/ai/src/session     # one folder
```

No test calls a real model or a network: models are scripted, provider HTTP is replayed from fixtures. Anything that
calls a real model costs money; ask before running it.

## The harness cases

The 59 essential cases run on the in-memory harness (`createHarness`), as plain test files in `src/harness/_test/` (one
per row below, `helpers.ts` beside them). They are the harness's contract: a change to a contract (a record, an
operation, a decision, an ending) comes with a case here. Which cases were kept, and why:
`.v6/ai/simplification/conformance.md`.

| File (`src/harness/_test/`) | Ids | Cases | What                                                                             |
| --------------------------- | --- | ----: | -------------------------------------------------------------------------------- |
| `journal.test.ts`           | J   |     2 | stop during tools leaves `main` as it was; the message hash is stable            |
| `turns.test.ts`             | T   |     5 | begin while active, Stop (during a step, before a settlement), reset, catch-up   |
| `homes.test.ts`             | H   |     3 | negative cases: each concern stays in its home                                   |
| `steps.test.ts`             | E   |     9 | step endings, rounds, retries, park and resume                                   |
| `tools.test.ts`             | X   |     9 | denied, unbound, failed, skipped calls, keys, parked permission, sibling failure |
| `nested.test.ts`            | N   |     5 | nested runs: records, depth, pointers, usage                                     |
| `budgets.test.ts`           | B   |     4 | charging, token shortfall, `endBy`                                               |
| `projection.test.ts`        | C   |     5 | digest, prefix, window admission, unrepresentable, unknown tool                  |
| `observation.test.ts`       | O   |     1 | observers cannot affect a turn                                                   |
| `performance.test.ts`       | P   |     1 | P11: no round hashes more than the first                                         |
| `coverage.test.ts`          | V   |    12 | unbound calls, drive errors, reset, nested decisions, limits, defaults, filtered |
| `regressions.test.ts`       | R   |     3 | bugs found, each pinned                                                          |
| `state.test.ts`             | —   |     4 | not cases: the reads (empty state, one entry per message, paging, `watch`)       |

Case names start with their id and cite the design section (`'R4 stop-refuses-exit: … (harness.md §5, §6)'`).

`@molstar/ai/testing/threads`: `runThreadConformance({ name, create, durable? })`, the threads' suite
(`src/threads/_test/`, its own entry because it imports Vitest), run on memory (`src/threads/_test/memory.test.ts`) and
on the viewer's IndexedDB rows (`prism/viewer/src/journal/_test/`).

## Test doubles (`@molstar/ai/testing`, `src/testing/`)

| Export                                                                    | What                                                                                                                                                                                 |
| ------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `scriptedModel(steps, options?)`                                          | a `Provider` playing `ModelStep`s (`text`, `calls`, `reasoning`, `finish`, `usage`, `fail`, `waitMs`, `deltas`); `.provider`, `.bindings`, `.projections`, `.invocations`, `.push()` |
| `scriptedTools(specs)`                                                    | `ToolBindings` from functions or specs; records calls                                                                                                                                |
| `SimulatedClock`, `settle()`                                              | virtual time (`advance`, `runUntil`); `settle` lets pending promises run                                                                                                             |
| `simulatedIds()`                                                          | deterministic ids                                                                                                                                                                    |
| `user`, `system`, `assistant`, `toolDefinition`, `declareTools`, `textOf` | message builders                                                                                                                                                                     |

`createSession(setup)` (from the same entry) makes a harness over a scripted model and tools on a simulated clock for
harness-level tests: `harness`, `observations`, `send(prompt, { configuration?, wait? })` (a turn, run to its outcome
unless `wait: false`), `run(handle)`, `resume()`, `events()`, `main()`, `turns()`; `setup.models` and
`setup.toolBindings` replace the scripted bindings.

## Testing at each level

- **Harness change**: a case in `src/harness/_test/` with `createSession`, scripted model and tools, and assertions on
  events, `main()` and turns.
- **Tools**: `src/tools/_test/`: `tool.test.ts` checks tools directly (`directRun`, handlers), `turn.test.ts` runs them
  in harness turns (`createSession`).
- **Session and drive**: `src/session/_test/helpers.ts` has a toy host (a counter with a version), its toolkit and an
  `AgentSession` over `scriptedModel`; tests drive `turn`, `resume`, `stop` and check `TurnResult` and hooks.
- **Providers**: SSE fixtures and a recording `fetch` (`src/providers/_test/support.ts`); the tests record the HTTP
  bodies and check what they say (models.md § Tests). A change to projection or replay rules raises the provider's
  `revision`.
- **Trace**: `src/trace/_test/helpers.ts` runs turns on `createSession` (with segmenting model bindings) and reads them
  back with `readTrace`.

## Pitfalls

- Hashing is synchronous (core's MurmurHash3); on a simulated clock, wait for a condition (`until` in
  `harness/_test/helpers.ts`, `settle()`) or `runUntil`, not a fixed number of ticks.
- Some cases still wait a fixed number of event-loop rounds and retry under load (`.v6/ai/plan.md`, "The suite waits on
  conditions"); prefer waiting on the condition in new cases.
- Tests inside the package import by relative path (`../../testing/index.js`); other packages by name
  (`@molstar/ai/testing`), which Vitest resolves to sources through the `molstar-src` condition.
