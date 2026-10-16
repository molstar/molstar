# Dead and speculative code

Method: every exported symbol of every entry point (and the deep imports found in the repository) was classified by
who imports it: outside `packages/ai`, another AI package, tests only, its own package only, or nobody (a TypeScript
checker scan over import declarations; the components package's dynamic imports were then checked by hand). Counts
(symbols: outside / other AI package / tests only / own package / nobody):

| Package        | Symbols | Outside | AI | Tests | Own | Nobody |
| -------------- | ------: | ------: | -: | ----: | --: | -----: |
| agent          |      31 |       3 | 25 |     0 |   3 |      0 |
| harness        |      99 |       1 | 52 |    24 |  21 |      1 |
| components     |      10 |       1 |  2 |     7 |   0 |      0 |
| model-adapters |      20 |       7 |  0 |     5 |   8 |      0 |
| conformance    |      45 |       2 |  0 |    12 |  28 |      3 |
| tools          |      60 |      11 | 13 |    13 |  22 |      1 |
| session        |      73 |       8 |  0 |    37 |  27 |      1 |
| trace          |      54 |      30 |  0 |     7 |  17 |      0 |
| fibers         |      79 |      15 |  4 |    31 |  28 |      1 |
| narration      |      25 |       0 |  0 |    17 |   8 |      0 |

(components by hand: outside, only `createIndexedDbRows`; in the AI packages, `createMemoryStorage` and
`createInProcessSupervisor` from the session; the other seven tests only or nobody.)

"Own package only" mostly means a symbol is live inside and its export is surplus; one package with explicit entries
(packaging.md §2) drops those exports for free. The table below is what is dead or speculative as **code**.

## 1. Whole pieces with no consumer

| What                                         | Where                                            | Src | Test | Evidence                                                                 | Recommendation |
| -------------------------------------------- | ------------------------------------------------ | --: | ---: | ------------------------------------------------------------------------ | -------------- |
| **Narration** (TTS providers, WAV, Mediabunny) | `packages/ai/narration/`                       | 651 |  636 | no import outside `packages/ai` (only `prism/ui/src/_test/opt-in-boundary.test.ts:32` names it); built for a story mode that is not built (`.v6/ai/plan.md` "Also") | **kept** (decided 2026-10-16): a `narration/` folder of `@molstar/ai` (packaging.md §2) |
| **`runOnce`**, `terminalToolFor`, `Settlement` | `session/src/once.ts`                          | 161 |  122 | no caller; meant for "the command menu's assistant proposing a filled command" (`_spec/session.md` § One-off runs), not built | delete (decided 2026-10-16) |
| **Session prices**: `PriceTable`, `dollarsOf`, `creditsOf`, `ModelPrice` | `session/src/prices.ts`             |  48 |   34 | no caller; collides with the catalog's names (packaging.md §3)           | delete         |
| **`callTool`** (run a handler without a harness) | `harness/src/tools/call.ts`                  |  56 |    0 | no importer and no internal use                                          | delete         |
| **Fiber deltas**: `appendDelta`, `deltas`, delta pages, gap detection | `fibers/src/kernel.ts:278–353, 585–650`, `host.ts:597+`, `domain.ts` | ~250 | ~250 | the viewer never streams deltas (threads.md §1) | **kept** (decided 2026-10-16: streaming may need them) |
| **Fiber cleanup pages** (`nextFiberCleanupPage`, `cleanupPageSize`, `deltaFlush*`, `heartbeatIntervalMs` validated only) | `fibers/src/kernel.ts:651–663`, `limits.ts` | ~40 | — | "validated but the host does not use them" (`_spec/fibers.md` Gotchas) | delete |
| **`@molstar/ai-tools/zod`** (`zodInput`, `verbatim`) | `tools/src/zod.ts`                         | 118 |  205 | no import outside `tools/src/_test/zod.test.ts`                          | drop the entry; zod itself per the TODO note (open question 7) |
| **Runner retention inspection** (`RetentionSnapshot`, `retention()`) | `harness/src/contracts/inspect.ts` (20), `runner/state.ts:259`, `runner/index.ts:728` | ~30 | — | read only by the P10 cases (`conformance/src/cases/performance.ts:103`) | delete with P |
| **Storage meter** (`meterStorage`, `StorageCounts`) | `conformance/src/meter.ts`                  | 158 | — | the P cases only                                                         | delete with P  |
| **Conformance generators** (`FIBER_INVALID_ENVELOPES`, `fiberCleanupCursorPages`, `fiberEnvelopeIsValid`) | `conformance/src/fibers/generators.ts:23–47` | 25 | — | no importer | delete |
| Small: `ceilingsOf` (`session/src/drive.ts:38`), `definitionOf` (`tools/src/tool.ts:133`), `FIBER_TABLES` (`fibers/src/storage.ts:63`) | | ~25 | — | no importer | delete or unexport |

## 2. Speculative: built for a host that does not exist yet

| What                                  | Where                                                      | Note                                                                                       |
| ------------------------------------- | ---------------------------------------------------------- | ------------------------------------------------------------------------------------------ |
| Pricing end to end                    | catalog `price` fields and `dollarsOf` (`model-catalog.ts:241`), the runner's `price` option and `cost` limit, the session's `price` (`session.ts:144,442`), the trace's credits hook (`trace/src/read.ts:119,719`) | nothing in Prism passes a price (`git grep price` in `prism/viewer/src/agent`, `prism/assistant/src/viewer`: none). The `cost` limit is a harness contract (V6 price-cost); keep the option, delete the two unused `dollarsOf` (`.v6/ai/plan.md` §2.5 says pricing is the host's) |
| Remote deployments                    | `HarnessClient` / `RunnerPort` split, `unknown-outcome`, serializable operations (`.v6/ai/harness.md` §15) | no remote host exists; the split costs little and is the protocol. Keep; only the tests of backend faults go (conformance.md) |
| Durable takeover, leases, supervisor  | persistence.md §4, harness.md                              | removed (decided 2026-10-16)                                                                                   |
| `ModelChoice.contextWindow` / `maxOutput`, `AgentSessionOptions.contextWindow`, adapter `modelLimits` | duplicates.md §4 | three ways to state a model's limits                                         |
| Native tool additions                 | providers.md §3                                            | two catalog models                                                                         |

## 3. Legacy paths

Reading entries "written before labels" (persistence.md §5): `session/src/compaction/domain.ts:50`,
`compaction/plan.ts:129`, `session/src/session.ts:361`, `tools/src/toolkit.ts:330`, `trace/src/read.ts:56`; test
`tools/src/_test/toolkit.test.ts:374`. With nothing persisted, ~40–80 lines (inferred) and one test go.

## 4. Totals

| Item                                       | Src    | Test   |
| ------------------------------------------ | -----: | -----: |
| `runOnce`                                  |    161 |    122 |
| session prices                             |     48 |     34 |
| `callTool`                                 |     56 |      0 |
| fiber cleanup pages                        |    ~40 |      — |
| `./zod`                                    |    118 |    205 |
| retention, meter, generators, small        |   ~240 |      — |
| pre-label fallbacks                        |    ~60 |    ~15 |
| **total**                                  | **~700**   | **~350**   |

The meter, retention and generators are also counted in conformance.md's cut; the README's totals count them once.
