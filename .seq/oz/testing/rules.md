# Testing rules {#testing-rules}

Test tiers, application-level checks, CI, and the conformance library.

## In brief

**Tests** (§tests).
- A test costs about what the code it exercises costs; a slow test relative to its subject is a defect.
- Tiers by file suffix: **quick** (`*.quick.test.ts`, unit, in-memory, seconds), **medium** (compositions, still
  in-memory), **guard** (allocation, heap and scaling budgets), **smoke** (the built artifact, no external services),
  **full** (everything, including real services; nightly), and **eval** (a real model or metered service: never in the
  gate or nightly by default; only with explicit acknowledgement and a stated budget: an `evals` task taking the path, `--ack`, `--runs N` and `--budget <usd>`).
- Targeted first (the module's own tests), then widen. The gate = build + quick, medium, guard, smoke; nightly = full.
  Run everything through the entry point's `test` and `gate` tasks; CI runs the same ones.
- The final composition is tested whole: composition tests drive the app through commands on `local()` and check
  representations; smoke tests boot the built artifact, every target it ships (worker bundles included).
- A new test is mutation-checked once (break the code, see it fail). Responses consumed by other modules or agents have
  pinned shapes.
- Wait on conditions, never fixed sleeps; size load budgets for a busy machine; a lane's "gate green" is re-run.
- Evals sit next to the code as their own tier, run N of M times with near-miss cases for anything that acts unasked,
  results in `.tmp/`, scored over time. Apps are headless and their state reconstructible, so an eval viewer can show the
  UI where an agent's run left it (§replay-in-ui; aimed at, not built).

**Conformance** (§conformance-suites).
- Every port ships a suite (plain data plus functions) next to it; every adapter, including in-memory ones, passes it.
  Suites never reach a runtime entry: the runner binding and runner-using suites sit behind a subpath, gate-checked.
- Applications rerun the same suites against their compositions. Optional capabilities are declared; a case needing one is skipped visibly, never silently passed.
- Conformance is the contract for agents: a second implementation is done when it passes the suite.
- Cross-language cases are data (fixtures), so another language's implementation runs the same cases.

## Testing {#tests}

1. **Tests are fast in proportion to what they test**: {#proportionate-tests} a test should cost about what the code it exercises costs. Slow tests cripple agent loops as the codebase grows. There is no fixed global time budget; a test that is slow relative to its subject is a defect to fix or move to a higher tier.
2. **Tiers** {#test-tier-suffixes} (conformance cases carry the same cost tags, §suites-are-data). → [§test-tiers](bootstrap-patterns/test-tiers/README.md)
   - **quick**: unit, in-memory, seconds. Runs constantly during work.
   - **medium**: compositions of modules, still in-memory.
   - **guard** ("heap guards"): allocation, heap and complexity-scaling budgets (§performance).
   - **smoke**: a tiny set run against the built artifact (§whole-composition-tested).
   - **full**: everything, including anything slow or touching external services, except evals.
   - **eval**: checks that call a real model or another metered service (§evals). A separate tier: never in the gate, never in nightly by default, run only on explicit acknowledgement with a stated budget, because every run costs real money.
3. **Targeted first.** {#targeted-first} Every module and submodule can run its own tests; run those first, then widen.
4. **The gate** {#the-gate} is what every change must pass before merge: build, the unfiltered quick, medium, guard and smoke tiers, and the repository-wide checks: the package's id, path and merge checks ([`_checks/`](../_checks/README.md)) and the project's own, such as the drift check (§mechanical-drift-checks). **Nightly** runs the full tier. There is no third trigger.
5. **The final app composition is tested as a whole**, {#whole-composition-tested} through its real composition root (§one-composition-root):
   Two cheap practices keep tests honest: a new test is **mutation-checked** once (break the code it guards and see it fail), and API responses consumed by other modules or agents have **pinned shapes** that fail on any unannounced change.
   - **Composition tests** (medium) run the real wiring on in-memory adapters (`local()`, §compositions-are-values), driving the app through commands (§command-reachable) and checking representations (§declarative-representations).
   - **Smoke tests** prove the built artifact starts and its critical paths work (boots, loads a representative state, runs one command per major area, renders a frame), using only adapters that need no external service. They run at the gate.
   - **Smoke tests build every target the app ships**, including worker bundles: package export conditions and globals differ there (a worker has no `document`), and only a test of the built artifact sees it.
   - **Provider-backed checks** (real adapters, real services, the full-cost cases of their conformance suites) are full tier: nightly, plus after each deploy. Checks that spend money per run (a real model) are evals (§evals).
6. **CI is the enforcement layer.** {#ci-enforces} CI runs the same scripts agents run locally (one entry point per tier, no CI-only logic): the gate (§the-gate, which includes the merge check and the drift check) on every change; the full tier nightly. A failing gate blocks merge. Flaky tests are fixed or quarantined with an issue, never ignored. Most flakes are timing under load: a test waits on the condition it needs (a value, an event, a settled promise), never on a fixed sleep, and a budget for loading code (a lazy chunk, a cold start) is sized for a busy machine. A lane's claim that the gate is green is checked by re-running it (§upward-authority), ideally while other work loads the machine.

7. **Evals measure behavior that tests cannot pin down** {#evals} (agent workflows through the project's agent tools; model-backed features; quality of results). They are their own tier (§test-tier-suffixes), kept next to the code they exercise like any test (a suffix marks the tier), and run on the in-memory composition with only the model real (§compositions-are-values). They run **only on explicit acknowledgement with a stated budget** (who asked, how many runs, the cost ceiling; one `evals` task with an acknowledgement, `--runs N` and `--budget <usd>`), never in the gate and never by default, since they spend real money; results are recorded (gitignored scratch, §machine-invariant) and scored over time rather than pass/fail. Each case runs several times (N of M), and every behavior that acts without being asked has **near-miss cases**, inputs that must not trigger it: a model offered a ready action applies it where it should not, which no scripted test shows. A deterministic eval case that guards a fixed behavior is promoted to a test tier.
8. **Agent runs can be replayed in the UI.** {#replay-in-ui} Every app has a headless mode (§headless-mode) and its model and view state are reconstructible from data (§reconstructible-state), so a viewer in `tooling/` can render the app where an agent's run left it, from the run's recorded state, beside its transcript and scores.

## Conformance {#conformance-suites}

A small, composable **conformance library**: reusable checks for "does this correctly implement or compose X". It is a testing helper, not a framework. → [§conformance](bootstrap-patterns/conformance/README.md)

1. **Suites are plain data plus functions.** {#suites-are-data} A suite is a named list of cases, each tagged with a **cost** (its tier, §test-tier-suffixes) and, if any, the optional capability it requires (§optional-capabilities). Suites compose freely (combine, adapt to a larger subject, filter by cost), so the same checks run cheaply at the gate and exhaustively nightly.
2. **Libraries ship suites next to their ports.** {#suites-next-to-ports} Every port or abstract capability in the common libraries (§compose-through-ports) has a suite that every implementation (in-memory, provider adapters, executors, encodings) runs. The behavioral contract lives in the suite, not in prose. A new adapter is done when it passes; a bug found in one adapter becomes a new case for all. A suite never reaches a runtime entry point: the suite core is runner-free, the binding to a test runner and any suite that needs one sit behind a separate subpath, and a gate check follows each package's runtime entries to prove no test runner is imported.
3. **Applications reuse the same suites** {#applications-reuse-suites} against their bindings as composed, at the cost level that fits the target (local composition at the gate, deployed one nightly). Enumerate bindings from the wiring where cheap (§enumerable-bindings).
4. **Optional capabilities are explicit.** {#optional-capabilities} A case that needs an optional feature is skipped visibly when an implementation does not claim it, never silently passed.
5. **Conformance is the contract for agents.** {#conformance-contract} When an agent adds an adapter, executor, encoding or language, the first step is to run the relevant suite; its failures define the remaining work.
6. **Cross-language cases are data.** {#cross-language-cases} Conformance cases for anything with (or likely to get) implementations in several languages, such as codecs, parsers, query semantics, operation application and the wire protocol, are stored as fixtures (input and expected output in the project's text spelling) next to the module. Each language's runner executes the same fixtures, so every implementation is held to one contract (oz's `language-neutral-contracts` (not adopted here)).
