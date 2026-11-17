# Test tiers and targeted runs {#test-tiers}

- **Implementations:** [TypeScript](typescript.md)

Rules §tests, §suites-are-data: tests cost about what their subject costs; the tier (quick / medium / guard / smoke / full) is in the filename; targeted runs go first; the gate runs the build, the unfiltered quick, medium, guard and smoke tiers and the repository-wide checks; nightly runs everything.

## Problem

When every test runs on every change, slow tests cripple the agent loop as the codebase grows; when tiers are marked
by tags or config, nobody can tell a file's cost at a glance and a slow test hides in a fast tier. Tests that wait on
fixed sleeps flake under load, and checks that live only in CI differ from what people and agents run locally.

## Shape

**The tier is the file suffix**: `<name>.<tier>.test.<ext>`, kept next to the code it tests.

| Tier | What it runs | When |
|---|---|---|
| `quick` | unit, in-memory | constantly; at the gate |
| `medium` | module compositions, in-memory | at the gate |
| `guard` | allocation, heap and scaling budgets (§performance) | at the gate |
| `smoke` | the built artifact on the local composition, no external services; in `{apps,services}/<product>/<name>/smoke/` | at the gate |
| `full` | slow or external (network, real provider) | nightly, and after each deploy |
| `eval` | a real model or other metered service (§evals) | only with an acknowledgement and a budget |

**One runner**, `test <tier | set> [dir]`, behind the repository's task runner:

- a tier name runs that tier; a directory makes it a *targeted* run, scoped to that directory;
- sets: `gate` = quick, medium, guard, smoke; `nightly` = gate + full; sets always run unfiltered;
- evals are not a tier of this runner: `evals <path> --ack --runs N --budget <usd>` refuses to start without all three,
  stops when the budget is spent, and writes results to gitignored scratch (§machine-invariant).

Typical loop: targeted first (`test quick lib/storage`), then wider (`test quick`), then `gate` before merge.

**Conformance suites by tier**: a file registers a suite (§conformance) up to its own tier, so the same suite runs its
quick cases in a quick file and all cases in a full file.

**The whole app** (§whole-composition-tested):

- *Composition tests* (medium): the real wiring on in-memory adapters (`local()`, §ports-and-adapters), driven through
  commands (§commands) and observed through read commands.
- *Smoke tests*: start the built artifact with the local composition and an ephemeral port, then call one real read
  command per major area; every shipped target (worker bundles too) is built and tested.
- *Provider-backed checks*: the full-cost cases of the same conformance suites.

**Guard tier**: allocation is measured after a warm-up, with a collection before measuring and none after (garbage made
in the loop must still be counted); scaling compares the best of several runs at two input sizes. Guards run alone,
after the other tiers. A timing case runs for several milliseconds at n (a sub-millisecond case measures the
scheduler); a memory case takes its baseline after the collector has settled (collect, yield a macrotask, collect
again), since stores freed late read as growth; a positive control plants growth that stands clear of the case's
noise on a loaded machine, several times the bound. A scaling guard drives the path as the app calls it (its
operation, its event, its frame), at sizes where the per-item part dominates the fixed cost, and its positive control
injects a scan into that same path, not into a separate function. A guard is never quieted by a raised bound or a
dropped control. Why: guards that flake under a loaded gate get ignored, and a guard timed on a function the app
does not call, or at sizes where the fixed cost dominates, stays green on a quadratic path.

## CI layer

CI calls the same tasks people and agents run; there is no CI-only logic (§ci-enforces, §one-entry-point). Toolchains
come from the entry point's pinned configuration (mise in the examples), so adding a language changes the tasks, not
the workflow.

```yaml
# .github/workflows/ci.yml
name: ci
on:
  pull_request: {}
  push: { branches: [main] }
  schedule: [{ cron: "0 2 * * *" }]
jobs:
  gate:                                              # every change
    if: github.event_name != 'schedule'
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
        with: { fetch-depth: 0 }                     # the drift check diffs against the merge base
      - uses: jdx/mise-action@v2                     # installs the pinned toolchains from mise.toml
      - run: mise run setup
      - run: mise run gate                           # build; quick + medium + guard + smoke; id, path and merge checks; drift check (§spec-mirror)
  nightly:
    if: github.event_name == 'schedule'
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: jdx/mise-action@v2
      - run: mise run setup
      - run: mise run nightly                        # all tiers, including full (real providers)
```

## Invariants

- A test file's tier is its suffix and nothing else.
- The gate is build, the unfiltered quick, medium, guard and smoke tiers, and the repository-wide checks (the package's
  three and the project's own, §the-gate); nightly adds full; evals run in neither.
- `gate` and `nightly` cannot be scoped to a directory.
- There is no global time limit; a test slow relative to the code it exercises is a defect, fixed or moved up a tier.
- A file registering a conformance suite runs its cases only up to the file's own tier.
- Smokes use no external services and stay tiny.
- Tests wait for the condition they need (a value, an event, a settled promise), never a fixed sleep; budgets for
  loading code are sized for a busy machine, seconds rather than one.
- The guard tier fails rather than passes when it cannot measure (for example, collection is not exposed).
- CI runs the same tasks as local runs.

## Conformance cases

- `test quick` → only files with the `quick` suffix; `test medium lib/x` → only `medium` files under `lib/x`.
- `test gate` → quick, medium, guard, smoke; `test nightly` → those plus full.
- `test gate lib/x` → refused, usage error, nothing runs.
- `test bogus` → "unknown tier or set", usage error, nothing runs.
- `evals <path>` missing `--ack`, `--runs` or `--budget` → refuses to start.
- Guard: after 100 warm-up frames, 1,000 frames grow the heap by less than 16 KiB.
- Guard: building an index for 40,000 items takes less than 8× the time for 10,000 (best of 5); quadratic would be ~16×.
- Composition: `doc.create {title: "a"}`, `doc.rename {id, title: "b"}`, `doc.get {id}` → `{title: "b"}`.
- Smoke: the built server on the local composition answers one of the project's read commands per major area (e.g.
  `doc.list`, §commands) with success.
