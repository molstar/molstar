# Conformance library {#conformance}

- **Implementations:** [TypeScript](typescript.md)

Rules §conformance-suites: a small, composable library for conformance checks; suites are plain data, cost is a tag, composition is ordinary function composition; no registry, no mandatory manifest.

## Problem

A port's behavioral contract written only in prose is read differently by each adapter, and a bug found in one adapter
is fixed there alone. Shared checks run only at one cost either slow the gate or never run their expensive cases. A
case that needs an optional feature, silently passed by an implementation that lacks it, claims conformance that is not there.

## Shape

| Name | Shape |
|---|---|
| `Cost` | `quick` < `medium` < `guard` < `full`, the test tiers (§test-tiers) |
| `Case<S>` | `title`, `cost`, optional `needs` (an optional capability, §optional-capabilities), `run(subject) -> nothing \| async` |
| `Suite<S>` | `id` and an ordered list of cases |
| `kase(title, cost, run, needs?)` | builds a case |
| `suite(id, cases...)` | builds a suite |
| `combine(id, suites...)` | one suite with all cases, in order |
| `adapt(suite, id, pick)` | the same cases on a part of a bigger subject: each runs on `pick(subject)` |
| `upTo(suite, maxCost)` | the cases whose cost is at most `maxCost` |
| `plan(suite, {upTo?, caps?})` | per case: run, or skip with a reason (`cost`, or `needs <cap>`); pure |
| `register(suite, make, {label?, upTo?, caps?})` | the binding to a test runner: one test per case, a fresh subject from `make` per case |

A suite is named `<subject>Suite` (`documentStoreSuite`, `commandCatalogSuite`). That is the whole library; everything
else is ordinary code using these pieces.

**Packaging.** The core (everything but `register`) is runner-free and is the library's main entry. `register`, and any
suite that uses the runner's assertions, sit behind a separate entry, so a runtime entry point never reaches a test
runner (one that does breaks the app that imports it; a story crashed this way); a gate check follows each package's runtime entries to prove it (§suites-next-to-ports). A suite written with
plain checks that raise may live in the main entry.

**Where suites run.** A library ships its suite beside its port; every implementation registers it at the cost that
fits (the in-memory adapter in a quick file, a provider adapter in a full file against a real deployment). An
application lists, in an ordinary test file, which suites run on which parts of its composition
(§applications-reuse-suites): the local composition at the gate, the deployed one nightly and after deploy; it derives
the list from the wiring when that is cheap (§enumerable-bindings).

## Invariants

- Suites and cases are plain values; there is no global registry and no manifest.
- `upTo` keeps the suite's id and the order of the cases it keeps.
- `combine` keeps every case of every suite, in argument order.
- `adapt` keeps titles and costs; only the subject changes.
- `plan` and `register` agree: a case is skipped for cost when its cost exceeds the limit (default `full`), otherwise
  for a capability when it `needs` one not in `caps`.
- A skipped case shows as skipped, never as passed (§optional-capabilities).
- A case that needs a capability carries it in its displayed title.
- Every case gets a fresh subject; no state is shared between cases (an application uses a fresh composition per case).
- No runtime entry point of a package imports a test runner.
- Every port in the common libraries has a suite that every implementation passes; a bug found in one adapter becomes
  a case for all (§suites-next-to-ports).
- Cases for anything with implementations in several languages are stored as data fixtures next to the module, and
  every language's runner executes the same fixtures (§cross-language-cases).

## Conformance cases

- Cases with costs `[quick, medium, guard, full]`, `plan(upTo: medium)` → `[run, run, skip cost, skip cost]`.
- No `upTo` → no case is skipped for cost.
- A `quick` case that `needs: "transactions"`, `caps: []` → skip `needs transactions`; `caps: ["transactions"]` → run.
- A `full` case that needs a capability not claimed, `upTo: quick` → skip `cost`.
- `upTo(s, quick)` → the same id, only the quick cases.
- `combine("X", a, b)` → id `X`, `a`'s cases then `b`'s.
- `adapt(s, "Y", t -> t.store)` run on `t` → each case of `s` ran on `t.store`.
- `register` with three runnable cases → `make` called three times.
- The `DocumentStore` suite: `put {id: "a", title: "x"}` then `get "a"` → the same; `get "nope"` → absent; `watch`
  emits the current value first, then changes (quick); 10k puts scale ~linearly (guard); writes in a transaction are
  atomic (medium, needs `transactions`); `watch` delivers in write order under load (full, needs `ordered-watch`).

## Other suites follow the same shape

`observableSuite` (§api-surfaces), `executorSuite` (reconciles, ignores unknown nodes, disposes cleanly; §representations), `encodingSuite` (every atom kind including `-0`, NaN and bigint round-trips; [§representations](../../../representations/bootstrap-patterns/representations/README.md)), `commandCatalogSuite` (§commands), `vocabularySuite`.

## Open

- `Cost` has no `smoke` (or `eval`) level, though §test-tier-suffixes says conformance cases carry the same cost tags as
  the tiers; what a smoke file registers up to is not said.
- §cross-language-cases asks for shared fixtures, but the library shows no way to load fixture data into cases.
