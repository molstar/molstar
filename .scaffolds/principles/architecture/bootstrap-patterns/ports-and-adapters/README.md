# Ports and adapters, composition root {#ports-and-adapters}

- **Implementations:** [TypeScript](typescript.md)

Rules §structure, §composition-root: common code depends on ports; provider implementations live in their own `*-<provider>-runtime` packages (§provider-runtime-packages) so their dependencies never leak; each app and service has one composition root that wires everything with plain typed code. An in-memory adapter makes everything testable with nothing deployed.

## Problem

When business code imports a provider SDK, every test needs that provider deployed, swapping the provider means
rewriting callers, and an SDK upgrade touches every package. When code constructs its own collaborators or reaches for
globals, nobody can substitute an in-memory one, and nobody owns disposing what was created.

## Shape

- **Port**: an interface in common code, named for the capability (`DocumentStore`), never for the provider. The
  example port, for documents that carry an `id`:
  - `get(id) -> async doc | absent`, `put(doc) -> async`, `delete(id) -> async`;
  - `watch(id) -> Observable<doc | absent>`: the current value first, then every change (streaming surface, §api-surfaces).
- **In-memory adapter**: ships in the same package as the port, with the port's conformance suite (§conformance).
- **Provider adapter**: in its own `<module>-<provider>-runtime` package (§provider-runtime-packages); it wraps the
  provider's reactive stream into an `Observable` inside the adapter.
- **Business code** receives its ports as constructor or function parameters, narrowed to what it uses
  (`NoteService(deps: { store, clock })`).
- **`Scope`**, the one wiring helper (no container):
  - `lazy(make, dispose?) -> thunk`: creates the instance on the first call, memoized within the scope; `dispose` is
    registered when the instance is created.
  - `child() -> Scope`: a per-session or per-request scope on top of the app-wide one.
  - `dispose() -> async`: children first, then this scope's instances in reverse creation order.
- **Compositions**, as values over one dependency shape:
  - `Deps`: a record of thunks (`store`, `docs`, `sync`, `clock`, `llm`, and services built on them such as `notes`, `doc`).
  - `Bindings`: the provider-facing subset (`store`, `docs`, `sync`, `clock`, `llm`).
  - `wire(scope, bindings) -> Deps`: builds the services lazily from the bindings.
  - `local(scope)`: in-memory bindings (memory stores, loopback sync, fake clock, scripted model); tests, stories, evals.
  - `test(overrides, scope)`: `local` with single bindings replaced by value.
  - `production(scope)`: in each app's or service's own root; binds the runtime packages and reads the environment.
  - `createApp(deps, scope) -> App { invoke, state, scope }`: the whole app as one object for tests, smokes and agents
    (`invoke` is the `Invoker` of §commands).

Where the pieces live:

| Piece | Place |
|---|---|
| Port, in-memory adapter, conformance suite | `lib/<module>/` |
| Provider adapter | `lib/<module>-<provider>-runtime/` (the only package that depends on the provider SDK) |
| `Deps`, `wire`, `local`, `test`, `createApp` | the product's compose package, `packages/<product>/compose/`; no provider dependencies |
| `production` and the entry point | each `apps/<product>/<name>/` and `services/<product>/<name>/` root |

An entry point builds a scope, calls `production(scope)`, starts, and on shutdown stops its server and then disposes the
scope.

## Invariants

- Common code (`lib/`, product packages) depends on ports only; only `*-<provider>-runtime` packages depend on a provider SDK.
- Imports point one way: apps and services → product packages → lib (an import-boundary lint, §common-vs-product).
- Every port has an in-memory adapter in the port's package and a conformance suite that every adapter passes.
- No code outside a composition root constructs an adapter, reads the environment, or imports a provider runtime package.
- One composition root per deployable; `lib/` contains none.
- No reflection, decorators, string-token lookup or service locator; no module-level singleton holds a dependency.
- A provider's stream or connection type never appears outside its adapter (§no-transport-primitives).
- `watch` reads the value at subscribe time, not when `watch` was called.
- `lazy`: `make` runs at most once per scope, and not at all if the thunk is never called.
- `dispose`: children are disposed before the scope's own instances; instances in reverse creation order; a failing
  disposer does not stop the others; all failures are raised together at the end; a disposed child is no longer held
  by its parent.
- `test` overrides are memoized like every other binding.
- A list of bindings (for conformance runs or docs) is derived from `wire`, never kept by hand (§enumerable-bindings).

## Conformance cases

- Store: `put {id: "a", title: "x"}`, then `get "a"` → `{id: "a", title: "x"}`; `get "nope"` → absent.
- Store: `w = watch("a")`, then `put a`, then subscribe to `w` → the first value is the stored `a`.
- Store: subscribe to `watch("a")`, then `put a`, `delete a` → absent, `a`, absent.
- `lazy` thunk called twice → `make` ran once; both calls return the same instance.
- `lazy` thunk never called, then `dispose` → neither `make` nor its disposer ran.
- Instances created in order A, B → disposed B, then A.
- Parent instance P and child instance C → C is disposed before P.
- Disposers of A and B both fail → both ran; `dispose` raises one aggregate error with both failures.
- Child disposed on its own, then parent disposed → the child is not disposed a second time.
- `test({ clock: c })` → `clock()` is `c`'s instance on every call; the other bindings are the in-memory ones.

## Open

- Calling a `lazy` thunk after its scope is disposed returns the already-disposed instance; neither the code nor the
  rules say whether that should be refused.
