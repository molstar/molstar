# Architecture rules {#architecture-rules}

Structure, boundaries, wiring, errors, commands, performance and access.

## In brief

**Structure** (§structure). `lib/` is common and knows no product; `packages/<product>/*` are a product's building blocks;
`apps/` and `services/` are thin compositions; `tooling/` never ships (oz's `top-level-layout` (not adopted here)). Dependencies point lib ← packages ← apps/services.
Every boundary is a port (an interface) with adapters; a provider SDK lives only in its own `<module>-<provider>-runtime`
package. Promote to `lib/` on the 2nd-3rd real use, not before. `lib/core/` depends on nothing else.

**API surfaces** (§api-boundaries). Sync for cheap and pure, a promise for one result, an observable for many. State is a *behavior*
(current `.value` plus changes). Never expose a transport (WebSocket, SSE) through an API.

**Composition** (§composition-root). Dependencies are passed in, typed as narrow ports; no globals or singletons. `production()` lives
only in each deployable's `compose.ts` (the only reader of env and importer of provider runtimes); `local()` and
`test(overrides)` live in a product package so stories, tests and evals share them; `lib/` has no root. No heavy
container or service locator: typed factories plus a small `Scope` helper. Whoever creates a resource disposes it.

**Errors** (§stable-error-codes). Every error a caller can handle has a stable namespaced code (`realm.not_found`), defined with
`defineErrors`, cataloged, with an actionable hint. Codes are API. Cancellation and limits have codes. Details carry no
secrets. Programmer errors stay plain exceptions.

**Commands** (§command-reachable). One path to behavior: UI, agents, tests and scripts call the same data-described commands; queries are
commands too. A command runs inside the transaction it is given and never commits itself. Mutating commands support dry-run, which
the transaction gives as a preview, so a command needs no dry-run branch of its own; idempotency keys where retries are
plausible; destructive commands require confirmation; effects outside the documents are outbox entries performed after
commit; a write that changes nothing says so, never plain
success (§no-op-writes-say-so). Commands are safe to automate, least-privilege (one Actor type; the policy sees the parsed input; agents act within
their verified requester's authority, no allow-list), and every invocation is recorded as data. The registry is
the machine-readable spec; a feature with no command is a gap.

**Performance** (§performance). No accidental quadratics; pick structures for the real access pattern; zero allocations in
steady-state hot loops; budgets are guard-tier tests, not intentions.

**Access and identity** (§access). One policy port `(actor, action, locator) → allow | deny | filter`, asked at every read,
write, subscription and tool call; policies belong to the application, deny by default. Denials are coded and do not leak
existence (uniform unless the decision is public). Sync and catch-up are boundaries too. Every change records its actor
and, for an agent, whom it acts for. Agents act only on *requests* from actors allowed to direct them; everything
else is *material*. Effects nobody requested (inferred edits, agent-initiated view changes, watches) follow one
consent model: apply / propose / off, per workspace, person and kind of effect; an unrequested effect quotes its reason
from the material, checked mechanically, and runs with the authority of the person whose consent allows it. An agent's
reads and writes are checked as both the agent and the person it acts for. Identity is a port (local picker, password, OIDC, service credentials as adapters); the session is
reactive state; sign-in and sign-out are commands; credentials never enter model state, views or agent context; the
server derives the actor from verified credentials only; a served picker is not authentication. The base ships no
roles: roles are application policy. A space can be declared unrecorded: its writes leave nothing behind (never
logged, persisted, derived from or indexed), a guarantee only the base can give (§unrecorded-spaces).

Patterns: §ports-and-adapters, §api-surfaces, §hot-loop-allocations, §commands, §error-codes (bootstrap).

## Structure: common and product {#structure}

1. **Split common from product, and compositions from building blocks.** {#common-vs-product} `lib/*` holds common modules: abstract, configurable, no product knowledge. `packages/<product>/<module>` (package `@scope/<product>-<module>`, oz's `family-folders` (not adopted here)) holds a product's building blocks (models, commands, adapters, runtimes, UI). `apps/*` (user-facing) and `services/*` (deployable backends) hold only top-level compositions: a composition root, entry points, configuration and deployment, no building blocks. Dependencies point one way (apps and services → product packages → lib, never the reverse), enforced by an import-boundary lint rather than by convention.
2. **Compose through ports.** {#compose-through-ports} Common code depends on interfaces (ports/clients), never on concrete providers. The application layer implements the ports and wires them in a single composition root (§composition-root). → [§ports-and-adapters](bootstrap-patterns/ports-and-adapters/README.md)
3. **Every dependency boundary is an interface**: {#interface-boundaries} cloud services and storage (a hosted database, a cloud provider, a sandbox host, ..., on client and server alike), transport, and rendering/compute backends. Provider code never leaks past its adapter; swapping a provider means writing one adapter.
4. **Everything is testable in isolation.** {#testable-in-isolation} Each module and each app composition runs on in-memory implementations with nothing deployed, which is what makes local tests and evals cheap.
5. **Not everything must be reusable.** {#not-everything-reusable} Product-specific logic stays in the product's packages.
6. **Promote on the 2nd-3rd real usage**: {#promote-on-reuse} a product package module moves to `lib/` when a second product needs it. Exception: foundational utilities (e.g. the tree representation and its parsers, `lib/repr/core`) start in `lib/`.
7. **Provider implementations are packages of their own, to isolate dependencies.** {#provider-runtime-packages} A module never depends on a provider SDK; its implementation on a provider is a separate package named `<module>-<provider>-runtime`: `lib/storage-<provider>-runtime` for a common module, `packages/<product>/<provider>-runtime` for product-specific code. When the provider hosts a product's backend as its own deployable (server functions and their schema, a worker image), that deployable is a service, `services/<product>/<provider>-runtime`, composing the `lib/*-<provider>-runtime` pieces with the product's own. Installing or upgrading an SDK then touches only the packages that use it.
8. **`lib/core/` is the base layer.** {#lib-core-base} What nearly every other package uses (streams and behaviors, DI scopes, coded errors, the conformance library) lives in the `lib/core/` family. Core packages depend only on each other and on no provider; everything else in the workspace may depend on them. Something joins core when most packages need it and it has no dependency outside core; foundational modules with their own domain (the representation core, `lib/repr/core`) stay in their family.

## API boundaries {#api-boundaries}

Every API uses the surface that fits what it does, and the three are never mixed up. → [§api-surfaces](bootstrap-patterns/api-surfaces/README.md)

| Surface | Shape | Use for |
|---|---|---|
| Synchronous | plain return value | cheap, pure, in-memory |
| Asynchronous | promise-like | one request, one result |
| Streaming | observable-like | many values over time |

1. {#streaming-library-irrelevant} The streaming library is irrelevant (rxjs, Effect, own); the pattern is not.
2. {#return-observables} Return an observable rather than accepting a callback.
3. {#events-vs-state} Distinguish **events** (plain observable) from **state**: state is exposed as a *behavior*, an observable that always has a current value, readable synchronously via `.value`. Models expose read-only behaviors and keep the writable one private. → [§api-surfaces](bootstrap-patterns/api-surfaces/typescript.md#behavior-subject-observe-the-current-value)
4. {#no-transport-primitives} Never expose a transport primitive (WebSocket, a provider's reactive query, SSE) through an API. Wrap it in an observable inside an adapter.

## Composition root and dependency injection {#composition-root}

How the pieces from §structure get wired together. → [§ports-and-adapters](bootstrap-patterns/ports-and-adapters/typescript.md#composition-root-and-lightweight-di)

1. {#dependencies-passed-in} Dependencies are passed in (constructor/function parameters, typed as ports). No code constructs its own collaborators or reaches for globals or singletons.
2. **One composition root per deployable**: {#one-composition-root} each app and each service (and each test/eval composition) has its own, and it is the only place that knows concrete implementations and wires them together. Shared, provider-free wiring (the dependency shape, `local()`, `test()`) lives in a product package so stories, evals and tests reuse it; `production()` lives in each root, which is the only code that imports provider runtime packages (§provider-runtime-packages). `lib/` never contains a composition root.
3. **No heavy container**: {#no-heavy-container} no reflection, decorators, string-token lookup or auto-wiring. Wiring is plain typed code (factories and a small helper for lazy, memoized, disposable instances). A service locator (resolving dependencies from inside business code) is forbidden.
4. {#narrow-dependencies} Dependencies are explicit and narrow: depend on the port(s) actually used, not on a bag of everything.
5. {#explicit-lifetimes} Lifetimes are explicit: app-wide, per-session, per-request; whatever creates a resource owns disposing it (scopes dispose their instances in reverse order).
6. **Alternative compositions are values**: {#compositions-are-values} `production()`, `local()` (in-memory), `test(overrides)`, so tests, stories (§stories-on-real-architecture) and evals swap bindings by passing different ones, not by patching modules.
7. {#enumerable-bindings} Where useful, the composition can enumerate its bindings (for conformance runs, §applications-reuse-suites, or docs); derive that list from the wiring rather than maintaining it by hand.

## Errors have stable codes {#stable-error-codes}

Errors that cross a boundary (parsing, repair, validation, queries, commands, ports, APIs) are structured and machine-readable. → [§error-codes](bootstrap-patterns/error-codes/README.md)

1. **Every such error has a stable, namespaced code** {#namespaced-codes} (`repr.parse.unexpected_end`, `query.unknown_function`, `command.not_found`, `store.conflict`), a human message, and optional `hint`, position (a locator: a reference, as data, to where in the input or document) and `details`. Callers (UI, agents, tests, retries) branch on the **code**, never on message text.
2. **Codes are API.** {#codes-are-api} They are never reused or silently renamed; a retired code stays reserved. Messages and hints may improve freely.
3. **One shared shape** {#one-error-shape} (`CodedError` in the common library) for all layers, so a parse error from a command input, a query error and a port error look alike and serialize alike. Errors are plain data and can be represented as a node (§declarative-representations), streamed (§api-boundaries) and logged.
4. **Hints are actionable**: {#actionable-hints} they name what to try next (the known function names, the command to list candidates, the expected arity), which is what lets an agent self-correct (§safe-to-automate, §parser-repairs-slips).
5. **Codes are cataloged** {#cataloged-codes} in one place per module (code, description, retryable?, example). The catalog drives docs (§rendered-docs) and one cheap check: every thrown code is declared. Example tests per code are welcome where cheap, not required.
6. **Cancellation and limits have codes too** {#cancellation-codes} (`cancelled`, `deadline_exceeded`, `limit_exceeded`); aborts are not reported as failures.
7. {#programmer-errors} Programmer errors (broken invariants, impossible states) stay plain exceptions/assertions; coded errors are for conditions callers can reasonably handle.
8. **Error details carry no secrets.** {#no-secrets-in-errors} Values that may be sensitive (input fragments, row contents, tokens) are left out of messages and details or redacted where errors are created, so logs, wire frames and agent contexts can carry errors freely.

## Business logic is reachable by commands {#command-reachable}

Everything (or almost everything) the application can do is reachable through a **command**: a named, typed, self-describing operation. UI, tests, scripts and agents all drive the system through the same commands. → [§commands](bootstrap-patterns/commands/README.md)

1. **One path to behavior.** {#one-path-to-behavior} UI actions, API endpoints, agent tools and CLI all call the same command; none has private logic. A model action (§actions-are-methods) is a thin, typed wrapper over a command, or is the command's handler.
2. **Commands are data-described.** {#data-described-commands} Each has a name, an input schema, an output schema, a short description, and a declared effect (read, mutate, destructive). From the registry you generate docs, a CLI, an agent tool list (e.g. MCP) and validation. Schemas evolve like codes (§codes-are-api): additive changes freely, breaking changes get a new name.
3. **Queries are reachable too.** {#reachable-queries} State inspection (what is selected, what the model contains, the status of job X) is exposed the same way, as read-only commands or streaming queries (§api-boundaries; through the project's query language, if it has one), so an agent can observe as well as act.
4. **Near-native agent access.** {#near-native-agent-access} Through skills and a small read/write tool surface on top of the registry, agents operate the real system rather than a simulation of it, and can test and debug it the way a user would. The registry, the dispatch and a transport-agnostic invoker belong in common libraries; skill text, transport, auth and app-specific safety rules are app-specific.
5. **Commands are safe to automate.** {#safe-to-automate} Dry-run for mutating commands, idempotency keys where retries are plausible, structured errors with stable codes (§stable-error-codes). Several commands that must succeed together run in one transaction (§optimistic-transactions); a command runs inside the transaction it is given and never commits on its own.
6. **Exceptions are explicit.** {#explicit-exceptions} Purely visual or ephemeral interactions (hover, drag in progress, scroll) need not be commands; their *committed* result (the final move) does (§committed-vs-transient). If something is not command-reachable, say why in the spec.
7. **The registry is the machine-readable part of the spec.** {#registry-is-spec} `_spec/` holds rationale and links to the registry instead of copying it (§spec-folders); a feature with no command is a gap (§recursive-loop).
8. **Callers act with least privilege.** {#least-privilege} A policy layer decides per actor (the actor of §recorded-actors, one type shared by commands, access and transactions) and per target (it sees the parsed input) which commands are allowed; an agent acts within its verified requester's authority, not through a separate allow-list (§requester-authority); destructive commands require confirmation; secrets never travel through command inputs, outputs or representations.
9. **Invocations are recorded as data** {#recorded-invocations} (name, actor, input hash, result code, duration). One record gives audit, replay, debugging and basic observability.
10. **A write that changes nothing says so.** {#no-op-writes-say-so} A mutating command whose input changes nothing (a field that does not apply to the target, a value equal to the current one, an id that addresses a part the command does not edit) returns a coded "unchanged" result or is refused with a hint; it never reports plain success. A caller that loops on success, an agent above all, otherwise believes it did the work and moves on.

## Performance {#performance}

1. {#no-accidental-quadratics} No accidentally quadratic behavior (unless N is provably small; say so in a comment).
2. **Choose efficient algorithms and data structures for the context.** {#efficient-structures} Pick them for the actual access pattern, sizes and lifetime (lookup vs scan vs ordered iteration, read-heavy vs write-heavy, bulk vs incremental, cache locality), not by habit. When the choice is non-obvious, or a simple one is knowingly chosen, state the assumption in a comment (expected N, hot or cold path).
3. {#scratch-arenas} Minimize allocations generally. For temporary working memory use **scratch arenas** (bump/linear allocators with mark/release or per-frame reset) instead of per-call allocations; allocate long-lived buffers once and reuse them. Arena memory never outlives its scope. → [§hot-loop-allocations](bootstrap-patterns/hot-loop-allocations/typescript.md#scratch-arenas)
4. **Zero allocations in steady-state render/hot loops.** {#zero-alloc-hot-loops} Preallocate, reuse, pool, use arenas. → [§hot-loop-allocations](bootstrap-patterns/hot-loop-allocations/README.md)
5. {#performance-tested} Performance rules are enforced by tests (the guard tier, §test-tier-suffixes), not by good intentions.

## Access: enforcement in the base, policy in the application {#access}

Who may see and change what is decided by the application; that every read and write is checked is guaranteed by the base.

1. **One policy port, asked everywhere.** {#one-policy-port} Every boundary that reads or writes asks the same port: `(actor, action, locator)` → allow, deny, or allow with a filter. Query results are filtered, subscriptions and observations are scoped, transaction steps are authorized, agent tool calls pass through it (wherever the project runs agent tools, and no agent rule can override it). Sync is a boundary too: each connection is scoped per document at connect, at catch-up and on every change, so a membership change takes effect live, and catch-up answers nothing about documents the actor was never given. There is no path that skips it.
2. **Policies are the application's.** {#application-policies} Equal visibility, private branches, delegated capabilities the delegate has accepted, separation of duties ("the author of a change never approves it"), consent before data is exported: all are policies over actors, actions, locators, data classes (labels the application attaches to values, carried to what is derived from them) and recorded provenance. The base ships no policy beyond deny-by-default for undeclared actions.
3. **Decisions are explainable and recorded, without leaking.** {#explainable-decisions} A denial is a coded error; the rule that denied it is named in the record and to the operators, but to the caller only when naming it reveals nothing the caller may not know (refusals that would reveal whether something exists are uniform; a denial carries whether it may be shown, and the default is not). Decisions that matter (who saw what, who changed what) are recorded like invocations (§recorded-invocations).
4. **Unrecorded spaces are a base mechanism.** {#unrecorded-spaces} A space can be declared unrecorded: its writes go to an ephemeral store, are never logged, persisted, derived from or indexed, and leave nothing behind when the space closes. Policies decide where such spaces are allowed; only the base can guarantee that nothing leaks into history.
5. **Actors are recorded on every change.** {#recorded-actors} Every operation, transaction, log entry and commit carries its actor (person, team, agent, service; and, for an agent, the person it acts for), so policies about who did what (separation of duties, reputation, review) can be expressed and checked by the application.
6. **Requests and material are different things.** {#requester-authority} An agent acts on *requests* from actors the policy allows to direct it; everything else it reads (documents, tool output, web pages, other agents' messages, chat from people who may not direct it) is *material*. The distinction is policy at the boundary (who may direct which agent, §one-policy-port) and is carried into the agent's context as provenance, not left to prompt wording. Authors of messages an agent treats as requests are verified by the server, never taken from the content. Provenance reaches the model as structure (each message one record with its verified author, standing and an escaped text), and every read, dry-run and write the agent makes for a request is checked as both the agent and the verified requester (`Actor.on`, the requester's actor as verified, §recorded-actors), recorded on the change; an agent acting for no one reads only what it may itself and writes nothing.
7. **Identity is a port; the session is reactive state.** {#identity-port} Who an actor is comes from an identity port with adapters (a local picker for development and tests, passwords, OAuth/OIDC, service and agent credentials), each a runtime package (§provider-runtime-packages); the base ships the port, a session model and the local adapter. The session is a reactive model (signed out, signing in, signed in as an actor, expiring, expired) that other models read (§nested-models); sign-in, sign-out and switching actors are commands (§command-reachable). Credentials stay inside the adapter and the transport: never in model state, views, logs or an agent's context. The server derives the actor from verified credentials on every request and never trusts an actor named by the client (§requester-authority), and refuses a message that names one; a local picker served over the network is named for what it is (it is not authentication) and is refused off loopback without a signing secret; agents sign in through the same port as any actor (§near-native-agent-access).
8. **One consent model for what nobody asked for.** {#one-consent-model} An agent acts on requests (§requester-authority); a requested effect needs nothing beyond the command policy (§least-privilege). An effect nobody requested (an edit the agent infers from a document, such as ticking a task a message satisfied; a change to a person's view the agent decides on, §command-steered-views; the effect of a watch, work an agent starts on its own from a document event) runs under a **consent mode**: `apply` (done, shown, undoable in one step), `propose` (a proposal the person accepts with one action) or `off`. Modes are one setting model: a default per workspace, overridden per person, keyed by the kind of effect (view changes, edits per concept). Defaults are application policy; destructive effects are never applied without a request. An unrequested effect runs with the authority of the person whose consent setting authorizes it (the agent acts `on` that person, §recorded-actors), checked like a request, never with the agent's authority alone. An unrequested effect states its reason, quoted from the material that implied it, and the quote is checked mechanically against that material before the effect runs: a model offered a ready write will otherwise apply it to the wrong target.

