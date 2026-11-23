# Ports and adapters, composition root: TypeScript

The concept, invariants and conformance cases are in [README.md](README.md); this file is the TypeScript implementation.

## Layout

```
lib/
  storage/                   # port + in-memory adapter + conformance suite (§conformance)
  storage-cloud-runtime/    # the provider adapter ("cloud" stands for whichever provider the project uses): the only lib package that depends on its SDK (§provider-runtime-packages)
  agents/                    # uses ports only
packages/
  product-x/compose/         # @scope/product-x-compose: Deps, wire, local(), test(), createApp: no provider dependencies
apps/
  product-x/web/
    compose.ts               # a composition root: production() binds *-runtime packages, reads env
services/
  product-x/api/compose.ts   # another root, same shape
```

## Port (in common)

```ts
// lib/storage/port.ts
export interface DocumentStore<T extends { id: string }> {
  get(id: string): Promise<T | undefined>;
  put(doc: T): Promise<void>;
  delete(id: string): Promise<void>;
  /** Streaming surface (§api-surfaces): the current value first, then every change. */
  watch(id: string): Observable<T | undefined>;
}
```

Port names describe the capability (`DocumentStore`), never the provider (`CloudClient`).

## In-memory adapter (ships with the port)

```ts
// lib/storage/memory.ts
export function memoryStore<T extends { id: string }>(): DocumentStore<T> {
  const data = new Map<string, T>();
  const changed = createSubject<string>();                     // ids that changed
  return {
    async get(id) { return data.get(id); },
    async put(doc) { data.set(doc.id, doc); changed.next(doc.id); },
    async delete(id) { data.delete(id); changed.next(id); },
    watch: (id) => observable((emit) => {
      emit(data.get(id));                                      // read at SUBSCRIBE time, not when watch() was called
      return changed.subscribe((k) => { if (k === id) emit(data.get(id)); });
    }),
  };
}
```

## Provider adapter (in its own runtime package)

```ts
// lib/storage-cloud-runtime/cloud-store.ts
import type { CloudClient } from "cloud-sdk";      // the provider's SDK, imported only here
import { api } from "./functions";               // its server functions ship with it; services/<product>/cloud-runtime deploys them

export function cloudStore<T extends { id: string }>(client: CloudClient, table: string): DocumentStore<T> {
  return {
    get: (id) => client.query(api.docs.get, { table, id }),
    put: (doc) => client.mutation(api.docs.put, { table, doc }),
    delete: (id) => client.mutation(api.docs.remove, { table, id }),
    // the provider's reactive stream is wrapped here and never escapes
    watch: (id) => observable((emit, fail) => client.onUpdate(api.docs.get, { table, id }, emit, fail)),
  };
}
```

## Composition root and lightweight DI

Dependencies are parameters. One composition root wires them; there is no container.

```ts
// Business code asks for ports; it never constructs collaborators or looks them up.
export class NoteService {
  constructor(private readonly deps: { store: DocumentStore<Note>; clock: Clock }) {}
  add = async (text: string) => { await this.deps.store.put({ id: newId("n"), text, at: this.deps.clock.now() }); };
}
// BAD: `new CloudClient(...)` inside, `container.get("store")` (a service locator), or a module-level singleton.
```

A tiny helper gives lazy, memoized, disposable instances and nested lifetimes (no reflection):

```ts
// lib/core/di/scope.ts
export class Scope {
  private disposers: Array<() => void | Promise<void>> = [];
  private children = new Set<Scope>();
  private parent?: Scope;

  /** Lazy singleton within this scope; `dispose` is registered when the instance is first created. */
  lazy<T>(make: () => T, dispose?: (t: T) => void | Promise<void>): () => T {
    let made = false, value!: T;
    return () => {
      if (!made) { value = make(); made = true; if (dispose) this.disposers.push(() => dispose(value)); }
      return value;
    };
  }

  /** Per-session/per-request scope. Children are disposed before their parent's instances, which they may depend on. */
  child(): Scope { const c = new Scope(); c.parent = this; this.children.add(c); return c; }

  async dispose(): Promise<void> {
    this.parent?.children.delete(this);                        // a disposed request scope is not retained by its parent
    const errors: unknown[] = [];
    for (const c of [...this.children]) { try { await c.dispose(); } catch (e) { errors.push(e); } }
    this.children.clear();
    for (const d of this.disposers.splice(0).reverse()) { try { await d(); } catch (e) { errors.push(e); } }   // reverse creation order
    if (errors.length) throw new AggregateError(errors, "scope dispose failed");
  }
}
```

```ts
// packages/product-x/compose/compose.ts: the shape and the provider-free compositions
export interface Deps {
  store: () => DocumentStore<Note>;          // thunks keep creation lazy and cycle-free
  docs: () => DocumentStore<DocRecord>;      // used by the doc.* commands (§commands)
  sync: () => SyncPort;                      // collaboration transport (§committed-transient-identity): server log in production, loopback locally
  clock: () => Clock;
  llm: () => LlmClient;
  notes: () => NoteService;
  doc: () => Doc;                            // the committed/transient document (§committed-transient-identity); its `effective` is the app's TreeState
}

export type Bindings = Pick<Deps, "store" | "docs" | "sync" | "clock" | "llm">;

export function wire(scope: Scope, b: Bindings): Deps {
  return {
    ...b,
    notes: scope.lazy(() => new NoteService({ store: b.store(), clock: b.clock() })),
    doc: scope.lazy(() => createDoc(emptyDocument(), createHistory(), b.sync())),
  };
}

export function local(scope = new Scope()): Deps {             // in-memory: tests, stories, evals; no cloud account
  return wire(scope, {
    store: scope.lazy(() => memoryStore<Note>()),
    docs:  scope.lazy(() => memoryStore<DocRecord>()),
    sync:  scope.lazy(() => loopbackSync()),
    clock: scope.lazy(() => fakeClock()),
    llm:   scope.lazy(() => scriptedLlm([])),
  });
}

/** Tests override single bindings by value, not by patching modules. Overrides are memoized like everything else. */
export function test(over: Partial<Bindings> = {}, scope = new Scope()): Deps {
  const base = local(scope);
  const pick = <K extends keyof Bindings>(k: K): Bindings[K] => (over[k] ? (scope.lazy(over[k]!) as Bindings[K]) : base[k]);
  return wire(scope, { store: pick("store"), docs: pick("docs"), sync: pick("sync"), clock: pick("clock"), llm: pick("llm") });
}
```

Each deployable has its own root, the only code that imports provider runtime packages and reads the environment (§one-composition-root, §provider-runtime-packages):

```ts
// apps/product-x/web/compose.ts: the ONLY place concrete implementations meet, for this app
import { cloudStore, serverSync } from "@scope/storage-cloud-runtime";
import { hostedLlm } from "@scope/llm-hosted-runtime";   // "hosted" stands for whichever model provider the project uses
import { wire, type Deps } from "@scope/product-x-compose";

export function production(scope = new Scope()): Deps {
  const cloud = scope.lazy(() => cloudClient(env.CLOUD_URL), (c) => c.close());   // the client owns the connection
  return wire(scope, {
    store: scope.lazy(() => cloudStore<Note>(cloud(), "notes")),
    docs:  scope.lazy(() => cloudStore<DocRecord>(cloud(), "docs")),
    sync:  scope.lazy(() => serverSync(cloud())),                     // the client owns the connection
    clock: scope.lazy(() => systemClock),
    llm:   scope.lazy(() => hostedLlm(env.LLM_API_KEY)),
  });
}
```

When tests, smokes and agents need the whole app as one object, define the facade once, next to the composition (used by [§test-tiers](../../../testing/bootstrap-patterns/test-tiers/README.md) and [§commands](../commands/README.md)):

```ts
export interface App { invoke: Invoker; state: TreeState; scope: Scope }       // Invoker: §commands; TreeState: §representations
export const createApp = (deps: Deps, scope: Scope): App =>
  ({ invoke: bindInvoker(registryFor(deps)), state: deps.doc().effective, scope });
```

```ts
// services/product-x/api/main.ts: build, run, dispose (its compose.ts has its own production())
const scope = new Scope();
const deps = production(scope);
const server = startServer(deps);
process.on("SIGTERM", async () => { await server.stop(); await scope.dispose(); });
```

Per-session or per-request lifetimes use `scope.child()`: the child creates its own instances (an authenticated actor, a request cache) on top of app-wide ones and is disposed when the request ends.

If useful, the composition can enumerate its bindings for conformance runs ([§conformance](../../../testing/bootstrap-patterns/conformance/README.md)) or docs; derive that list from `wire`, never maintain it by hand (§enumerable-bindings).

## Checks an agent can run

Invariants: see [README.md](README.md#invariants). In a TypeScript workspace:

- `grep -rl '"cloud-sdk' --include=package.json lib packages` lists only `*-cloud-runtime` packages (§provider-runtime-packages); an import-boundary lint enforces apps/services → packages → lib (§common-vs-product).
- No code outside a composition root (`apps/*/*/compose.ts`, `services/*/*/compose.ts`) calls `new` on an adapter, reads `process.env`, or imports a provider runtime package.

Helpers assumed: `newId`, `cloudClient`, `systemClock`, `fakeClock`, `hostedLlm`, `scriptedLlm`, `serverSync`, `loopbackSync`, `emptyDocument`, `createHistory` (§committed-transient-identity), `createDoc` (§committed-transient-identity), `registryFor`/`bindInvoker`/`Invoker` (§commands), `startServer`, `env`.
