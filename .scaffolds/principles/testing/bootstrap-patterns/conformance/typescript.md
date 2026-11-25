# Conformance library: TypeScript

The concept, invariants and conformance cases are in [README.md](README.md); this file is the TypeScript implementation.

## The whole library

```ts
export type Cost = "quick" | "medium" | "guard" | "full";       // the test tiers (§test-tiers); guard = scaling/allocation cases
const order: Record<Cost, number> = { quick: 0, medium: 1, guard: 2, full: 3 };

export interface Case<S> {
  title: string;
  cost: Cost;
  needs?: string;                                                // optional capability this case requires (§optional-capabilities)
  run(subject: S): void | Promise<void>;
}
export interface Suite<S> { id: string; cases: readonly Case<S>[] }

// build
export const kase = <S>(title: string, cost: Cost, run: Case<S>["run"], needs?: string): Case<S> => ({ title, cost, run, needs });
export const suite = <S>(id: string, ...cases: Case<S>[]): Suite<S> => ({ id, cases });

// compose
export const combine = <S>(id: string, ...ss: Suite<S>[]): Suite<S> => ({ id, cases: ss.flatMap((s) => s.cases) });
export const adapt = <S, T>(s: Suite<S>, id: string, pick: (t: T) => S): Suite<T> =>      // reuse a suite on a part of a bigger subject
  ({ id, cases: s.cases.map((c) => ({ ...c, run: (t: T) => c.run(pick(t)) })) });
export const upTo = <S>(s: Suite<S>, max: Cost): Suite<S> =>                               // filter by cost
  ({ id: s.id, cases: s.cases.filter((c) => order[c.cost] <= order[max]) });

// plan: which cases run and which are skipped, with the reason; pure, so any runner (or a CLI) can use it
export const plan = <S>(s: Suite<S>, opts: { upTo?: Cost; caps?: readonly string[] } = {}) =>
  s.cases.map((c) => ({ c, skip: order[c.cost] > order[opts.upTo ?? "full"] ? "cost"
                                : c.needs !== undefined && !opts.caps?.includes(c.needs) ? `needs ${c.needs}` : undefined }));
```

Packaging: see [README.md](README.md#shape). Here the runner binding is the subpath `@scope/core-conformance/vitest`,
and a package's suite that uses the runner's assertions is its own subpath (`@scope/storage/conformance`).

```ts
// @scope/core-conformance/vitest: the only part that touches the test runner
export function register<S>(
  s: Suite<S>,
  make: () => S | Promise<S>,
  opts: { label?: string; upTo?: Cost; caps?: readonly string[] } = {},
) {
  const max = order[opts.upTo ?? "full"];
  describe(`${s.id}${opts.label ? ` (${opts.label})` : ""}`, () => {
    for (const c of s.cases) {
      const skip = order[c.cost] > max || (c.needs !== undefined && !opts.caps?.includes(c.needs));   // skipped visibly, never passed
      (skip ? it.skip : it)(c.needs ? `${c.title} [${c.needs}]` : c.title, async () => c.run(await make()));  // fresh subject per case
    }
  });
}
```

That is all of it; everything else is ordinary code using these pieces.

## A library suite (ships beside its port)

```ts
// lib/storage/src/conformance.ts, exported as "./conformance" (it uses the runner's expect)
export const documentStoreSuite = suite<DocumentStore<DocRecord>>("DocumentStore",
  kase("round-trips a document", "quick", async (s) => {
    await s.put({ id: "a", title: "x" });
    expect(await s.get("a")).toEqual({ id: "a", title: "x" });
  }),
  kase("get of missing id is undefined", "quick", async (s) => expect(await s.get("nope")).toBeUndefined()),
  kase("watch emits current value first, then changes", "quick", async (s) => { /* ... */ }),
  kase("10k puts scale ~linearly", "guard", async (s) => { /* ratio test, see §test-tiers */ }),
  kase("writes in a transaction are atomic", "medium", async (s) => { /* ... */ }, "transactions"),
  kase("watch delivers in write order under load", "full", async (s) => { /* ... */ }, "ordered-watch"),
);
```

## Every implementation runs it, at the cost that fits

```ts
// lib/storage/memory.quick.test.ts
register(documentStoreSuite, () => memoryStore<DocRecord>(), { label: "memory", upTo: "quick", caps: ["ordered-watch"] });
// guard cases run from memory.guard.test.ts with { upTo: "guard" }

// lib/storage-cloud-runtime/cloud.full.test.ts (real deployment, nightly)
register(documentStoreSuite, () => cloudStore<DocRecord>(testClient(), uniqueTable()), { label: "cloud", caps: ["transactions"] });
```

Where suites run, and at what cost: see [README.md](README.md#shape).

## Composing for larger subjects

```ts
// an app's "notes" subsystem exposes a store; reuse the store suite on it, plus its own cases
const notesSuite = combine<NotesSubsystem>("Notes",
  adapt(documentStoreSuite, "Notes.store", (n) => n.store),
  suite("Notes.behavior", kase("add then list", "quick", async (n) => { /* ... */ })),
);
register(notesSuite, () => createNotes(local()), { label: "local", upTo: "medium", caps: ["ordered-watch"] });
```

## Application conformance: a test file, not a framework

See [README.md](README.md#shape) (§applications-reuse-suites).

```ts
// packages/product-x/compose/conformance.medium.test.ts   (local composition, gate)
// a fresh composition per case: local()'s bindings are memoized, so a shared one would leak state between cases
register(documentStoreSuite, () => local().docs(), { label: "composed store", upTo: "medium", caps: ["ordered-watch"] });
register(commandCatalogSuite, () => registryFor(local()), { upTo: "quick" });   // every command described, valid effect (§commands)
register(vocabularySuite, () => vocabularySamples(), { upTo: "quick" });      // each node kind's sample validates and round-trips (§representations)

// services/product-x/api/conformance.full.test.ts        (deployed composition, nightly + after deploy)
register(documentStoreSuite, () => production().docs(), { label: "deployed store" });
```
