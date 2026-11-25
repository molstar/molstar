# Query DSLs over representations: TypeScript

The concept, invariants and conformance cases are in [README.md](README.md); this file is the TypeScript implementation.

## Language core (`repr-core/query`)

```ts
/** Everything a run needs besides the data and the query. Created in the composition root (§composition-root) and passed down. */
export interface QueryContext {
  readonly cache: QueryCache;                                  // compiled queries, per-snapshot indexes
  readonly signal?: AbortSignal;                               // cooperative cancellation
  readonly limits?: { maxResults?: number; deadlineMs?: number; maxDepth?: number };   // maxDepth: embedded evaluation (§core-concepts-and-embedding)
  readonly trace?: (e: { kind: string; ms?: number; detail?: unknown }) => void;
  readonly derived?: Materializer;                             // derived forms (meso's `derived-representations` (not adopted here))
}
export type QueryEval = QueryContext & { readonly root: Node };

export type Pred = (n: Node, path: Path, ctx: QueryEval) => boolean;
/** A function compiles to a synchronous predicate, or to a *staged* one that needs async preparation first. */
export type Staged = { prepare(ctx: QueryEval): Promise<Pred> };
export type Compiled = Pred | Staged;
export const isStaged = (c: Compiled): c is Staged => typeof c !== "function";
export type Compile = (q: Node) => Compiled;
/** `langs` lets a function compile an argument in another language, late-bound by id (§core-concepts-and-embedding). */
export type FnCompiler = (args: readonly Node[], compile: Compile, langs: LanguageRegistry) => Compiled;

export class QueryLang {
  private fns = new Map<string, FnCompiler>();
  constructor(readonly id: string, private readonly langs: LanguageRegistry = new LanguageRegistry()) {}

  register(name: string, c: FnCompiler): this {
    if (this.fns.has(name)) throw new Error(`duplicate query fn ${name}`);      // programmer error (§programmer-errors)
    this.fns.set(name, c); return this;
  }

  compile: Compile = (q) => {
    const name = tag(q);
    if (name === undefined) throw queryError("query.not_a_query", `not a query: ${print(q)}`, { hint: "write (fn args...)" });
    const c = this.fns.get(name);
    if (!c) throw queryError("query.unknown_function", `unknown function ${name}`, { hint: `known: ${[...this.fns.keys()].sort().join(" ")}` });
    return c((q as readonly Node[]).slice(1), this.compile, this.langs);
  };

  /** (select <pred>) -> compiled body, cached per language + query text. */
  compileSelect(q: Node, ctx: QueryContext): Compiled {
    const parts = q as readonly Node[];
    if (tag(q) !== "select" || parts.length !== 2) throw queryError("query.expected_select", undefined, { hint: "(select <pred>)" });
    return ctx.cache.compiled(this, q, () => this.compile(parts[1]!));
  }

  /** The one walk: signal checked before and periodically, deadline and maxResults enforced. */
  evaluate(pred: Pred, root: Node, ctx: QueryContext): Path[] {
    const ev: QueryEval = { ...ctx, root }, out: Path[] = [];
    const max = ctx.limits?.maxResults ?? Infinity;
    const deadline = ctx.limits?.deadlineMs !== undefined ? performance.now() + ctx.limits.deadlineMs : Infinity;
    const check = () => {
      if (ctx.signal?.aborted) throw queryError("query.cancelled");
      if (performance.now() > deadline) throw queryError("query.deadline_exceeded");
    };
    check();
    let visited = 0;
    walk(root, [], (n, p) => {
      if ((++visited & 1023) === 0) check();                     // cheap periodic checks
      if (!pred(n, p, ev)) return;
      if (out.length === max) throw queryError("query.limit_exceeded", `more than ${max} results`, { hint: "narrow the query or raise limits.maxResults" });
      out.push(p);
    });
    return out;
  }

  /** Synchronous run; staged queries need runAsync. */
  select(q: Node, root: Node, ctx: QueryContext): Path[] {
    const c = this.compileSelect(q, ctx);
    if (isStaged(c)) throw queryError("query.needs_async", undefined, { hint: "use runAsync" });
    return this.evaluate(c, root, ctx);
  }
}

export const queryErrors = defineErrors("query", {                                     // 14
  "query.not_a_query":       { description: "Expected a function form (fn args...)." },
  "query.unknown_function":  { description: "No such function in this language." },
  "query.bad_argument":      { description: "Argument of the wrong kind." },
  "query.expected_select":   { description: "Top-level form must be (select <pred>)." },
  "query.needs_async":       { description: "Query has staged parts; use runAsync." },
  "query.cancelled":         { description: "Query was cancelled." },
  "query.deadline_exceeded": { description: "Query exceeded its deadline.", retryable: true },
  "query.limit_exceeded":    { description: "Query matched more results than allowed." },
});
export const queryError = queryErrors.error;
```

`walk` allocates a path per node; fine for the scan engine on moderate trees. Large representations use engines with indexes (below), guided by §performance.

Core functions every language has (shared, so domains stay small):

```ts
/** Tag and attribute names are symbols; strings are accepted for convenience. */
export function nameOf(a: Node | undefined): string {
  if (a !== undefined && isSym(a)) return a.sym;
  if (typeof a === "string") return a;
  throw queryError("query.bad_argument", `expected a name, got ${a === undefined ? "nothing" : print(a)}`);
}
export const equalNode = (a: Node | undefined, b: Node | undefined): boolean =>
  Object.is(a, b) || (a !== undefined && b !== undefined && hashNode(a) === hashNode(b));   // hashNode: 15; 0 and -0 differ, as their hashes do

/** Combine compiled children; if any is staged, the result is staged and prepares them in parallel. */
export const lift = (xs: Compiled[], combine: (ps: Pred[]) => Pred): Compiled =>
  xs.every((x) => !isStaged(x))
    ? combine(xs as Pred[])
    : { prepare: async (ctx) => combine(await Promise.all(xs.map((x) => (isStaged(x) ? x.prepare(ctx) : x)))) };

export const withCore = (L: QueryLang) => L
  .register("and", (a, c) => lift(a.map(c), (ps) => (n, p, x) => ps.every((f) => f(n, p, x))))
  .register("or",  (a, c) => lift(a.map(c), (ps) => (n, p, x) => ps.some((f) => f(n, p, x))))
  .register("not", ([a], c) => lift([c(a!)], ([f]) => (n, p, x) => !f!(n, p, x)))
  .register("tag", ([t]) => { const name = nameOf(t); return (n) => tag(n) === name; })
  .register("attr", ([k, v]) => { const key = nameOf(k); return (n) => equalNode(attrOf(n, key), v); })
  .register("has-child", ([a], c) => lift([c(a!)], ([f]) => (n, p, x) =>
    Array.isArray(n) && (n as readonly Node[]).some((ch, i) => f!(ch, [...p, i], x))));
```

## A domain adds functions, not an engine

```ts
// lib/repr/molecule/src/query.ts
// spatial predicates (within, inside, intersects) are core functions (§core-concepts-and-embedding); the domain says how to read a position
export const molLang = (langs: LanguageRegistry) =>
  withSpatial(withCore(new QueryLang("mol.select", langs)), { position: atomPosition })
    .register("residue-name", ([name]) => (n) => tag(n) === "residue" && equalNode(attrOf(n, "name"), name));
```

The language is built in the composition root with the shared `LanguageRegistry`, so functions that embed another language resolve it late (§core-concepts-and-embedding). Domain functions document their semantics in `_spec/` and are covered by the language's conformance suite.

## Composable by default: generators, modifiers, filters, combinators (§composable-queries)

Selections, operator roles and an example: see [README.md](README.md#shape). The predicates above are the tests used inside generators and filters.

```ts
/** A selection: groups of nodes (a residue = the paths of its atoms). */
export type Group = readonly Path[];
export type Selection = readonly Group[];
export type SelQuery = (ev: QueryEval) => Selection;                   // staged variants work as for predicates (lift)
export type CompileSel = (q: Node) => SelQuery;

export type Role = "generator" | "modifier" | "filter" | "combinator";
export interface SelOp { role: Role; compile(args: readonly Node[], c: { pred: Compile; sel: CompileSel }): SelQuery }

/** Core operators every selection language gets; domains add their own with the same signature. */
export const selectionOps = (lang: QueryLang): Record<string, SelOp> => ({
  // generator: one walk, a predicate test, groups by an ancestor level (none: one group per node)
  atoms:  { role: "generator", compile: ([test, by], c) => { const p = asPred(c.pred(test!)); const level = by && nameOf(by);
            return (ev) => groupBy(lang.evaluate(p, ev.root, ev), level, ev.root); } },
  // modifier: applied per group, results flattened (the monadic bind: sel.flatMap(f))
  expand: { role: "modifier", compile: ([s, level], c) => { const inner = c.sel(s!); const lv = nameOf(level);
            return (ev) => inner(ev).flatMap((g) => [wholeAncestor(g, lv, ev.root)]); } },
  // filters: keep groups that pass a group-level test
  pick:   { role: "filter", compile: ([s, test], c) => { const inner = c.sel(s!); const t = groupTest(test!, c);
            return (ev) => inner(ev).filter((g) => t(g, ev)); } },
  within: { role: "filter", compile: ([s, r, target], c) => { const inner = c.sel(s!), other = c.sel(target!);
            return (ev) => { const near = nearIndex(other(ev), ev); return inner(ev).filter((g) => near(g, r as number)); }; } },
  // combinators: several selections -> one
  union:     { role: "combinator", compile: (ss, c) => { const qs = ss.map(c.sel); return (ev) => dedupeGroups(qs.flatMap((q) => q(ev))); } },
  intersect: { role: "combinator", compile: (ss, c) => { const qs = ss.map(c.sel); return (ev) => intersectGroups(qs.map((q) => q(ev))); } },
});
// helpers assumed: asPred (rejects staged predicates in sync paths), groupBy, wholeAncestor, groupTest, nearIndex (spatial, §core-concepts-and-embedding), dedupeGroups, intersectGroups
```

Sugar, laws and opting out: see [README.md](README.md#shape) and [conformance cases](README.md#conformance-cases); the `within` naming question is in [Open](README.md#open).

## Several languages and dialects per representation (§queryable-representations, §dialects-lower)

Languages vs dialects, with the molecular example: see [README.md](README.md#shape).

```ts
/** Another syntax for an existing query language: parsed, lowered, never executed directly. */
export interface Dialect {
  readonly id: string;                         // "mol.short"
  readonly lowersTo: string;                   // canonical language id: "mol.select"
  parse(text: string): Node;                   // its own syntax -> its own tree; coded errors with text positions (§error-codes)
  readonly lower: Pass;                        // 19: dialect tree -> canonical query, with a SourceMap back to the text
  print?(canonical: Node): string;             // optional: show a canonical query in this notation
}

/** One entry point for every caller (UI search box, agents, commands): text + optional language/dialect id -> canonical query. */
export async function resolveQuery(langs: LanguageRegistry, repr: string, text: string, id?: string, ctx?: DeriveCtx) {
  const dialect = id ? langs.dialect(repr, id) : undefined;
  if (dialect) {
    const { out, sourceMap, diagnostics } = await dialect.lower.run(dialect.parse(text), [], ctx!);
    return { lang: langs.query(repr, dialect.lowersTo), canonical: out, sourceMap, diagnostics, repairs: [] as Repair[] };
  }
  const { node, repairs } = reprArg(text);                                    // canonical text: lenient at the edge (§parser-repairs-slips)
  return { lang: langs.query(repr, id), canonical: node, repairs };           // id omitted: the representation's selection language
}
```

Rules for dialects: see [README.md](README.md#invariants).

## Query context: caching, cancellation, limits

What the context carries and who creates it: see [README.md](README.md#shape).

```ts
export interface QueryCache extends RootMemo {                                        // RootMemo: forRoot, owned by 15
  /** Compiled predicate per language + query text. Bounded (LRU) in practice. */
  compiled(lang: QueryLang, q: Node, make: () => Compiled): Compiled;
}

export function createQueryCache(): QueryCache {
  const compiled = new Map<string, Compiled>();
  const perRoot = new WeakMap<object, Map<string | symbol, unknown>>();
  return {
    compiled(lang, q, make) {
      const k = `${lang.id}:${print(q)}`;
      let c = compiled.get(k); if (!c) { c = make(); compiled.set(k, c); } return c;
    },
    forRoot<T>(root: Node, key: string | symbol, build: () => T): T {
      if (typeof root !== "object" || root === null) return build();          // atoms are not cacheable by identity
      let m = perRoot.get(root); if (!m) perRoot.set(root, (m = new Map()));
      if (m.has(key)) return m.get(key) as T;
      const v = build(), memo = m;
      memo.set(key, v);
      if (v instanceof Promise) v.catch(() => memo.delete(key));              // only successes stay cached
      return v;
    },
  };
}

export const createQueryContext = (over: Partial<QueryContext> = {}): QueryContext => ({ cache: createQueryCache(), ...over });

// policy varies by caller without touching the engine (all in the composition root)
const cache    = createQueryCache();
const uiCtx    = createQueryContext({ cache });                                   // trusted, unbounded
const agentCtx = createQueryContext({ cache, limits: { maxResults: 1000, deadlineMs: 200 } });
```

Why per-root caching is safe: see [README.md](README.md#shape).

## Engines: same semantics, different execution

The three surfaces: see [README.md](README.md#shape).

```ts
export interface QueryEngine {
  /** `ctx` changes cost and bounds, never results (§query-context). It is required: there is no default context. */
  run(root: Node, q: Node, ctx: QueryContext): Path[];
  runAsync(root: Node, q: Node, ctx: QueryContext): Promise<Path[]>;
  stream(root: Node, q: Node, ctx: QueryContext): Observable<readonly Path[]>;    // batches; completes when done
}

export function scanEngine(lang: QueryLang): QueryEngine {
  const runAsync = async (root: Node, q: Node, ctx: QueryContext): Promise<Path[]> => {
    const c = lang.compileSelect(q, ctx);
    let pred: Pred;
    try { pred = isStaged(c) ? await c.prepare({ ...ctx, root }) : c; }       // may await remote work
    catch (e) { throw ctx.signal?.aborted ? queryError("query.cancelled") : e; }
    return lang.evaluate(pred, root, ctx);                                     // re-checks the signal before walking
  };
  return {
    run: (root, q, ctx) => lang.select(q, root, ctx),
    runAsync,
    stream: (root, q, ctx) => observable((emit, fail, done) => {
      const ac = new AbortController();
      const signal = ctx.signal ? AbortSignal.any([ctx.signal, ac.signal]) : ac.signal;
      runAsync(root, q, { ...ctx, signal }).then((ps) => { emit(ps); done(); }, (e) => { if (!isCancelled(e)) fail(e); });
      return () => ac.abort();                                                 // unsubscribe cancels
    }),
  };
}
```

Indexed and remote engines: see [README.md](README.md#shape). An indexed engine keeps its tag index in `ctx.cache.forRoot(root, TAG_INDEX, ...)`. Long main-thread runs use `runChunked` ([§reactive-models](../../../../../principles/ui/bootstrap-patterns/reactive-models/typescript.md#keeping-the-main-thread-responsive)) inside `runAsync`/`stream`.

## Async in queries

Staging: see [README.md](README.md#shape); invariants: see [README.md](README.md#invariants).

```ts
// the service is a port passed in by the composition root, never a global
interface Embeddings { nearest(text: string, k: number, opts?: { signal?: AbortSignal }): Promise<readonly string[]> }

export const withSimilarity = (L: QueryLang, embeddings: Embeddings) => L.register("similar-to", ([text, k]) => {
  if (typeof text !== "string" || typeof k !== "number") throw queryError("query.bad_argument", "(similar-to \"text\" k)");
  return {
    async prepare(ctx) {
      // shared per snapshot and per ALL arguments; the fetch is not tied to any one caller's signal,
      // so one caller aborting does not cancel the others (a data port would also abort it when no waiter is left, meso's `data-references` (not adopted here))
      const hits = ctx.cache.forRoot(ctx.root, `similar-to:${print(text)}:${k}`, () => embeddings.nearest(text, k));
      const ids = new Set(await untilAborted(hits, ctx.signal));            // this caller stops waiting on abort
      return (n) => ids.has(attrOf(n, "id") as string);                     // evaluation stays sync and pure
    },
  };
});
declare function untilAborted<T>(p: Promise<T>, signal?: AbortSignal): Promise<T>;   // rejects with query.cancelled on abort
```

Properties: see [README.md](README.md#invariants).

### Queries are cancelable

How each surface cancels and the rules: see [README.md](README.md#shape) and [README.md](README.md#invariants). A cancelled query rejects with `query.cancelled` (`isCancelled(e)`, §error-codes).

## Live queries

Invariants: see [README.md](README.md#invariants).

```ts
export function observeQuery(state: TreeState, engine: QueryEngine, q: Node, ctx: QueryContext): Observable<Path[]> {
  return observable((emit, fail) => {
    let ac: AbortController | undefined, seq = 0;
    const run = (root: Node) => {
      ac?.abort(); ac = new AbortController();
      const mine = ++seq;
      const signal = ctx.signal ? AbortSignal.any([ctx.signal, ac.signal]) : ac.signal;
      engine.runAsync(root, q, { ...ctx, signal }).then(
        (r) => { if (mine === seq) emit(r); },                                // drop superseded results
        (e) => { if (mine === seq && !isCancelled(e)) fail(e); });
    };
    const off = state.subscribe(run);                                          // behavior: current root, then each change
    return () => { seq++; off(); ac?.abort(); };                               // nothing is emitted after teardown
  });
}
```

## Reachable by command (§command-reachable, §queries-by-command)

```ts
registry.register({
  name: "state.query", effect: "read",
  description: "Run a representation query; returns matching nodes and their paths.",
  input: schema<{ query: string; lang?: string; limit?: number }>(),        // lang: a query language or dialect id
  output: schema<{ matches: { path: Path; node: string }[]; normalized: string; repairs: Repair[] }>(),
  async run({ query, lang, limit = 100 }, cmd) {
    const { lang: ql, canonical: q, repairs } = await resolveQuery(langs, "doc", query, lang, deriveCtx());   // dialects lowered
    const normalized = print(q);
    const root = state.value;                                                  // ONE snapshot for the whole run
    const max = Math.min(limit, agentCtx.limits?.maxResults ?? Infinity);      // the caller may lower the policy cap, never raise it
    const ctx = { ...agentCtx, signal: cmd.signal, limits: { ...agentCtx.limits, maxResults: max } };
    const paths = await langs.engine("doc", ql.id).runAsync(root, q, ctx);     // the engine of the resolved language (§core-concepts-and-embedding)
    return { matches: paths.map((p) => ({ path: p, node: print(nodeAt(root, p)!) })), normalized, repairs };
  },
});
// helpers assumed: schema (§commands), state: TreeState, langs, agentCtx, deriveCtx (a DeriveCtx, meso's `derived-representations` (not adopted here)) from the composition root
```

Usage: see [README.md](README.md#shape).

## Conformance and tests

```ts
// lib/repr/core/query/conformance.ts: every engine passes the same semantic cases, on every surface (§conformance)
const tree = parse('(scene (mesh "a") (cam))');
export const queryEngineSuite = suite<QueryEngine>("QueryEngine",
  kase("tag selects by head symbol", "quick", (e) =>
    expect(e.run(tree, parse("(select (tag mesh))"), createQueryContext())).toEqual([[1]])),
  kase("run, runAsync and stream agree", "quick", async (e) => { /* same cases through all three */ }),
  kase("unknown function: query.unknown_function, known names in the hint", "quick", (e) => { /* ... */ }),
  kase("too many results: query.limit_exceeded", "quick", (e) => { /* limits: { maxResults: 0 } */ }),
  kase("abort: rejects with query.cancelled, port calls aborted, cache unchanged", "quick", async (e) => { /* ... */ }),
  kase("staged query equals its pre-resolved equivalent", "medium", async (e) => { /* ... */ }),
  kase("results equal the scan engine on generated trees and queries", "medium", async (e) => { /* property test */ }),
);

// composable languages (§composable-queries): the laws, on generated selections
export const selectionLawsSuite = suite<Record<string, SelOp>>("SelectionLaws",
  kase("chained modifiers are associative", "quick", (ops) => { /* expand(expand(s, a), b) vs grouped */ }),
  kase("(pick s true) equals s; union is commutative and idempotent", "quick", (ops) => { /* ... */ }),
  kase("a generator equals the union of its groups", "medium", (ops) => { /* property test */ }),
);

// dialects (§dialects-lower): written form -> expected canonical form, errors point into the written text
export const dialectSuite = (cases: ReadonlyArray<[written: string, canonical: string]>) => suite<Dialect>("Dialect",
  ...cases.map(([w, c]) => kase(`lowers ${w}`, "quick", async (d) => expect(print((await d.lower.run(d.parse(w), [], deriveCtx())).out)).toBe(c))),
);
register(queryEngineSuite, () => scanEngine(withCore(new QueryLang("test"))), { label: "scan", upTo: "medium" });

// queries are data
test.prop([arbQuery])("parse(print(q)) == q", (q) => expect(print(parse(print(q)))).toBe(print(q)));
```

Conformance cases: see [README.md](README.md#conformance-cases).
