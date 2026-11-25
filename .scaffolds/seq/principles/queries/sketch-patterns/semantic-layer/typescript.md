# Semantic layer: concepts, spaces, mappings: TypeScript

The concept, invariants and conformance cases are in [README.md](README.md); this file is the TypeScript implementation.

## 1. Concepts: what things are

```ts
export interface ConceptSpec {
  id: string;                                                   // "Residue", "Shape", "Selection"
  description: string;
  isa?: readonly string[];                                      // "Residue isa Substructure"
  attrs?: Record<string, { concept?: string; dimension?: Dimension }>;   // names, optionally typed by a core concept (§core-concepts-and-embedding)
  relations?: readonly { rel: string; to: string }[];           // part-of, refers-to, derived-from, depicts, ...
  owner: string;                                                // the module that owns the domain (§representation-modules)
}

export const concepts: ConceptSpec[] = [
  { id: "Shape",   description: "Something with geometry in a scene.", owner: "repr-scene",
    attrs: { name: {}, position: { concept: "Vec3", dimension: "length" }, bounds: { concept: "Box" } } },
  { id: "Residue", description: "A residue of a polymer chain.", isa: ["Substructure"], owner: "repr-molecule",
    attrs: { name: {}, index: {}, position: { concept: "Vec3", dimension: "length" } },
    relations: [{ rel: "part-of", to: "Chain" }, { rel: "depicts", to: "Shape" }] },
];
```

Where shared concepts live and how much to declare: see [README.md](README.md#shape).

## 2. Realizations come from the language registry

Realizations: see [README.md](README.md#shape). On the form itself:

```ts
// in repr-scene's language: { tag: "mesh", realizes: "Shape", attrs: { name: "name", position: "pos" }, slots: {...}, doc: "..." }
```

## 3. Spaces and locators: how things are addressed

```ts
export type SpaceId = string;

export interface SpaceSpec {
  id: SpaceId;
  description: string;
  stability: "stable" | "positional" | "ephemeral";             // ids | positions that go stale | fleeting (screen)
  coords: Schema<readonly Node[]>;                              // shape of the coordinates (§commands's Schema)
  spatial?: { dims: 2 | 3; units: string; handedness: "right" | "left"; up?: "y" | "z"; parent?: SpaceId };
  owner: string;
}

export interface Locator {
  readonly space: SpaceId;
  readonly at: readonly Node[];                                 // coordinates in that space
  readonly concept?: string;                                    // what it denotes, when known: (is X)
  readonly anchor?: { hash?: string; rev?: number };            // snapshot it was made against (positional spaces)
}
declare function locatorToNode(l: Locator): Node;               // Locator <-> (at <space> <coords...> [(is X)] [(anchor ...)])
declare function toLocator(n: Node): Locator;                   // throws locator.bad_locator
```

Locator forms, units and the rule against bare coordinates: see [README.md](README.md#shape) and [invariants](README.md#invariants).

## 4. Mappings: runtimes connect spaces

```ts
export type Fidelity = "exact" | "approximate" | "lossy";
export interface MappingResult { results: readonly Locator[]; fidelity: Fidelity }   // results may be empty or many
export interface MappingCtx extends QueryContext {                                     // cache, signal, limits (§query-dsls)
  snapshotHash(space: SpaceId): string | undefined;                                    // current snapshot of a space's representation
}

export interface Mapping {
  readonly from: SpaceId;
  readonly to: SpaceId;
  readonly cost: number;
  readonly fidelity: Fidelity;                                  // declared upper bound; results may report worse
  translate(l: Locator, ctx: MappingCtx): MappingResult | Promise<MappingResult>;     // sync or async surface (§api-surfaces)
}
```

The runtime that knows each mapping registers it in the composition root ([§ports-and-adapters](../../../../../principles/architecture/bootstrap-patterns/ports-and-adapters/README.md)). Mappings close over their sources, with no globals (§dependencies-passed-in):

```ts
mappings.push(
  { from: "screen", to: "doc",    cost: 3, fidelity: "approximate", translate: (l, c) => renderer.pick(l, c) },   // picking is a mapping
  { from: "screen", to: "world",  cost: 2, fidelity: "approximate", translate: (l) => renderer.unproject(l) },    // a pixel is a ray
  { from: "world",  to: "screen", cost: 1, fidelity: "exact",       translate: (l) => renderer.project(l) },
  { from: "doc",    to: "world",  cost: 1, fidelity: "exact",       translate: (l) => worldTransformOf(l, doc.effective.value) },
  { from: "doc",    to: "mol",    cost: 1, fidelity: "exact",       translate: (l) => molRangeOf(l, doc.effective.value) },
  { from: "mol",    to: "doc",    cost: 1, fidelity: "exact",       translate: (l) => entitiesCovering(l, doc.effective.value) },
  { from: "mol",    to: "text",   cost: 1, fidelity: "exact",       translate: (l) => sequenceRangeOf(l) },
  { from: "data",   to: "doc",    cost: 2, fidelity: "exact",       translate: (l) => entitiesReferencing(l, doc.effective.value) },
);
// helpers assumed: renderer (WebGPU runtime), worldTransformOf, molRangeOf, entitiesCovering, sequenceRangeOf, entitiesReferencing
```

Spatial transforms and hot paths: see [README.md](README.md#shape).

## 4b. Spaces and translations as a language

The `space` language, with the molecular numbering example: see [README.md](README.md#shape). Its mapping forms compile as follows:

```ts
/** Mapping forms; each compiles to a range-aware, cacheable Mapping (section 4). */
export type MappingSpec =
  | { kind: "fields"; pairs: readonly [to: string, from: string][] }                       // rename / reorder
  | { kind: "offset"; field: string; by: number }
  | { kind: "table"; source: Node; match: readonly string[]; yield: readonly string[] }     // lookup; source is a data ref or (data-of <domain>)
  | { kind: "segments"; source: Node; match: readonly string[]; yield: readonly string[] }  // range -> ranges (alignments)
  | { kind: "compose"; steps: readonly MappingSpec[] }
  | { kind: "inverse"; of: string };                                                       // only for (bijective true)

export function compileMapping(m: { from: SpaceId; to: SpaceId; fidelity: Fidelity; spec: MappingSpec }, data: DataPort): Mapping {
  const indexes = new LruMap<string, Promise<Index>>(4);                                     // one index per source snapshot
  return {
    from: m.from, to: m.to, cost: costOf(m.spec), fidelity: m.fidelity,
    async translate(l, ctx) {
      const snap = ctx.snapshotHash(m.from) ?? "";
      let p = indexes.get(snap);
      if (!p) {                                                                             // shared: no caller's signal, failures not cached (§query-dsls)
        p = buildIndex(m.spec, data, snap); indexes.set(snap, p); p.catch(() => indexes.delete(snap));
      }
      const index = await untilAborted(p, ctx.signal);
      const out = translateRanges(index, coordKey(l), ctx);                                 // a range may become several ranges, or none
      return { results: out.map((k) => toSpaceLocator(m.to, k, l.concept)), fidelity: out.length ? m.fidelity : "lossy" };
    },
  };
}
// buildIndex: hash map for `table` (O(1) per key), sorted segments with binary search for `segments` (O(log n + k)) (§efficient-structures)
// helpers assumed: LruMap, Index, costOf, coordKey (a locator's key fields), translateRanges, toSpaceLocator; untilAborted: 13
```

The agent conversion example: see [README.md](README.md#shape); cases generated from declarations: see [README.md](README.md#conformance-cases); the plain `mol` space: see [Open](README.md#open).

## 5. The router

```ts
export const { error: locatorError } = defineErrors("locator", LOCATOR_CODES);       // codes listed in section 9 (§error-codes)
export const { error: semError } = defineErrors("sem", SEM_CODES);
const rank: Record<Fidelity, number> = { exact: 0, approximate: 1, lossy: 2 };

export class SemanticGraph {
  private chains = new Map<string, readonly Mapping[]>();
  constructor(
    readonly concepts: ReadonlyMap<string, ConceptSpec>,
    readonly spaces: ReadonlyMap<SpaceId, SpaceSpec>,
    readonly mappings: readonly Mapping[],
    readonly langs: LanguageRegistry,                            // realizations come from here (§core-concepts-and-embedding)
  ) {}

  /** Best chain: lowest worst-fidelity first, then lowest total cost. Ephemeral spaces only as endpoints. Memoized. */
  path(from: SpaceId, to: SpaceId): readonly Mapping[] {
    const k = `${from}->${to}`;
    let chain = this.chains.get(k);
    if (!chain) {
      const via = (s: SpaceId) => s === from || s === to || this.spaces.get(s)?.stability !== "ephemeral";
      chain = bestPath(this.mappings, from, to, via, (m) => [rank[m.fidelity], m.cost]);
      if (!chain) throw locatorError("locator.no_path", `no mapping ${from} -> ${to}`, { hint: `reachable from ${from}: ${reachable(this.mappings, from).join(", ")}` });
      this.chains.set(k, chain);
    }
    return chain;
  }

  async convert(l: Locator, to: SpaceId, ctx: MappingCtx): Promise<MappingResult> {
    this.checkAnchor(l, ctx);                                     // only the input carries a meaningful anchor
    const max = ctx.limits?.maxResults ?? Infinity;
    let current: readonly Locator[] = [l], fidelity: Fidelity = "exact";
    for (const m of this.path(l.space, to)) {
      if (ctx.signal?.aborted) throw locatorError("locator.cancelled", "cancelled");
      const outs = await Promise.all(current.map((c) => m.translate(c, ctx)));
      const next = new Map<string, Locator>();                   // fan-out, deduped by locator hash
      for (const o of outs) for (const r of o.results) next.set(hashNode(locatorToNode(r)), r);
      if (next.size > max) throw locatorError("locator.limit_exceeded", `${next.size} results at ${m.to}`, { hint: "narrow the locator or raise limits.maxResults" });
      current = [...next.values()];
      for (const o of outs) if (rank[o.fidelity] > rank[fidelity]) fidelity = o.fidelity;   // a chain is as good as its weakest hop
    }
    return { results: current, fidelity };
  }

  private checkAnchor(l: Locator, ctx: MappingCtx) {
    const spec = this.spaces.get(l.space);
    if (!spec) throw locatorError("locator.bad_locator", `unknown space ${l.space}`, { hint: `spaces: ${[...this.spaces.keys()].join(" ")}` });
    if (spec.stability === "positional" && l.anchor?.hash && l.anchor.hash !== ctx.snapshotHash(l.space))
      throw locatorError("locator.stale", `${l.space} locator made against another snapshot`, { hint: "re-resolve from a stable (id) locator" });
  }

  /** Forms realizing a concept or any of its sub-concepts (via isa). */
  realizationsOf(concept: string): { lang: Language; form: FormSpec }[] {
    const is = (c: string | undefined): boolean => c === concept || (!!c && (this.concepts.get(c)?.isa ?? []).some(is));
    return this.langs.all().flatMap((lang) => lang.forms.filter((f) => is(f.realizes)).map((form) => ({ lang, form })));
  }
}
// helpers assumed: bestPath (for each fidelity rank in order, Dijkstra on cost over mappings of at most that rank; the first rank
//   that connects wins. A single Dijkstra on lexicographic [max rank, summed cost] is not order-preserving and can miss the best chain),
//   reachable, LOCATOR_CODES / SEM_CODES (ErrorSpec records)
```

## 6. Navigate by meaning

Navigate by meaning: see [README.md](README.md#shape). Queries are built as Nodes, never strings:

```ts
export async function semanticSelect(
  g: SemanticGraph, roots: ReadonlyMap<string, Node>, concept: string, where: readonly Node[], ctx: MappingCtx,
): Promise<Locator[]> {
  const out = await Promise.all(g.realizationsOf(concept).map(async ({ lang, form }) => {
    const root = roots.get(lang.id);
    const engine = lang.queries?.[0]?.engine;                        // the representation's selection language (§queryable-representations)
    if (!root || !engine) return [];
    const q: Node = [sym("select"), [sym("and"), [sym("tag"), sym(form.tag)], ...where.map((w) => localize(w, form))]];
    const paths = await engine.runAsync(root, q, ctx);
    return paths.map((p) => locatorAt(lang.id, root, p, concept));   // (at <space> (id ...)) where the node has an id, else a path + anchor
  }));
  return out.flat();
}
// helpers assumed: localize ((name "ALA") -> (attr name "ALA") through form.attrs), locatorAt
```

The `select-concept` example and relation traversal: see [README.md](README.md#shape); how language ids relate to spaces: see [Open](README.md#open).

## 7. The layer is data, reachable by command

The `sem.*` commands and what `sem.concepts` returns: see [README.md](README.md#shape).

```ts
// inside registryFor(deps): graph and semCtx come from deps
registry.register({ name: "sem.convert", effect: "read", description: "Convert a locator to another space.",
  input: schema<{ locator: string; to: SpaceId }>(), output: schema<{ results: string[]; fidelity: Fidelity; normalized: string; repairs: Repair[] }>(),
  async run({ locator, to }, cmd) {
    const { node, normalized, repairs } = reprArg(locator);                          // lenient at the agent edge (§parser-repairs-slips)
    const r = await graph.convert(toLocator(node), to, { ...semCtx, signal: cmd.signal });
    return { results: r.results.map((l) => print(locatorToNode(l))), fidelity: r.fidelity, normalized, repairs };
  } });
```

## 8. Everything that points uses locators

What points with locators: see [README.md](README.md#shape). A shared selection and a view's projection of it:

```ts
export class SelectionModel {
  private readonly _state = createState<{ items: readonly Locator[] }>({ items: [] });
  readonly state: ReadonlyBehavior<{ items: readonly Locator[] }> = this._state;
  select(items: readonly Locator[]) { this._state.update({ items }); }
}

/** A view's projection of the selection; locators with no path to `space` are skipped, not fatal. */
export async function inSpace(items: readonly Locator[], g: SemanticGraph, space: SpaceId, ctx: MappingCtx): Promise<Locator[]> {
  const rs = await Promise.allSettled(items.map((l) => g.convert(l, space, ctx)));
  const out: Locator[] = [];
  for (const r of rs) {
    if (r.status === "fulfilled") out.push(...r.value.results);
    else if (!(isCoded(r.reason) && r.reason.code === "locator.no_path")) throw r.reason;   // cancelled, stale, limits still fail
  }
  return out;
}
// 3D view: pick (screen -> doc) -> SelectionModel; sequence view shows it via doc -> mol -> text; tree view via doc
```

## 9. Conformance and tests

```ts
export const mappingSuite = suite<Mapping>("Mapping",
  kase("results validate against the target space's coords schema", "quick", async (m) => { /* ... */ }),
  kase("reported fidelity is never better than declared", "quick", async (m) => { /* ... */ }),
  kase("honors the abort signal", "quick", async (m) => { /* async mappings */ }),
);

/** Generated from declarations (4b): no hand-written cases per space or mapping. */
export const declaredMappingSuite = suite<{ spec: MappingSpec; from: SpaceSpec; to: SpaceSpec; sample: readonly Locator[] }>("DeclaredMapping",
  kase("bijective mappings round-trip every sample key and range", "quick", async (d) => { /* translate, then inverse */ }),
  kase("a range maps to ranges that cover exactly the translated points", "quick", async (d) => { /* ... */ }),
  kase("keys outside the table or in alignment gaps map to nothing, never to a neighbor", "quick", async (d) => { /* ... */ }),
);

export const semanticGraphSuite = suite<SemanticGraph>("SemanticGraph",
  kase("realized concepts, isa and relations name existing concepts", "quick", async (g) => { /* no dangling ids */ }),
  kase("no_path error lists reachable spaces", "quick", async (g) => { /* ... */ }),
  kase("prefers an exact chain over a cheaper lossy one; never routes through ephemeral spaces", "quick", async (g) => { /* ... */ }),
  kase("stale positional locator is rejected; an id locator survives edits", "quick", async (g) => { /* ... */ }),
  kase("exact round trips: convert(convert(x, B), A) contains x", "medium", async (g) => { /* where both directions exist */ }),
  kase("semanticSelect(concept) equals the union of per-language queries", "medium", async (g) => { /* ... */ }),
);
```

Error codes and the remaining cases: see [README.md](README.md#conformance-cases). Codes are declared in this module with `defineErrors` ([§error-codes](../../../../../principles/architecture/bootstrap-patterns/error-codes/README.md)).
