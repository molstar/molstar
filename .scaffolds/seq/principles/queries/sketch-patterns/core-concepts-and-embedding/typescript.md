# Core concepts and embedded languages: TypeScript

The concept, invariants and conformance cases are in [README.md](README.md); this file is the TypeScript implementation.

## Core forms and types

The core forms: see [README.md](README.md#shape).

```ts
export type Dimension = "length" | "angle" | "time" | "mass" | "ratio" | (string & {});
export interface Quantity { value: number; unit: string }
export type Vec<N extends 1 | 2 | 3 | 4> = readonly number[] & { readonly length: N };
export interface Box<N extends 1 | 2 | 3> { min: Vec<N>; max: Vec<N> }
export interface Interval<T> { lo: T; hi: T }                                          // closed
export interface Signal<T> { keys: Table<{ t: number; v: T }>; interp: "step" | "linear" | "cubic" }   // t in seconds; Table: 25
export interface Tensor { shape: readonly number[]; dtype: "f32" | "f64" | "i32" | "u8"; data: DataRef }

declare function convertQuantity(q: Quantity, to: string): Quantity;   // factor table per dimension; sem.unit_mismatch across dimensions
declare function toSeconds(n: Node): number;                            // (t 12.5 "s"), (duration 2 "ms"), (q 3 "min")
declare function valueAt<T>(s: Signal<T>, t: number, lerp: (a: T, b: T, u: number, out?: T) => T, out?: T): T;   // columns + binary search, writes into out (§hot-loop-allocations, meso's `tables` (not adopted here))
```

Vocabularies validate these with schema combinators from `repr-concepts` (`S.vec(3)`, `S.quantity("length")`, `S.box(3)`, `S.signal(S.vec(3))`, `S.table({ ... })`).

## Core concepts, declared once (§semantic-layer)

```ts
const core = (id: string, description: string, isa?: string[]): ConceptSpec => ({ id, description, isa, owner: "repr-concepts" });

export const coreConcepts: ConceptSpec[] = [
  core("Scalar", "A single number."), core("Quantity", "A scalar with a unit."),
  core("Vector", "A fixed-size numeric vector."), ...[1, 2, 3, 4].map((n) => core(`Vec${n}`, `${n}D vector.`, ["Vector"])),
  core("Interval", "Closed range over an ordered type."), core("Box", "Axis-aligned region in 1/2/3D."),
  core("Transform", "Affine transform between spaces."), core("Time", "Instant on a timeline."), core("Duration", "Time difference."),
  core("Signal", "A value that varies over Time."), core("Table", "Records of one schema (a dataframe), in row or column layout."), core("Tensor", "N-dimensional typed array; bytes are referenced data."),
  core("Id", "Stable entity id."), core("Span", "A range in a text file."),
];
// domains type their attributes by core concept: Residue.position = { concept: "Vec3", dimension: "length" } (§semantic-layer)
```

## Generic operations, defined once on core types

Generic operations: see [README.md](README.md#shape).

```ts
// lib/repr/concepts/src/query.ts
export interface SpatialAccess { position(n: Node): Vec<3> | undefined; bounds?(n: Node): Box<3> | undefined }

export function withSpatial(L: QueryLang, acc: SpatialAccess): QueryLang {
  const GRID = Symbol("spatial-grid");                                        // per accessor: languages never share a grid
  return L
    .register("within", ([radius, target], compile) => lift([compile(target)], ([t]) => (n, p, ctx) => {
      const pos = acc.position(n); if (!pos) return false;
      const grid = ctx.cache.forRoot(ctx.root, GRID, () => buildSpatialGrid(ctx.root, acc));   // per snapshot
      return grid.anyWithin(pos, radius as number, t, ctx);                   // radius in the space's units; candidates carry paths
    }))
    .register("inside", ([box]) => { const b = toBox3(box); return (n) => { const p = acc.position(n); return !!p && contains(b, p); }; })
    .register("intersects", ([box]) => { const b = toBox3(box); return (n) => { const bb = acc.bounds?.(n); return !!bb && intersects(bb, b); }; });
}

export interface TemporalAccess { time(n: Node): number | Interval<number> | undefined }   // seconds
export const withTemporal = (L: QueryLang, acc: TemporalAccess): QueryLang => L
  .register("before", ([t]) => { const s = toSeconds(t); return (n) => { const v = acc.time(n); return v !== undefined && startOf(v) < s; }; })
  .register("during", ([iv]) => { const w = toInterval(iv); return (n) => { const v = acc.time(n); return v !== undefined && overlaps(v, w); }; });
// helpers assumed: buildSpatialGrid, toBox3, contains, intersects, startOf, toInterval, overlaps

// a domain reuses them: molecular language = core + spatial + its own functions
export const molLang = (langs: LanguageRegistry) => withSpatial(withCore(new QueryLang("mol.select", langs)), { position: atomPosition })
  .register("residue-name", ([name]) => (n) => tag(n) === "residue" && equalNode(attrOf(n, "name"), name));
```

## Embedded languages: slots, late binding, mutual recursion

Slots and late binding: see [README.md](README.md#shape).

```ts
export type Slot =
  | { kind: "core"; concept: string }                          // Vec3, Interval, ...
  | { kind: "atom"; type: "string" | "number" | "boolean" }
  | { kind: "optional"; of: Slot }
  | { kind: "list"; of: Slot }                                 // all values of the attribute: (points (vec3 ..) (vec3 ..))
  | { kind: "any-of"; options: readonly Slot[] }
  | { kind: "lang"; lang: string; form?: string }              // a DATA form of another language
  | { kind: "query"; lang: string }                            // a (select ...) in another language's query DSL
  | { kind: "link"; to: string; edge: string }                 // (link "<id>") to an entity of form `to`, edge kind `edge` (meso's `graphs` (not adopted here))
  | { kind: "table"; of: RecordSpec }                         // (rows ...) or (columns ...) of one record schema (meso's `tables` (not adopted here))
  | { kind: "any" };                                          // any well-formed node; validated by its head's language if it has one

export interface FormSpec { tag: string; slots: Readonly<Record<string, Slot>>; realizes?: string; attrs?: Record<string, string>; doc: string }
export interface Language {
  id: string;
  forms: readonly FormSpec[];
  queries?: readonly { lang: QueryLang; engine?: QueryEngine }[];   // one or more query languages (§queryable-representations); the first is the selection language
  dialects?: readonly Dialect[];                                    // surface syntaxes lowered to one of them (§query-dsls, §dialects-lower)
}

export const embedErrors = defineErrors("repr", {
  "repr.embed.unknown_language": { description: "A slot names a language that is not registered.", retryable: false },
  "repr.embed.too_deep":         { description: "Embedding deeper than the limit.",                retryable: false },
  "repr.validate.type":          { description: "Value of the wrong kind for its slot.",           retryable: false },
  "repr.validate.missing":       { description: "Required slot is missing.",                       retryable: false },
  "repr.validate.unknown_tag":   { description: "Not a form of the slot's language.",              retryable: false },
  "repr.ref.cycle":              { description: "Named definitions or data refs form a cycle.",    retryable: false },
});

// repr-core: built in the composition root, passed to QueryLang and the semantic graph (§semantic-layer)
export class LanguageRegistry {
  private langs = new Map<string, Language>();
  register(l: Language) {
    if (this.langs.has(l.id)) throw new Error(`duplicate language ${l.id}`);           // programmer error (§programmer-errors), as in QueryLang
    this.langs.set(l.id, l); return this;
  }
  get(id: string): Language {
    const l = this.langs.get(id);
    if (!l) throw embedErrors.error("repr.embed.unknown_language", `no language ${id}`, { hint: `registered: ${[...this.langs.keys()].join(" ")}` });
    return l;
  }
  all(): readonly Language[] { return [...this.langs.values()]; }
  /** A representation's query language by id; no id = its first (selection) language. */
  query(repr: string, id?: string): QueryLang { return this.queryEntry(repr, id).lang; }
  engine(repr: string, id?: string): QueryEngine {
    const e = this.queryEntry(repr, id).engine;
    if (!e) throw embedErrors.error("repr.embed.unknown_language", `no engine for ${id ?? repr}`);
    return e;
  }
  dialect(repr: string, id: string): Dialect | undefined { return this.get(repr).dialects?.find((d) => d.id === id); }
  private queryEntry(repr: string, id?: string) {
    const qs = this.get(repr).queries ?? [];
    const e = id ? qs.find((q) => q.lang.id === id) : qs[0];
    if (!e) throw embedErrors.error("repr.embed.unknown_language", `no query language ${id ?? "(default)"} for ${repr}`, { hint: qs.map((q) => q.lang.id).join(" ") });
    return e;
  }
  /** Composition check: every language id named by a slot resolves. */
  unresolved(): string[] {
    const ids = (s: Slot): string[] => s.kind === "lang" || s.kind === "query" ? [s.lang]
      : s.kind === "optional" || s.kind === "list" ? ids(s.of) : s.kind === "any-of" ? s.options.flatMap(ids) : [];
    return this.all().flatMap((l) => l.forms.flatMap((f) => Object.values(f.slots).flatMap(ids))).filter((id) => !this.langs.has(id));
  }
}
```

Mutually recursive languages, by id:

```ts
langs.register({ id: "scene", queries: [{ lang: sceneLang, engine: scanEngine(sceneLang) }], forms: [
  { tag: "mesh", realizes: "Shape", attrs: { position: "pos" }, doc: "...", slots: {
      id: { kind: "atom", type: "string" }, pos: { kind: "core", concept: "Vec3" },
      style: { kind: "optional", of: { kind: "lang", lang: "style", form: "rule" } } } },
]});
langs.register({ id: "style", forms: [
  { tag: "rule", doc: "Paint mol targets while the time window holds.", slots: {
      target: { kind: "query", lang: "mol" }, paint: { kind: "atom", type: "string" },
      when: { kind: "optional", of: { kind: "core", concept: "Interval" } } } },
]});
langs.register({ id: "mol", queries: [{ lang: molQuery, engine: scanEngine(molQuery) }], forms: [
  { tag: "residue", realizes: "Residue", attrs: { name: "name", index: "seq" }, doc: "...", slots: {
      name: { kind: "atom", type: "string" }, pos: { kind: "core", concept: "Vec3" } } },
]});
// sceneLang, molQuery = molLang(langs): query languages built with the same registry
// scene -> style -> (mol query) -> visible-in-scene -> (scene query): a cycle in DATA; the MODULES stay acyclic.
```

The sample tree these registrations validate: see [README.md](README.md#shape).

## Validation and queries dispatch on the slot

```ts
export function validate(n: Node | undefined, slot: Slot, langs: LanguageRegistry, path: Path = [], depth = 0, max = 64): ErrorInfo[] {
  if (depth > max) return [{ code: "repr.embed.too_deep", message: "embedding too deep", path }];
  if (n === undefined) return slot.kind === "optional" ? [] : [{ code: "repr.validate.missing", message: "missing", path }];
  switch (slot.kind) {
    case "optional": return validate(n, slot.of, langs, path, depth, max);
    case "core":     return coreValidator(slot.concept)(n, path);                       // Vec3 arity, Quantity unit, Box ordering
    case "atom":     return typeof n === slot.type ? [] : [{ code: "repr.validate.type", message: `expected ${slot.type}`, path }];
    case "list":     return validate(n, slot.of, langs, path, depth, max);              // called once per value by the form case
    case "any-of": {
      const tries = slot.options.map((o) => validate(n, o, langs, path, depth, max));
      return tries.some((e) => e.length === 0) ? [] : tries[0];
    }
    case "query": {                                                                      // compile with that language's query DSL
      if (tag(n) !== "select") return [{ code: "query.expected_select", message: "expected (select ...)", path }];
      try { langs.query(slot.lang).compile((n as readonly Node[])[1]); return []; }
      catch (e) { if (isCoded(e)) return [{ code: e.code, message: e.message, hint: e.hint, path }]; throw e; }
    }
    case "link":                                                                         // shape only; the target is checked at commit (meso's `graphs` (not adopted here))
      return tag(n) === "link" && typeof (n as readonly Node[])[1] === "string" ? [] : [{ code: "repr.validate.type", message: "expected (link \"<id>\")", path }];
    case "any":                                                                          // any well-formed node
      return [];                                         // Open: also validate by the language that owns the head tag, if any
    case "table":                                                                        // (rows ...) or (columns ...) of one record schema (meso's `tables` (not adopted here))
      return validateTable(n, slot.of, langs, path, depth + 1, max);
    case "lang": {                                                                       // late binding by id
      const lang = langs.get(slot.lang);
      const form = lang.forms.find((f) => f.tag === tag(n) && (!slot.form || f.tag === slot.form));
      if (!form) return [{ code: "repr.validate.unknown_tag", message: `not a ${slot.lang} form`, path, hint: lang.forms.map((f) => f.tag).join(" ") }];
      const kids = n as readonly Node[];
      return Object.entries(form.slots).flatMap(([k, s]) => {
        const i = kids.findIndex((c) => tag(c) === k);
        if (i < 0) return s.kind === "optional" || s.kind === "list" ? [] : [{ code: "repr.validate.missing", message: k, path }];
        return s.kind === "list"
          ? attrsOf(n, k).flatMap((v, j) => validate(v, s, langs, [...path, i, j + 1], depth + 1, max))
          : validate(attrOf(n, k), s, langs, [...path, i, 1], depth + 1, max);
      });
    }
  }
}
// helpers assumed: coreValidator (per core concept), validateTable (fields per RecordSpec, packing, repr.table.* codes; meso's `tables` (not adopted here))
```

Cross-language query functions: see [README.md](README.md#shape). Registered in the composition root:

```ts
// composition root
langs.query("mol").register("visible-in-scene", ([q], _compile, langs) => {
  const sceneEngine = langs.engine("scene");                                        // late-bound by id
  return { async prepare(ctx) {
    const root = sceneState.value;                                                  // TreeState of the scene (§representations)
    const shared = ctx.cache.forRoot(ctx.root, `visible-in-scene:${hashNode(root)}:${hashNode(q)}`,   // keyed by BOTH roots
      async () => new Set((await sceneEngine.runAsync(root, q, { ...ctx, signal: undefined })).flatMap((p) => idAt(root, p) ?? [])));
    const ids = await untilAborted(shared, ctx.signal);                             // shared run, as in 13: callers stop waiting on abort
    return (n) => ids.has(depictedShapeIdOf(n));                                    // ids, not Path arrays: Set compares by value
  } };
});
// helpers assumed: idAt (idOf(nodeAt(root, p))), depictedShapeIdOf; untilAborted: 13
```

Termination and the registry check: see [README.md](README.md#invariants); where depth is counted: see [Open](README.md#open).

## Conformance and tests

```ts
export const languageRegistrySuite = suite<LanguageRegistry>("LanguageRegistry",
  kase("every embedded language id resolves", "quick", async (r) => expect(r.unresolved()).toEqual([])),
  kase("the mutually recursive sample validates and round-trips", "quick", async (r) => { /* scene -> style -> mol -> scene */ }),
  kase("depth limit yields repr.embed.too_deep, not a stack overflow", "quick", async (r) => { /* pathological nesting */ }),
  kase("errors carry the path through the embedding", "quick", async (r) => { /* ... */ }),
  kase("within agrees with brute force in every language that uses it", "medium", async (r) => { /* shared cases via accessors */ }),
);
```

Conformance cases: see [README.md](README.md#conformance-cases).
