# Core concepts and embedded languages {#core-concepts-and-embedding}

- **Implementations:** [TypeScript](typescript.md)

Rules §core-concepts: a small core vocabulary (quantities, vectors, intervals and boxes, transforms, time, signals, tables, tensors, ids, spans) that every domain reuses, with generic operations defined once on it. DSLs may embed each other, even mutually recursively, through slots resolved late by id, while modules stay acyclic.

## Problem

Without a shared core, every domain reinvents vectors, units, time and intervals, and rewrites spatial and temporal
predicates with slightly different semantics. Languages that need each other (a scene's style embedding a molecular
query that refers back to the scene) then force either cyclic module dependencies or one merged language that nobody
owns.

## Shape

- **Core forms**, one spelling for every domain:
  ```
  (q 1.5 "nm")                                    ; Quantity: value + unit (dimension from the unit)
  (vec3 1.0 2.0 3.0)  (vec2 3 4)  (vec1 7)        ; vectors (vec4 too); units come from the owning attribute
  (box3 (vec3 0 0 0) (vec3 1 1 1))                ; axis-aligned box
  (interval 0 2.5)                                ; closed interval over any ordered core type
  (t 12.5 "s")  (duration 2 "s")                  ; instant, duration: time quantities named for readability
  (transform (mat4 1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1))
  (signal (interp linear) (keys (columns (t 0 2) (v (f64 (width 3) 0 0 0  1 0 0)))))   ; keys are a Table; Vec3 column packed (meso's `tables` (not adopted here))
  (columns (id "a1" "a2") (element "C" "N") (b-factor (f32 12.1 9.8)))           ; Table (dataframe), or (rows (atom ...) ...)
  (tensor (shape 1024 3) (dtype f32) (data (blob (hash "sha256:...") (size 12288))))   ; bytes behind a data ref (meso's `data-references` (not adopted here))
  (id "m_8f3a")  (span (file "a.ts") (range 120 168))
  ```
  Types: a quantity is `value` + `unit`, its dimension (`length`, `angle`, `time`, `mass`, `ratio`, ...) from the unit;
  an interval is closed (`lo`, `hi`); a signal is a table of keys (`t` in seconds, `v`) with `interp` `step` | `linear`
  | `cubic`; a tensor has `shape`, `dtype` (`f32`, `f64`, `i32`, `u8`) and a data ref. Value functions:
  `convertQuantity(q, unit)` (a factor table per dimension), `toSeconds(node)`, and `valueAt(signal, t, lerp, out)`
  (binary search over the columns, writing into `out`, §hot-loop-allocations).
- **Core concepts**, declared once (§semantic-layer), owner `repr-concepts`: `Scalar`, `Quantity`, `Vector` (`Vec1` to
  `Vec4` isa `Vector`), `Interval`, `Box`, `Transform`, `Time`, `Duration`, `Signal`, `Table`, `Tensor`, `Id`, `Span`.
  Domains type their attributes by core concept (`Residue.position`: `Vec3`, dimension `length`).
- **Generic operations** (§generic-operations), written once on core types; a domain supplies only an accessor:
  - spatial access `position(node) -> Vec3?`, optional `bounds(node) -> Box3?`; adds `(within r <pred>)`,
    `(inside <box>)`, `(intersects <box>)`. `within` uses a spatial grid built per snapshot and per accessor.
  - temporal access `time(node) -> seconds or interval?`; adds `(before <t>)`, `(during <interval>)`.
  - The molecular language is core + spatial (with the atom position) + its own `residue-name`.
- **Slots**: a form declares where core forms or other languages may appear. Slots name languages by id; the registry
  resolves them at use time.

  | Slot kind | Holds |
  |---|---|
  | `core <concept>` | a core form (`Vec3`, `Interval`, ...) |
  | `atom <type>` | a string, number or boolean |
  | `optional <slot>`, `list <slot>`, `any-of <slots>` | combinators; `list` is every value of the attribute |
  | `lang <id> [form]` | a data form of another language |
  | `query <id>` | a `(select ...)` in another language's query DSL |
  | `link <form> <edge>` | `(link "<id>")` to an entity of that form, with that edge kind (meso's `graphs` (not adopted here)) |
  | `table <record>` | `(rows ...)` or `(columns ...)` of one record schema (meso's `tables` (not adopted here)) |
  | `any` | any well-formed node |

- **`FormSpec`**: `tag`, `slots`, optional `realizes` (a concept) and `attrs` (concept attribute → local attribute),
  `doc`. **`Language`**: `id`, `forms`, optional `queries` (one or more query languages with engines; the first is the
  selection language, §queryable-representations) and `dialects` (§query-dsls).
- **`LanguageRegistry`**, built in the composition root and passed to query languages and the semantic graph:
  `register(language)`, `get(id)`, `all()`, `query(repr, id?)` and `engine(repr, id?)` (no id: the selection
  language), `dialect(repr, id)`, and `unresolved()` (language ids named by slots that are not registered).
- **Mutual recursion, by id**: `scene` forms hold an optional `style` rule; a `style` rule holds a `mol` query; that
  query may call `visible-in-scene` with a `scene` query. A cycle in data; the modules stay acyclic.
  ```
  (mesh (id "m_1") (pos (vec3 0 1 0))
    (style (rule (target (select (and (residue-name "ALA") (visible-in-scene (select (inside (box3 (vec3 -1 -1 -1) (vec3 1 1 1))))))))
                 (paint "red")
                 (when (interval (t 0 "s") (t 2 "s"))))))
  ```
  Three languages and the core forms sit in one uniform tree, so it parses, prints, hashes and streams as one
  (§representations, §committed-transient-identity).
- **Validation** `validate(node, slot, langs, path, depth, max) -> errors` dispatches on the slot kind; a `lang` slot
  finds the form by tag and validates each of its slots, a `query` slot compiles with that language's query DSL.
- **Cross-language query functions**: a function compiler receives the registry (§query-dsls), so it can compile an
  argument in another language at compile time. A function that needs another representation's live state (like
  `visible-in-scene`) is registered in the composition root, the only place that knows that state
  (§dependencies-passed-in); it is staged and its shared result is keyed by both roots and the embedded query.
- **Errors**: `repr.embed.unknown_language`, `repr.embed.too_deep`, `repr.validate.type`, `repr.validate.missing`,
  `repr.validate.unknown_tag`, `repr.ref.cycle`; none retryable.

## Layout

```
lib/
  repr/core/              # tree, print/parse, ops, hashing, query core (§representations, §query-dsls, §committed-transient-identity) + LanguageRegistry (below)
  repr/concepts/          # THIS (@scope/repr-concepts): concepts, forms, generic query functions, unit/time conversions
  repr/molecule/          # depends on repr-core + repr-concepts; registers its language
  repr/scene/             # same; style is a language of its own, not a dependency of mol (or vice versa)
```

Languages refer to each other by id, in data. Only the composition root sees them all.

## Invariants

- Modules stay acyclic; cycles between languages exist only in data, through ids resolved late.
- The composition root registers every language, then asserts `unresolved()` is empty. Registering a language id twice
  is a programmer error; an unknown id is `repr.embed.unknown_language` with the registered ids in the hint.
- Termination (§bounded-recursion): trees are finite, so syntactic embedding terminates; evaluation and validation depth
  is bounded (`repr.embed.too_deep`, never a stack overflow); only references (named definitions, data refs) are
  cycle-checked (`repr.ref.cycle`).
- Validation errors carry the path through the embedding.
- Validation by slot: a missing slot that is neither optional nor a list is `repr.validate.missing`; an atom of the
  wrong type is `repr.validate.type`; a tag that is not a form of the slot's language is `repr.validate.unknown_tag`,
  with the language's tags in the hint; a `query` slot that is not `(select ...)` is `query.expected_select`; `any-of`
  passes when any option passes; a `link` slot checks shape only, the target is checked at commit (meso's `graphs` (not adopted here)).
- Spatial and temporal predicates are defined once; a domain supplies only how its nodes expose a position or a time.
  Each accessor has its own spatial grid per snapshot: languages never share a grid.
- A cross-language function's shared result is keyed by both roots, and callers stop waiting on abort without cancelling
  the shared run (as in §query-dsls).
- Unit conversions are exact value functions; converting across dimensions is `sem.unit_mismatch`.
- Core forms round-trip through print/parse with stable hashes. `valueAt` writing into `out` allocates nothing.

## Conformance cases

- The registry → `unresolved()` is `[]`.
- The mutually recursive sample (scene → style → mol → scene) → validates and round-trips.
- Pathological nesting → `repr.embed.too_deep`, not a stack overflow.
- An invalid node deep in an embedding → an error whose path runs through the embedding.
- `within` → agrees with brute force in every language that uses it (shared cases via accessors).
- Core forms → print/parse round-trip with stable hashes.
- Unit conversions → golden values, invertible within tolerance. To be decided: the values.
- `valueAt` with `out` → no allocation (guard tier).

## Open

- Where the evaluation depth is counted: the query evaluation context carries no depth, and `validate` uses its own
  `max`; whether a depth field rides on the context passed to embedded runs (and `validate` reads `limits.maxDepth`) is
  not settled.
- Whether an `any` slot also validates its node by the language that owns the head tag.
