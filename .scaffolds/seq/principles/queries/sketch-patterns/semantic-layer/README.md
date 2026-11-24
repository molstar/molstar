# Semantic layer: concepts, spaces, mappings {#semantic-layer}

- **Implementations:** [TypeScript](typescript.md)

Rules §semantic-concepts: a light semantic layer above all representations. Concepts say what things are, language forms realize them, spaces say how things are addressed, locators point, mappings (in runtimes) connect spaces, and a router composes the mappings. It is all data, queryable and reachable by command. Declare what is useful and skip the rest.

```
              ┌──────────── semantic layer (data) ────────────┐
              │   concepts    spaces + locators    mappings   │
              └──────┬──────────────┬──────────────────┬──────┘
   representations:  doc tree   UI state   scene   mol structure   DB schema   text   data blobs
   runtimes provide the mappings: renderer (project/pick), document (ids), mol module, data port ...
```

## Problem

Without it, each module addresses things its own way: bare indexes and `(x, y)` pairs cross module boundaries, numberings
are converted by shifting numbers by hand, each view invents its own way to show a shared selection, and a question by
meaning ("every residue named ALA", across the document and the structure) needs glue per pair of modules. Positions
that went stale are used silently.

## Shape

- **Concept** (`ConceptSpec`): `id`, `description`, `isa` (super-concepts), `attrs` (name → optional core concept and
  dimension, §core-concepts-and-embedding), `relations` (`rel`, `to`: part-of, refers-to, derived-from, depicts),
  `owner` (the module that owns the domain). A few attributes and relations are enough; there is no ontology machinery.
  Shared concepts live in `repr-concepts`, and products extend them.
- **Realization**: a language form says which concept it realizes (`realizes` plus an `attrs` map from concept
  attribute to local attribute, §core-concepts-and-embedding); there is no second registry.
- **Space** (`SpaceSpec`): `id`, `description`, `stability` (`stable`: ids; `positional`: positions that go stale;
  `ephemeral`: fleeting, like screen pixels), `coords` (the shape of the coordinates), optional `spatial` (`dims` 2 or 3,
  `units`, `handedness`, `up`, `parent`), `owner`.
- **`Locator`**: `space`, `at` (coordinates), optional `concept` (what it denotes), optional `anchor` (`hash`, `rev`: the
  snapshot a positional locator was made against). Printed as `(at <space> <coords...> [(is X)] [(anchor ...)])`:
  ```
  (at doc    (id "m_8f3a"))                                       ; stable entity
  (at screen (viewport "v1") (xy 412 188))                        ; ephemeral pixel
  (at world  (xyz 1.5 0.0 -3.2))                                  ; m, right-handed, y-up
  (at mol.auth (chain "A") (num 42 45) (is Residue))             ; a range in one of several residue numberings
  (at doc (realm "studio") (page "welcome") (block "w3") (is Task)) ; a block of a document, as a concept instance (meso's `documents` (not adopted here))
  (at doc.text (page "welcome") (block "w3") (range 4 9) (anchor (hash "v:3c1a")))   ; UTF-16 offsets, positional
  (at thread "th_1" (msg "m1"))                                   ; a message in a thread (meso's `documents` (not adopted here), agent-harness's `anchored-threads` (not adopted here))
  (at text   (file "src/foo.ts") (range 120 168))
  (at data   (blob (hash "sha256:9f2c...")) (bytes 4096 8192))    ; the coordinate is the data ref node (meso's `data-references` (not adopted here))
  ```
  Quantities carry their unit (`(q 1.5 "nm")`); conversions are exact value functions (§core-concepts-and-embedding),
  and spatial spaces declare their units, so px/world/Å conversions are explicit.
- **`Mapping`**: `from`, `to`, `cost`, `fidelity` (`exact` | `approximate` | `lossy`, a declared upper bound) and
  `translate(locator, ctx) -> {results: locators (none, one or many), fidelity}`, sync or async (§api-surfaces). The
  context is a `QueryContext` (§query-dsls) plus `snapshotHash(space)`. The runtime that knows a mapping registers it in
  the composition root (§ports-and-adapters), closing over its sources: picking is `screen -> doc` (approximate), a
  pixel is a ray `screen -> world`, projection is `world -> screen` (exact). Spatial transforms are ordinary mappings
  with matrices, units and handedness checked; hot paths write into caller-provided buffers (§hot-loop-allocations).
- **The `space` language** (§addressing-language): most spaces and many mappings are data, compiled to mappings. Any
  domain with several addressing schemes has this shape (text offsets vs line/column; frames vs timecode; row index vs
  primary key):
  ```
  (space (id mol.label) (of Residue) (stability stable)              ; numbering assigned by the format
    (key (entity string) (asym string) (seq int))  (range seq))
  (space (id mol.auth) (of Residue) (stability stable)               ; numbering chosen by the authors
    (key (chain string) (num int) (ins (optional string)))  (range num ins))
  (space (id seq) (of Residue) (stability stable)                    ; position in a sequence
    (key (entity string) (pos int))  (range pos))
  (space (id ext.uniprot) (of Residue) (stability stable)            ; position in an external database
    (key (acc string) (pos int))  (range pos))

  (mapping (from mol.label) (to mol.auth) (fidelity exact) (bijective true)
    (via (table (source (data-of mol)) (match asym seq) (yield chain num ins))))     ; the correspondence the structure carries
  (mapping (from mol.label) (to seq) (fidelity exact) (bijective true)
    (fields (entity entity) (pos seq)))                                              ; a rename
  (mapping (from seq) (to ext.uniprot) (fidelity approximate)
    (via (segments (source (ref (uri "...alignment..."))) (match entity pos) (yield acc pos))))  ; aligned ranges; gaps map to nothing
  (mapping (from mol.auth) (to mol.auth-v1) (fidelity exact) (bijective true)
    (offset num -1000))                                                              ; arithmetic, when it really is arithmetic
  ```
  Mapping forms: `fields` (rename, reorder), `offset`, `table` (lookup in a data ref or `(data-of <domain>)`),
  `segments` (range to ranges: alignments), `compose`, `inverse` (only of a bijective mapping). A compiled mapping keeps
  one index per source snapshot: a hash map for `table`, sorted segments with binary search for `segments`
  (§efficient-structures).
- **The router** (`SemanticGraph`): built from concepts, spaces, mappings and the `LanguageRegistry`.
  `path(from, to)` picks the best chain, memoized: for each fidelity rank in order (exact, approximate, lossy), the
  cheapest chain over mappings of at most that rank; the first rank that connects wins. `convert(locator, to, ctx)`
  runs the chain hop by hop, fanning out and deduplicating. `realizationsOf(concept)` lists the forms realizing a concept
  or any sub-concept (via `isa`). Relation traversal (`depicts`, `part-of`), when needed, is mappings tagged with `rel`.
- **Navigate by meaning**: `semanticSelect(graph, roots, concept, where, ctx) -> locators` fans out over every form that
  realizes the concept and runs that language's selection query (§query-dsls), built as a tree:
  ```
  (select-concept Residue (where (name "ALA")))
    => (at doc (id "r_17") (is Residue))  (at mol (chain "A") (res 42) (is Residue))
  ```
- **By command** (§commands): `sem.concepts`, `sem.path`, `sem.convert`, `sem.select`, all read-only. `sem.concepts`
  returns the layer as one representation (concepts, realizing forms, spaces, mappings with cost and fidelity), so the
  query DSL answers meta questions ("which spaces reach `mol`?") and agents discover routes. `_spec/` links to it, and the
  gate checks those links (§spec-mirror, §drift-triage). Agents convert instead of computing:
  ```
  → sem_read  (sem.convert (at mol.auth (chain "A") (num 42 45)) mol.label)
  ← (results (at mol.label (entity "1") (asym "A") (seq 38 41))) (fidelity exact)
  → sem_read  (sem.convert (at mol.label (entity "1") (asym "A") (seq 38 41)) ext.uniprot)
  ← (results (at ext.uniprot (acc "P69905") (pos 87 89))) (fidelity approximate)      ; one residue fell in an alignment gap
  ```
- **Everything that points uses locators**: query results, command arguments, error positions (§error-codes), data
  ranges (meso's `data-references` (not adopted here)), presence cursors (§committed-transient-identity). A selection shared between views
  (`SelectionModel`) is a set of locators; each view projects it into its own space (`inSpace`).
- **Errors**: `locator.no_path`, `locator.stale`, `locator.unresolved`, `locator.bad_locator`, `locator.cancelled`,
  `locator.limit_exceeded`, `sem.unknown_concept`, `sem.unit_mismatch`.

## Invariants

- Locators are representations: printable, hashable and writable by agents. A bare `(x, y)` or index crossing a module
  boundary is a bug (§spaces).
- A mapping's results validate against the target space's coordinates, and the fidelity it reports is never better than
  the one it declares. A compiled mapping that finds nothing reports `lossy`.
- A chain is as good as its weakest hop: `convert` reports the worst fidelity met.
- Routing prefers better fidelity over lower cost, and never routes through an ephemeral space (only as an endpoint).
- No chain is `locator.no_path`, with the spaces reachable from the source in the hint. An unknown space is
  `locator.bad_locator`.
- A positional locator whose anchor differs from the space's current snapshot is `locator.stale`; only the input's
  anchor is checked. An id locator survives edits.
- Each hop deduplicates results by locator hash; more than `maxResults` is `locator.limit_exceeded`; an aborted signal
  between hops is `locator.cancelled`.
- Keys outside a table or in an alignment gap map to nothing, never to a neighbor; a range maps to ranges that cover
  exactly the translated points; a bijective mapping round-trips.
- A shared index build is not tied to any caller's signal, and a failed build is not cached (as in §query-dsls).
- Realized concepts, `isa` and relations name existing concepts.
- Concept queries are built as trees, never strings, and equal the union of the per-language queries.
- Projecting a shared selection skips locators with no path to the view's space; cancelled, stale and limit errors still
  fail.

## Conformance cases

- Every mapping: results validate against the target coordinates; reported fidelity ≤ declared; an async mapping honors
  the abort signal.
- Generated from declarations, with no hand-written cases: every bijective mapping round-trips every sample key and
  range; ranges cover exactly the translated points; keys outside the table or in gaps map to nothing; every declared
  space validates its keys.
- The graph: no dangling concept ids; `locator.no_path` lists reachable spaces; an exact chain is chosen over a cheaper
  lossy one; no route passes through an ephemeral space; a stale positional locator is rejected and an id locator
  survives edits.
- Where both directions exist: `convert(convert(x, B), A)` contains `x` for exact chains.
- `semanticSelect(concept)` equals the union of the per-language queries.
- Spatial conversions: golden values for units and handedness. To be decided: the values.
- A unit conversion across dimensions → `sem.unit_mismatch`.
- To be decided: cases for `locator.unresolved` and `sem.unknown_concept`.

## Hardened by a prototype

Until a space language exists, an application declares its spaces as data (name, stable or positional, the locator
forms) and parses its in-app links as a dialect of those locators (dorothy's `doc`, `doc.text`, `thread`).

## Open

- The plain `mol` space (used in the mappings and the concept select above, and in §skills-and-tools and
  §agent-protocol, with `(res n)`) predates the numberings; whether it is an alias of `mol.auth` (`(num n)`) or a space
  of its own is not settled.
- How a language id relates to the spaces of its representation (`semanticSelect` keys roots by language id, its results
  are in `doc` and `mol`), and which space a result locator is in, is not settled.
- `locator.unresolved` and `sem.unknown_concept` are declared, but nothing in the sketch raises them yet
  (`sem.unit_mismatch` is raised by unit conversion across dimensions, §core-concepts-and-embedding).
