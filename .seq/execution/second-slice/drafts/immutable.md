# Persistent collections in `@molstar/core/data/immutable` (draft, for the author's review)

A draft for Gate 0 of the second slice, beside `patches.md`. Once approved, it is agreed into
`.seq/design/immutable.md` and base-library.md lists the entry. Lane R builds exactly the surface in
[The API](#the-api); anything past it stops the lane.

The author, 2026-11-18:
- "shall we do a little yak shaving and implement packages/core/src/immutable? based on immutable.js and possibly other
  similar (perhaps non js) libraries as well? that way we can control all allocations, hashing, etc, make it
  configurable and extend as we need (e.g. we could do clojure hashtries etc when we need)";
- "we would move these to our immutable as well (in a later separate PR, here we would just implement core/immutable
  and use it in repr)".
- "the collections should line in core/src/data/immutable, next to the other collections" (`repr` stays where it is: it
  is a language framework, a peer of `state`, and builds on `data`).

## Decisions for the author first

1. **What `repr` uses it for** ([In `repr`](#in-repr)). Recommendation: **A**, persistent indexes beside the tree,
   with the node shape unchanged. The alternative is **B**, nodes holding persistent collections. B is kept as a
   measured trigger, not built.
2. **The map is CHAMP**, not Bagwell's HAMT as in immutable.js and Clojure. Recommendation: yes. Its canonical form
   gives structural equality and a stable iteration order for free, and it iterates faster.
3. **The API is our own, close to the subset of immutable.js Mol\* uses.** Recommendation: yes. The state tree's later
   move is then mostly renames, and we carry none of immutable.js's Seq, Record or lazy operations.
4. **Configurable per collection** (hash, equality, the small-collection threshold), not globally. Recommendation:
   yes. Why: one global setting would change every other user's maps.

## Sources

Read for the design; the lane reads the code where it needs to:

| Source | What we take | What we leave |
|---|---|---|
| **immutable.js 5** (the repository's dependency, 5.1.9) | the API shape Mol\* uses (`get`, `set`, `has`, `delete`, `asMutable`, `asImmutable`, `withMutations`); `ownerID` transients (one token per batch, nodes it made edited in place); `ArrayMapNode` (a plain array of entries up to 8); `OrderedMap` as a map plus an entry list, compacted when holes pile up | `Seq` and lazy operations, `Record`, `hashCode`/`equals` on arbitrary objects, a global string-hash cache |
| **Clojure** (`PersistentHashMap`, `PersistentVector`, transients) | the vector: a 32-way trie with a tail of up to 32, so an append touches one array almost always; transients that fail when used after `persistent!` | thread-owned edits (one thread here) |
| **CHAMP** (Steindorfer and Vinju, *Optimizing Hash-Array Mapped Tries for Fast and Lean Immutable JVM Collections*, OOPSLA 2015; Scala 2.13's `HashMap`) | two bitmaps a node (inline entries, sub-nodes), entries first; deletion that compacts back to one canonical shape, so equal maps are equal trees; iteration over contiguous entries | the JVM-specific node specialisations |
| **Bagwell** (*Ideal Hash Trees*, 2001) | the 32-way, 5-bits-a-level trie; hash-collision nodes at the bottom | — |
| **Rust `im` and `rpds`** | `im`'s RRB vector (Bagwell and Rompf 2011; Stucki et al. 2015) as the later step for concatenation and slicing; `rpds`'s small, explicit API | `im`'s B-tree `OrdMap` (no sorted map is needed yet) |

## Decision

### The collections

- **`HashMap<K, V>` and `HashSet<K>`**: CHAMP tries with 5 bits a level. Up to `smallMax` entries (default 8), a
  collection is one flat array scanned linearly, as immutable.js's `ArrayMapNode` is. A hash collision at the bottom
  is a collision node.
- **`OrderedMap<K, V>` and `OrderedSet<K>`**: insertion order, as Mol\*'s state tree and custom properties need. A
  `HashMap` from key to slot, plus a `Vector` of entries. A removal leaves a hole, and the vector is compacted when
  holes pass half its length (immutable.js's rule).
- **`Vector<T>`**: a 32-way trie with a tail (Clojure). Append, update, index and pop are effectively O(1) (log32).
  Concatenation and slicing copy for now: the RRB vector is the trigger-gated next step.

### Allocation

- **Reads allocate nothing**: `get`, `has`, `size`, `forEach`, and the cursors (`entries(cursor)` with a reusable
  cursor object). No iterator object is made on a hot path. `[Symbol.iterator]` exists for convenience and is
  documented as allocating.
- **A persistent edit** copies the path: at most ⌈32/5⌉ = 7 small arrays for a map, far fewer in practice. An edit
  that changes nothing (setting an equal value, deleting an absent key) returns the same object.
- **Transients**: `asMutable()` gives a collection that edits its own nodes in place under an owner token.
  `asImmutable()` freezes it. Using a transient after it was frozen throws (Clojure's rule): a programmer error, not an
  issue. `withMutations(fn)` wraps the two. Why: a batch (building an index, a patch's ops) pays one path copy a node,
  not one an edit.
- **No global caches.** Nothing like immutable.js's string-hash cache: a collection's hash function is its own.

### Hashing and equality

- **Per collection**, given at creation (`config`): `hash(key) => int32` and `equals(a, b) => boolean`.
  - **Defaults**: strings by core's `hashString`, numbers by their bits (−0 and 0 apart, one NaN), `Object.is` for
    equality, and objects by identity (a `WeakMap` id).
  - **Preset configs**: `ValueKeys` uses `repr`'s `valuesEqual` and the cached `hashNode` for nodes and columns, so a
    node or a column can be a key by content.
- **Collections are values**: `equals(a, b)` is structural, short-circuited by identity. CHAMP's canonical shape makes
  it a walk of two equal trees, never a lookup per entry. `hashOf(collection)` is order-independent for maps and sets,
  ordered for vectors, cached on the collection object.
- **This is not `repr`'s Merkle hash.** A trie's 32-bit key hash places keys; `hashNode`'s 128-bit content hash names
  revisions. The two never mix.

### Configuration

Each collection takes an optional config with these items and defaults. The trie's 5 bits a level is fixed: making it
configurable costs a branch on every access, for nothing measured.

| Item | Default | Why |
|---|---|---|
| `hash` | `defaultHash` | as above |
| `equals` | `Object.is` | as above |
| `smallMax` | 8 | below it, a flat array beats a trie (immutable.js's measured choice; ours re-measured in the lane) |

### In `repr`

**A (recommended): persistent indexes beside the tree; nodes unchanged.** `repr`'s nodes stay plain JSON-shaped
objects. They are the public contract the codecs, the builder, Seq\* and later meso\* read, and the shape just agreed.
The collections back the indexes `repr` needs, so each one is derived from the previous revision's in O(changes):

- **A child index a node**: id to position, built when a node with more than `smallMax` children is first looked up,
  and cached on the node (a `WeakMap`). `resolveAddress` uses it instead of its linear `find` at each step. A patch
  derives the new node's index from the old one by the ops' edits, not by a rebuild.
- **The address index of a revision**: node to address. `addressOf` uses it instead of a walk of the whole tree.
  Derived per patch from the touched paths.
- **The ref index of a revision** (`patches.md`, Applying 3): the ref values by the address they name, built at the
  first move. Each later revision derives it from its parent's by the touched nodes, so a move costs O(refs it touches)
  plus the path, never a walk of the spec.
- **`diff` and `merge3`**: their per-parent id maps (transients, frozen once built).

**B (not built; a measured trigger): nodes holding persistent collections**, `children` as a `Vector` and `params` as
a `HashMap`. An edit would copy log32 of a container instead of its whole children array. The cost is the public node
shape: every reader of `node.children`, the codecs, the builder, the hash's canonical text and Seq\*'s whole core. The
trigger: a measured spec in which copying a children array dominates an edit. A guess: a container with about 10,000
children, such as a refs list of bookmarks. Then a change to the node shape, its own design.

## The API

The agreed surface for lane R. The folder `packages/core/src/data/immutable/`, next to `data/int`, `data/generic` and
`data/db`; its entry `data/immutable.ts` (as `data/int.ts` is for `int/`), so `@molstar/core/data/immutable` through
the package's wildcard export, and `Immutable` beside `Int` and `Generic` in `data/index.ts`. One module a collection,
plus `config` and `hash`.

```ts
interface CollectionConfig<K> {
  readonly hash?: (key: K) => number;
  readonly equals?: (a: K, b: K) => boolean;
  readonly smallMax?: number;
}
const DefaultConfig: CollectionConfig<unknown>;
function defaultHash(key: unknown): number;

interface HashMap<K, V> {
  readonly size: number;
  get(key: K): V | undefined;
  has(key: K): boolean;
  set(key: K, value: V): HashMap<K, V>;
  delete(key: K): HashMap<K, V>;
  update(key: K, fn: (value: V | undefined) => V | undefined): HashMap<K, V>;   // undefined removes
  forEach(fn: (value: V, key: K) => void): void;
  entries(cursor?: MapCursor<K, V>): MapCursor<K, V>;     // a reusable cursor: next() fills .key/.value
  asMutable(): MutableHashMap<K, V>;
  withMutations(fn: (m: MutableHashMap<K, V>) => void): HashMap<K, V>;
  equals(other: HashMap<K, V>): boolean;
  hashCode(): number;
  [Symbol.iterator](): Iterator<[K, V]>;                  // allocates: not for hot paths
}
interface MutableHashMap<K, V> { /* the same edits, in place, returning this */ asImmutable(): HashMap<K, V> }
namespace HashMap {
  function empty<K, V>(config?: CollectionConfig<K>): HashMap<K, V>;
  function of<K, V>(entries: Iterable<readonly [K, V]>, config?: CollectionConfig<K>): HashMap<K, V>;
}
// HashSet<K>, OrderedMap<K, V>, OrderedSet<K>: the same shape (sets: add, delete, has; ordered: first, last, in order)
// Vector<T>: get, set, push, pop, size, forEach, cursor, asMutable, withMutations, equals, hashCode; Vector.of, Vector.empty

// repr's preset (in @molstar/core/repr/tree, so immutable needs no repr)
const ValueKeys: CollectionConfig<Value>;
```

`repr`'s own surface does not change: `resolveAddress` and `addressOf` keep their signatures and get faster. The
indexes are internal, except what `patches.md` exposes (`touched`).

## Tests

- **Model-based property tests**: random sequences of operations on each collection, checked against a JS `Map`, `Set`
  or `Array` after each step. Seeded, short in the quick lane, long in the heavy lane.
- **Persistence**: every earlier version is unchanged after any later edit, transients included.
- **Transients**: in-place edits under one token, isolation from other versions, a throw after `asImmutable`.
- **Collisions**: a deliberately bad hash (every key 0) gives correct results.
- **Canonical shape**: equal contents reached by different edit orders give equal trees (CHAMP), and `equals` and
  `hashCode` agree.
- **Allocation** (the opt-in tier): `get`, `has` and cursor iteration allocate 0 bytes warm. An edit stays within a
  stated bound.
- **Against immutable.js** (already a dev dependency here): get, set, delete, iterate and build at 10, 1,000 and
  100,000 entries, as heavy-lane ratios. The target is no slower than immutable.js on any of them.
- **In `repr`**: `resolveAddress` on a 1,000-child container, before and after; the ref index derived across 100
  patches equals one built from scratch at every step.
- **Bundle**: the entry's minified size, recorded. The state tree's later move then shows what replacing immutable.js
  saves.

## Not in this slice

- **Moving Mol\*'s state tree and `@molstar/model`'s custom properties off immutable.js**: a later, separate PR. The
  API's closeness to the subset they use (decision 3) is what keeps that small.
- **The RRB vector, sorted maps, and B**: each when its trigger is met.
