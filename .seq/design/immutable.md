# Persistent collections in `@molstar/core/data/immutable`

Agreed by the author at Gate 0 of the second slice (2026-11-18, "rest as recommended"), after two independent reviews
(Fable 5.2). base-library.md lists the entry. Lane R builds exactly the surface in [The API](#the-api); anything past
it stops the lane.

The author, 2026-11-18:
- "shall we do a little yak shaving and implement packages/core/src/immutable? based on immutable.js and possibly other
  similar (perhaps non js) libraries as well? that way we can control all allocations, hashing, etc, make it
  configurable and extend as we need (e.g. we could do clojure hashtries etc when we need)";
- "we would move these to our immutable as well (in a later separate PR, here we would just implement core/immutable
  and use it in repr)";
- "the collections should line in core/src/data/immutable, next to the other collections". `repr` stays where it is:
  it is a language framework, a peer of `state`, and builds on `data`.

**Draft 2** follows the second independent review (Fable 5.2), whose points 1 to 5 the author took ("let's go with
1-5"), then its focused re-review ("apply N-1 to N-14 as recommended"). Findings are cited where they land: R2-S12 to
R2-S16, R2-C7 to R2-C10, N-9 to N-11.

## Decisions (agreed)

1. **`repr` uses the collections behind an interface, and only once measured** (R2-S15, R2-S16; the author's point 5).
   Nodes keep their plain shape. Lane K builds `repr`'s ref index on plain `Map`s behind an internal interface, and that
   index moves onto this library's persistent trie in the step after lanes R and K, when a measured case asks for it
   ([In `repr`](#in-repr)). Recommendation: yes.
2. **The map is CHAMP**, not Bagwell's HAMT as in immutable.js and Clojure. Recommendation: yes. Deletion compacts to
   one shape for the trie part, so it iterates faster and compares by a walk where shapes can match. The flat and
   collision parts are not canonical (R2-S13).
3. **The API is our own, close to the subset of immutable.js Mol\* uses.** Recommendation: yes. A type-check-only
   spike of `StateTree` on it, inside lane R, validates this before the API is frozen (R2-C7).
4. **Configurable per collection**, for hash and equality only; the small-collection threshold is a measured constant
   (R2-S13). Recommendation: yes.

## Sources

Read for the design; the lane reads the code where it needs to:

| Source | What we take | What we leave |
|---|---|---|
| **immutable.js 5** (the repository's dependency, 5.1.9) | the API shape Mol\* uses (`get`, `set`, `has`, `delete`, `asMutable`, `asImmutable`, `withMutations`, `forEach` with a context, `first`, `toArray`, `values`); `ownerID` transients that read and iterate like the persistent ones, `asMutable()` on a transient returning itself; `ArrayMapNode` (a plain array of entries up to 8); `OrderedMap` as a map plus an entry list, compacted only past 32 entries when holes pass half | `Seq` and lazy operations, `Record`, `hashCode`/`equals` on arbitrary objects, a global string-hash cache |
| **Clojure** (`PersistentHashMap`, `PersistentVector`, transients) | the vector: a 32-way trie with a tail of up to 32 | thread-owned edits; transients that cannot be read or iterated |
| **CHAMP** (Steindorfer and Vinju, *Optimizing Hash-Array Mapped Tries for Fast and Lean Immutable JVM Collections*, OOPSLA 2015; Scala 2.13's `HashMap`) | two bitmaps a node (inline entries, sub-nodes), entries first; deletion that compacts back to one shape; iteration over contiguous entries | the JVM-specific node specialisations |
| **Bagwell** (*Ideal Hash Trees*, 2001) | the 32-way, 5-bits-a-level trie; hash-collision nodes at the bottom | — |
| **Rust `im` and `rpds`** | `im`'s RRB vector (Bagwell and Rompf 2011; Stucki et al. 2015) as the later step for concatenation and slicing; `rpds`'s small, explicit API | `im`'s B-tree `OrdMap` |

## Decision

### The collections

- **`HashMap<K, V>` and `HashSet<K>`**: CHAMP tries with 5 bits a level.
  - Up to a measured constant (about 8) entries, a collection is one flat array in insertion order, scanned linearly,
    as immutable.js's `ArrayMapNode` is.
  - A deletion demotes a trie back to flat at the threshold, so equal contents reach one of two shapes, not many.
  - A hash collision at the bottom is a collision node (an unordered bucket).
- **`OrderedHashMap<K, V>` and `OrderedHashSet<K>`**: insertion order, as Mol\*'s state tree and custom properties need
  (R2-C8). The names keep clear of `Int.OrderedSet`, beside them in `data/index.ts`.
  - A `HashMap` from key to slot, plus a `Vector` of entries.
  - Setting a key already present keeps its place.
  - A removal leaves a hole. The vector is compacted when holes pass half its length and it holds at least 32 entries
    (R2-C10).
- **`Vector<T>`**: a 32-way trie with a tail (Clojure). Append, update, index and pop are effectively O(1). It backs the
  ordered collections, and is public with what Mol\*'s uses need (R2-C7): `get`, `set`, `push`, `pop`, `shift`, `size`,
  `first`, `last`, `indexOf`, `find`, `remove(index)`, `forEach`, `toArray`. `shift` and `remove(index)` copy (O(n)) until the
  RRB vector, its trigger a measured use. `shift` is in the list too (N-11).

### Transients (R2-S12)

The contract `StateTree`'s `TransientTree` relies on:

- **A transient reads like the persistent collection**: `get`, `has`, `size`, `forEach`, `values`, `first`, `toArray`
  and cursors. `MutableHashMap extends` the read interface.
- **Edits happen in place** under the transient's owner token. Every node checks the owner before it is written, collision
  nodes and compacted entries included. A node not made under the token is copied first. Regression case: deleting
  from colliding keys must not write into a node another version shares (the bug `TransientTree.renewOwner` works
  around in immutable.js).
- **`asMutable()` on a transient returns itself.**
- **`asImmutable()` returns a new persistent collection sharing the transient's nodes, and the transient handle is dead**
  (N-9): any edit or `asImmutable()` through it afterwards throws, a programmer error, not an issue. One small
  allocation, and no write is ever lost to a handle that silently turned persistent. This departs from immutable.js,
  which freezes in place; `TransientTree` works with it, since it always reassigns from the result and calls
  `asMutable()` before its next write. `withMutations(fn)` follows the same rule.
- **Iterating a transient while editing it is undefined.**

### Allocation

- **Reads allocate nothing**: `get`, `has`, `size`, `forEach(fn, ctx)` (the context saves a closure, as `StateTree`
  passes one), and the cursors (`cursor(c)` with a reusable cursor object). `[Symbol.iterator]`, `keys()`, `values()`
  and `toArray()` exist and are documented as allocating. One exception: a bigint key's hash makes a string (N-11).
- **A persistent edit** copies the path: at most 7 small arrays for a map, far fewer in practice. An edit that changes
  nothing returns the same object: deleting an absent key, or setting a value equal to the stored one by the
  collection's `valueEquals`.
- **No global caches.** A collection's hash function is its own. A bounded string-hash memo is opt-in: a documented
  `hash` a config may pass, never the default (the author, 2026-11-18).

### Hashing and equality (R2-S14)

- **Per collection**, given at creation: `hash(key) => int32`, `equals(a, b)` for keys, `valueEquals(a, b)` for values.
- **Default key equality is SameValueZero**, as JS `Map`, `Set` and immutable.js treat keys: 0 and −0 are one key, and
  NaN equals NaN. A collection that needs them apart says so with `equals: Object.is`.
- **Default value equality** is `Object.is`.
- **`defaultHash`**, for every key type:

  | Key | Hash |
  |---|---|
  | string | core's `hashString`, then a 32-bit finalizer (MurmurHash3's `fmix32`), since `hashString`'s low bits are weak and CHAMP reads the low bits first |
  | number | an int32 by its value through `fmix32` (as immutable.js hashes small integers by value: no collisions, measured 7–10 ns faster); any other number its float64 bits mixed, both words through `fmix32`; −0 hashed as 0 and every NaN as one |
  | bigint | its decimal string's hash |
  | boolean, `null`, `undefined` | fixed constants |
  | symbol, object, function | an identity id from a `WeakMap`, assigned when the key is first inserted |

  A lookup of an object that has no id yet returns absent without assigning one, so reads allocate nothing. A symbol
  that cannot be a `WeakMap` key falls back to a `Map` id.
- **Collections are values**:
  - `equals(a, b)` is structural, short-circuited by identity.
  - Over the trie parts with the same config, it is a walk of the two trees. Otherwise, and for flat and collision
    nodes, it is size and lookup.
  - Collections with different configs are unequal (N-10). The same config means the same resolved `hash`, `equals`,
    `valueEquals` and `valueHash` functions, so `empty()` and `empty({})` match. Then `a.equals(b)` is `b.equals(a)`,
    and equal collections hash alike.
  - For the ordered collections, order is part of equality (as in immutable.js): an ordered walk.
  - `hashCode()` is order-independent for maps and sets, ordered for vectors and the ordered collections, and cached on
    the persistent collection object. It covers keys by `hash` and values by `valueHash` (default: `defaultHash`).
  - **The contracts**: a custom `equals` agrees with `hash` (equal keys hash alike), and a custom `valueEquals` with
    `valueHash`. A config that sets `valueEquals` without `valueHash` gets a `hashCode` over keys and size only.
- **This is not `repr`'s Merkle hash.** A trie's 32-bit key hash places keys; `hashNode`'s 128-bit content hash names
  revisions.

### Configuration

| Item | Default | Why |
|---|---|---|
| `hash` | `defaultHash` | as above |
| `equals` | SameValueZero | JS's own key equality |
| `valueEquals` | `Object.is` | no-op edits and `equals` |
| `valueHash` | `defaultHash` | `hashCode`, agreeing with `valueEquals` |

Fixed, not configurable:
- the trie's 5 bits a level;
- the flat threshold, a measured constant: about 8, re-measured in the lane against immutable.js.

Why: a configurable value costs a branch or a lookup on every access, for nothing measured.

### `update` and absent values (R2-C9)

- `set(k, undefined)` stores `undefined`.
- `update(k, fn)` stores what `fn` returns, `undefined` included.
- Removal is only ever `delete`.

### In `repr`

**Nodes keep their plain JSON shape.** They are the public contract the codecs, the builder, Seq\* and later meso\*
read. What the review found (R2-S15), and what follows:

- **The child index** (id to position) is a lazily built native `Map`, cached per children array in a `WeakMap`. It
  needs no persistent collection, since an insert shifts every later position. A path copy makes a new children array,
  so its index is rebuilt on its first lookup: O(fanout), what the edit's array copy already cost, and accepted (N-11).
  `resolveAddress` uses it above the flat threshold. It is built in `repr` by lane K, which already edits
  `tree/address.ts`, not by lane R (N-11).
- **The address index** (node to address) is dropped. One node object can sit at two addresses, a move changes every
  descendant's address, and `addressOf` has no caller outside tests.
- **The ref index** (`patches.md`, Applying 4), behind the internal `RefIndex` interface that answers by address prefix:
  - **Lane K** builds it on plain `Map`s, rebuilt per apply when a patch has a structural op. It is correct and simple.
  - **Measured**: a trie by id segment, its levels this library's `HashMap`s, derived from the parent revision's index
    by the changed nodes and keyed by root and language. It replaces the plain one when a measured case (a spec with
    10,000 bookmark refs, edited by moves) shows the rebuild matters.
  - That step comes after lanes R and K merge, behind the same interface, with nothing else in `repr` changing.
- **`diff` and `merge3`'s per-parent maps** are native `Map`s: built once, thrown away.
- **`ValueKeys`** (a preset keying by `repr` content) is dropped until a consumer needs it.
- **Nodes holding persistent collections** (draft 1's B) stay unbuilt. Their trigger is a measured spec in which copying
  a children array dominates an edit, and then a change to the node shape, with its own design.

## The API

The agreed surface for lane R:
- **The folder** is `packages/core/src/data/immutable/`, next to `data/int`, `data/generic` and `data/db`, with its
  `_test` folder's `null` export entry.
- **The entry** is `data/immutable.ts` (as `data/int.ts` is for `int/`), so `@molstar/core/data/immutable` through the
  package's wildcard export.
- **`Immutable`** sits beside `Int` and `Generic` in `data/index.ts`.
- One module a collection, plus `config` and `hash`.

```ts
interface CollectionConfig<K, V = unknown> {
  readonly hash?: (key: K) => number;
  readonly equals?: (a: K, b: K) => boolean;          // must agree with hash
  readonly valueEquals?: (a: V, b: V) => boolean;
  readonly valueHash?: (value: V) => number;          // must agree with valueEquals
}
interface MapCursor<K, V> { next(): boolean; readonly key: K; readonly value: V }   // reused; next() fills key and value
interface Cursor<T> { next(): boolean; readonly value: T }
function defaultHash(key: unknown): number;
function sameValueZero(a: unknown, b: unknown): boolean;

interface ReadonlyHashMap<K, V> {
  readonly size: number;
  get(key: K): V | undefined;
  has(key: K): boolean;
  forEach<C>(fn: (this: C, value: V, key: K) => void, ctx?: C): void;
  cursor(c?: MapCursor<K, V>): MapCursor<K, V>;           // reset and return c (or a new one)
  keys(): IterableIterator<K>;                            // allocates
  values(): IterableIterator<V>;                          // allocates
  toArray(): [K, V][];                                    // allocates
  equals(other: ReadonlyHashMap<K, V>): boolean;
  hashCode(): number;
  [Symbol.iterator](): Iterator<[K, V]>;                  // allocates
}
interface HashMap<K, V> extends ReadonlyHashMap<K, V> {
  set(key: K, value: V): HashMap<K, V>;
  delete(key: K): HashMap<K, V>;
  update(key: K, fn: (value: V | undefined) => V): HashMap<K, V>;
  asMutable(): MutableHashMap<K, V>;
  withMutations(fn: (m: MutableHashMap<K, V>) => void): HashMap<K, V>;
}
interface MutableHashMap<K, V> extends ReadonlyHashMap<K, V> {
  set(key: K, value: V): this;
  delete(key: K): this;
  update(key: K, fn: (value: V | undefined) => V): this;
  asMutable(): this;
  asImmutable(): HashMap<K, V>;                           // a new persistent map; this handle is dead
}
namespace HashMap {
  function empty<K, V>(config?: CollectionConfig<K, V>): HashMap<K, V>;
  function of<K, V>(entries: Iterable<readonly [K, V]>, config?: CollectionConfig<K, V>): HashMap<K, V>;
}
// HashSet<K>: add, delete, has, and the same reads, transients and statics.
// OrderedHashMap<K, V>, OrderedHashSet<K>: the same, with first(), last(), and iteration in insertion order.
// Vector<T>: get, set, push, pop, shift, size, first, last, indexOf, find, remove(index), forEach(fn, ctx), toArray,
//   cursor(c?: Cursor<T>), asMutable, asImmutable, withMutations, equals, hashCode; Vector.of, Vector.empty.
```

`repr`'s own surface does not change in this slice.

## Tests

- **Model-based property tests**: random sequences of operations on each collection, checked against a JS `Map`, `Set`
  or `Array` after each step. They include −0 and NaN keys, which the default equality treats as JS does. Seeded,
  short in the quick lane, long in the heavy lane.
- **Persistence**: every earlier version is unchanged after any later edit, transients included.
- **Transients**:
  - reads and iteration on a transient;
  - `asMutable()` returns itself;
  - `asImmutable()` returns a new persistent collection, after which the transient handle throws;
  - isolation from other versions;
  - the colliding-keys delete regression.
- **Collisions**: a deliberately bad hash (every key 0) gives correct results.
- **Equality**:
  - equal contents reached by different edit orders compare equal, and their `hashCode`s agree;
  - the trie parts have one shape;
  - a deletion demotes to flat at the threshold;
  - different configs are unequal; `empty()` and `empty({})` are the same config;
  - a structural `valueEquals` with its `valueHash` gives equal collections equal `hashCode`s.
- **`StateTree` spike**: `StateTree` and `TransientTree` type-check against the API (no code moved).
- **Allocation** (the opt-in tier): `get`, `has`, `forEach` with a context and cursor iteration allocate 0 bytes warm,
  and a lookup of an unseen object assigns no id. An edit stays within a stated bound.
- **Against immutable.js** (already a dev dependency here): get, set, delete, iterate and build at 10, 1,000 and
  100,000 entries, as heavy-lane ratios. The key sets include 22-character string ids, Mol\*'s refs. The target is no
  slower than immutable.js on collections of 1,000 or more, and within 2× below (the author, 2026-11-18: small maps pay
  for CHAMP's bitmaps, and Seq\*'s indexes are large); a miss is reported with its case. Integer keys at 1,000 miss (up to 1.9×: `fmix32` against immutable.js's identity hash) and are accepted, each case with its own bound (the author, 2026-11-18); a config wanting dense integer keys passes an identity `hash`.
- **In `repr`** (lane K): `resolveAddress` on a 1,000-child container, before and after the child index.
- **Bundle**: the entry's minified size, recorded.

## Not in this slice

- **Moving Mol\*'s state tree, `@molstar/model`'s custom properties and the plugin's uses off immutable.js**: a later,
  separate PR, with `ChildSet.map`'s change recorded as a v6 break then. The custom properties' `OrderedMap` used
  mutably is a native `Map` there.
- **The persistent ref index, the RRB vector, sorted maps, and nodes holding persistent collections**: each when its
  trigger is met.
