# Mol\* persistent collections: the specification

`@molstar/core/data/immutable` (`packages/core/src/data/immutable/`, the entry `data/immutable.ts`, `Immutable` in
`data/index.ts`): persistent hash maps, sets, their insertion-ordered kinds and a vector, with transients and cursors.
The rules the code relies on and why; code comments cite the sections (`README §3`).

| File                  | What                                                                              |
| --------------------- | --------------------------------------------------------------------------------- |
| `config.ts`           | `CollectionConfig`, `MapCursor`, `Cursor`, `sameValueZero`, resolved configs (§4) |
| `hash.ts`             | `defaultHash` and its lookup form, `memoHash`, `fmix32` (§4)                      |
| `champ.ts`            | the flat array and the CHAMP trie behind maps and sets (§2)                       |
| `hash-map.ts`         | `HashMap`, `MutableHashMap`                                                       |
| `hash-set.ts`         | `HashSet`, `MutableHashSet`                                                       |
| `ordered-hash-map.ts` | `OrderedHashMap`: a `HashMap` to slots and a `Vector` (§2)                        |
| `ordered-hash-set.ts` | `OrderedHashSet`                                                                  |
| `vector.ts`           | `Vector`: a 32-way trie with a tail (§2)                                          |
| `cursor.ts`           | the one reusable cursor (§5)                                                      |

## 1. Surface

The agreed API is the design's (`get`, `set`, `has`, `delete`, `update`, `forEach(fn, ctx)`, `cursor`, `keys`, `values`,
`toArray`, `equals`, `hashCode`, `asMutable`, `asImmutable`, `withMutations`, `of`, `empty`; sets `add`; ordered kinds
`first`, `last`; `Vector` `push`, `pop`, `shift`, `remove(index)`, `indexOf`, `find`). Nothing past it without a design
change: no lazy sequences, no `map`, no records. `Vector` has no iterator (its `toArray` and cursor cover it).

## 2. The collections

- **Maps and sets**: one flat array `[k0, v0, k1, v1, …]` in insertion order up to `FlatThreshold` (16) entries, scanned
  linearly with no hashing, as immutable.js's `ArrayMapNode`; a CHAMP trie above it, 5 bits a level. A set stores each
  key as its own value.
- **The threshold is one constant both ways**: an insert past it builds the trie; a delete down to it flattens. So a
  size has one representation. Measured (§8).
- **A CHAMP node** holds a data bitmap and a node bitmap, its entries from the front of one array in bit order and its
  sub-nodes from the back. Why one array: an edit copies one array a level.
- **The trie part is canonical** for equal contents: an entry sits at the shallowest level its hash prefix allows, and
  deletion compacts back (a sub-node left with one entry is inlined in its parent; a node left with nothing but one
  sub-node of one entry goes up whole).
- **Full collisions** are a collision node (an unordered bucket) at the first level the two keys meet, not a chain to
  the bottom. Its canonical place is kept on deletion: a node left with nothing but one collision node is replaced by
  it, up to the root, which keeps a node above it. The flat and collision parts are not canonical in order: they compare
  by lookup.
- **`OrderedHashMap`/`OrderedHashSet`**: a `HashMap` from key to slot, plus a `Vector` of keys (a set) or of keys and
  values interleaved (a map, slot `s` at `2s` and `2s + 1`: one vector, one path copy an update). A removal leaves a
  hole; trailing holes are trimmed, so `last()` reads the end. Compacted when holes pass half the slots and there are at
  least `CompactMin` (32) slots, as immutable.js's `OrderedMap`. Setting a present key keeps its slot. The vector holds
  keys and the `Hole` marker, so it uses the default config; the map's `valueEquals` is applied by the map.
- **`Vector`**: Clojure's trie of 32-way nodes with a tail of up to 32. `get`, `set`, `push`, `pop` are effectively
  O(1); `shift` and `remove(index)` move the later elements (O(n − index)) until the RRB vector. `set` outside
  `[0, size)` throws a `RangeError`; `get` outside returns `undefined`.

## 3. Transients

- A transient has an owner token; a node made under it is written in place, every other node is copied first: CHAMP
  nodes, collision nodes, flat roots and vector leaves and tails alike. So deleting colliding keys never writes into a
  node another version shares (the bug `TransientTree.renewOwner` works around in immutable.js, state `tree.md §2`).
- A transient reads like the persistent collection: every read, cursors and iteration.
- `asMutable()` on a transient returns itself. `asImmutable()` returns a new persistent collection sharing the
  transient's nodes, and the handle is dead: an edit, `asMutable()` or `asImmutable()` through it throws. Reads on a
  dead handle still answer. Why: no write is ever lost to a handle that silently turned persistent; one small object.
- A persistent edit that changes nothing returns the same object; so does a transient's (itself).
- Iterating a transient while editing it is undefined.

## 4. Hashing and equality

- **Per collection**: `hash`, `equals` for keys, `valueEquals` and `valueHash` for values, resolved at creation.
  Defaults: `defaultHash`, `sameValueZero` (0 and −0 one key, NaN equal to NaN), `Object.is`, `defaultHash`. A config
  with `valueEquals` and no `valueHash` has none: `hashCode` covers keys and size.
- **Two configs are the same** when the four resolved functions are; an empty config resolves to the default one, so
  `empty()` and `empty({})` are the same collection. Different configs are unequal.
- **`defaultHash`**: strings by core's `hashString` then MurmurHash3's `fmix32` (its low bits, which the trie reads
  first, are weak alone); an int32 by its value through `fmix32` (as immutable.js hashes small integers by value: no
  collisions, 7 to 10 ns faster than the bits), any other number by its float64 bits, both words through `fmix32`; −0 as
  0 and every NaN alike; bigints by their decimal string; constants for booleans, `null`, `undefined`; an identity id
  for objects, functions and symbols, from a `WeakMap` (a registered symbol, which cannot be a `WeakMap` key, from a
  `Map`).
- **A lookup assigns no id**: reads use the lookup form, which answers NaN for an object without an id, and the map
  answers absent. The hash and equality functions of an edit must not edit collections (an edit reports what it did
  through one module-level record).
- **A custom hash is coerced to int32** (`| 0`) wherever it places or finds a key: a fraction would recurse past the
  trie's 32 bits, NaN would file a key no lookup reaches.
- **`memoHash(capacity = 4096)`**, opt-in, never the default: `defaultHash`'s values with a bounded string memo of its
  own, at most `capacity` strings in two generations of half (the older dropped when the newer fills, an older hit moved
  up). One memo shared by the collections that must compare equal (§4: a config is its functions). A miss allocates; a
  lookup assigns no object id (the memo carries a lookup form, `LookupHash`). Why opt-in: a global cache is shared state
  with no owner; the memo pays only for strings hashed again and again.
- **`equals`** is structural and short-circuits on identity, on a cached `hashCode` mismatch, on size and config. Two
  tries compare by a walk of both (§2's canonical shape); flat and collision parts by lookup. Ordered kinds and vectors
  compare in order.
- **`hashCode`**: maps and sets a sum of their entries' hashes (order-independent), plus the size mixed; vectors and
  ordered kinds an ordered fold. Cached on persistent collections, computed afresh on transients.

## 5. Allocation

- Reads allocate nothing: `get`, `has`, `size`, `first`, `last`, `forEach(fn, ctx)` and a reused cursor (`cursor(c)`
  resets `c`; one cursor class serves every collection). A trie cursor walks on a fixed stack of 9 levels.
  `[Symbol.iterator]`, `keys()`, `values()` and `toArray()` allocate. A bigint key's hash makes a string.
- A persistent edit copies its path: a node object and its array a level (about 840 bytes for a set in a 1,000-entry
  map, §8). No `splice` on a copied array (its result is an allocation of its own).
- No global cache: a collection's hash is its config's; a memo is opt-in (§4).

## 6. Configuration

Fixed, not configurable: 5 bits a level, `FlatThreshold` (16) and `CompactMin` (32). Why: a configurable value is a
branch or a lookup on every access, for nothing measured.

## 7. `update` and absent values

`set(k, undefined)` stores `undefined`; `update(k, fn)` stores what `fn` returns, `undefined` included; only `delete`
removes.

## 8. Tests and measurements

- Model-based runs against JS `Map`, `Set` and `Array` (−0 and NaN keys among strings, numbers, objects, symbols and
  22-character ids), every earlier version checked at the end; seeded and short in the quick lane, long in the heavy
  lane. Shapes across edit orders, demotion at the threshold, all-zero and four-valued hashes, configs, transients.
- The allocation cases: `*.quick.ts`, `npx vitest run -c packages/core/vitest.perf.config.ts`, by hand; each read and
  walk of each kind its own case, measured every step (0 bytes, 30 cases).
- Against immutable.js (`_test/immutable-js.test.ts`, heavy): ratios for get, set, delete, iterate and build at 10,
  1,000 and 100,000 entries. The target (the author, 2026-11-18): no slower than immutable.js at 1,000 entries or more,
  within 2× below (small maps pay for CHAMP's bitmaps; Seq\*'s indexes are large). A case over the target prints as a
  miss; it fails past 1.3 (1,000 or more) or 2.5 (below), the target with a noisy nightly's headroom, or a known miss's
  own bound (its worst of three runs and a quarter). Three runs (2026-11-18): `Vector` 0.1 to 0.9 on every case; maps
  and sets of string ids 0.35 to 0.93 at 1,000 and more. The misses: integer keys at 1,000 (a lookup 1.5, a build 1.9,
  the ordered map 1.1 to 1.3: `fmix32` scatters what immutable.js's identity hash packs densely; an identity hash reads
  as fast as immutable.js), and inserts across the flat threshold (3.4 for string ids, 2.2 ordered: every key hashed
  then; `memoHash` is the opt-in remedy).
- The flat threshold (same file): a lookup and a persistent set on a flat array and on a trie of 1 to 32 entries. The
  flat array won both at every size to 16 (string ids: a lookup 14 to 20 ns flat, 52 to 60 ns as a trie; integers: 17 to
  24 against 33 to 38), integers crossing at about 25 and string ids past 32; so 16.
- The bundle: the entry minified, 28.2 kB (7.9 kB gzip; immutable.js 5.1.9's `immutable.min.js` 68.5 kB, 19.1 kB gzip),
  2026-11-18:
  `npx esbuild packages/core/src/data/immutable.ts --bundle --minify --format=esm --outfile=/tmp/immutable.min.js`.

## 9. The `StateTree` spike

`_test/state-tree-spike.ts` type-checks `StateTree`'s and `TransientTree`'s collection code against this API without
changing the state tree. What does not fit, each a `@ts-expect-error` or a note there:

- the interfaces: `ChildSet.map`; `forEach` returning immutable.js's count; `remove` (here `delete`); `_Map.get`
  returning `T` (here `T | undefined`); `asImmutable()` on a persistent map;
- behaviour: immutable.js's `forEach` passes the collection third and stops on `false`, ours neither; a transient has no
  `withMutations` (immutable.js keeps a mutable collection mutable through it; ours would freeze the handle);
- code: `dump` reads `keySeq()`/`valueSeq()` (here `keys()`/`values()`); `transient.ts` casts to `OrderedSet`.

Restated on this API, the rest fits; `renewOwner` is not needed (§3), and the tree reassigns from each `asImmutable()`
as it already does. Type-checked with the tests (`tsc -b` leaves `_test/` out).
