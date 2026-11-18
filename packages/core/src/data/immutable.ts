/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Persistent collections with transients and cursors (immutable/_spec/README.md).

export { type CollectionConfig, type Cursor, type MapCursor, sameValueZero } from './immutable/config.js';
export { defaultHash, type MemoHash, memoHash } from './immutable/hash.js';
export { HashMap, type MutableHashMap, type ReadonlyHashMap } from './immutable/hash-map.js';
export { HashSet, type MutableHashSet, type ReadonlyHashSet, type SetConfig } from './immutable/hash-set.js';
export {
  type MutableOrderedHashMap,
  OrderedHashMap,
  type ReadonlyOrderedHashMap,
} from './immutable/ordered-hash-map.js';
export {
  type MutableOrderedHashSet,
  OrderedHashSet,
  type ReadonlyOrderedHashSet,
} from './immutable/ordered-hash-set.js';
export { type MutableVector, type ReadonlyVector, Vector, type VectorConfig } from './immutable/vector.js';
