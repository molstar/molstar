## Core: persistent collections

- `@molstar/core/data/immutable`: persistent `HashMap`, `HashSet` (CHAMP, flat below a threshold), `OrderedHashMap`,
  `OrderedHashSet` and `Vector`, with transients (`asMutable`/`asImmutable`, a dead handle after `asImmutable`), a hash
  and equality per collection (SameValueZero keys and `defaultHash` by default) and structural `equals`/`hashCode`.
  `memoHash` is an opt-in bounded string-hash memo. The entry is 28.2 kB minified (7.9 kB gzip). Mol\*'s state tree
  still uses immutable.js.
