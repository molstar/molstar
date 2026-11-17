# The tree core in @molstar/core

## Intent

The author, 2026-11-17, on bringing meso\*'s tree core (the tree, schema-as-data, patches, addresses, the S-expression
reader and printer, the builder engine) here: "yes i think we can bring this to packages/core here. standard headers,
since both are mol*".

Seq\* and meso\* need the same foundation: one node shape, languages described as data, validation with every issue,
addresses with query segments, patches with inverses, codecs and hashes, a typed builder. meso\* built and hardened it
(`@molstar/meso-core` in the molstar-ludo repository, about 7,000 lines of source and 5,700 of tests). It moves into
`@molstar/core`, generic and documented, so Seq\* builds on it and meso\* can later use it from here instead of its own
copy.

## Decision

- **What moves**, from meso-core's folders: `tree` (the node, JSON and S-expression codecs, Merkle hashes, addresses
  with query segments, issues), `schema` (languages as data: kinds, params, unions, refs, definitions and uses,
  validation, docs), `patch` (operations, apply with inverses, history, refs rewritten on move, both spellings),
  `builder` (the typed builder's engine) and `descriptors` (languages, IRs and runtimes described as data, gathered by
  the caller).
- **What stays in meso\* for now**: `locators`, `context` and `units` (Settled, below).
- **Headers**: the standard Mol\* header with the year and `@author` lines (both are Mol\*'s; no provenance line
  needed beyond the authors, who are the same).
- **MessagePack from `@molstar/io`**: its encoder and decoder, extended with typed columns as extension types
  `0x11`–`0x18` (Seq\*'s default encoding, [format.md](format.md), [base-library.md](base-library.md)); JSON and MessagePack as the stored encodings, S-expressions as the printing, per language
  user's choice (meso\* keeps S-expressions as its spelling, Seq\* stores MessagePack or JSON).
- **Spec and tests**: a `_spec` beside the module, from meso-core's (`schema.md`, `patch.md`, `builder.md` and its
  README), cut to rules and reasons; the tests brought over within the fast-test budgets ([testing.md](testing.md)):
  property tests at a small fixed number of seeded cases, the scaling tests moved to the heavy lane.
- **Its own work**: a work id of its own (`tree-core`), done before the Seq\* language work that builds on it; meso\*
  is not changed by it (it lives in another repository).

## Settled

- **The entry is `@molstar/core/repr`** (the author, 2026-11-17: "ok" to the recommendation): `repr/tree`,
  `repr/schema`, `repr/patch`, `repr/builder`, `repr/descriptors`, after the principles' representations module;
  `state` is Mol\*'s state tree, and `spec` would read as a `_spec` folder.
- **Units wait** (the author: "ok let's wait"): quantities with units stay in meso\* until something here needs them;
  the schema's tie to them becomes a hook, so the move does not bring them along.

## Consequences

- `@molstar/seqx-core` imports the tree core from `@molstar/core`; its own language levels are descriptions on it.
- A public API is added to `@molstar/core`; the change goes into the v6 changelog entry of its branch.
