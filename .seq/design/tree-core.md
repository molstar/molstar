# The tree core in @molstar/core

## Intent

The author, 2026-11-17, on bringing meso\*'s tree core here: "yes i think we can bring this to packages/core here.
standard headers, since both are mol*"; then, after the review found MessagePack could not live in core: "make the
repr "codec indepedent" and about abstract patterns then?"

Seq\* and meso\* need the same foundation: one node shape, values with identity, languages described as data,
validation with every issue, addresses with query segments, patches with inverses, a typed builder. meso\* built it
(`@molstar/meso-core` in the molstar-ludo repository). Its abstract part moves into `@molstar/core/repr`; codecs stay at
the edge, beside the libraries they need.

## Decision

- **`@molstar/core/repr` is codec-independent**: the abstract patterns only.
  - `repr/tree`: the node, values and their identity (typed columns by element type and bytes, hashed once and cached on
    the node), Merkle hashes, addresses with query segments, issues.
  - `repr/schema`: languages as data (kinds, params with types including column types, unions, refs), validation with
    every issue, defaults, docs.
  - `repr/patch`: git's model ([patches.md](patches.md)): revisions are immutable roots named by their hash, commits a
    root and its parents; six operations by address of ids, applied all or nothing with inverses, normalized and
    validated for what they change, refs following a move, link integrity; `diff`; `merge3` (its own agreement).
  - `repr/builder`: the typed builder's engine.
  - `repr/codec`: the codec interface (a tree to bytes or text and back, with the value model's guarantees) and its
    conformance suite; no codec in it.
- **Codecs live where their libraries do**:
  - JSON and the S-expression reader and printer: `@molstar/core/repr/codecs/json` and `…/sexpr` (no dependencies).
  - MessagePack: in `@molstar/core/repr/codecs/msgpack`, on `@molstar/core/data/msgpack` (moved down from io, which
    re-exports both at their old paths), as a repr codec with typed columns `0x11`–`0x18`; `@molstar/seqx-core` imports
    it from core, with no io dependency (amended 2026-11-23: the code's way, the author as recommended).
- **The first move is minimal** (the review's first slice): `tree`, `schema`, `builder`, the codec interface and the
  JSON and S-expression codecs, enough for the first Seq\* slice; `patch` and `descriptors` follow when the slice needs
  them. Definitions and uses (`def`, arguments, overrides) wait ([composition.md](composition.md)).
- **What stays in meso\* for now**: `locators`, `context`, `units`; the schema's tie to units becomes a hook (the
  author: units "ok let's wait").
- **Issues**: one shape, `{ code, severity, path, message, hint? }` with a stable namespaced code; `path` the step path
  meso\* uses, with the address beside it where the place has one; severity means one thing everywhere (error: invalid
  input; warning: something is lost; info: never changes what is drawn). meso\*'s issue gains the severity in the move.
- **History**: `repr/patch` gives inverses; the core keeps no undo stack. The stack is the host application's, reached
  through the consumer's own port (Seq\*'s history port, [changes.md](changes.md)); `@molstar/ui/history` is one adapter.
- **Headers**: the standard Mol\* header with the year and `@author` lines.
- **Spec and tests**: a `_spec` beside the module with what is built; tests within the fast budgets
  ([testing.md](testing.md)).
- **The entry's name is `@molstar/core/repr`** (the author: "ok"). The review noted `@molstar/graphics/repr` means
  visual representations; the name stands, and the module's README says which it is.

## Consequences

- `@molstar/seqx-core` describes its one language on `repr` and picks its codecs (JSON by default, MessagePack from core's `repr/codecs/msgpack`).
- A public API is added to `@molstar/core`; the change goes into the v6 changelog entry of its branch.

## In the specs

| Decision | Where |
|---|---|
| `@molstar/core/repr` codec-independent; its folders and import direction | repr › Imports go one way |
| JSON and S-expression codecs in core | repr › Codecs |
| MessagePack in core's `repr/codecs/msgpack`, on `data/msgpack` (io re-exports) | repr › Codecs |
| The first move minimal; `patch` later | design-only: history (the slice's plan) |
| `locators`, `context`, `units` stay in meso\*; units a hook | repr › Units are a hook |
| Issues: one shape, namespaced codes, one meaning of severity | repr › Issues |
| No undo stack in the core; inverses only | repr › A revision |
| The history port and its adapters | viewer › History |
| Headers | design-only: the repository's rule (root `AGENTS.md`) |
| A `_spec` and tests within the quick budgets | repr › Tests |
| The entry's name, and what `repr` means here | repr › repr: the tree core |
