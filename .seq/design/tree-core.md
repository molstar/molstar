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
  - `repr/patch`: operations by address, apply with inverses, refs rewritten on move, `diff`.
  - `repr/builder`: the typed builder's engine.
  - `repr/codec`: the codec interface (a tree to bytes or text and back, with the value model's guarantees) and its
    conformance suite; no codec in it.
- **Codecs live where their libraries do**:
  - JSON and the S-expression reader and printer: `@molstar/core/repr/codecs/json` and `…/sexpr` (no dependencies).
  - MessagePack: in `@molstar/io`, on its existing `encodeMsgPack` and `decodeMsgPack` and their extension hooks (no
    change to io's codec), as a repr codec with typed columns `0x11`–`0x18`; `@molstar/seqx-core` imports it through
    its allowlist.
- **The first move is minimal** (the review's first slice): `tree`, `schema`, `builder`, the codec interface and the
  JSON and S-expression codecs, enough for the first Seq\* slice; `patch` and `descriptors` follow when the slice needs
  them. Definitions and uses (`def`, arguments, overrides) wait ([composition.md](composition.md)).
- **What stays in meso\* for now**: `locators`, `context`, `units`; the schema's tie to units becomes a hook (the
  author: units "ok let's wait").
- **Issues**: one shape, `{ code, severity, path, message, hint? }` with a stable namespaced code; `path` the step path
  meso\* uses, with the address beside it where the place has one; severity means one thing everywhere (error: invalid
  input; warning: something is lost; info: never changes what is drawn). meso\*'s issue gains the severity in the move.
- **History**: `repr/patch` gives inverses; the undo stack is `@molstar/ui/history`'s (`HistoryPort`), not a second one.
- **Headers**: the standard Mol\* header with the year and `@author` lines.
- **Spec and tests**: a `_spec` beside the module with what is built; tests within the fast budgets
  ([testing.md](testing.md)).
- **The entry's name is `@molstar/core/repr`** (the author: "ok"). The review noted `@molstar/graphics/repr` means
  visual representations; the name stands, and the module's README says which it is.

## Consequences

- `@molstar/seqx-core` describes its one language on `repr` and picks its codecs (JSON by default, MessagePack from io).
- A public API is added to `@molstar/core`; the change goes into the v6 changelog entry of its branch.
