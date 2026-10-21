# Nodes

Node selectors: CSS-shaped text naming cells of the plugin's state tree, compiled to Mol\*'s `StateSelection.Builder`.
The code is the truth for how (`src/api/nodes.ts`, `features/state-tree/operations.ts`).

## 1. Selectors

- **A query, not a name**: an operation acts on every match
  (`structure[label="1CBS"] component[type=ligand] > representation`, `#ref`, `representation:hidden, volume`). A bare
  ref works too.
- **Kinds as the Structure panel reads cells** (`NodeKinds`, `nodeKind`): a particle list is `particles`, a shape
  provider (a mesh file, MVS primitives) `shape`, their representations `representation`; a structure inside a structure
  is a component; decorators have no kind and fold into what they decorate (`>` matches through them). Ghost cells never
  match; the root only by its ref.
- Attributes, operators, states: `NodeAttributes`, the parser. Text compares without case, except refs.

## 2. In operations

- **`node: NodeParam()`**: a many-cell operation acts on every match of its kind as one step (`nodesFor`), skipping
  other kinds; a one-cell one (`state.get`, `state.action`) fails `ambiguous` listing matches (`selectOneNode`).
- **Errors**: nothing of the kind is `not-found` saying what matched; a plain word that is no kind nor ref is
  `not-found`; other unparsable text is `invalid`.
- `state.select` lists matches (read-only, to try a selector); `state.apply` applies a transformer under every match as
  one step, `tags` for later selectors, never on the root (`blocked`: loading is the loaders').

## 4. In scripts

`node(kind, attrs)` (`@molstar/prism-viewer/api`) writes the same text, chainable; its JSON is the text, so the history
records a selector that replays.

## 5. Limits

No `:has()`, no positions; the structure-scoped operations (`representation.set`, `.color`, `component.*`) take a
component and structure, not `node`.
