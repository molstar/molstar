# Nodes: addressing the state tree

An operation that acts on cells of the plugin's state tree (a representation, a component, a volume's isosurface, any
cell) says which with a **node selector**: text in the shape of a CSS selector, read the same way by the API, the
commands, the history and the assistant. It compiles to Mol\*'s `StateSelection.Builder`, the monadic selection the
state already has (`state.select(...)`), so a selector is a query over the tree, not a name for one thing: an operation
that takes a `node` acts on every cell it matches. Status: **accepted** (2026-10-16). Code:
`packages/prism/viewer/src/api/nodes.ts`.

Before, each area had its own address: a ref (`state.*`, `representation.update`, `.visibility`), a component and a type
(`representation.set`, `.remove`), a component (`material`), an owner and a ref or a type (`volume.*`, `particles.*`).
Now the raw-ref fields are `node`; the friendly scopes stay beside it (a component by label, a volume by name, a type).

## 1. Selectors

```text
structure[label="1CBS"] component[type=ligand] > representation[type=ball-and-stick]
#Ag3k1mXq…              a cell by its ref (a ref alone, without #, works too)
representation:hidden, volume representation
```

- **Kinds**: `data`, `trajectory`, `model`, `structure`, `component`, `representation`, `volume`, and `*` for any. They
  are read as the Structure panel's hierarchy reads cells (`nodeKind`): a structure inside a structure is a component;
  a representation is any 3D representation (of a structure, a volume, particles, a shape). Decorators (a structure's
  transform, custom properties) have no kind: they fold into the cell they decorate, as the State tree tab folds them.
- **Attributes** `[name op value]`: `label` (the cell's label), `type` (a representation's type; a component's static
  kind: `ligand`, `water`, `polymer`), `ref`, `tag` (`.tag` too), `transformer` (its id). Operators `=`, `*=`
  (contains), `^=` (starts), `$=` (ends). Text compares without case, except refs. A value with spaces is quoted.
- **States**: `:hidden`, `:visible`, `:error`, `:not(compound)`.
- **Combinators**: a space for anything below, `>` for a child (a decorator between does not count:
  `structure > component` matches through a transform). A comma for a list; each cell once, in tree order.
- **Ghosts and the root**: ghost cells (the hierarchy's internal ones) never match; the root matches only by its ref.

## 2. In operations

- `node: NodeParam()` is a text param. One that acts on several (`visibility`, `update`, `remove`, `state.apply`,
  `volume.isovalue`) acts on every match of its kind as one step; one that acts on one (`state.get`, `state.action`)
  fails `ambiguous` listing the matches (label and ref), so the caller narrows it or picks by `#ref`.
- A match of another kind is skipped: `representation.visibility({ node: 'component[type=ligand] *' })` hides the
  ligand's representations. Nothing of the kind matched fails `not-found`, saying what did match.
- A plain word that is neither a kind nor a ref is `not-found` ("There is no node 'x'"); other text that does not parse
  is `invalid`, with where and what was expected.
- The areas: `state.get`, `.select`, `.visibility`, `.update`, `.apply`, `.action`, `.remove`; `representation.update`,
  `.visibility`, `.remove`; `volume.isovalue`, `.update`, `.visibility`, `.remove`; `particles.update`, `.visibility`,
  `.remove`. The scene's rows pass `#ref`-like refs, so a row acts on its one cell.

## 3. `state.select` and `state.apply`

- `state.select({ node })` lists what a selector matches (`{ ref, label, kind, transformer, parent, visible, status }`):
  read only, the way to try a selector before acting on it.
- `state.apply({ node, transformer, params, tags })` applies a transformer under every match (Mol\*'s
  `build().to(ref).apply(t, params)`), one history step; the new refs are the value, `tags` let a later selector find
  them (`.mine`). The transformer is an id or its name without the namespace (`structure-representation-3d`); its params
  are checked against the transformer's definition for each cell, defaults filling the rest. A transformer that does
  not apply to a match fails `invalid`, naming those that do; the root is the loaders' (no downloads here).

## 4. In scripts

`node(kind, attrs)` writes the same text, chainable, so a script never builds selector strings by hand:

```ts
import { node } from '@molstar/prism-viewer/api';

await viewer.api.representation.visibility({
  node: node('structure', { label: '1CBS' }).find('component', { type: 'ligand' }).child('representation'),
  visible: false,
});
node.ref(ref); // '#ref'
node('structure').compile().subtree(); // Mol*'s StateSelection.Builder, to keep composing
```

An operation given a `node(...)` takes its text (`checkValues`), and its JSON is the text, so the history and the
session record show the selector, which replays.

## 5. Not yet

- `:has(...)`, positions (`:first`, `:nth(n)`), and a kind for shapes and particle lists.
- The structure-scoped operations (`representation.set`, `.color`, `component.*`) keep their component and structure
  fields; a `node` there (several components at once) waits for a use.
