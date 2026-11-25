# meso\*

meso\* (the molstar-ludo repository: `.meso/`, `packages/meso/*`, `apps/meso`) is a hierarchy of declarative languages
for mesoscale biology, with its app a projection of one spec. It is built on the same principles package, and it is
where the lessons on building such languages were learned. Its `.meso/` is the precedent for `.scaffolds/seq/`.

## What we take

- **One node shape**: `kind`, a local `id`, params, children. Ids are unique among siblings; an address is the path of
  ids from the root (`/data/p53`).
- **Languages as data**: a language declares its kinds (allowed parents, id policy), typed params with defaults,
  docs and bounds, refs and slots for embedded languages. Validation, defaults, the builder's types, the docs and
  (here, later, a JSON Schema) are all derived from that one description.
- **Validation reports every issue**, each `{code, message, path, hint?}` with a stable namespaced code.
- **A typed builder** derived from the language: a kind's method exists only where the kind may sit; misuse fails
  where it is written; the builder produces the tree the reader would.
- **Patches**: `set`, `unset`, `insert`, `remove`, `move`, each naming one node by address; a batch is all or nothing
  and one undo step; inverses are computed when applied; the spec changes only by patches (UI, commands and agents
  alike).
- **Pure passes**: spec → resolved model → render form, each a pure function cached by hash; computed data never goes
  back into the spec (§derived-data-outside; pure passes in its representations module).
- **Value or reference**: content is written inline or referenced and resolved by a port.
- **Versioning**: the language has its own version; migrations are explicit tree-to-tree functions.
- **Descriptors, not registration on import**: a caller gathers what it needs.
- **Storybook as documentation**: pages generated from the language description, recipes that run, tests that keep
  them true.
- **Its design note on sequences** (`.meso/design/meso/sequences.md`): a position is a locator in a declared space,
  spaces meet through declared mappings composed by a router, the sequence view is a projection.

## What we leave

- **S-expressions as the canonical spelling and `.sexpr` files**: here JSON is the default encoding, MessagePack an
  option, S-expressions a printing (the author's decisions, [author](../author/README.md)).
- Its time, animation, illustration and packing languages, the Prism adapter, the converters.
- Its own id checker (`.meso/tooling`): here the package's checks run as they are.

## Reconciliation

- Its principles copy is richer than the principles package's (queries, agent-interface, a larger representations module).
  queries and agent-interface are taken as this instance's own modules; citations of its representations ids the principles
  package lacks are reworded ("not adopted here"). If the language work needs one of them, it is brought in through
  the self-check.
