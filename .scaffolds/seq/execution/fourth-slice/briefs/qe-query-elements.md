# Lane QE: the query language by role, property access, `elements` and `groups`

**Task**: build the language half of [locations-and-queries](../../../design/locations-and-queries.md) (Decision: Three
result types, The language: forms by role), [properties](../../../design/properties.md) (all of it), the schema
selections of [addressing](../../../design/addressing.md) (Schema selections; the grammar's positional values) and the
`query-form` and `query-property` kinds of [extension-openings](../../../design/extension-openings.md) and
[viewer-composition](../../../design/viewer-composition.md); the queries of [molecules](../../../design/molecules.md)
follow. The agreed draft, with the examples, is `drafts/elements-form.md` beside this folder.

- **Forms by role**: generators `elements` (a Location; its items one record, columns, or `:items [(fields …) …]`,
  normalized to typed columns, printed in the shorter layout; the finest property's level decides residues or atoms;
  the index walk, rows plus matches) and `groups` (Keys; `:molecule-test`, `:unit-test`, `:residue-test`,
  `:atom-test`, `:by`); `filter :test`; modifiers `expand :to` and `linked :via :type`; combinators `merge`,
  `intersect`, `except` (renamed from `or`, `and`, `not`: every use in `packages/seq`, the viewer's `dataFilter` and
  the stories included; `except` with one operand the complement within the prefix); operators `=`, `!=`, `<`, `<=`,
  `>`, `>=`, `in-range`, `in`, `and`, `or`, `not`, `prop`, `value`, with positional operands, accepted by expression
  forms only (the address grammar reads positional values before keyword params; any other form with one is
  `seq.query.form.param`). Types checked when the query is read (`seq.query.form.type`). Every form monadic over a
  common prefix: nodes, a Location, Keys.
- **Property access**: the built-ins at four levels with their short (printed) and long (MolViewSpec) names; the
  `QueryModel.property` hook (`column | 'pending' | undefined`); `QueryUnit.entity` (replacing the key-of-`of` lookup),
  `QueryAtoms.element` and `residue`; values channels from the store, pending staged; `seq.query.property.level`,
  `seq.query.property.absent` (a warning, the element yields nothing), `seq.query.extension.missing`.
- **Registry**: `query-form` (with its role) and `query-property` entries in `SeqRegistryEntry`, namespaced and checked
  at composition with the descriptor; `queryResolvers` composed from the built-ins and the registry.
- **Rename and remove**: the result type *Elements* to `Keys` (`ElementKey` to `Key`, `elements.ts` to `keys.ts`); ML2
  built the forms this replaces in `packages/seq/core/src/query/molecule-forms.ts` (`atom`, `entity`, `component`):
  remove them (0.5.0 is unreleased: no migration); `chain`, `molecule` (gaining `:entity`) and `atoms` stay.

Stay on the language and the query model: no drawing, no mmCIF reading (MC writes `b-factor` and `occupancy` channels;
until then they are absent), no structure-viewer wiring (ST2).

**Output**: the forms, properties and registry kinds with quick-lane tests:
- the examples on 4HHB (`(elements :component HEM :atom FE)`: four irons; `(elements :auth-chain [A C] :auth-seq 87)`),
  1FIN (`(elements :chain A :beg 10 :end 20 :atom CA)`: eleven atoms; a `filter` over a property), 3SGJ's glycan
  (`(elements :atom C1)` over its residues: one atom a residue; `linked` from the Asn to the glycan; record a small
  fixture of its units and links if none exists), a modified residue (1EHZ: `(groups :residue-test (= (prop modified)
  true))`, eleven keys; an MSE's selenium);
- properties: monadic composition (a form over a union of prefixes equals the union of the form over each; an element
  lacking a property yields nothing and the rest stays); `elements` equal to its `filter` sugar for random items;
  SoA ≡ AoS ≡ records after normalization; aliases round trip (short and long parse equal, print short;
  `parse(print(a))` equal as a value); combinators over one type;
- an extension's property (by `{ values, field }` and by a function) and its filter: registered, missing, namespace
  refused; a pending values channel staged;
- `_spec/README.md` "Queries" (Language, Properties, Keys, Extensions) and "Context" (registry entries); a quick
  representative of the index walk (rows × units not scanned).

**Model**: strong (opus).

**Where**: A worktree on a branch `seq/qe-lane` from the base the orchestrator names (after ML2 merges). Hot files
(plan.md) are the orchestrator's: list the changes you need in the report. Every modified file's header updated; a
quick-lane test beside the code, failing on the parent for a fix.

**Report**: Under 300 words: commits; checks run and their result; tests shown failing on the parent; the public surface
added and removed (names, for the author's acceptance); drift notes; the `.v6/changelogs/` entry's text.
