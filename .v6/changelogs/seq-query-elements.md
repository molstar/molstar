## Seq* (next): the query language

- [Breaking] Queries by role: the generators `elements` (rows of typed fields, as records or as columns) and `groups`;
  a generic `filter :test` over typed properties with positional predicates (`(> (prop b-factor) 60)`), `and`, `or`,
  `not`; the modifiers `expand` and `linked`; the combinators `merge`, `intersect`, `except` (were `or`, `and`, `not`);
  `atom`, `entity` and `component` removed. The result type `Elements` is `Keys`.
- Built-in properties by level with short and long (mmCIF) names, the short printed; `b-factor` and `occupancy` as
  residue-level values channels; an atom field on residues with no atoms held reads their standard component's atoms.
- Extensions add `query-form` and `query-property` registry entries; a pending property stages the result as pending.
