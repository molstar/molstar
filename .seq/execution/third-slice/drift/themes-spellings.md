# Drift: themes' spellings (TH)

What the code spells differently from [themes](../../../../design/themes.md), and a recommendation each:

- **Clause kinds are words** (`eq`, `lt`, `one-of`, …), not symbols: a kind must be a word in the language's
  description. Recommend: keep the words; the design's examples follow at the close.
- **`non_positive` with an underscore**, as Mol\*'s palette spec spells it (MolViewSpec parity). Recommend: keep.
- **A data column named `kind`** collides with a form's `kind` key in the builder, and the S-expression codec refuses
  it as an inline column (reserved). The secondary layer reads `kind`, else `category`. Recommend: the builder's docs
  say so; a `kind` column in a source table still works.

**Decided** (the author, 2026-11-19: "ok to have words. we could support optional aliases with symbols for = etc if
needed later?"): clause kinds stay words; the design's examples follow at finalize; optional symbol aliases (`=`,
`<`, `<=`, …) read as the words are a later addition if a user needs them (followups). The other two items stand as
recommended unless the author says otherwise.
