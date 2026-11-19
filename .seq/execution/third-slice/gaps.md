# Gaps: the third slice

Questions for the author, each yes or no with a recommendation. Gaps the loop decided as recommended while the author
was away are listed with "(decided as recommended)"; a public API change is never decided that way.

## Open

- **A lists-free palette module in `@molstar/core`** (TH's review, 2026-11-19). `util/color/value-palette.ts` imports
  `lists.ts` statically, so the lazy `ValuePaletteLists` saves work but not bundle size: every reader of value
  palettes ships all of Mol\*'s color lists (about 13.6 kB in Seq\*'s headless rendering and viewer). The themes design
  reserves the module split for you. Split the evaluation into a lists-free module that `value-palette.ts` re-exports
  (with `ValuePaletteLists`), non-breaking for Mol\*'s and MolViewSpec's code? Recommend **yes**.
- **The default protein scheme is not color-blind safe** (TH's drift 6): `classes` stays the default and its safe
  substitute is drawn under the reader's `theme.cvdSafe`. Keep `classes` as the default? Recommend **yes** (the
  reader's setting covers it; changing the default changes every figure).
