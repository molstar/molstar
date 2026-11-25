## Mol* core and Seq* (next): colour lists and the gaps decided

- core: `util/color/value-palette-eval`, the value-palette evaluation without Mol*'s colour lists; `value-palette`
  re-exports it with `ValuePaletteLists` (`ValuePalettePropsOptions.lists` is new). Each colour list is a named export of
  `util/color/lists` (`ColorListViridis`, …); `ColorLists` is unchanged and tree-shakable, so importing one list ships
  only that list.
- Seq* (next): `classes-safe` (residue classes on Okabe–Ito hues), the colour-blind substitute of every unsafe protein
  scheme; in dark, letter cells draw a hairline edge (`seq-cell-edge`); one motion item, `view.motion: platform | cut |
  animate`, for jumps and stories (`story.reducedMotion` removed; a legacy `auto` reads as `platform`); clades named by
  unique node ids that are not support values, else their first and last leaves, and a `collapsed` name no clade has
  reported; Seq* ships only the eight colour lists it names (about −13 kB).
