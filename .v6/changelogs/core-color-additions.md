## Core: color additions (non-breaking)

- Oklab conversions and interpolation (`Color.interpolateOklab`), color-vision-deficiency simulation (Machado 2009) and
  a separation check, residue color schemes as data (`util/color/residue-schemes`).
- The hydrophobicity scale moves to `util/color/residue-schemes`; `@molstar/model` re-exports it under its name.
- Value palettes gain an optional `scale` (`linear`, `log`, `sqrt`) with `mid` and `non_positive`; categorical palettes
  gain `unlisted: 'hash'` (FNV-1a over UTF-8); `ValuePaletteLists` is lazy. Defaults are unchanged.
