## Graphics: labels

- Labels are anti-aliased: glyph and outline edges are smoothed over one pixel on WebGL and WebGPU (they were cut
  hard at the distance field's middle: stair-stepped edges, broken strokes in small text)
- Labels are kerned and spaced by the font's measured advances (they advanced by each glyph's atlas cell, rounded up:
  `AV`, `To`, `LT` had wide gaps); glyphs reaching left of the pen or past their advance (`j`, `f`, italics) are no
  longer clipped; emoji, accented letters with combining marks and astral characters are one glyph each
- `FontFamilyFonts` maps each label font family (`sans-serif`, `monospace`, `serif`, `cursive`) to the fonts it draws
  with (Prism's theme stacks by default); `setFontFamilyFonts` points one at other fonts, web fonts included (CSS
  `@font-face`, Google Fonts, a `FontFace`), which download only when a label needs them: text visuals wait for them
  (`loadFontAtlasFonts`, `withFontAtlasFonts`)
- `TextMeasure` and `splitGraphemes` in `@molstar/core/util/text-measure`: graphemes and kerned pen positions of a text
  from a measure function, adapted from pretext (MIT)
