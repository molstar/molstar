## Graphics: labels

- Labels are anti-aliased: glyph and outline edges are smoothed over one pixel on WebGL and WebGPU (they were cut
  hard at the distance field's middle: stair-stepped edges, broken strokes in small text) (#2030)
- Labels are kerned and spaced by the font's measured advances (they advanced by each glyph's atlas cell, rounded up:
  `AV`, `To`, `LT` had wide gaps); glyphs reaching left of the pen or past their advance (`j`, `f`, italics) are no
  longer clipped (up to 1.5 lineheights of ink); emoji, accented letters with combining marks and astral characters are one glyph each (#2030)
- `FontFamilyFonts` maps each label font family (`sans-serif`, `monospace`, `serif`, `cursive`) to the fonts it draws
  with (Prism's theme stacks by default); `setFontFamilyFonts` points one at other fonts, web fonts included (CSS
  `@font-face`, Google Fonts, a `FontFace`), which download only when a label needs them: text visuals wait for them
  (`loadFontAtlasFonts`, `withFontAtlasFonts`; `clearFontAtlasCache` drops cached atlases) (#2030)
- `TextMeasure` and `splitGraphemes` in `@molstar/core/util/text-measure`: graphemes and kerned pen positions of a text
  from a measure function, adapted from pretext (MIT) (#2030)
- Labels, measurements and MVS labels share one GPU texture per font atlas instead of one per text object (WebGL and
  WebGPU), and upload it only when glyphs were added: 10 distance measurements 983.7 → 613.8 MiB of textures, an MVS
  scene with labels on 6 structures 819.3 → 613.8 MiB (41 MiB per text object saved). `TextureSpecOptions.shared` opts a
  texture of a renderable schema in; `TextureImage.version` marks in-place writes (#2030)
