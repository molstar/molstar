# Labels: fonts, layout and rendering

Text geometry (`packages/graphics/src/geo/geometry/text/`): structure labels, measurements, MVS labels and primitives,
the camera helper's axes. A label is one quad per glyph from a signed-distance-field (SDF) font atlas, drawn by
`text.frag` / `text.wgsl`.

## 1. One atlas per font

- `getFontAtlas(props)` keys by the CSS font the atlas draws with (`getFontAtlasFont`: style, variant, weight, size from
  `fontQuality`, the family's fonts) and gives the atlas only those five props. Why: keyed by every Text prop, each
  label size step or MVS label set (its `items` are props) kept a new 43 MB atlas (perf-audit-2 08/B6).
- The cache is an LRU of 4. An evicted atlas lives on in the texts built with it (each `Text` holds its texture);
  texts built later get a new one. More than 4 fonts in use at once rebuild atlases when all their texts rebuild: time,
  not memory.

## 2. Fonts

- `fontFamily` stays a `PD.Select` of the generic families (`sans-serif`, `monospace`, `serif`, `cursive`): states, MVS
  and the UIs are unchanged. What a family draws with is `FontFamilyFonts[family]`, a CSS font-family list, generic
  last. Defaults are Prism's theme stacks (`--msp-font-sans`, `--msp-font-mono`). Why not a free font param: the
  maintainers chose a table over extending `PD.Select` (2026-10-22).
- `setFontFamilyFonts(family, fonts)` changes an entry, at startup: texts already built keep their glyphs until they
  rebuild.
- A family string not in the table (an old state, a value set through the API) is used as is.

## 3. Web fonts load when an atlas needs them

- An entry may name a web font: CSS `@font-face`, a Google Fonts stylesheet, a `FontFace` added to `document.fonts`.
  None downloads until used; adding a `FontFace` is how a font by URL is offered.
- Text visuals wait before building (`withFontAtlasFonts(props, build)`): `document.fonts.check(font)`, else
  `document.fonts.load(font)`, one download per font however many wait. Wired into the structure label visual, the
  loci shapes (label, distance, angle, dihedral), MVS custom and annotation labels and MVS primitives; their
  `createGeometry` / `getShape` may return a promise.
- A text built without waiting while its font loads (the camera helper) gets an atlas of fallback glyphs; that atlas
  is dropped from the cache when the font arrives, so the next build draws the font.
- No `document` (Node, workers): nothing to wait for. node-canvas resolves the family list itself (Pango); register a
  font file with its `registerFont` before drawing.
- A failed download resolves: the fallback fonts stay. `load` is asked with the default text, so a Google Fonts subset
  by `unicode-range` other than the first loads only when a glyph in it is drawn, which then draws a fallback once.

## 4. Layout

- A label's text is split into graphemes (`splitGraphemes`, `Intl.Segmenter`; ASCII by code unit): an emoji with its
  modifiers, a letter with combining marks, an astral character are one glyph each.
- Each grapheme's pen position comes from `TextMeasure` (`@molstar/core/util/text-measure`, adapted from pretext's
  `pair-context` fit, MIT, attribution in the file): the advance of a grapheme after another is the two measured
  together less the first, so the pair's kerning is in it. Widths are cached per atlas. Why: advancing by the atlas
  cell width (rounded up, clamped) put no kerning and uneven gaps in every label.
- A glyph's cell holds its ink (`actualBoundingBoxLeft/Right`) plus the field's buffer on both sides; `FontAtlasMap.left`
  is where the pen starts in it. Why: ink left of the pen (`j`) or past the advance (`f`, italics) was clipped.
- Spaces advance but get no quad.
- Not done: shaping across graphemes (ligatures, Arabic joining, Indic clusters), bidi, line breaks.

## 5. Rendering (WebGL and WebGPU alike)

- Coverage is anti-aliased over one pixel: `sdfPerPixel` is the screen-space change of the field
  (`length(dFdx, dFdy)`, sampled before any branch: derivatives need uniform control flow), coverage
  `clamp((sdf - 0.5) / sdfPerPixel + 0.75, 0, 1)`. Why: a hard `sdf < 0.5` discard gave stair-stepped edges and broken
  thin strokes at small sizes.
- The +0.75 (not +0.5) puts the edge a quarter pixel out: one-pixel stems of small text stay solid, large text is
  unchanged to the eye.
- The fractional edge has alpha < 1 and is drawn in the transparent pass (`check_transparency` splits by alpha; text
  counts as half transparent in `opacityAverage`, so the pass runs whenever text is shown). Interior fragments stay in
  the opaque pass with depth.
- The border mixes into the fill by the fill's coverage; fill fragments at coverage ≥ 0.5 are pushed forward in depth
  as before, so a neighbouring glyph's border never covers them.
- Pick, marking, emissive and depth variants keep the hard edge at the field's middle.

## Open

- Atlas memory: 350 × lineHeight × 0.75 lineHeight bytes allocated up front (43 MB at `fontQuality` 3, 256 px glyphs),
  and each text renderable uploads its own copy of `tFont`. A growing atlas and a shared texture would cut both.
- Minification: no mipmaps on `tFont`; text a few pixels tall samples the 256 px field sparsely.
- Multi-channel SDF (MSDF) would keep corners sharp on very large labels; not needed at today's sizes.
- Multi-line labels (`\n`) and shaping (§4).
