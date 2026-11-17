# Layout and draw lists

## Intent

The second review found no layer-type contract, no layout interface and a draw-list format named only by its primitive
kinds: no units, passes, clips, text model, hit-test index, recording format or buffer layout, and nothing about colors
and fonts without CSS. The first slice's step 3 needs all of them. This design gives the contracts; the numbers are
settled on the slice.

## Decision

### The layer-type contract (a handful of members)

```ts
interface LayerType<P, Prep, Lay> {
  readonly name: string;                         // 'letters', 'ranges', 'values'
  readonly fallback?: string;                    // the layer type a target without this one draws
  prepare(layer: ResolvedLayer<P>, model: StaticModel): Prep;                       // tier 1, cached by hash
  layout(prep: Prep, ctx: LayoutContext, data: StoreReader, overlay: Overlay): Lay;  // per viewport, reused buffers
  draw(lay: Lay, pass: Pass, out: DrawLists): void;
  hit(lay: Lay, x: number, y: number): Hit | undefined;
  legend?(lay: Lay): LegendSpec | undefined;
}
```

- Built-ins are written to it, as extensions will be ([extensibility.md](extensibility.md)); the previous
  implementation's more than 25 optional members are not carried over.
- `LayoutContext` holds the track's box, the axis transform (positions to x, with pieces and layout maps), the device
  pixel ratio, the density thresholds, the token table and the text measurer.

### The layout interface (internal)

A section layout (tracks, text, rings, dot plots) takes a section of the resolved view, the configuration and a size,
and gives the axis transform per block, each track's box, and the order of its layers; it calls the layer types. It is
internal until an outside application needs one ([ux-foundations.md](ux-foundations.md)).

### Draw lists

- **Units**: CSS pixels, origin top left. The layout snaps cell edges to the device grid it is given (the device pixel
  ratio is in the context), so backends draw what they get.
- **Passes**, each its own list, drawn in order: `background`, `data`, `highlight`, `selection`, `hover`, `cursor`,
  `labels`. The Canvas2D backend draws `background` to `data` on one canvas and the emphasis passes on an overlay
  canvas, so a hover redraws only the overlay (the 2 ms target, [performance.md](performance.md)). SVG writes them as
  groups in order.
- **Primitives**, each a struct of arrays in growable typed arrays reused frame to frame (nothing allocated once warm):
  - `rect`: x, y, w, h, fill, stroke, corner (`square`, or `squircle` with a radius);
  - `path`: commands and points in a shared `Float32Array` (move, line, cubic, close), fill, stroke, dash;
  - `glyphs`: a run of single letters in one font, an x per letter, one y, a color per letter or one for the run;
  - `text`: a string from an interned table, x, y, anchor, font, color, a maximum width (faded past it, never
    ellipsized, as `@molstar/ui`'s rule);
  - `image`: RGBA bytes with a width and height, for dense matrices and grids;
  - `clip` push and pop, by rect.
- **Colors** are 32-bit RGBA in the lists; a theme resolves to them through the token table.
- **The recording backend** writes one line per primitive, pass by pass, numbers rounded to two decimals
  (`rect 10 20 5.5 16 #ff0000ff squircle 3`): the format of the few committed records ([testing.md](testing.md)).

### Hit-testing

- Each layer's layout keeps a hit index: cells by arithmetic on the axis transform; items by an interval index of their
  drawn extents per lane, with lanes' y bands; marks larger than their cell by their drawn box (EC-123).
- A hit is `{ layer, element?, position?, system }`: an element key and its Location, or a position; the viewer turns
  it into hover and selection ([viewer-models.md](viewer-models.md)).

### Stable lanes while data streams

- Lanes are assigned by item key and kept while the item is in view (sticky lanes); a new item takes the first lane free
  over its extent; lanes are repacked only when the view settles and the lanes have become sparse. Why: items must not
  jump when a neighboring region arrives.

### Colors without CSS

- **A token table** (`SeqTokens`: every `--seq-*` data token and the `--msui-*` tokens the draw lists use, resolved to
  colors) is part of the layout context. In the browser it is read once from the computed style of the viewer's root;
  in Node it comes from a built-in table for the light and dark themes. Its version is part of the draw lists' key, so a
  theme switch redraws.
- **Colors in a spec** are literals (CSS color strings) or token names (`token:accent`); an export resolves tokens with
  the table of its theme and writes literals.

### Fonts

- **Letters are centered in their cells** and need no measuring; their size follows the cell
  ([ux-foundations.md](ux-foundations.md)).
- **Labels are measured** through the text metrics port. The defaults are `@molstar/ui`'s font tokens (a sans stack and
  a mono stack); the browser measures with a canvas (`TextMeasure`, cleared when a web font loads, then relaid), Node
  with a bundled metrics table for an Arial-compatible sans and a fixed-advance mono. An export from Node may place a
  label's fade a few pixels differently from the browser; letters and data do not move.
- The port's `measure(font, text)` keeps one `TextMeasure` per font behind it.

## Consequences

- The first slice builds the contract with three layer types (`letters`, `ranges`, `values`), the tracks layout, the
  passes, the recording and SVG backends and a Canvas2D backend with its overlay.
- Tests: draw lists bounded by the screen (device columns × lanes), hits where marks are, lanes stable across an
  arrival, the token table's version invalidating, the recording format.
