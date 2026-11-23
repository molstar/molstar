# Rendering: ports, SVG, and forms other viewers read

## Intent

A spec must draw the same picture on screen, in an exported SVG or PNG, headless in Node (tests, a CLI, a server),
and in other viewers. The JBrowse reply showed the last one working: its author drew Seq\* documents in JBrowse, beside
Seq\*'s own SVG renders, and compared the two. That must stay possible, and get easier.

## Precedents

- The previous implementation: a headless context laid out and hit-tested a view; a composer turned it into batches
  for an abstract backend (Canvas2D with a glyph atlas, SVG and PNG export, a recording backend for tests); a
  `renderSvg` and a render CLI; a resolved view (`resolveView`) that adapters read (Nightingale) and that the JBrowse
  spike compared against ([sources/seqstar](../sources/seqstar/README.md), [sources/jbrowse-reply](../sources/jbrowse-reply/README.md)).
- The reply's cost: its converter read the raw JSON tree instead of the core, because the core weighed 3.5 MB and
  needed host stream manifests before anything resolved.

## Decision

The author, 2026-11-17: "make sure we keep the render ports, the ability render to svg etc. because in the "review
reply" Colin was able to render it in Jbrowse".

Four forms, each public, documented in the core's `_spec`, and each a stable place another renderer can attach:

1. **The spec** (JSON by default, documented by the language's generated reference pages): a viewer that maps the
   language onto its own model reads this, as the JBrowse converter did.
2. **The resolved view** (core, pure; [composition.md](composition.md)): views, sections, tracks and layers with refs
   followed, defaults filled, maps resolved, themes as values or value functions, each layer's fallback chain, and how
   each layer's data reaches its section. No data, no pixels, no fonts, no target: the same for every renderer. Other
   viewers (JBrowse, Nightingale) read it, and read data through the same source protocol and store.
3. **Draw lists** (core): the resolved view laid out for a size, a viewport, a target (which picks each layer's first
   drawable fallback) and the data store's entries ([data-store.md](data-store.md)), as typed lists of primitives
   (rects, paths, glyph runs, text, images) in passes, with clips. Built into reused buffers: pure by contract (same
   inputs, same lists). Hit-testing reads the same layout.
4. **Backends** (the render port): a backend draws draw lists. Canvas2D on screen, SVG (a string, no DOM), PNG (from
   the SVG or a canvas), and a recording backend for tests. An application may add its own (WebGL, a PDF writer).

- **Headless**: the core plus the SVG backend run in Node with no DOM. Layout needs text metrics, so measurement is a
  port: the browser measures with a canvas, Node with a bundled metrics table for the default font; both give the same
  layout for the default font.
- **Text metrics build on what the base library has**: `@molstar/core/util/text-measure` (`TextMeasure`, adapted from
  pretext: grapheme advances with kerning, a width cache per font, used by the graphics font atlas) takes any measure
  function, so the port's adapters are its measure function (a canvas, or the bundled table). What Seq\* needs beyond it
  (pretext's prepare and layout split, lines for the text view, inline chips) is proposed as an extension of that
  module in the core, not a private copy (root `AGENTS.md`, reuse before writing helpers). Pretext itself throws without
  a canvas, so it is a candidate for the browser adapter of the text view's inline layout, not the port's only one.
- **One picture**: the Canvas2D and SVG backends receive the same draw lists, so they differ only in rasterization; the
  quick lane tests the draw lists through the recording backend, the heavy lane compares the two backends' pixels.
- **A render entry and a CLI**: `renderSvg(spec, { size, viewport?, runtimes })`, asynchronous (it waits for the data
  the export needs to settle, [data-store.md](data-store.md)), and a command that renders a spec file
  to SVG or PNG, so a reply like the JBrowse one can put both pictures side by side without writing code.
- **`renderSvg` answers `{ svg, settled, issues, lists }`**: the SVG carries `data-*` attributes for whether it settled
  and the open domains it computed; the lists are returned for tests. A section whose length is unknown takes the length
  its source's coverage gives once attached, else the export reports `seq.render.axis.pending` and waits.
- **Light to depend on**: the resolved view and draw lists are importable without the viewer, the UI or any format
  runtime; reading a spec needs only the core's language entry.
- **A second consumer, in tests only**: a subset of the resolved view (letters, ranges, values, variants) is also
  drawn by a test-only consumer, plain SVG with no dependency, that reads every field it needs from the resolved view
  and the data store alone; its test fails if it needs anything private. Why: it keeps the resolved view general, a
  form another renderer (JBrowse, Nightingale) can draw without Seq\*'s layout (the author, 2026-11-17: "a sanity check
  for the renderer being general"); Nightingale itself dropped (the author, 2026-11-19: "drop it as the public
  feature … only keep it as test/demonstration", "our own test consumer only"): a pinned outside renderer tests the
  translation more than the form.

## Rejected alternatives

- Canvas only, export by screenshot: no headless rendering, no vector output, no side-by-side comparison.
- Exposing only the spec: every other viewer would reimplement maps, themes and fallbacks, and draw something subtly
  different.

## In the specs

| Decision | Where |
|---|---|
| Four public forms, each a place another renderer attaches | core › Render |
| The spec, read by a viewer mapping the language | core › Language |
| The resolved view: no data, pixels or target; other viewers read it | core › Resolved view |
| Draw lists: laid out for a size, viewport, target and the store; pure by contract; hits on the same layout | core › Draw lists |
| Backends: Canvas2D, SVG, PNG, recording; an application may add its own | core › Backends |
| Headless: the core and the SVG backend in Node; text measured through a port, Node's bundled table | core › Text |
| Seq\*'s needs past `TextMeasure` proposed as an extension of core's module | design-only: a proposal for `@molstar/core`, not this package's |
| One picture: both backends draw the same lists; the recording in the quick lane | core › Backends |
| `renderSvg` asynchronous, waiting to settle; the render CLI | core › The render CLI |
| `renderSvg` answers `{ svg, settled, issues, lists }`; the SVG's `data-*` | core › Call options |
| A length unknown waits (`seq.render.axis.pending`) | core › Tracks layout |
| Light to depend on: the resolved view and draw lists without the viewer, UI or runtimes | core › Entries carry what their readers need |
| A second consumer, in tests only | core › A second consumer |
| Canvas only; the spec only | design-only: rejected alternatives |

## Consequences

- The resolved view's and draw lists' shapes are contracts with their own section in the core's `_spec`; a change to
  them is a language-level change (versioned, listed in the changelog).
- Every example renders to SVG in a test; a few representative views keep a draw-list record ([testing.md](testing.md));
  Storybook shows the canvas and the SVG side by side for any example.
