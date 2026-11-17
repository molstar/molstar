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

1. **The spec** (JSON, documented by the language's generated reference pages): a viewer that maps the language onto its own model reads this, as the
   JBrowse converter did.
2. **The resolved view** (core, pure; the view IR of [composition.md](composition.md)): the spec with refs followed, defaults filled, maps resolved to segment tables,
   themes resolved to values, fallbacks chosen per layer type the target declares it can draw, and data reached through
   the source protocol by region. Renderer-independent: no pixels, no fonts. Other viewers (JBrowse, Nightingale) read
   this to get Seq\*'s meaning without reimplementing it.
3. **Draw lists** (core, pure): the resolved view laid out for a size and a viewport, as plain typed lists of
   primitives (rects, paths, glyph runs, text, images) in passes (background, data, overlay), with clips. Hit-testing
   reads the same layout.
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
- **One picture**: the Canvas2D and SVG backends receive the same draw lists, so they differ only in rasterization. A
  test compares their recorded draw calls; a visual test (heavy lane) compares pixels.
- **A render entry and a CLI**: `renderSvg(spec, { size, viewport?, runtimes })` and a command that renders a spec file
  to SVG or PNG, so a reply like the JBrowse one can put both pictures side by side without writing code.
- **Light to depend on**: the resolved view and draw lists are importable without the viewer, the UI or any format
  runtime; reading a spec needs only the core's language entry.
- **A reference renderer that is not ours**: a subset of the resolved view (letters, ranges, values, variants) is also
  drawn by Nightingale, through a small private adapter, beside Seq\*'s own picture in Storybook and in a test that
  checks the adapter reads every field it needs from the resolved view alone. Why: it keeps the resolved view general,
  a form another renderer can draw without Seq\*'s layout (the author, 2026-11-17: "we could also still keep the
  reference nightingale rendering of a subset as a sanity check for the renderer being general").

## Rejected alternatives

- Canvas only, export by screenshot: no headless rendering, no vector output, no side-by-side comparison.
- Exposing only the spec: every other viewer would reimplement maps, themes and fallbacks, and draw something subtly
  different.

## Consequences

- The resolved view's and draw lists' shapes are contracts with their own section in the core's `_spec`; a change to
  them is a language-level change (versioned, listed in the changelog).
- Every example renders to SVG in a quick-lane test, with its draw-list record kept; Storybook shows the canvas and the
  SVG side by side for any example.
