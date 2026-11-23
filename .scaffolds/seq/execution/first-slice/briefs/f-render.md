# Lane F: layout and draw lists

**Goal**: in `packages/seq/core/src/render/`: tier 3 for the first slice ([layout-and-draw-lists.md](../../../design/layout-and-draw-lists.md)):
the layer-type contract (`name`, `fallback`, `prepare`, `layout`, `draw`, `hit`, `legend`), the layer types `letters`,
`ranges` and `values`; the tracks layout (sections stacked, each with its axis transform from its system and range,
tracks with heights in px, one label column); draw lists in CSS pixels with the passes (background, data, highlight,
selection, hover, cursor, labels) and primitives as structs of arrays in reused growable typed arrays (rect with
squircle corners, path, glyph runs, text with a max width, image, clip); hit indexes (cells by arithmetic, items by an
interval index over drawn extents per lane, hits as `{ layer, element?, position?, system }`); stable lanes (sticky by
item key); the density steps (letters, cells, dots, bars per run, density) with thresholds from a config object; the
token table (`SeqTokens`, a built-in light and dark table for Node) and its version in the draw lists' key; the text
metrics port (`measure(font, text)` over `@molstar/core/util/text-measure`, with a bundled metrics table for Node); the
recording backend (one line per primitive, pass by pass, numbers to two decimals) and the SVG backend (a string, text
escaped); `renderSvg(spec, { size, viewport?, runtimes })`, asynchronous, waiting for its data through a store wanter
(`seqx-core/data`). Data comes from the store's entries (inline data as `seq` sources); ranges sliced across tiles are
deduped by `id` or `row`. Draw lists are bounded by the screen (device columns × lanes).

**Read first**: `.seq/design/layout-and-draw-lists.md`, `rendering.md`, `composition.md` (tier 3), `ux-foundations.md`
(density, shapes, emphasis order), `performance.md` (cost follows pixels), `data-store.md`, `protocol.md`, `security.md`
(SVG escaping), the core `_spec` (Resolve, Resolved view, Data, Queries), `.seq/design/testing.md`.

**Model**: strong.

**Base**: `packages-seq` at the commit named in the launch instructions. Branch `seq/f-render`; commit on it only.

**Touches**: `packages/seq/core/src/render/**` (new) and its `_test/`; the core `_spec`'s new section "Render" only.
Read-only: `lang/`, `query/`, `data/`, `resolve/`, `view/` (ask in the report for changes there). Not in scope: the
Canvas2D backend and the viewer (lanes G and H), rings, text layout, dot plots.

**Tests**: two draw-list records (p53 tracks; p53's protein on its gene, gene section following); properties: draw lists
bounded by the screen at every density step (a 3M-position synthetic system and p53 give the same primitive counts per
device column), hits where marks are, lanes stable when a neighboring tile arrives, the token version invalidating,
SVG escaping; `renderSvg` on the p53 example waits and renders.

**Hot files, never staged**: `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`, `biome.json`,
`vitest.config.ts`, `storybook/packages.ts`, `.v6/changelogs/`, `packages/seq/*/package.json` (ask in the report for
the entry exports and `_test` nulls you need), anything in `.seq/` and in `packages/core/`. Ask in the report for
what you need there.

**Headers**: the Mol\* header, `Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more
info.`, with `@author David Sehnal <david.sehnal@gmail.com>`. No mention of `packages/sequences` or `.seq/` in code,
comments or tests.

**Checks**, each its own command, green only on exit 0 with a non-zero test count:
- `pnpm test:unit packages/seq/core` (your folder under 2 s, a test under 100 ms)
- `npx tsc -b packages/seq/core`
- `pnpm lint`
- `pnpm format:check`

**Stop when** the checks are green, or a check has failed twice after a fix, or a permission is refused (report the
action and why; do not work around it).

**Report**: the commits; each check with its exit code and test count; the quick lane's time for your folder; the
exports and `_test` nulls you need in `packages/seq/core/package.json`; any place the designs were unclear or wrong (as
drift: what, where, options, your recommendation); followups.

**Configurable by default**: where a choice could reasonably go more than one way, make it a config item (or a param) with
the recommended choice as its default, and name it in the report.
