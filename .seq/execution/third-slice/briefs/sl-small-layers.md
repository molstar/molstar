# Lane SL: the small layers

**Goal**: build the small layers of [themes](../../../design/themes.md), section "The small layers: what each reads
from the theme", as agreed:
- `ruler`: its `numbering` param (the numbering it labels by, the section's default when unset; switching schemes is a
  patch of it), insertion labels (100A, 111.1) as minor ticks labeled only with room, steps of 1, 2, 5 × 10ⁿ, a label
  skipped when it does not fit (never shrunk), on a reversed axis counting down (the axis, built by AX); fallback
  omitted;
- `secondary`: helices and strands as blocks, the strand's arrowhead on the piece holding its C-terminal end whatever
  the direction or the origin (EC-130), cut ends square, 200 copies drawing what one does (EC-67); fallback `ranges`;
- `variants`: a stem and a head a site; at `bars` and `density` a count strip; fallback `ranges`;
- `values`: `encoding` `line` and `area` beside `bars` (EC-57: color changing halfway, held over half cells at the ends
  and gaps, missing breaking the line, the area under it);
- `ranges`: line styles anchored to the content (the dash phase moves with the data), an invalid style drawn solid and
  said once (`seq.render.line-style`), an arrowhead on the piece holding the end it points to, never on a cut end.

Each layer type's kind in the language with its default fixed by the language (never the host or the data), its
fallback, its hits (`Hit` records as the existing layers), its SVG export, and its tests drawing against the data.

**Coordination with TH** (the theme lane, beside you): build against today's color path (a layer's existing `color`
param and the tokens); TH owns the theme slot, the channel accessor and `render/tokens.ts` and wires your layers to the
channels at its merge, which lands after yours. Where you need a token the table lacks, ask in the report (TH adds
it). Do not edit `lang/` theme files or `render/tokens.ts`.

Not in scope: themes, legends, the theme slot; layouts other than tracks; the `bands` layer (NV).

**Read first**: [themes](../../../design/themes.md) (the small layers section, Decision's author's decisions, the
Configuration rows for these layers); [axis](../../../design/axis.md) (what reads the axis: a reversed and collapsed
axis); [layout-and-draw-lists](../../../design/layout-and-draw-lists.md); [edge-cases](../../../design/edge-cases.md)
(EC-57, EC-67, EC-130, EC-234); `packages/seq/AGENTS.md`; the core `_spec`; `packages/seq/core/src/render/` (the
existing `ranges`, `values`, `letters` layers as models).

**Model**: fast (sonnet).

**Base**: `packages-seq` at `<set at start>`. Your worktree may be cut from master: with a clean tree, first
`git checkout -B seq/sl-small-layers <base>`, and check `git merge-base --is-ancestor <base> HEAD`.

**Touches**: `packages/seq/core/src/{lang,render}/**` (the new kinds and layers; not the theme files or tokens),
`packages/seq/examples/src/**` (the p53 views gain a ruler, secondary structure and variants where examples.md says),
`packages/seq/storybook/src/_stories/**` (a story per layer, light and dark), the core `_spec` (Layer types, Issues,
Configuration), `packages/seq/TODO.md`.

**Two consumers**: each layer drawn by the canvas and by `renderSvg` from one spec, their hits agreeing; on a reversed
collapsed axis as on a plain one.

**Tests**: per layer, drawing against the data (every datum in view has its marks, every mark is hit) at the fit view,
across a tile edge, and on a reversed collapsed axis; property tests for the ruler's steps (labels never overlap, every
label a position on the numbering) and the arrowhead rule; the fallbacks.

**Performance representative**: a frame of each layer at the fit view under the frame target (`itWithin`, best of
three); draw and hover warm at 0 bytes (the opt-in tier: run it and report).

**Whose each choice is**: the design's Configuration rows for these layers; defaults fixed by the language; nothing
configurable beyond it.

**Minimal code**: reuse the existing layers' machinery (tiles, merge, hits, the axis transform); no private copies; sizes
reported.

**Checks**: each its own command, green only on exit 0 with a non-zero test count.
- `pnpm test:unit packages/seq/`
- `pnpm test:heavy packages/seq/` (never while another heavy run is going; if one is, wait)
- `npx tsc -b packages/seq/core packages/seq/viewer packages/seq/examples`
- `pnpm -C packages/seq/core typecheck:test`, the same for `viewer` and `examples`, and `pnpm -C packages/seq/storybook typecheck`
- `pnpm lint`, `pnpm format:check`
- `pnpm -C packages/seq/core size` (report the numbers)

**Headers**: a new file has the Mol\* header `Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for
more info.` with `@author David Sehnal <david.sehnal@gmail.com>`; a modified file keeps its header, its year range
extended to 2026 and that author line added if missing. No mention of `packages/sequences` or `.seq/` in code, comments,
tests or `_spec`s; code cites its own `_spec` sections.

**Fixes land with a test shown failing on the fix's parent**: the report names the command.

**Hot files, never staged**:
- `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`, `biome.json`, `vitest.config.ts`;
- `storybook/packages.ts`, `.v6/changelogs/` (write the changelog entry's text in your report);
- every `package.json` and `tsconfig*.json` (ask in the report);
- anything in `.seq/`.

**Stop when** the checks are green, or a check has failed twice after a fix, or the work needs a public API past the
design, or a target cannot be met without changing a design (report the measurement and the options), or a permission
is refused (report the action and why; do not work around it).

**Report**: the commits; each check with its exit code and test count; the sizes and the performance numbers; the
hot-file changes you need; the tokens you need from TH; the changelog entry's text; config items with their defaults;
drift (what, where, options, your recommendation); followups. You cannot see a browser: name the stories to look at,
and what to look for in each (light and dark).
