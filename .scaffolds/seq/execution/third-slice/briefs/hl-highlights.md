# Lane HL: overlapping highlight channels

**Goal**: [ux-foundations](../../../design/ux-foundations.md) section 6, "Overlapping highlight channels stay apart"
(the author, 2026-11-19: "the packages/sequences is able to select/highlight/mark multiple overlapping sections with
different colors. can we have this functionality too?"). Today `drawHighlights` (`packages/seq/core/src/render/
tracks.ts`) paints every channel as a 0.3 tint over the whole band, so two overlapping channels blend into a third
color. Build:
1. Where two or more channels light one position of a track, the band splits into equal horizontal stripes, one a
   channel present there, in channel order; up to `emphasis.highlightStripes` (4, a `RenderConfig` value), the rest
   merged into the last stripe; one channel alone takes the whole band. In tracks, text lines and rings (a ring's band
   radially), canvas and SVG alike.
2. A channel's `style`: `tint` (the default, today's look), `band` (a bar along the track's foot), `outline` (its
   color's outline, no fill), set with the channel (`viewer.highlight.set`'s value) and kept through its carry onto
   linked sections.
3. The legend part lists the active channels with their colors, each removable (`highlight.clear` with its name), as
   [themes](../../../design/themes.md) plans ("Highlight channels in the legend part").
4. `Features/Highlights` (examples-and-stories: a feature's story under `Seq* (next)/Features/<Feature>`): on p53,
   three overlapping channels (a domain, a motif search result, a host's own) set and cleared from buttons, each
   style shown; its headless twin checks the stripes' geometry and colors in canvas draw lists and SVG.

The `_spec`s (core's Render, Emphasis; the viewer's highlight model) say the rule. Minimal code; allocation once warm
unchanged (the emphasis guard).

**Tests** (fail on the parent): two and five overlapping channels in tracks, text and a ring, canvas and SVG, both
themes (stripes, order, the merged last); each style; a carried channel keeps its style; the legend lists and clears.

**Not in scope**: anything else. IZ (navigators) and FS-A1 (data) run beside you.

**Model**: fast (sonnet). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B seq/hl-highlights
<that commit>`.

**Checks**: `pnpm test:unit packages/seq/`; `pnpm test:medium packages/seq/`; `npx tsc -b packages/seq/core
packages/seq/viewer packages/seq/examples`; `pnpm -C packages/seq/<core|viewer|examples> typecheck:test`; `pnpm -C
packages/seq/storybook typecheck`; `pnpm lint`; `pnpm format:check`; `pnpm -C packages/seq/core test:perf`. Exit codes
explicit (`; echo exit=$?`), non-zero test counts.

**Headers**: a modified file keeps its header, year range to 2026, `@author David Sehnal <david.sehnal@gmail.com>` added
if missing; a new file has the Mol\* header with that author. No mention of `.scaffolds/seq/` in code.

**Hot files, never staged**: `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`,
`biome.json`, `vitest.config.ts`, `storybook/packages.ts`, `.v6/changelogs/` (write the entry's text in the report),
every `package.json` and `tsconfig*.json`, anything in `.scaffolds/seq/`.

**Stop when** green, or a check fails twice after a fix.

**Report**: commits; checks; the story to look at in light and dark; drift.
