# Lane IZ: the zoom strip

**Goal**: [navigators](../../../design/navigators.md), "The zoom strip" (the author, 2026-11-19: "the window on the gene
strip doesn't work. cant interact with it and it isn't synced", then "Yes, the zoom strip"). Features › Navigators ›
Gene Region On Chr 17 (`packages/seq/examples/src/views.ts`, `gene-region-ideogram`; p53's 6 kb region on chr17's
83 Mb): the window's true extent is 1/13,800 of the ideogram, drawn at the floor width, so it never visibly moves and a
drag moves the view by 78 kb a pixel. Build: where a navigator's reachable span is narrower than
`navigator.zoomStrip.minPx` (48 CSS px), the navigator draws the span as a mark (an outline at least
`navigator.window.minPx` wide, centred, no grips, the default cursor, not a drag target; hovering it names the region)
and adds a zoom strip under it: the same navigator kind's form over the reachable span alone (an ideogram's bands
there; a strip's axis or tracks), with the real window, its grips and its gestures; two connector hairlines
(`seq-window`) from the mark's ends to the zoom strip's ends. `navigator.zoomStrip.heightPx` (16) more height. The pair
is one navigator: one mount, one hit target, one navigator record, one draw list set (`window` pass as now). Config in
`ViewerConfig` with the defaults above. The viewer `_spec`'s navigator section and Configuration say it. Minimal code:
the zoom strip is the existing navigator model over a narrower axis, not a second implementation.

**Tests** (fail on the parent): on the gene region with the ideogram: a drag on the zoom strip's body pans the view by
the zoom strip's scale; a grip zooms; a view zoom or pan moves the zoom strip's window; a press on the mark is not a
drag; the connector lines and both strips drawn in canvas draw lists and SVG, both themes; a navigator whose reachable
span is wider than the floor keeps one strip (the plasmid's, p53's protein strip). Allocation once warm unchanged
(the navigator's guard).

**Not in scope**: anything else. FS-A1 (under review, not merged) touches data and sources only.

**Model**: strong (opus). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B seq/iz-zoom-strip
<that commit>`.

**Checks**: `pnpm test:unit packages/seq/`; `pnpm test:medium packages/seq/`; `npx tsc -b packages/seq/core
packages/seq/viewer packages/seq/examples`; `pnpm -C packages/seq/<core|viewer|examples> typecheck:test`; `pnpm -C
packages/seq/storybook typecheck`; `pnpm lint`; `pnpm format:check`; `pnpm -C packages/seq/core test:perf`; `pnpm -C
packages/seq/viewer test:perf` if it exists. Exit codes explicit (`; echo exit=$?`), non-zero test counts.

**Look**: run the storybook from your worktree and look at the story in light and dark before reporting (a drag on
the zoom strip, a zoom in the view); say what you saw.

**Headers**: a modified file keeps its header, year range to 2026, `@author David Sehnal <david.sehnal@gmail.com>` added
if missing; a new file has the Mol\* header with that author. No mention of `.scaffolds/seq/` in code.

**Hot files, never staged**: `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`,
`biome.json`, `vitest.config.ts`, `storybook/packages.ts`, `.v6/changelogs/` (write the entry's text in the report),
every `package.json` and `tsconfig*.json`, anything in `.scaffolds/seq/`.

**Stop when** green, or a check fails twice after a fix, or the work needs a public API past the config above.

**Report**: commits; checks; what you saw in the browser; drift.
