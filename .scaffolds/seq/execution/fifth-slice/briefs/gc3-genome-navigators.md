# Lane GC3: genome navigators and the viewer

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/gc3`. Model:
strong. BE runs beside you in `runtimes/bed` (and `seqx-examples`' recording script): no overlap expected; shared files
(`data/protocol.ts`) only grow, listed.

**Task**: build [genome-coordinates](../../../design/genome-coordinates.md)'s navigator and viewer part (Decision: "A
navigator over several systems"; Consequences: "navigators (GC3)" under `seqx-core`, and `seqx-viewer` (GC3)),
reading GC1's params and GC2's blocks (`blockPositionAt`, `blockUnitOf`, `BlockPlacement.runAt`, `placesOf`; core's
`_spec` Render: Blocks):

- core navigators: `navigator.system` an assembly (or a list) drawn as a section of its blocks (`bands`, `axis`,
  `tracks` by block); the window over blocks (one window across a gap; a section of one chromosome washes the others);
  the zoom strip of the block holding the window's centre (`navigator.zoomStrip.minPx`); `navigator.doubleClick:
  toggle-fit`; two axes over blocks × blocks (the square: its lattice and window, `navigator.square.maxCells`), one
  section on both axes; `navigator.switch` (`none`, `include`); the readout and slider text naming the block
  ("chr17:7,668,402"; `aria-valuetext` "chr17: 7,668,402–7,687,550 of GRCh38", across blocks both ends). **The
  square's reads by block and block pair are RG's** (region reads, stage 3): draw the lattice and what the static
  data gives; say in core's `_spec` what RG adds;
- viewer: `navigate.acrossBlocks`; the cursor across blocks (EC-154 amended); `block.next`, `block.previous`,
  `block.reverse`, `section.include` commands (change kinds `navigate`, `navigate`, `patch`, `patch`) with keys in
  `view.keys`; saved and history viewports by `{ system, position }`; the label menu's entries ("Reverse", "Show
  only", "Show all"; lazy with the menu chunk); a run's hover; the locus box's fallback to `findNames` (BZ's
  protocol types; test it on a fake source that answers names — GG builds the GFF index later);
- **GC2's leftovers from the orchestrator's look** (`fixes.md`, Open), in your area: a double-click on a block in the
  view (the 25-block row story) zooms to bases instead of fitting the block (`toggle-fit`: fit the block, then back to
  the whole); the reversed middle block of chr2 (B in A's order) shows no label though "chr2 ←" fits; the hover card
  calls a chromosome a "Chain" (name the row by the block's system kind); the ring's group labels ("peach", "grape")
  straddle the group arc and sit among the block labels (draw them in their own band outside the block labels, never
  over an arc);
- GC1's declared params this lane reads lose their `seq.resolve.param.unread` info (`navigator.switch` and the rest
  naming GC3); the design's Configuration rows marked GC3 built; the "Not built" line (genome-coordinates.md) is a
  hot-file need: give its new text;
- tests: switching blocks keeps the selection's addresses (property); a drag across blocks selects exactly the shown
  positions (property); history round-trips by `{ system, position }`; EC-154 (amended) and EC-184 named in titles;
  the window's move on a navigator over 25 blocks in 0.5 ms (perf tier, B6) and 0 bytes warm;
- stories: on the Design page "Genome coordinates": a chromosome at a time with the ideogram navigator (switching by a
  click), the square over blocks × blocks, a ring with `navigator.switch: include`; light and dark, canvas and SVG,
  each with its headless twin. Toy assemblies generated in the story; **no downloads**.

**Where**: a worktree on `seq/gc3`. Headers updated (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); hot files
never staged (`package.json`s, lockfile, `biome.json`, `.seq/`, `.v6/changelogs/`, `storybook/packages.ts`,
`packages/seq/TODO.md`): list them. Public API changes allowed in `packages/seq`, core's `repr` and `@molstar/ui`;
anything else in `@molstar/core` or in `@molstar/io` stops the lane and comes back as a question.

**Checks**: as `fx1-checks-and-carried.md` lists them, both perf tiers (`pnpm -s -C packages/seq/{core,storybook}
test:perf`), `pnpm -s -C packages/seq/core size` (bounds: reading 240, headless 870, default 1,400, minimal 1,000 kB;
the menu entries lazy), `pnpm test:medium`, `pnpm test:heavy packages/seq/`, the four `.seq/_checks` and
`features-planned.ts --plan .seq/execution/fifth-slice/plan.md`. A fix's test fails on its parent. Do not claim a
failure is the base's without running it on the base.

**Report**: under 300 words: commits; checks with numbers; the story ids for the orchestrator's look (light and dark,
canvas and SVG); the public surface added; what RG must add to the square; sizes; hot-file needs; the
`.v6/changelogs/` text.
