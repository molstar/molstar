# Lane LA3: grids over two sections, categorical axes, the liftover

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/la3`. Model:
strong (the grid's two axes, a permutation through the pyramid, projection through a map). Alone in resolve, the store
and the grid layer. **SX** (section scroll) is in its fix round in `render/layout/`, the stacked layout and the rail:
touch them only where a grid section must, and list each touch.

**Why**: [across-sections](../../../design/across-sections.md) Consequences steps 5–6, LA's rest (LA did 1–4, LA2 the
per-side reads): "Contact maps and dot plots across sections", "Categorical grid axes", "A liftover" (EC-209, EC-212,
EC-223, EC-35).

**Task**, a green commit at each step; **stop after step 2** if the lane runs long and report what is left:

1. **A grid over two sections of blocks**: `follow_y` naming a second section of blocks gives the grid its y blocks;
   each block pair a strip pair (cells clipped at both blocks' edges, labels on both axes, gaps and runs per axis, a
   reversed block flipping its side, EC-145); reads by block pair at the finer axis's level, the upper triangle only
   when one assembly and `symmetric`; a dot plot of two genomes (the `pairs` layer in the grid: dots, diagonal
   segments by strand, bins by count). Still refused: a grid following blocks in a stack or on a ring.
2. **Categorical grid axes**: `of`/`of_y` a rows system, `member_order`/`member_order_y` (declared in `language.ts`:
   **no language bump**; anything undeclared stops the lane with a report); the axis transform holding the order's
   permutation and its inverse (`Int32Array`), the pyramid keyed by the order's version; chrome (no ruler, member
   labels thinned past `rows.labelMinPx`, x labels at 45°, a box fixed under zoom); hover names both members and the
   value; a label click selects the member.
3. **The liftover** (EC-212), from existing nodes only: an old section, a bridge over the alignment's map, a new
   section with the old annotations drawn through the map; runs read by either side (`RegionRequest.side`); a
   projection the union of its images, the card listing them ("2 places"); a part that maps nowhere not drawn.
   `mapSummary` and the builder's `liftover()` stay deferred.
4. **Examples and stories**: peach × grape as a dot plot (an inversion an anti-diagonal), chr1 × chr4 contacts, a
   kinship grid beside a genotype matrix (one rows system), a liftover over a generated chain (truths computed in
   code, `PCG`); `Features/Across sections/*` beside LA's, light and dark, canvas and SVG, headless twins on a fake
   clock; the Design page "Across sections" gains them.

**Tests**: a block pair's cells equal the whole grid's clipped to both blocks; the upper triangle mirrors exactly; a
reorder rebuilds the pyramid once and hits read members through the permutation (a property over random orders); a
lifted item's drawn parts equal its images through the map (many to one, one to many, a deletion); EC-35, EC-145,
EC-209, EC-212, EC-223 in titles. **Performance**: a grid frame over two sections 0 B warm (an allocation case), a
reorder over 2,000 members under its bound (heavy `itWithin`), the liftover's hover at 0 B warm.

**Where**: a worktree on `seq/la3`. Headers (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); commit at each
green step; hot files never staged (`package.json`s, `tsconfig*.json`, lockfile, `biome.json`, `.seq/`,
`.v6/changelogs/`, `storybook/packages.ts`, `packages/seq/TODO.md`): list them. The `_spec` sections that move change in
the same commits; designs to amend listed for the orchestrator. Anything in `@molstar/core` or `@molstar/io` stops the
lane. A file past 800 lines is a question, past 1500 a failure. **Sizes** (headless 940 — at 937.9 —, default 1,495 —
at 1,481.3 —, minimal 1,060, the bridge chunk 28 kB): soft. Trim first (duplicates, dead code, code on the wrong side
of an existing lazy door); keep the liftover and categorical chrome behind the grid's or the bridge's door only if they
are rare and load without jank; report any growth with its cause. Do not edit bounds. Timing and allocation bounds are
the gates: never loosen them.

**Checks**: the quick lane (under 2 s a touched folder), `npx tsc -b packages/seq/core packages/seq/viewer
packages/seq/examples`, each `typecheck:test`, the storybook typecheck, `pnpm lint`, `pnpm format:check`, both perf
tiers to completion, `pnpm -s -C packages/seq/core size` (read the rows for FAIL), `pnpm test:medium packages/seq/`,
the heavy files you add when load is under 6, the `.seq/_checks` (features-planned with
`--plan .seq/execution/fifth-slice/plan.md`). A fix's test fails on its parent.

**Report**: under 400 words: commits; steps done and left; checks with numbers; story ids; sizes; the public surface;
SX-area files touched; designs to amend; hot-file needs; the `.v6/changelogs/` text.
