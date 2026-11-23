# Lane LA4: the liftover

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/la4`. Model:
strong (projection through a map by either side, many-to-one images). HC is changing hover, emphasis and the card
across layers: the "2 places" card is the one card change here; build it on the card's existing rows and list it, so
HC's rule absorbs it at merge. PF is in perf test files only.

**Why**: [across-sections](../../../design/across-sections.md), "A liftover (EC-212; cut last, Gate 0, 1)", LA3's step 3
(LA3 merged steps 1–2 and 4). Stage 3's last lane before its gate.

**Task**, a green commit at each step:

1. **Runs read by either side**: a map's runs asked by region of `from` or of `to` (`RegionRequest.side`, built by LA2
   for pairs): forward projection reads by `from`, a layer drawn on the `to` section by `to`; the `seq` runtime reads
   its pairs whole and indexes both ends (no second read).
2. **The figure from existing nodes only**: a section of the old assembly, a bridge whose `pairs` layer's data is the
   alignment's map (`map.source` over the `seq` runtime's runs), a section of the new assembly with the old assembly's
   annotations drawn through the same map (a layer's data along a path, by region). Unaligned stretches are gaps in the
   ribbons; a part of an item that maps nowhere is not drawn on the new side (EC-21).
3. **Many to one, one to many**: a projection is the union of every image (EC-20); a projected hover's card lists the
   images ("2 places"); hover and selection projected through the map both ways.
4. **Example and stories**: a liftover over a generated chain (an old and a new toy assembly with a duplication, a
   deletion, an inversion and a translocation; truths computed in code, `PCG`, in `seqx-examples`);
   `Features/Across sections/Liftover` light and dark, canvas and SVG, its headless twin on a fake clock; the Design
   page "Across sections" gains it.

Deferred (the design's): `mapSummary` and the builder's `liftover(old, map, new)`.

**Tests**: a lifted item's drawn parts equal its images through the map (many to one, one to many, a deletion, an
inversion's strand); a run read by `to` equals the whole read clipped to the window; the card's images; EC-20, EC-21,
EC-212 in titles. **Performance**: a lifted layer's warm frame 0 B, its hover 0 B warm (allocation cases), a pan of the
new section relaying that section and the bridge only (counted).

**Where**: a worktree on `seq/la4`. Headers (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); commit at each green
step; hot files never staged (`package.json`s, `tsconfig*.json`, lockfile, `biome.json`, `.scaffolds/seq/`, `.v6/changelogs/`,
`storybook/packages.ts`, `packages/seq/TODO.md`): list them. No language change (a needed param stops the lane with a
report). Anything in `@molstar/core` or `@molstar/io` stops the lane. A file past 800 lines is a question. Sizes soft
(headless 960 at ~951, default 1,510 at ~1,500.5): trim first, keep the liftover behind the bridge's lazy door where it
is rare, report growth with its cause. Timing and allocation bounds are the gates. Allocation readings in a few cases
vary with file order (fixes.md, PF is on it): if one flips, rerun the file alone and report both.

**Checks**: the quick lane, `npx tsc -b packages/seq/core packages/seq/viewer packages/seq/examples`, each
`typecheck:test`, the storybook typecheck, `pnpm lint`, `pnpm format:check`, both perf tiers to completion,
`pnpm -s -C packages/seq/core size` (read the rows for FAIL), `pnpm test:medium packages/seq/`, the `.scaffolds/seq/_checks`
(features-planned with `--plan .scaffolds/seq/execution/fifth-slice/plan.md`). A fix's test fails on its parent.

**Report**: under 350 words: commits; checks with numbers; story ids; sizes; the card change (for HC); designs to amend;
hot-file needs; the `.v6/changelogs/` text.
