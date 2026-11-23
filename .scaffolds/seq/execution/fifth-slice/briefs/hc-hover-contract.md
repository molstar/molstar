# Lane HC: one hover contract for every mark, and a sweep that holds it

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/hc`. Model:
strong (a rule across every layer type and layout). LA3 is in the grid layer, resolve and the store: where a grid's
hover needs a change, list it rather than reshape LA3's code.

**Why**: `fixes.md` Open, in its words: "One hover contract for every mark, and a sweep that holds it" (the author,
2026-11-21: "we should not be fixing this one-by-one btw ... we need a general rule for this and global fix"); "A
legend's hover works one way only"; "A VCF allele's hover card names the record, not the allele"; and the items they
close: the density row's missing highlight (in "The hover card: placement, wrapping, ..."), "Hovering a ligand's score
box lights its residues but not its arcs", the strands' density bands on a ring and a strip. Read each item whole.

**Task**, a green commit at each step:

1. **The rule** in `ux-foundations.md`'s words (list the text for the orchestrator, it is in `.seq/`) and in the core
   `_spec`'s Emphasis: every drawn mark of every layer type, in every layout (row, ring, blocks, rows, grid, bridge) and
   at every level of detail (item, lane, density bin or band, summary, chip, reduced tile), answers a hover with (1) a
   hit naming what is under the pointer (the item, else the bin's items or stretch, else the position), (2) an
   emphasis mark in the hovered view and every linked view, reaching every layer the covered positions reach (arcs,
   columns, the ring's position box), (3) a card with the hit's fields. A layer type that cannot is refused at
   registration (a check issue, catalogued), never silent. A drag suspends the card.
2. **The sweep**: one medium-lane test, canvas and SVG, over every story's spec and every layer instance, probing points
   on drawn marks taken from the draw list (not guessed), asserting the hit, at least one emphasis primitive in the
   hovered view and each linked view, and a non-empty card. A new layer type or layout is covered by adding nothing. A
   quick-lane representative (one story a layout).
3. **Fix what the sweep finds**, by rule, not per story: the items above close with it (each named in a test title).
4. **Legends both ways**: an entry that names items lights them in view when hovered (every layer type, canvas and
   SVG); an entry with no action shows the default cursor. Swept: a legend entry's hover gives an emphasis over its
   items in view, or no pointer.
5. **The VCF allele's card**: `pos:ref>alt`, its category and its INFO fields (the design's matrix row "a variant
   allele").

**Tests**: the sweep; the closed items each a test failing on the parent; the registration refusal. **Performance**: a
hover step at 0 B warm per layer type the sweep touches (the allocation tier's cases; the card's 2,168 B a step today:
do not raise it), a sweep run inside the medium lane's budget (seconds, not minutes: sample marks per layer).

**Where**: a worktree on `seq/hc`. Headers (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); commit at each green
step; hot files never staged (`package.json`s, `tsconfig*.json`, lockfile, `biome.json`, `.seq/`, `.v6/changelogs/`,
`storybook/packages.ts`, `packages/seq/TODO.md`): list them. No language change (a needed param stops the lane with a
report). Anything in `@molstar/core` or `@molstar/io` stops the lane. A file past 800 lines is a question. Sizes soft:
trim first, report growth with its cause; timing and allocation bounds are the gates.

**Checks**: the quick lane, `npx tsc -b packages/seq/core packages/seq/viewer packages/seq/examples`, each
`typecheck:test`, the storybook typecheck, `pnpm lint`, `pnpm format:check`, both perf tiers to completion,
`pnpm -s -C packages/seq/core size` (read the rows for FAIL), `pnpm test:medium packages/seq/`, the `.seq/_checks`
(features-planned with `--plan .seq/execution/fifth-slice/plan.md`). A fix's test fails on its parent.

**Report**: under 400 words: commits; the rule's text; what the sweep found and how each was fixed; checks with
numbers; story ids; sizes; hot-file needs; the `.v6/changelogs/` text.
