# Lane RX2: the ring's square as a window, its cells, the stem menu and keys

**Base**: `packages-seq` at the orchestrator's head after RX merges (`git reset --hard packages-seq` first). Branch
`seq/rx2`. Model: strong. **The language, alone** if a param is needed (a minor bump as RX's; say so first in the report
if one is).

**Why**: [ring-extension](../../../design/ring-extension.md), the parts RX could not ship (its report: core's square
draws no picture; no element menus; the `.` and `,` keys unbound; layer types cannot read the extension's config; a
links source that fails to open keeps its hatch), and **Decided 16** (the author, 2026-11-23: the square over a ring is
a real window). Read "The square navigator", "Hover, selection, emphasis, the keyboard", the matrix's square row, and
Decided 16. The author's words, on trying the square: "i dont understand how this control works ... its impossible to
drag or resize the box. what works is clicking once to select a grid point, and then the 2nd click enlarges that
section... very confusing".

**Task**, a green commit at each step:

1. **Measure first**: why a click inside the square's window selects a point and a second click switches
   (`viewer/src/navigator.ts`'s press and click around the `switch` branch, line ~642). A test that fails on the parent:
   one click on any cell (inside the window or out) shows that pair at once.
2. **The square over a ring is a window** (Decided 16): its box is the blocks shown (`section.include` a run on the
   diagonal, or a pair off it); dragging the box moves the run along the diagonal, an edge or corner grows or shrinks it,
   snapped to whole blocks; the cursor says which (move, resize); each gesture one patch, one undo step; the title or
   Home shows all; keyboard on the focused square too (arrows move, Shift+arrows resize, Enter on a cell). Core's square
   (GC3's), not the extension's: a square over a grid section keeps today's behaviour. Tests per gesture, its patch and
   its undo; EC-110 kept (a ring's window followed or not).
2b. **A non-contiguous selection** (the author, 2026-11-23: "would it be possible to have "non-continous" window too?
    and e.g. show chr4 and 9 only for example?"; then "yes, fold it into RX2"): Cmd/Ctrl+click a cell adds or removes
    its two chromosomes, Cmd/Ctrl+click a row or column label adds or removes that chromosome, a click on a label shows
    only it. A contiguous run keeps the one draggable box; a scattered set (chr4, chr9) is outlined cell by cell (each
    shown chromosome's diagonal cell and the cross cells between them), with no drag grips. The keyboard adds and removes
    too; Home shows all. Each change one patch, one undo step. Tests per gesture, the outline's cells, and the ring
    showing exactly the set.
3. **A navigator picture** (CN's contract, "Without CN" in the design): core's square calls a picture hook when one is
   registered; `rondo:square` returns one `rondo:cells` track over the bound ribbons' data, a cell a block pair by count
   on a square-root ramp, the diagonal included, mirrored for unsided data, sided data in both triangles; hover "chr4 ×
   chr7 · 37 links"; the window outline over the cells shown. Without the extension the lattice as today.
4. **Element menus** (a viewer opening: the label menu's shape for a layer type's hit) and the stem menu (show the
   partner pair, select partners, open in `opens`), and the keys `.` and `,` (`rondo:partner`, announced), through the
   scoped bindings RX made work for every gesture.
5. **The extension's config readable by its layer types** (the descriptor's `config` through the extension entry), so
   `rondoWith` is a convenience, not the only door.
6. **A failed source's status** on the store reader the extension reads, so the ring's hatch ends and the issue shows
   (core-wide: the smallest opening; stop with a report if it is not small).
7. **The one-chromosome look** (the author, 2026-11-23, after the meaning of the links was explained: "yes, fold the
   stem and stub idea into RX2"): on a ring of one chromosome (or a pair), the stems' trumpets dominate while the
   information is in the ribbons and stubs. Draw the stems lighter there (a `rondo.*` config value, default from a look
   in a visible pane at the fit and on chr4), so the ribbons read over them. Stubs of a count under a threshold
   (`rondo.stubMin`, default 2) fold into one "other chromosomes · n" stub per stem, its card listing the partners, and
   unfold on hover. Neck additivity holds (the folded stub's width the sum). A test per rule, a look before and after in
   the report.
8. **Stories**: the square beside the ring dragged, resized and clicked (twins replaying each gesture), the stem menu,
   the keys; light and dark, canvas and SVG.

**Where**: a worktree on `seq/rx2`. Headers (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); hot files never
staged (`package.json`s, `tsconfig*.json`, lockfile, `biome.json`, `.scaffolds/seq/`, `.v6/changelogs/`, `storybook/packages.ts`,
`packages/seq/TODO.md`, `packages/seq/_measure/sizes.mjs`): list them. Anything in `@molstar/core` or `@molstar/io`
stops the lane. A file past 800 lines is a question. Sizes soft (menus and the window's gestures lazy where they can be;
report growth with its cause); timing and allocation bounds are gates (a square drag frame 0 B warm).

**Checks**: `pnpm test:unit packages/seq/ cli/seqx-render/`, `npx tsc -b packages/seq/core packages/seq/viewer
packages/seq/examples packages/seq/ext-rondo`, each `typecheck:test`, the storybook typecheck, `pnpm lint`,
`pnpm format:check`, both perf tiers to completion, `pnpm -s -C packages/seq/core size` (read the rows), `pnpm
test:medium packages/seq/`, the `.scaffolds/seq/_checks`.

**Report**: under 400 words: commits; checks with numbers; per step what was built or stopped; story ids; sizes; the
`.v6/changelogs/` text; designs to amend (navigators, extension-openings, ring-extension).
