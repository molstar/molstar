# Lane CL11: the drag rule

**Model**: strong (a class: every drag kind's work per move, across the selection, emphasis, hover, the SVG twin and
linked panels). **Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch
`seq/cl11`. Beside you: CL7 (patches on streamed things: `core/src/{resolve,lang,data}/` refs and the viewer's
`models/{address,spec}.ts`). Never edit its files; if you must, stop and report.

**Why**: the fifth slice's open class, carried to this plan (`plan.md`, "The drag rule"): a drag's work is bounded by
the selection, not the screen. The author, 2026-11-21: selecting a rect on the tRNA dot plot "is really laggy", a range
on its navigator too, "similar lag to selecting on the genome ring with the contigs"; and on the PAE plot "showing the
hover when dragging … makes the view feel really laggy" (the card now hides during drags; the rest is open).

**Read first**: `.scaffolds/seq/execution/fifth-slice/fixes.md`, the item "A drag's work is bounded by the screen, not
by the selection" (what was read: every move runs a full hit test and writes a new selection value, and each value
re-runs every derived view) and "The hover card follows a drag on the grid"; `.scaffolds/seq/design/performance.md`
(draw lists bounded by the screen, allocation-free warm frames); `ux-foundations.md` (selection, drags); the viewer
spec's Selection, Input and History; the code in `packages/seq/viewer/src/{select,input}.ts`,
`models/selection.ts` (`rectDragTo`), the emphasis readers in `core/src/render/`, the navigator's strip
(`navigator-blocks`), the ring's arcs in `ext-rondo`; the perf tier's cases (`*.quick.ts`) and
`navigator-blocks.quick.ts`'s 256 kB bound.

**Items**, a commit at each:

0. **Profile first**: a perf-tier scenario a drag kind (a rect on the 4,000 × 4,000 dot plot; a range on its ruler; a
   sweep across 12 chromosomes on the 3,000-contig ring; a box on the PAE grid); report each kind's work per move
   (time, bytes, what re-ran) before changing code. A finding against the design stops the lane: a drift note.
1. **The selection written at most once a frame**: pointer moves coalesced to the frame, not the event; the full
   selection resolved once on release (one history step, as now). Proof: a test counting writes over many moves in one
   frame (fake clock).
2. **Derived work in proportion to the screen**: during a drag, emphasis derived for the visible part only and
   incrementally where the shape allows (a rect's growth, not the whole rect again); no hover lookup or card; the SVG
   twin and linked panels updated on the frame. Proof: tests per drag kind on what re-ran.
3. **The 3,000-contig strip**: the navigator's strip not laid out again per move (FX8's 110 kB a frame). Proof: its
   perf case's bound lowered to what the fix measures.
4. **The bounds**: each drag kind's warm frame at 0 bytes in the perf tier; a wall-clock bound per kind in the heavy
   lane under `itWithin`. Proof: the cases; the tier green.
5. **The specs and designs**: the rule in `performance.md` and `ux-foundations.md` (its Storybook Design pages), the
   viewer spec's Selection; the fifth slice's two items closed by pointer. Proof: `designs-in-specs` green.

**How**: host neutral (the frame is the viewer's draw scheduling, never a global rAF the host can't stop); reuse the
viewer's frame and the store's version; nothing per frame; module size (§module-size); the behavior on release
unchanged (tests that pin a drag's final selection stay green unchanged; a changed expectation says why).

**Where**: a worktree on `seq/cl11`. Out of scope: CL7's files. Headers per the root `AGENTS.md`. Hot files never
staged: `package.json`s, `tsconfig*.json`, the lockfile, `biome.json`, `.v6/changelogs/`, `storybook/packages.ts`,
`packages/seq/TODO.md`, `packages/seq/_measure/sizes.mjs`, the work folder. Never a bare `git stash`. No network.

**Checks**, green only on exit 0 with a non-zero test count: `pnpm test:unit packages/seq/ cli/seqx-render/`; `pnpm
test:medium packages/seq/`; `pnpm -s -C packages/seq/core test:perf` and `pnpm -s -C packages/seq/storybook test:perf`;
the heavy cases you add, alone; `pnpm -s -C packages/seq/core size` (every row: the default viewer at 1,551.7 of
1,560, the minimal at 1,092.1 of 1,095); `npx tsc -b`; `pnpm lint`; `pnpm format:check`; the scaffold checks in
`.scaffolds/seq/AGENTS.md` (stage first). Read timings only with the machine quiet; say the load.

**Stop when** green, or a check red twice after a fix, or a step is refused by the permission system (report the
action and why; don't work around it, not even by another tool).

**Report** (under 450 words): commits; item 0's profile per drag kind, before and after; each item with its proof; the
checks with exit codes and counts; the size rows; what to look at (the four drags, light and dark); drift, gaps and
followups; the changelog text.
