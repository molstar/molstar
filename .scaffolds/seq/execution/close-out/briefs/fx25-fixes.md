# Lane FX25: the close-out's first fix point

**Model**: strong (one item is a class: the canvas and the SVG of one spec laying out differently, already in the TODO;
classify before fixing). **Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first).
Branch `seq/fx25`. Beside you: CL9 (`packages/seq/viewer/src/{selection,input,navigator,keys}.ts` and the models'
selection and cursor, `core/src/lang/` validation, the Antibody example) and CL13 (`core/src/theme/provenance.ts`, a
folder-graph test in `core`, the lanes band's history, `rowsValue` on navigators, `seq.render.pairs.filter`). Never
edit their files; if a fix needs one, stop and report it.

**Why**: the loop's fix point after GATE 1: the open items of `.scaffolds/seq/execution/close-out/fixes.md`, each found
by the author's use, the orchestrator's look or a review, and the overnight decision on the zoom hatch
(`drift/zoom-hatch-empty.md`).

**Read first**: `fixes.md` (every open item: what is seen, where, the kind, the author's words);
`drift/zoom-hatch-empty.md`; `packages/seq/TODO.md`'s line "The canvas and the SVG of one spec lay out differently";
the specs of what you touch (core Render, viewer Header parts, ext-rondo Stubs, Hits, emphasis, words).

**Items**, a commit at each:

1. **TP53-201 first on the canvas too** (Features › Genome formats › GFF3 › Tp53): the canonical transcript is first in
   its gene's stack in the SVG and third on the canvas. Find why the two backends' layouts differ; if it is the TODO's
   class (one spec, two layouts), fix the class, not the instance, with a twin test over every story's settled spec
   comparing the two copies' lane order (medium lane), and a quick representative. Proof: the twin test; Tp53's lane
   order in both.
2. **A variant's card title names the variant** (Examples › Genomes › Views › Samples: "20001:601:f95ba8e6>C"): its
   kind and position, or its VCF ID when it has one. Proof: a card test.
3. **Stub labels never overprint** (the author's look: the ring zoomed to chr2 from the square, "chr11 · 335" and
   "chr1 · 3…" in one smear): stub labels follow the ring's label collision rule (stagger, or drop the weaker behind a
   "+n"), as the chromosome labels do. Proof: a test on the zoomed ring's label boxes not overlapping; a look.
4. **F9's pressed slot on canvas headers** (Features › Structure link › Color By This Track): the canvas header draws
   a pressed slot at rest, one at a time, as the DOM parts do. Proof: a test on the canvas draw list.
5. **The square's failure text by its own viewer** (`ext-rondo/src/viewer.ts` `HoleIssue.of` with no section): pass
   the square's ring's section or viewer. Proof: the two-rings test extended to their squares.
6. **The zoom hatch only where data exists** (the overnight decision): the hatch where this level does not draw data
   the source has; a region with no items gets the empty-well row, never the hatch. Proof: tests on both cases.

**How**: each fix's test fails on its parent (say so); a class fixed by its rule and a sweep; host neutral; module
size (§module-size); nothing per frame (the allocation tiers stay green: `pnpm -s -C packages/seq/core test:perf`);
an icon in an `<svg>` gets its icon class or a size; lazy parts awaited.

**Where**: a worktree on `seq/fx25`. Touches the files the items name and their tests and specs. Out of scope: CL9's
and CL13's files (above). Headers per the root `AGENTS.md`. Hot files never staged: `package.json`s, `tsconfig*.json`,
the lockfile, `biome.json`, `.v6/changelogs/`, `storybook/packages.ts`, `packages/seq/TODO.md`,
`packages/seq/_measure/sizes.mjs`, the work folder. Never a bare `git stash` (the stack is shared). No network.

**Checks**, green only on exit 0 with a non-zero test count: `pnpm test:unit packages/seq/ cli/seqx-render/`; `pnpm
test:medium packages/seq/`; `pnpm -s -C packages/seq/core test:perf`; `pnpm -s -C packages/seq/core size` (every row:
the ring at 45.0 of 48, the minimal at 1,090.1 of 1,095); `npx tsc -b`; `pnpm lint`; `pnpm format:check`; the
scaffold checks in `.scaffolds/seq/AGENTS.md` (stage first).

**Stop when** the checks are green, or the same check failed twice after a fix, or a step is refused by the permission
system (report it; don't work around it).

**Report** (under 450 words): commits; each item with its proof and whether its test fails on the parent; the checks
with exit codes and counts; the size rows; what to look at; drift, gaps and followups; the changelog text.
