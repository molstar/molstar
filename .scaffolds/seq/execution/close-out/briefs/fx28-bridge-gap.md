# Lane FX28: a bridge's ribbons missing on the canvas after its side moves

**Model**: Sonnet (one bug, bisected, its class fixed). **Base**: `packages-seq` at the orchestrator's head
(`git reset --hard packages-seq` first). Branch `seq/fx28`. Beside you: FX25 and CL8 wait to merge (FX25: `lanes.ts`,
`ranges.ts`, headers' draw, `ext-rondo/` stub labels, `store-wants.ts`; CL8: the values layer). Never edit their
files; one you need: stop and report.

**Why**: the author, 2026-11-24 (a screenshot of Features › Across sections › One Side Zoomed In): "there is a draw
error on the canvas version of this". The canvas leaves the bridge's ribbons out over part of the zoomed side (Pp01's
in the author's view at the whole genome; about 1.55–1.58 Mb of Pp02 after a Ctrl+wheel on the orchestrator's), while
the story's SVG twin of the same view draws them. Reproduced without FX25 and CL8; it stays (not a pending tile).

**Read first**: `.scaffolds/seq/design/across-sections.md` (bridges: ribbons drawn once by ownership, the lock between
sides), `region-reads.md`; the core spec's Render › Bridges; `core/src/view/project.ts` (indexed by CL11, `cc4fab363`),
the bridge's lay and draw in `core/src/render/` (`bridge*.ts`), its region reads and caches; the story in
`storybook/src/_stories/features/across-sections/`.

**Items**, a commit at each:

0. **Reproduce and bisect**: a headless test of the story's composition that pans and zooms the zoomed side (and zooms
   it out to the whole genome) and compares the canvas draw list's ribbons with the SVG export's, item by item; then
   `git bisect` over `packages-seq`'s history to the commit that introduced it (say it; if it predates the close-out,
   say the oldest commit tested). Proof: the test failing at the head, the bisect's result.
1. **The fix, by its class**: whatever the canvas path keeps stale (an index, a cache keyed without the side's window, a
   ribbon owner chosen by the old view), fixed for every bridge kind, not this story's numbers. Proof: the test green; a
   sweep over the across-sections stories (canvas against SVG after a pan and a zoom each, medium lane), with a quick
   representative.
2. **The records**: the core spec's Bridges if a rule moves; the twin listed in `testing.md` if it is a new kind of
   check. Proof: `designs-in-specs` green.

**How**: the frame's work bounded by the screen; nothing per warm frame (the allocation tiers green); host neutral;
module size.

**Where**: a worktree on `seq/fx28`. Headers per the root `AGENTS.md`. Hot files never staged: `package.json`s,
`tsconfig*.json`, the lockfile, `biome.json`, `.v6/changelogs/`, `storybook/packages.ts`, `packages/seq/TODO.md`,
`packages/seq/_measure/sizes.mjs`, the work folder. Never a bare `git stash`. No network. Do not run the heavy lane.

**Checks**, green only on exit 0 with a non-zero test count: `pnpm test:unit packages/seq/ cli/seqx-render/`;
`pnpm test:medium packages/seq/`; `pnpm -s -C packages/seq/core test:perf` and `pnpm -s -C packages/seq/storybook
test:perf`; `pnpm -s -C packages/seq/core size` (every row); `npx tsc -b`; `pnpm lint`; `pnpm format:check`; the
scaffold checks in `.scaffolds/seq/AGENTS.md` (stage first).

**Stop when** green, or a check red twice after a fix, or a step is refused by the permission system (report it;
don't work around it).

**Report** (under 350 words): the bisect's result and the cause; commits; each item with its proof (the test failing on
its parent); the checks with exit codes and counts; the size rows; what to look at; followups; the changelog text.
