# Lane FX36: the measurements' four overs

**Model**: strong (profiling, causes across the viewer). **Base**: `packages-seq` at the orchestrator's head
(`git reset --hard packages-seq` first). Branch `seq/fx36`. No lane beside you.

**Why**: [measurements.md](../measurements.md), "Over target, real"; the author, 2026-11-24: "Yes, one perf lane". The
measurements' method (a queued `requestAnimationFrame` flushed after each input, then a microtask drain) found work
after the draw that the older method missed.

**Items**, each profiled first (Node's `--cpu-prof` over a headless twin of the story, or the perf tier's harness), its
cause named, then fixed, with a timing or allocation case failing on its parent (`itWithin` for a wall-clock bound, the
heavy tier; a quick representative), a commit each:

1. **The 200-copy split patch** (Examples › Complex › Views › Copies, then Split): the first frame 19.6 ms median
   (target: an edit under 16 ms to drawn), the main thread busy ~49 ms; the merge's task ~31 ms. Find what the patch's
   microtask tail does (a relayout twice, a re-read of the store, a history snapshot, listeners) and cut it.
2. **The square's release over the ring** (Features › Genome Ring › The Genome At The Fit: drag the square's window,
   release): ~36 ms to drawn, ~18 ms of microtasks then a ~17 ms draw (target: never over 16). The release's
   `section.include` and the ring's relayout: what is recomputed that the drag's frames already had.
3. **Search's Next and Previous** (Examples › p53 › Search, "PXXP"): ~16 ms median, 20 ms max, a 3.5 ms draw then 6–12
   ms of microtasks (target: never over 16; a step is a channel change and a reveal). Find what a step recomputes past
   the two channels.
4. **The ring's stem hover** (The Genome At The Fit): 1.8 ms median, 2.1 p95 with microtasks (target: under 2). What the
   hover's microtasks do (the card, the emphasis, listeners) and trim it under the bound with room.

**Records**: `performance.md`'s targets section and its Storybook page if a case or bound is added; the core or viewer
spec where a rule moves; `measurements.md` is the orchestrator's.

**How**: nothing per frame (the allocation tiers green); host neutral; module size (every bundle is within 1.3 to 4.8
kB of its bound: name a raise if one is needed, with its cause). Timing bounds are never loosened.

**Where**: a worktree on `seq/fx36`. Headers per the root `AGENTS.md`. Hot files never staged: `package.json`s,
`tsconfig*.json`, the lockfile, `biome.json`, `.v6/changelogs/`, `storybook/packages.ts`, `packages/seq/TODO.md`,
`packages/seq/_measure/sizes.mjs`, the work folder. Never a bare `git stash`. No network.

**Checks**, green only on exit 0 with a non-zero test count: `pnpm test:unit packages/seq/ cli/seqx-render/`;
`pnpm test:medium packages/seq/`; `pnpm test:heavy packages/seq/` (after `npx tsc -b`, alone, load under 3);
`pnpm -s -C packages/seq/core test:perf` and `pnpm -s -C packages/seq/storybook test:perf` (load under 3: `uptime`
first); `pnpm -s -C packages/seq/core size` (every row); `npx tsc -b`; `pnpm lint`; `pnpm format:check`; the scaffold
checks in `.scaffolds/seq/AGENTS.md` (stage first).

**Stop when** green, or a check red twice after a fix, or a step is refused by the permission system (report it;
don't work around it).

**Report** (under 400 words): each item's profile and cause, its fix and its numbers before and after (headless, and
how they map to the browser's); commits; the checks with exit codes, counts and the load; the size rows; the changelog
text.
