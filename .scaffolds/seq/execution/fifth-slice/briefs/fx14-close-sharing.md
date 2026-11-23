# Lane FX14: tracks of one kind share their height only when close

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/fx14`. Model:
Sonnet. SX is in `render/layout/` and the rail, LA3 in the grid layer, VP in `render/values.ts` and `lang/`: stay out.

**Why** (the author, 2026-11-22: "Share only when close"; fixes.md, "One tall track makes every track of its kind
tall"): `ContentSizing` (`packages/seq/core/src/render/sizing.ts`) gives every track of one kind the tallest lanes among
them, so one 40-lane Transcripts track draws thirty one-feature ranges tracks at `track.maxPx` (240 px each). The rule
(layouts.md, amended): tracks of a kind whose lanes (their content less the label lane) are within 2 lanes or 2x of each
other share the tallest of them; a track far from the rest fits its own. Clusters by sorted lane height, chained (a
track joins the cluster when it is close to the cluster's smallest), so the result does not depend on the tracks'
order. Two CDR tracks stay equal; the label lane stays each track's own.

**Task**: the rule in `ContentSizing` (its doc comment and the core `_spec`'s sizing line, README §3845 region); the
Features/Layouts story's tracks then fit (Gene tracks one lane, Transcripts at the cap). Nothing else: no config key,
no language change.

**Tests** (fail on the parent): one tall and many short of a kind (short ones fit their own); two close tracks share;
three in a chain (1, 2, 4 lanes: one cluster or two, as the rule says, stated in the title); the result independent of
order (a property over shuffles); the label lane kept per track; a track's `height` still overrides. Both perf tiers:
the sizing allocates nothing new per settle beyond today's (a Map today: keep it bounded by kinds, or reuse scratch).

**Where**: a worktree on `seq/fx14`. Headers (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); hot files never
staged (`package.json`s, `tsconfig*.json`, lockfile, `biome.json`, `.seq/`, `.v6/changelogs/`, `storybook/packages.ts`,
`packages/seq/TODO.md`). Anything in `@molstar/core` or `@molstar/io` stops the lane. Sizes soft: report growth.

**Checks**: the quick lane, `npx tsc -b packages/seq/core packages/seq/viewer packages/seq/examples`, each
`typecheck:test`, the storybook typecheck, `pnpm lint`, `pnpm format:check`, both perf tiers to completion,
`pnpm -s -C packages/seq/core size` (read the rows), `pnpm test:medium packages/seq/` (snapshots of heights that move:
list each with why), the `.seq/_checks`.

**Report**: under 250 words: commits; checks with numbers; tests whose heights moved; the `.v6/changelogs/` text.
