# Lane FX33: the Manhattan reference-line test, flaky in a full run

**Model**: Sonnet (one flake, its cause fixed). **Base**: `packages-seq` at the orchestrator's head (`git reset --hard
packages-seq` first). Branch `seq/fx33`. Beside you: none writing code (the novice pass reads the Storybook).

**Why**: `fixes.md`, "The Manhattan reference-line test once read [0, 5.78] for [0, 7.3]": seen in FX31's first full
quick run and twice in FX32's (load 7 to 22), each green alone and on a rerun. A second failure is a fix lane's
(`loop.md`, Gates). FX32 also saw `kinds-stories` fail once in a full run under load: look at it too, the same way.

**Read first**: `packages/seq/storybook/src/_stories/examples/genomes/_test/genomes-stories.test.ts` (the Manhattan
case), the Manhattan view in `packages/seq/examples/`, CL8's values layer with `reference` (`core/src/render/values-*.ts`,
the lazy reference-lines chunk), how the layer's open domain includes the lines (the union: 0 to 7.3 expected, 5.78 the
data alone, so the lines were not yet counted); the root `AGENTS.md` (Tests: `isolate: false`, files share a worker's
module registry; `vitest.setup.ts`).

**Items**, a commit at each:

0. **Reproduce**: make it fail on demand (the file after another that loads or resets the chunk, a fake clock that
   does not wait for the chunk, a domain cached before the lines land). Say the mechanism. Proof: the command that
   fails it every time on `packages-seq`.
1. **The cause, fixed**: in the code if the domain can be read before the lines are counted (a reader of the public
   surface would see the same: a settled view must include the lines), or in the test if it reads before settling.
   Never a retry or a longer timeout alone. Proof: the reproduction green ten times; the full quick lane twice beside a
   busy process (say how).
2. `kinds-stories`: the same, or "not reproduced" with what you tried.

**How**: nothing per frame; host neutral; module size.

**Where**: a worktree on `seq/fx33`. Headers per the root `AGENTS.md`. Hot files never staged: `package.json`s,
`tsconfig*.json`, the lockfile, `biome.json`, `.v6/changelogs/`, `storybook/packages.ts`, `packages/seq/TODO.md`,
`packages/seq/_measure/sizes.mjs`, the work folder. Never a bare `git stash`. No network. Do not run the heavy lane.

**Checks**, green only on exit 0 with a non-zero test count: `pnpm test:unit packages/seq/ cli/seqx-render/`;
`pnpm test:medium packages/seq/`; `pnpm -s -C packages/seq/core size` (if code moves); `npx tsc -b`; `pnpm lint`;
`pnpm format:check`; the scaffold checks in `.scaffolds/seq/AGENTS.md` (stage first).

**Stop when** green, or a check red twice after a fix, or a step is refused by the permission system (report it;
don't work around it).

**Report** (under 250 words): the mechanism; commits; each item with its proof; the checks with exit codes and counts;
the changelog text, if a reader could see it.
