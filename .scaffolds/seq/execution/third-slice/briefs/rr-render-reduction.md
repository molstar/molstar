# Lane RR: the render bundle reduced before the size targets move

**Task**: M3 (`../measurements.md`) measured every size row over its target: reading a spec 175.4 kB (< 160),
headless rendering 712.0 kB (< 300), the viewer, default, 1142.4 kB (< 600), viewer, minimal, 853.1 kB (< 600).
`seq/core/src/render` alone is about 370 kB in headless rendering and in the default viewer, and 237 kB in the minimal
one. The author (2026-11-20): "reduce render first, then move the size targets".

Find what each entry keeps that it does not use, and cut it, the way MV made the viewer's registries injectable
(viewer-composition.md, Creating a viewer):
- per layer type and per layout: is each a registry entry a composition can leave out, or does a barrel, a default
  table (`BuiltInLayerTypes`, the layouts table, `Renderer`'s default) or a static import keep them all?
- `"sideEffects": false` and top-level code that defeats tree shaking (module-level tables built eagerly, class
  statics);
- duplicated helpers between render modules, and generic code that belongs in `@molstar/core` (AGENTS.md, Reuse);
- reading a spec: the 15 kB over 160 (`lang`'s barrel re-exporting the editing validator, as the second slice found).

Measure with `pnpm -C packages/seq/core size` (its contributor lists) and a per-module breakdown from the bundler's
metafile. Stop and ask before a public API change beyond the injectable registries' pattern, or a default that changes
what a spec draws.

**Output**: the reductions, each a commit with the sizes before and after; tests unchanged and green (quick, medium,
perf for core and Storybook, tsc, lint, format); `_spec`s' Configuration tables kept. No target is edited by the lane.

**Model**: strong (opus). **Where**: a worktree, branch `seq/rr-render-reduction`. Do not stage `.seq/`, package.json
(but for `sideEffects`, reported for the orchestrator to apply), tsconfig, the lockfile, `.v6/`, `.claude/`.

**Report**: under 300 words: each entry's size before and after, what was cut and why it was kept before, what still
dominates, a proposed target per row (the measured size with room), any public API change, the checks, the branch.
