# Lane OP: PNG headless and the render CLI

**Goal**: the plan's OP and Gate 0's decision (the author, 2026-11-18): **PNG by rasterizing the SVG** (the same
picture as SVG, headless) and **the render CLI in `cli/` as `seqx-render`**, under the headless asset policy (EC-126),
with its cold-start target ([performance](../../../design/performance.md)). Read the "Done when" list in
`.seq/execution/third-slice/plan.md` ("PNG renders headless with the same picture as SVG; the render CLI turns a spec
file into SVG or PNG under the headless asset policy") and `packages/seq/core/src/render/png.ts` (a DOM-free PNG
encoder used for the SVG's images). Look at how `cli/mvs-render` rasterizes and packages before choosing the
rasterizer (reuse a dependency the workspace already has; a new dependency is a question in your report, not added).

Build:
1. `renderPng(spec | view, options)` in the headless entry (`@molstar/seqx-core/render` or where the headless SVG
   export lives): the SVG export rasterized at a scale and dpr, both themes, the images and the glyphs as SVG draws
   them. A test that the PNG's pixels match the SVG's picture (a few probe pixels per layer in light and dark) and that
   a story's export is deterministic.
2. `cli/seqx-render`: `seqx-render <spec.json|spec.seq> -o out.svg|out.png [--width --height --dpr --theme
   light|dark --view <id>]`, reading assets under the headless policy (local files under the spec's folder only; no
   network unless a flag allows it, EC-126), its usage text, exit codes, errors as the spec's issues printed. Tests run
   it built (a heavy-lane test, as the other CLIs: `pnpm test:heavy` after `npx tsc -b`) and in-process (quick).
3. The cold-start target: the time from process start to a written SVG and PNG for p53's tracks, measured, an
   `itWithin` bound in the heavy lane.
4. The core `_spec`'s Render › Backends and a `cli/seqx-render/README.md`.

**Not in scope**: anything else. Other lanes run beside you; you touch only the headless export path, `cli/`, and
their tests.

**Model**: fast (sonnet). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B seq/op-png-cli
<that commit>`.

**Checks**: `pnpm test:unit packages/seq/ cli/seqx-render`; `npx tsc -b`; `pnpm -C packages/seq/core typecheck:test`;
`pnpm lint`; `pnpm format:check`; `pnpm test:heavy cli/seqx-render` after `npx tsc -b`. Exit codes explicit (`; echo
exit=$?`), non-zero test counts.

**Headers**: a new file has the Mol\* header `Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE
file for more info.` with `@author David Sehnal <david.sehnal@gmail.com>`; a modified file keeps its header, year range
to 2026, that author added if missing. No mention of `.seq/` in code.

**Hot files, never staged**: `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`,
`biome.json`, `vitest.config.ts`, `storybook/packages.ts`, `.v6/changelogs/` (write the entry's text in the report),
every `package.json` and `tsconfig*.json` (the new `cli/seqx-render/package.json` and `tsconfig.json` included: write
them, leave them unstaged, list them and the inventory and root-tsconfig entries in your report), anything in `.seq/`.

**Stop when** green, or a check fails twice after a fix, or a new dependency is needed.

**Report**: commits; checks; the cold-start time; the hot-file entries; drift.
