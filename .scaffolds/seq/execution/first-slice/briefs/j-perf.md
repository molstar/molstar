# Lane J: the stage 7 performance fixes

**Goal**: close the two misses of `measurements.md` as the author decided (2026-11-18, "as recommended"):

1. **3M bars** (pan and zoom frames 9–15 ms against 8 ms):
   - **a per-position fast path in the values layer** (`packages/seq/core/src/render/values.ts`): when the layer's
     section is an identity map over its source (or the map is linear), place positions without `toSection` and
     `axis.spans` per position; aggregate per device column directly. Same draw lists as the slow path (a test compares
     both on the same inputs, including across a map's gap and at the edges of the view);
   - **a config item for the level choice**: `RenderConfig.values.coarser` (name as fits the existing config), the
     factor by which a layer may read a level coarser than a device pixel per bin; default **2**; 1 keeps today's
     choice. Documented in the core `_spec` (Render, Configuration).
   - Target: 3M bars at 100k and 200k positions in view under 8 ms median, never over 16 ms, in the existing heavy-lane
     bounds (`render/_test/frame-time.test.ts`); tighten those bounds to the targets.
2. **Allocation**:
   - **the draw and hover paths allocate nothing once warm**: `letters.draw` (about 12 KB a frame today), the
     renderer's frame key (`renderer.ts`, a `join` of nine parts on every call, hover included: make the reuse check
     compare the parts without building a string, or a numeric key), and anything else the hover path or a redraw with
     an equal key allocates;
   - **the layout's allocation is bounded per frame, not zero**: reusable scratch arrays where it is cheap (values'
     `pieces` and `spans`, ranges' layout, sorts and `Int32Array` copies), and a stated bound measured in the
     allocation tier (`render/_test/frame.quick.ts`, `vitest.perf.config.ts`): bytes a frame independent of the data's
     size at a fixed screen (3M bars no longer 44 MB);
   - the allocation cases assert: draw and hover 0 bytes warm (with the harness's tolerance), layout under a bound you
     state (per scenario, from the measurement, with headroom).

**Configurable over fixed**: any new knob (the level factor, scratch sizes) is a config item with the recommendation
as default; list each in the report.

**Read first**:
- `.scaffolds/seq/execution/first-slice/measurements.md` (the numbers and where the time goes);
- `.scaffolds/seq/design/performance.md`, `layout-and-draw-lists.md`, `rendering.md`, `testing.md`;
- the core `_spec` (`packages/seq/core/_spec/README.md`, Render);
- the code: `packages/seq/core/src/render/{values,letters,ranges,renderer,axis,layer-data,draw-lists,config}.ts`,
  `render/_test/{frame-time.test.ts,scenarios.ts,frame.quick.ts}`, `core/vitest.perf.config.ts`, and how the viewer
  calls the renderer on hover (`packages/seq/viewer/src/viewer.ts`, `canvas/`).

**Model**: strong (opus).

**Base**: `packages-seq` at the commit named in the launch instructions. Branch `seq/j-perf`; commit on it only, at
each small green step.

**Touches**: `packages/seq/core/src/render/**` and its `_test/`; `packages/seq/viewer/src/**` only where the hover or
draw path allocates; the core and viewer `_spec`s.

Not in scope: new features, the data store's summaries' factor, browser measurements.

**Tests**: the fast path's lists equal the slow path's (a property test over maps, views and levels); the level factor
(1 vs 2) choosing as stated; the heavy-lane bounds at the targets; the allocation cases as above. Draw lists stay pure
(same inputs, same lists). A test slow by nature goes in the medium or heavy lane with a quick representative.

**Hot files, never staged**:
- `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`, `biome.json`, `vitest.config.ts`;
- `storybook/packages.ts`, `.v6/changelogs/`;
- `packages/seq/*/package.json` and `tsconfig*.json` (ask in the report);
- anything in `.scaffolds/seq/` and in `packages/core/`.

**Headers**: the Mol\* header, `Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more
info.`, with `@author David Sehnal <david.sehnal@gmail.com>`. No mention of `packages/sequences` or `.scaffolds/seq/` in code,
comments or tests.

**Checks**: each its own command, green only on exit 0 with a non-zero test count.
- `pnpm test:unit packages/seq/` (each folder under 2 s, a test under 100 ms)
- `pnpm test:heavy packages/seq/` (the bounds; `--disableConsoleIntercept` prints the numbers)
- `pnpm -C packages/seq/core test:perf` (allocation; about 10 minutes, run it in the background)
- `npx tsc -b packages/seq/core packages/seq/viewer packages/seq/examples`
- `pnpm -C packages/seq/core typecheck:test`, and the same for `viewer` and `examples`
- `pnpm lint`
- `pnpm format:check`

**Stop when** the checks are green, or a check has failed twice after a fix, or a target cannot be met without
changing a design (report the measurement and the options), or a permission is refused (report the action and why; do
not work around it).

**Report**:
- the commits;
- each check with its exit code and test count; the before and after numbers per scenario (frame median and p95,
  bytes a frame for draw, hover and layout);
- the hot-file changes you need;
- the config items with their defaults;
- drift: what, where, options, your recommendation;
- followups.
