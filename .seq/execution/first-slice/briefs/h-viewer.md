# Lane H: the viewer and the views

**Goal**: the first views running end to end, in tests and in Storybook.
- **H1, the Canvas2D backend** (`packages/seq/viewer/src/canvas/`): it draws `seqx-core/render`'s draw lists onto a
  canvas, with an overlay canvas for the emphasis passes (hover, selection, cursor), so a hover redraws the overlay
  only. Squircles go through core's `squirclePath`. It has a text-metrics adapter over a canvas context (the text port's
  browser adapter). It is DPR-aware and reuses its buffers: nothing is allocated per frame once warm.
- **H2, composition** ([viewer-composition.md](../../../design/viewer-composition.md)):
  - `createSeqViewer(spec, target?)` and `DefaultSeqViewerSpec()` (`@molstar/seqx-viewer/default-spec`);
  - the three ports as factories with ownership: assets is `@molstar/core/util/assets`'s `AssetManager` behind a
    policy wrapper (allowed hosts, a size limit, no credentials, [security.md](../../../design/security.md)), adapted to
    core's `AssetReader`; text metrics is canvas in the browser and core's table in Node; structure is unset;
  - `createSeqContext(spec)` in `seqx-core` (headless: the registry, the config and the ports a render needs; the
    resolve passes and `renderSvg`);
  - the viewer composes lane G's models (`SpecModel`, `ViewModel`, `HoverModel`, `StoreModel`) and behaviors (`hover`,
    `navigate`) with F's `Renderer` and H1's backend. Frames and time come from `@molstar/ui/dom`'s `DomHost`;
  - pointer, wheel and key events go through the binding table; a minimal mounted UI is the canvas stack plus the
    section labels the layout already draws (no chrome beyond that in this lane);
  - `retire()` disposes what the viewer made, its owned ports included, and leaves a host's.
  - Registration is not built in this lane: the registry is a static list in the spec (atomic, undoable registration
    is a followup).
- **H3, the examples** (`packages/seq/examples`, created empty for you: `@molstar/seqx-examples`, depending on
  `@molstar/core` and `@molstar/seqx-core`):
  - the p53 dataset's fixtures. Start from the geometry and data in `packages/seq/core/src/resolve/_test/` and
    `packages/seq/core/src/render/_test/fixtures.ts`; copy them, don't import across test folders. Small committed
    files go in `packages/seq/examples/fixtures/p53/`; generated ones are made in code from a seed;
  - two views' builder modules, each exporting `build()`: `p53/tracks` (protein letters, domains, conservation) and
    `p53/protein-on-gene` (the protein following its gene through the codon map);
  - per-view tests: build, validate, resolve and lay out with no unexpected issues, and render to SVG with
    `renderSvg`; two draw-list records (one per view);
  - an assets adapter for Node that reads the package folder and refuses the network.
- **H4, the stories** (`packages/seq/storybook/src/_stories/examples/p53/`): `Seq* (next)/Examples/p53/Tracks` and
  `…/Protein on its gene`. Each composes a viewer from `DefaultSeqViewerSpec()` and the view's `build()`, with the
  canvas and the SVG side by side and the printed spec (the S-expression printing) beside them. Add a p53 index story
  listing its views. For fixtures in Storybook, write a Vite plugin in `packages/seq/storybook/vite.ts` that serves
  `packages/seq/examples/fixtures/` (the orchestrator registers it in `storybook/packages.ts` at the merge: name its
  export in the report).
- **Your two edge cases**:
  - EC-143: re-hit a resting pointer once per frame while the view moves under it.
  - Trackpad pinch: map it (a ctrl+wheel) onto the binding table's `pinch`/`wheel`.

**Configurable over fixed**: where a behavior could reasonably go more than one way, make it a config item, or an
option, with your recommendation as the default, and list each in the report.

**Read first**:
- `.seq/design/viewer-composition.md`, `ports.md`, `rendering.md` (the backends), `examples-and-stories.md`,
  `examples.md` (p53), `ux-foundations.md` (sections on hover, bindings and the grammar), `security.md` (assets),
  `testing.md`, `performance.md`;
- the code you compose:
  - the core `_spec` (`packages/seq/core/_spec/README.md`: Render, Resolve, Data);
  - the viewer `_spec` (`packages/seq/viewer/_spec/README.md`);
  - `packages/seq/core/src/render/{renderer,export,draw-lists,layer-type,text,squircle}.ts`;
  - `packages/seq/viewer/src/**`;
  - `packages/core/src/util/{assets,text-measure}.ts`;
  - `packages/ui/src/dom/{host,frame}.ts`;
- for how stories are written: one existing Storybook story under `packages/ui` and its use of `@molstar/ui/dom`, plus
  `packages/seq/storybook/src/_stories/preview.ts`.

**Model**: fast (sonnet).

**Base**: `packages-seq` at the commit named in the launch instructions. Branch `seq/h-viewer`; commit on it only, at
each small green step.

**Touches**:
- `packages/seq/viewer/src/**` and its `_test/`;
- `packages/seq/core/src/context/` (new: `createSeqContext`) and its `_test/`;
- `packages/seq/examples/src/**`, `fixtures/` and `_test/`;
- `packages/seq/storybook/src/_stories/examples/**` and `packages/seq/storybook/vite.ts`;
- the viewer and core `_spec`s (sections for what you build).

Read-only: the rest of `packages/seq/core` (report bugs you find there; fix one only if it blocks you, in a commit of
its own, and say so).

Not in scope: the selection model, chrome beyond the canvas, the registry's run-time registration, dynamic nodes,
region reads (stage 6).

**Tests**: Node, no DOM, except where a DOM test opts in with `// @vitest-environment jsdom`. The Canvas2D backend is
tested against a recording 2D context (a fake `CanvasRenderingContext2D` that logs calls). Cover:
- a hover redraws only the overlay;
- a viewer made, driven through bindings (wheel zoom, drag pan, hover) and retired leaks no listener or timer;
- `createSeqContext` + `renderSvg` equals the viewer's SVG of the same view.

A test that is slow by nature (a whole app session) goes in the medium lane with a quick representative.

**Hot files, never staged**:
- `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`, `biome.json`, `vitest.config.ts`;
- `storybook/packages.ts`, `.v6/changelogs/`;
- `packages/seq/*/package.json` and `tsconfig*.json` (ask in the report for exports, dependencies, references and
  `_test` nulls);
- anything in `.seq/` and in `packages/core/`.

If a dependency is missing (e.g. the examples package needs the viewer's Node assets adapter, or the storybook package
needs a reference), ask in the report rather than adding it.

**Headers**: the Mol\* header, `Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more
info.`, with `@author David Sehnal <david.sehnal@gmail.com>`. No mention of `packages/sequences` or `.seq/` in code,
comments or tests.

**Checks**: each its own command, green only on exit 0 with a non-zero test count.
- `pnpm test:unit packages/seq/` (each folder under 2 s, a test under 100 ms)
- `npx tsc -b packages/seq/core packages/seq/viewer packages/seq/examples`
- `pnpm -C packages/seq/core typecheck:test`, `pnpm -C packages/seq/viewer typecheck:test`,
  `pnpm -C packages/seq/examples typecheck:test`
- `pnpm -C packages/seq/storybook typecheck`
- `pnpm lint`
- `pnpm format:check`

**Stop when** the checks are green, or a check has failed twice after a fix, or a permission is refused (report the
action and why; do not work around it).

**Report**:
- the commits;
- each check with its exit code and test count, and the quick lane's time per folder;
- the hot-file changes you need: exports, dependencies, references, `_test` nulls, the Vite plugin's export name;
- the config items with their defaults;
- any place the designs were unclear or wrong, as drift: what, where, options, your recommendation;
- followups.
