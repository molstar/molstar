# Lane I: by region (tables, dynamic nodes, nested sources)

**Goal**: data read region by region, end to end, as stage 6 of `plan.md` lists it.
- **I1, tables by region:**
  - a synthetic asynchronous `seq` source with an injectable clock (the store's `StoreClock`): coverage, levels, tiles,
    pending vs empty, errors, cancellation, eviction, and summaries (value pyramids built per level, so a values layer
    zoomed out reads a coarse level, not level 1);
  - planning through the static codon map (a protein view wanting the gene's region);
  - FASTA with its index (`.fai`), read whole (Gate 0: range reads wait for bigWig), as a `FormatRuntime`.
- **I2, dynamic nodes** ([dynamic-nodes.md](../../../design/dynamic-nodes.md)):
  - an answer carries tables and nodes; identity is (source, key);
  - nodes are validated by the language description;
  - a transcript arriving with its rows' `map` gets an item coordinate system and a codon map on demand;
  - the overlay as tier 2 ([composition.md](../../../design/composition.md)): the static model and resolved view stay
    pure, and the overlay adds the store's nodes with a version that is part of the draw lists' key;
  - planning through a dynamic map;
  - node lifetime against eviction: a node stays while something holds it.
- **I3, nested sources and templates:**
  - the transcript's domains as a nested source, requested only when needed, never past a depth limit (EC-179);
  - shown by a track template (`for_each` over a source's nodes, by kind, with order and filter).
- **I4, the story**: `Seq* (next)/Examples/p53/Gene by region`. p53's gene region read region by region from the
  synthetic source with a delay, transcripts arriving as dynamic nodes, each with its domains via the template; canvas
  and SVG side by side, as the other p53 stories do (`example-root.ts`). Add the view's builder module and fixtures to
  `packages/seq/examples` with its tests (build, validate, resolve, lay out, SVG), and list it in the p53 index story.

**Configurable over fixed**: where a behavior could reasonably go more than one way (summary levels, nesting depth,
node retention, template defaults), make it a config item or an option with your recommendation as the default, and
list each in the report.

**Read first**:
- `.scaffolds/seq/design/dynamic-nodes.md`, `data-store.md`, `protocol.md`, `sources.md`, `composition.md`,
  `layout-and-draw-lists.md` (the overlay's version), `maps.md` (codon maps, paths), `edge-cases.md` (search for the
  EC numbers the designs cite), `testing.md`, `performance.md`;
- the code: the core `_spec` (`packages/seq/core/_spec/README.md`: Data, Resolve, Resolved view, Render, Context),
  `packages/seq/core/src/{data,resolve,view,render,context,lang}`, `packages/seq/examples/src/p53/`, and
  `packages/seq/storybook/src/_stories/examples/`;
- `.scaffolds/seq/execution/first-slice/followups.md` lists what earlier lanes deferred to you (the summary pyramid, the overlay
  type).

**Model**: strong (opus).

**Base**: `packages-seq` at the commit named in the launch instructions. Branch `seq/i-regions`; commit on it only, at
each small green step.

**Touches**:
- `packages/seq/core/src/**` (data, resolve, view, render, lang, context) and their `_test/`;
- `packages/seq/viewer/src/**` where the viewer must hold the overlay;
- `packages/seq/examples/**`;
- `packages/seq/storybook/src/_stories/examples/**`;
- the core and viewer `_spec`s.

Not in scope: bigWig and range reads, services, structure links, the selection model, the agent tools.

**Tests**:
- Node, with a fake clock, never real waiting;
- the synthetic source's states (pending, empty, error, cancelled, evicted);
- a pyramid level chosen by zoom;
- a protein view planning the gene's region through the codon map;
- FASTA with `.fai` round trip and a malformed index refused;
- a dynamic transcript's item system and codon map;
- a node surviving eviction while held and dropped after;
- the overlay's version changing the key, the static tiers unchanged;
- a nested source not requested until needed, and the depth limit;
- the template's order and filter.

A test slow by nature goes in the medium lane with a quick representative.

**Hot files, never staged**:
- `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`, `biome.json`, `vitest.config.ts`;
- `storybook/packages.ts`, `.v6/changelogs/`;
- `packages/seq/*/package.json` and `tsconfig*.json` (ask in the report for exports and `_test` nulls);
- anything in `.scaffolds/seq/` and in `packages/core/`.

**Headers**: the Mol\* header, `Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more
info.`, with `@author David Sehnal <david.sehnal@gmail.com>`. No mention of `packages/sequences` or `.scaffolds/seq/` in code,
comments or tests.

**Checks**: each its own command, green only on exit 0 with a non-zero test count.
- `pnpm test:unit packages/seq/` (each folder under 2 s, a test under 100 ms)
- `npx tsc -b packages/seq/core packages/seq/viewer packages/seq/examples`
- `pnpm -C packages/seq/core typecheck:test`, and the same for `viewer` and `examples`
- `pnpm -C packages/seq/storybook typecheck`
- `pnpm lint`
- `pnpm format:check`

**Stop when** the checks are green, or a check has failed twice after a fix, or a permission is refused (report the
action and why; do not work around it).

**Report**:
- the commits;
- each check with its exit code and test count, and the time per folder;
- the hot-file changes you need;
- the config items with their defaults;
- drift: what, where, options, your recommendation;
- followups.
