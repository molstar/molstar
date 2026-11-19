# Lane TH: themes and legends

**Goal**: build [themes](../../../design/themes.md) as agreed, except its small layers' drawing (lane SL, beside you):
its Decision (the author's decisions first: the writer says what is colored by what, a layer's default theme is the
language's, plain defaults, a field continuous unless `:mid` or `diverging`), its matrix cells, its Configuration table
(the "where each choice lives" column copied into the `_spec`s), its Performance targets, and its Consequences, which
are this lane's checklist. Tensions stay tensions.

In short:
- `@molstar/core`, non-breaking (the agreed surface: exactly the design's list): `ValuePaletteLists` a lazy getter;
  `scale`, `mid` and `unlisted: 'hash'` as optional fields; `util/color/spaces/oklab.ts`; `util/color/cvd.ts`;
  `util/color/residue-schemes.ts`; the hydrophobicity scale moved from `@molstar/model` with a re-export under its name;
  each with its test. Nothing kept privately in Seq\*; anything past this list stops the lane (ask in the report);
- the language: the `seq.theme` sublanguage and the channel slots on `layer` (`color`, `ink`, `size`, `opacity`,
  `line`, `arrow`, `shape`), the `scale` kind, `values.encoding` and `scale`, the builder's `theme.*` and
  `b.views.scale`, version `0.4.0`, the issues;
- `theme/` as its own entry `@molstar/seqx-core/theme` (a package `exports` change: ask in the report), the evaluator,
  palettes over Mol\*'s, the scheme registry and built-in schemes, measures, legends' data and converters;
- resolve and view: themes in the resolved view, shared open domains unioned at settle;
- render: channels read by the existing layer types (ranges, letters, values bars) and a channel accessor the small
  layers read (SL builds against today's color path; you wire them at your merge, which lands second), LUTs, the new
  tokens in both tables (`seq-window`, `seq-grip`, `seq-wash`, `seq-centromere`, the break ink, those the design
  names), `layoutLegends`, the export's `legend` option and the view's optional export params;
- the viewer: the legends state, the `Legend` part (`mountLegend(element, viewer, { of })`, a part the host places),
  `SeqViewerSpec.ui.legends`, the legend behavior provider, the scheme picker, the commands (`layer.filter` with its
  "n filtered" chip, decided), `theme.cvdSafe` as a reader setting, the legend's highlight rows with a clear running
  `highlight.clear` (built by C0);
- the examples: the p53 themes view and its editing view, three ways, a headless twin.

Not in scope: drawing the `ruler`, `secondary`, `variants` layers, values' line and area, ranges' line styles and
arrowheads (SL); navigators' bands drawing (NV); layouts.

**Coordination with SL**: SL adds those layer types' kinds and render code; you own the theme slot, the channel
accessor, `render/tokens.ts`, and the `lang/` theme files. Where both touch one file (a layer type's description, the
render registry), keep your edits to the theme parts; the orchestrator merges SL first and you rebase onto it before
your report if it has merged (`git rebase` onto the new `packages-seq`), wiring SL's layers to the channels.

**Read first**: [themes](../../../design/themes.md) whole; [ux-foundations](../../../design/ux-foundations.md) section
8; [features](../../../design/features.md); [base-library](../../../design/base-library.md); `packages/seq/AGENTS.md`;
the core and viewer `_spec`s; `packages/core/src/util/color/`, `util/legend.ts`; the code they name.

**Model**: strong (opus).

**Base**: `packages-seq` at `1e7e87870`. Your worktree may be cut from master: with a clean tree, first
`git checkout -B seq/th-themes 1e7e87870`, and check `git merge-base --is-ancestor 1e7e87870 HEAD`.

**Touches**: `packages/core/src/util/color/**` (and the hydrophobicity move in `packages/model`, its re-export only),
`packages/seq/{core,viewer,examples}/src/**`, `packages/seq/storybook/src/_stories/**` (not the Design pages), the
`_spec`s, `packages/seq/TODO.md`.

**Two consumers**: one theme evaluated by the canvas and by `renderSvg` (each mark's color its datum's, each legend entry
with a mark); a core addition read by Seq\*'s scheme and by a Mol\* value palette test.

**Tests**: the design's Consequences, Tests, all of it (property tests, the conformance suite over schemes and tokens,
drawing checked against the data, the size channel never moving a hit box, the allocation tier and the 1.1× bound).

**Performance representative**: the design's Performance targets, each as a quick case and its `itWithin` bound (best of
three): a theme rebuild, a LUT, the draw with channels against without (1.1×).

**Whose each choice is**: the design's table; nothing configurable beyond it.

**Minimal code**: reuse `@molstar/core`'s color code (`Color.luminance`, `Color.contrast`, the value palettes, the
lists); `render/tokens.ts` drops its own `luminance`; no private copies; sizes reported (they are guides: the author,
"don't overdo it on the limits ... but obviously try to minimize the amount of code").

**Checks**: each its own command, green only on exit 0 with a non-zero test count.
- `pnpm test:unit packages/core/src/util/color packages/model packages/seq/`
- `pnpm test:heavy packages/seq/` (never while another heavy run is going; if one is, wait)
- `npx tsc -b packages/core packages/model packages/seq/core packages/seq/viewer packages/seq/examples`
- `pnpm -C packages/seq/core typecheck:test`, the same for `viewer` and `examples`, and `pnpm -C packages/seq/storybook typecheck`
- `pnpm lint`, `pnpm format:check`
- `pnpm -C packages/seq/core size` (report the numbers)

**Headers**: a new file has the Mol\* header `Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for
more info.` with `@author David Sehnal <david.sehnal@gmail.com>`; a modified file keeps its header, its year range
extended to 2026 and that author line added if missing. No mention of `packages/sequences` or `.seq/` in code, comments,
tests or `_spec`s; code cites its own `_spec` sections.

**Fixes land with a test shown failing on the fix's parent**: the report names the command.

**Hot files, never staged**:
- `pnpm-lock.yaml`, `scripts/workspace/inventory.json`, the root `tsconfig.json`, `biome.json`, `vitest.config.ts`;
- `storybook/packages.ts`, `.v6/changelogs/` (write the changelog entries' text in your report: one for `@molstar/core`,
  one for Seq\*);
- every `package.json` and `tsconfig*.json` (ask in the report: the `theme` entry, any `_test` export nulls);
- anything in `.seq/`.

**Stop when** the checks are green, or a check has failed twice after a fix, or the work needs a public API past the
design's list, or a target cannot be met without changing a design (report the measurement and the options), or a
permission is refused (report the action and why; do not work around it).

**Report**: the commits; each check with its exit code and test count; the sizes and the performance numbers; the
hot-file changes you need; the changelog entries' text; config items with their defaults; the matrix cells filled;
drift (what, where, options, your recommendation); followups. You cannot see a browser: name the stories to look at,
and what to look for in each (light and dark, the schemes gallery under simulated deficiencies).
