# Packages

## Decision

Accepted as recommended (the author, 2026-11-17: "others as recommended"; the alternative not taken was one
`@molstar/seq` package with subpath entries, rejected because the core must stay free of DOM and UI dependencies that
the viewer needs).

| Folder | Package | What | May import |
|---|---|---|---|
| `packages/seq/core` | `@molstar/seqx-core` | the language (description, validation, builder, printing, patches), the coordinate model and maps, the source protocol and the core runtimes, queries, the pure passes (resolve, layout, draw lists), the SVG and recording backends, `renderSvg` and the render CLI ([rendering.md](rendering.md)) | `@molstar/core`, and an allowlist of other Mol\* entries |
| `packages/seq/viewer` | `@molstar/seqx-viewer` | DOM-free models (viewer, viewport, hover, selection, behaviors, commands, editing), the Canvas2D backend, the default UI | `@molstar/seqx-core`, `@molstar/core`, `@molstar/ui` |
| `packages/seq/nightingale` (later) | `@molstar/seqx-nightingale` (private) | the reference adapter: a subset of the resolved view drawn by Nightingale, a check that the resolved view is general ([rendering.md](rendering.md)) | `@molstar/seqx-core`, Nightingale |
| `packages/seq/examples` | `@molstar/seqx-examples` (private) | the datasets' fixtures and one builder module per view ([examples-and-stories.md](examples-and-stories.md)) | `@molstar/seqx-core` |
| `packages/seq/storybook` | `@molstar/seqx-storybook` (private) | docs pages, examples that run, stories of the viewer and its parts; the application-layer showcases (services, structure files and the Mol\* link, addressing) | all of the above, and Mol\*'s packages for the showcases |

- **Imports from other Mol\* packages, carefully and selectively** (the author, 2026-11-17: "we can import from other
  mol* parts now, but very carefully and seletively"): beyond `@molstar/core` and `@molstar/ui`, a package imports only
  entries on its allowlist, each with its reason (first: `@molstar/io`'s MessagePack codec, for the MessagePack option),
  enforced by the lint's restricted imports (the rule for `packages/seq` is added with its first import). A whole package is never allowed, only entries; a generic piece that many
  would want goes to `@molstar/core` instead (the tree core, [tree-core.md](tree-core.md)).
- **Names**: `seqx` while the previous `@molstar/seq-*` packages exist; renamed to `@molstar/seq-*` in the change that
  deletes them.
- **The core's folders follow the pipeline** ([composition.md](composition.md)): `lang/` (the description, the builder),
  `resolve/`, `view/`, `data/` (the protocol, the store, the runtime contract), `render/` (layout, draw lists, the SVG and
  recording backends).
- **Inside a package**, each folder is an entry (`@molstar/seqx-core/model`); across folders a module names the entry,
  never a path into another folder; the folder graph is a test (`_test/layers.test.ts`), and it has no cycles.
- **Storybook** is registered in `storybook/packages.ts` under the group `Seq* (next)` until the rename.
- **Tests** run in the repository's one Vitest run, in its three lanes; no per-package Vitest configs until a tier
  needs one.
- **Headers**: `Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.` and the
  `@author` lines.
