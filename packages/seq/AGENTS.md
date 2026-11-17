# Seq\* (next): packages/seq

Seq\* is SeqViewSpec, a language for sequence views (a tree, JSON by default) over one coordinate model, and a viewer
for it. A spec says what data there is (sequences, alignments, annotations, the maps between them) and how to show it
(sections of tracks of layers); data placed on one coordinate system shows on every system it maps to. This folder is a
new implementation, built on `@molstar/core/reactive` and `@molstar/ui`.

| Package                           | What                                                                                                                                                                            |
| --------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `core` (`@molstar/seqx-core`)     | the language, the coordinate model, the source protocol and data store, queries, the resolve passes, draw lists, the SVG and recording backends ([_spec](core/_spec/README.md)) |
| `viewer` (`@molstar/seqx-viewer`) | DOM-free models, the Canvas2D backend, the default UI ([_spec](viewer/_spec/README.md))                                                                                         |
| `storybook` (private)             | the docs pages and the examples, shown under `Seq* (next)` in the repository's Storybook                                                                                        |

- **Names**: `seqx` is temporary; the packages become `@molstar/seq-*` when they replace the current Seq\* packages.
- **The spec** of each package is its `_spec/`: what, never how; code comments cite its sections. Where code and spec
  disagree, fix whichever is wrong in the same commit.
- **Open items** are [TODO.md](TODO.md); what this version contains is the Storybook's Features page.
- **Imports**: core imports `@molstar/core` and an allowlist of other Mol\* entries (each with its reason, enforced by
  the lint once the rule is added with the first such import); the viewer adds `@molstar/seqx-core` and `@molstar/ui`.
  Never a whole package by default. Inside a package each folder is an entry; across folders name the entry, never a
  path into another folder.
- **State** lives in DOM-free models (`ReactiveModel`, `State`, `derive`), tested in Node; views are thin. No rxjs.
- **Performance**: targets for time, memory and size are the Storybook's Performance page, measured on the first slice
  and then held by tests; the draw and hover loops are allocation-free, checked by hand at checkpoints; the rest is
  plain code.
- **Tests**: the repository's lanes (root `AGENTS.md`); `pnpm test:unit packages/seq/` in the edit loop; the quick lane
  under 30 s, each folder under 2 s. Properties for codecs, patches and addresses; suites only where a second
  implementation exists; every example validated, laid out and rendered to SVG, a few with a committed record; models
  tested in Node; scaling ratios in the heavy lane; allocation checked by hand for the draw and hover loops.
- **Headers**: the Mol\* header with the year and `@author` lines.

**Development scaffold**: while this is being built, its principles, sources, designs and work records are in `.seq/` at
the repository root (start with `.seq/AGENTS.md`). That folder is temporary: never cite it from code, comments, tests or
`_spec/`.
