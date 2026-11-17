# Status

Start here. What is built, what is open, what comes next, and where to read.

## Built

- The scaffold: `.seq/` (principles manifest, sources, designs, records), the checks wired ([AGENTS.md](AGENTS.md),
  Checks).
- `packages/seq/{core,viewer,storybook}`: package skeletons, their `AGENTS.md`, `_spec/` stubs and `TODO.md`; the
  Storybook registered as `Seq* (next)` with the docs pages (Introduction, Building blocks, Format, Sources). No
  language code yet.

## Open, by importance

- **Next work: the tree core** ([design/tree-core.md](design/tree-core.md)): meso\*'s tree core into `@molstar/core`, its
  own work id; the entry is `@molstar/core/repr`; units wait.

0. **The scope** ([design/features.md](design/features.md), Provisional): review the `deferred` rows first.

1. **The root nouns** ([design/building-blocks.md](design/building-blocks.md), Provisional): the author's review.
2. **The language description** (`defineLanguage`): its shape, validation with every issue, the builder's types, the
   printer; the first slice is the spec, data and one view kind.
3. **The source protocol** ([design/sources.md](design/sources.md), Provisional): settled by inline `seq` data and
   FASTA by region.
4. **The coordinate model**: systems, maps with direction on every segment, paths, projection.
5. **The first example end to end**: the p53 text view, through resolve, layout, draw lists, and both the SVG and
   Canvas2D backends.
6. **Performance targets** ([design/performance.md](design/performance.md), Provisional): the author confirms the
   numbers; each gets its test as its area is built.
7. **UX foundations** ([design/ux-foundations.md](design/ux-foundations.md), Provisional): settled by the foundation
   stories and the grammar's conformance tests.
8. **Edge-case outcomes** ([design/edge-cases.md](design/edge-cases.md)): the proposed changes and flaw fixes await
   the author's confirmation; every other item gets its outcome when its area is built.

## Next

The language description, design first: a draft in `execution/<work-id>/`, agreed into `design/language.md`.

## Language version

Not yet versioned (`0.0.0` until the first slice runs).
