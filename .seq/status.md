# Status

Start here. What is built, what is open, what comes next, and where to read.

## Built

- The scaffold: `.seq/` (principles manifest, sources, designs, records), the checks wired ([AGENTS.md](AGENTS.md),
  Checks).
- `packages/seq/{core,viewer,storybook}`: package skeletons, their `AGENTS.md`, `_spec/` stubs and `TODO.md`; the
  Storybook registered as `Seq* (next)`, its design pages under `Design`. No language code yet.

## Open, by importance

The architecture is reconciled after the independent review; read [design/version-1.md](design/version-1.md) first.

1. **The first slice** (version-1, Order of work): its designs are written (language description, locations and
   queries, protocol, layout and draw lists, viewer models, examples and stories, security); start with the tree core.
2. **Provisional, to confirm on the first slice**: the root nouns ([design/building-blocks.md](design/building-blocks.md)),
   the source protocol's shape ([design/sources.md](design/sources.md)), the performance numbers
   ([design/performance.md](design/performance.md)), the UX foundations' numbers ([design/ux-foundations.md](design/ux-foundations.md)).
3. **Waiting for the author**: the proposed edge-case outcomes ([design/edge-cases.md](design/edge-cases.md)); the
   base-library changes, each when its work starts ([design/base-library.md](design/base-library.md)).
4. **To decide as the work reaches them**: worker boundaries, the asset manager's range reads, structure and map
   runtimes, accessibility, versioning under 0.x (version-1, Open).
5. **Provisional in the new designs**: the fork-or-migrate choice for meso\* and the slice's param names
   ([design/language-description.md](design/language-description.md)), the result types and default numbering
   ([design/locations-and-queries.md](design/locations-and-queries.md)), the tile size and table columns
   ([design/protocol.md](design/protocol.md)).

## Language version

Not yet versioned (`0.0.0` until the first slice runs).
