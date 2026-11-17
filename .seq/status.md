# Status

Start here. What is built, what is open, what comes next, and where to read.

## Built

- The scaffold: `.seq/` (principles manifest, sources, designs, records), the checks wired ([AGENTS.md](AGENTS.md),
  Checks).
- `packages/seq/{core,viewer,storybook}`: package skeletons, their `AGENTS.md`, `_spec/` stubs and `TODO.md`; the
  Storybook registered as `Seq* (next)`, its design pages under `Design`. No language code yet.

## Open, by importance

The architecture is reconciled after the independent review; read [design/version-1.md](design/version-1.md) first.

1. **The first slice**: planned and approved, Gate 0 settled; the loop is `execution/first-slice/loop.md`.
2. **Provisional, to confirm on the first slice**: the root nouns ([design/building-blocks.md](design/building-blocks.md)),
   the source protocol's shape ([design/sources.md](design/sources.md)), the performance numbers
   ([design/performance.md](design/performance.md)), the UX foundations' numbers ([design/ux-foundations.md](design/ux-foundations.md)).
3. **Waiting for the author**: the proposed edge-case outcomes ([design/edge-cases.md](design/edge-cases.md)); the
   base-library changes, each when its work starts ([design/base-library.md](design/base-library.md)).
4. **To decide as the work reaches them**: worker boundaries, the asset manager's range reads, structure and map
   runtimes, accessibility, versioning under 0.x (version-1, Open).

## Language version

Not yet versioned (`0.0.0` until the first slice runs).
