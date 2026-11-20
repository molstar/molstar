# Gaps: for the author's review (the fourth slice, overnight)

Decisions made while the author is away (loop.md, "Running while the author is away"), and public API changes made
under the author's standing permission (`state.md`, Decided). Each a yes or no question with what was done.

## Public API changes (allowed: `packages/seq`, core's `repr`, `@molstar/ui`)

- **ST2**: `JumpKind` gains `zoom` and `move`; `view.jump(move, record?, kind = 'move', location?)`;
  `NavigateContext.location`; core `toUnits`/`fromUnits` over an empty overlay read the model alone. Keep? Recommended
  yes (the review: unnamed jumps were recorded as `focus`).
- **ST3**: `StructureSide.onHover`'s callback takes the event's time (`at?: number`); `linkStructure` drops a side hover
  older than the pointer's last own hover (the review: a late non-null 3D hover took the pointer's). Keep? Recommended
  yes (EC-141 covered for hovers, not only leaves).

- **HD**: `withAdapter` is generic over the registry's entry type (a viewer's registry passes through); the examples'
  chain tables gain `ChainResidue.b`, `ChainTable.model`, `HostModel`. Keep? Recommended yes.
- **HD's drift**: no section can show a source's complex (`drift/host-complex-sections.md`): keep static sections with
  an issue, sections over a source's complex, or the source bringing `values` nodes? Recommended: sections over a
  source's complex in the next slice's view-model work, the issue until then.

## Decisions as recommended

- **ST2**: the link sends every jump that names a Location (`focus`, `reveal`, and back/forward over them); `fit`,
  `zoom`, `move` never go out. Keep `reveal` going out? Recommended yes (a reveal is a navigation to a place the
  structure can show).

## For the author (no default taken)

- **T-DM1's CCD ids** (FX3, item 6, open): the repo holds no CCD data for DM1 or the SMCC/MCC linker. Fetch their
  chemical component pages from the PDB once (two small files, recorded with their licence), or give the ids? A
  fetch is asked before it is made.
- **The viewer's quick-lane folder still over 2 s** (FX3, item 1): 2.78 s, its time imports (28%) and jsdom setup
  (20%) across 65 files. Fewer jsdom files, a lighter DOM environment, or a viewer-folder budget of 3 s?
  Recommended: fewer jsdom files first (the next fix lane), the budget unchanged.

- **The link's hover allocation**: a 2OCJ hover allocates 11.3 kB through the link (9.2 kB in `toUnits`, 1.7 kB per
  lit unit; was 31.8 kB), against performance.md's 0 B for warm hover loops; the target and the hover's Location are
  new each hover. Set a target for the link's path (performance.md names none), or ask for pooled results? Recommended:
  a bound per lit unit (2 kB) in performance.md, pooling later if a large assembly shows it.
