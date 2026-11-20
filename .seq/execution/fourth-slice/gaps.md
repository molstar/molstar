# Gaps: for the author's review (the fourth slice, overnight)

Decisions made while the author is away (loop.md, "Running while the author is away"), and public API changes made
under the author's standing permission (`state.md`, Decided). Each a yes or no question with what was done.

## Public API changes (allowed: `packages/seq`, core's `repr`, `@molstar/ui`)

- **ST2**: `JumpKind` gains `zoom` and `move`; `view.jump(move, record?, kind = 'move', location?)`;
  `NavigateContext.location`; core `toUnits`/`fromUnits` over an empty overlay read the model alone. Keep? Recommended
  yes (the review: unnamed jumps were recorded as `focus`).

## Decisions as recommended

- **ST2**: the link sends every jump that names a Location (`focus`, `reveal`, and back/forward over them); `fit`,
  `zoom`, `move` never go out. Keep `reveal` going out? Recommended yes (a reveal is a navigation to a place the
  structure can show).

## For the author (no default taken)

- **The link's hover allocation**: a 2OCJ hover allocates 11.3 kB through the link (9.2 kB in `toUnits`, 1.7 kB per
  lit unit; was 31.8 kB), against performance.md's 0 B for warm hover loops; the target and the hover's Location are
  new each hover. Set a target for the link's path (performance.md names none), or ask for pooled results? Recommended:
  a bound per lit unit (2 kB) in performance.md, pooling later if a large assembly shows it.
