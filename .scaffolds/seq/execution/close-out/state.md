# State: close-out (updated 2026-11-24, by the orchestrator)

Branch: `packages-seq` Head: `cc4fab363` (CL11 merged; its records after)

## Where it stands

**Overnight run** (the author, 2026-11-24: "keep going with the loop, overnight run so we will reconcile things in
the morning"): the loop's "Running while the author is away". A decision a gate asks for is made as recommended and
listed in `gaps.md` under "Decided overnight"; the stop conditions still stop a lane; no public API change beyond a
design and no language change is decided; the use gate is deferred. A look the pane cannot show visibly is done in the
DOM and pixels and listed in `gaps.md` for the author's morning look.

Stage 0 and GATE 0 are closed; waves 1 and 2 are merged and GATE 1 passed. Merged since: CL13 (`27470de3e`), CL9
(`19733d6c3`, the default viewer's bound raised to 1,560 kB under the author's leave), CL7 (`1e499c655`), CL11 (`cc4fab363`, the minimal viewer's bound raised to 1,100 kB under the author's leave).
Sizes at the head: spec 260.9/265, headless 979.0/980, default 1,557.2/1,560, minimal 1,096.4/1,100, ring 45.0/48.

## In progress

- FX25 stopped on a permission refusal (`gaps.md`), its worktree kept with items 1 and 2 uncommitted; waits for the
  author.

## Next

1. CL12 (`briefs/cl12-m5-leftovers.md`), alone; then GATE 2 waits for the author's word on CL5, CL6, CL8 and FX25.
2. The legend's eye chunk flake (fixes.md) to the next fix lane.
3. CL5, CL6 and CL8 wait for the author's word on their drafts (`drafts/search.md`, `track-groups.md`,
   `values-marks.md`); FX25 for the refusal. GATE 2 after wave 4.

## Blockers and open questions

`gaps.md`: FX25's refusal and its allele-title question; the three drafts (two are language changes); `rowsValue` over
a lazy alignment; the overnight decisions to confirm (the zoom hatch; a 0.x minor read as a newer minor; the strip's layout per hull, left for version 1). The headless
and default bounds have 1.0 and 2.8 kB left: CL12 keeps new code lazy, or the orchestrator trims or raises (the
author's leave). At the close: remove the locked, clean worktrees (CL4's, CL9's, CL13's, the
reviewers') and their branches.

## Verify

The loop's Checks; the scaffold checks were last green at `cc4fab363`.
