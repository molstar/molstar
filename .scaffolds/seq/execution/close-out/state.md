# State: close-out (updated 2026-11-24, by the orchestrator)

Branch: `packages-seq` Head: `1e499c655` (CL7 merged; its records after)

## Where it stands

**Overnight run** (the author, 2026-11-24: "keep going with the loop, overnight run so we will reconcile things in
the morning"): the loop's "Running while the author is away". A decision a gate asks for is made as recommended and
listed in `gaps.md` under "Decided overnight"; the stop conditions still stop a lane; no public API change beyond a
design and no language change is decided; the use gate is deferred. A look the pane cannot show visibly is done in the
DOM and pixels and listed in `gaps.md` for the author's morning look.

Stage 0 and GATE 0 are closed; waves 1 and 2 are merged and GATE 1 passed. Merged since: CL13 (`27470de3e`), CL9
(`19733d6c3`, the default viewer's bound raised to 1,560 kB under the author's leave), CL7 (`1e499c655`). Sizes at
the head: spec 260.9/265, headless 978.5/980, default 1,554.6/1,560, minimal 1,094.7/1,095, ring 45.0/48.

## In progress

- 2026-11-24 night: CL11, strong, the drag rule (`briefs/cl11-drag-rule.md`), `seq/cl11`; started at load 6.3, its
  timings read only under load 3.
- FX25 stopped on a permission refusal (`gaps.md`), its worktree kept with items 1 and 2 uncommitted; waits for the
  author.

## Next

1. On CL11's report: the loop's step 5. Then CL12 (its brief to write), alone.
2. The legend's eye chunk flake (fixes.md) to the next fix lane.
3. CL5, CL6 and CL8 wait for the author's word on their drafts (`drafts/search.md`, `track-groups.md`,
   `values-marks.md`); FX25 for the refusal. GATE 2 after wave 4.

## Blockers and open questions

`gaps.md`: FX25's refusal and its allele-title question; the three drafts (two are language changes); `rowsValue` over
a lazy alignment; the overnight decisions to confirm (the zoom hatch; a 0.x minor read as a newer minor). The headless
and minimal bounds have 1.5 and 0.3 kB left: CL11 and CL12 keep new code lazy, or the orchestrator trims or raises
(the author's leave). At the close: remove the locked, clean worktrees (CL4's, CL9's, CL13's, the
reviewers') and their branches.

## Verify

The loop's Checks; the scaffold checks were last green at `1e499c655`.
