# State: close-out (updated 2026-11-24, by the orchestrator)

Branch: `packages-seq` Head: `d9b53e580` (FX31 merged; GATE 2 passed) and its records

## Where it stands

The overnight run is over; the author reconciled `gaps.md` on 2026-11-24 ("Agree all": every recommendation as
written). Search, track groups and values marks are agreed (the drafts in `drafts/`, each moved to `design/` by its
lane); the language bump to `0.11.0` is CL6's, with CL8's `reference` after it in the same unreleased version.

Stage 0 and GATE 0 are closed; waves 1 and 2 are merged and GATE 1 passed. Merged since: CL13 (`27470de3e`), CL9
(`19733d6c3`), CL7 (`1e499c655`), CL11 (`cc4fab363`), CL12 (`245c82baf`), FX26 (`08a6ffe44`), FX27 (`938396a33`), CL5 (`a4668800a`), CL6 (`a1097142a`), FX25 (`9040adc8d`), CL8 (`a95dbf373`), FX29 (`271ab74ed`), FX30 (`60170c5b5`), FX28 (`d7857a8f0`), FX31 (`d9b53e580`). Sizes at the head: spec 265.0/268, headless 995.0/1,000, default 1,587.5/1,592, minimal 1,119.6/1,124, ring 45.0/48,
bridge chunk 26.0/28.

## In progress

GATE 2 passed (`d9b53e580`); the use gate ticked (the author, after FX32 `608d87ac2`); the novice pass ran (three
groups, `novice-findings.md`: 43 findings, 7 reproduced or read against the designs before placing). The old Seq\*
viewer's label churn over its baseline is in `fixes.md`, outside the close-out's packages.

## Next

1. FX34 (`293155bb1`) and FX35 (`90cdfe934`) merged: the novice pass's bugs and its four gaps.
2. The measurements (the bundles against their budgets, the representatives at full size), the last gate, FINALIZE.

## Blockers and open questions

The author: a name for a folded group's summary layer (its raw id
"ranges#summary" shows in its legend row and eye: TODO). Sizes raised at FX30's merge (headless 1,000). At the close: remove the locked, clean
worktrees (CL4's, CL7's, CL9's, CL11's, CL12's, CL13's, FX26's, the reviewers') and their branches.

## Verify

The loop's Checks; the scaffold checks green at the records commit after `72012efcf`.
