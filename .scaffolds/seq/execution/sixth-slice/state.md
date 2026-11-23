# State: the sixth slice (scaffolds)

## Decided

- Stage 0 (the author, 2026-11-23): `.scaffolds/`; the execution folders kept in the tree; meso\*'s modules kept in
  Seq\*'s instance for now; `.v6/` moves later; a closed scaffold kept, editable, never cited by code (plan, "Decided
  at stage 0").
- meso\* removed from `packages/core` (the author, 2026-11-23: "this is something that exists at a different "level"
  and should not leak into the common libraries"): `351d9e9fc`.
- Started (the author, 2026-11-23: "Yes, start it").

- Stage 1's fifteen drift notes (the author, 2026-11-23: "All as recommended"): designs amended to the code's way by
  S1F; the code items to the close-out plan.

## In progress

- 2026-11-23: stage 1 started. `.scaffolds/seq/_checks/designs-in-specs.ts` written (46 designs pending). Running: **S1A**
  (strong, `briefs/s1a-specs.md`) and **S1B** (strong, `briefs/s1b-specs.md`, merged `81a82d453`); S1A merged (`936242251`); **S1D** (Sonnet, `briefs/s1d-specs.md`) started beside it; S1C merged; stage 1's gate passed (15 drift notes to the author). S1E and S1F merged; stage 1b done. Stages 2 and 3 done by the orchestrator (`0be663880`, `88f55a78c`, `e0f554e49`), reviewed (no blockers; six stale lines fixed); stage 3's gate passed (quick 4,949, lint, format, check:workspace, every scaffold check from its new place, 612 links in the changed files resolve). Running: **S4** (strong, `briefs/s4-method.md`).

## Findings by lane

| Lane | Merged | Spec rows | Design-only | Drift | Model | Items ticked | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| S1B | 2026-11-23 (`81a82d453`) | about 245 rows across 13 designs | 41 | 2 (viewer-models-behaviors, header-color-slot) | strong | 13 of 13 designs out of PENDING | review: the clipboard rule wrong (a behavior's param, fixed), the clock's port, a drift list, the tooltips row split, 3 rows missing (added); the examples package has no spec (S1D writes it); a HitTarget rule owed to core (S1A) |
| S1A | 2026-11-23 (`936242251`) | 15 designs, about 310 rows (core, repr, viewer) | about 50 | 5 (msgpack-codec, unknown-kinds, language-description, locations-and-queries, molecules) | strong | 15 of 15 out of PENDING | review: drift notes against existing TODO entries, security's bytes in flight (TODO line added), ambiguous pointers, an HLA example in core's repr spec, a wrong path; fixed; Targets are closed added to core for S1B's row |
| S1D | 2026-11-23 | 7 designs, about 95 rows (ext-rondo, examples, core) | about 30 | 2 (ring-extension, packages-and-examples) | Sonnet | 7 of 7 out of PENDING; the examples spec written | review: the stubs row wrong (hit keys), draw order missing (added), lazy chunks and requires, unmarked drifts, five rows missing, the examples spec's io and entries; fixed; the orchestrator added ext-rondo and examples to packages/seq/AGENTS.md |
| S1C | 2026-11-23 | 12 designs, about 350 rows (core, viewer) | about 50 | 5 (store-over-budget, region-reads, genome-formats, grids-and-rows, layouts) + blosum62-provenance | strong | 12 of 12; PENDING empty | review: the hover outline's 3:1 a missed rule (added), five pointers on the wrong section, int32 positions' rule, reasons missing, a provenance string contradicting the code (drift); fixed; the Tables' Bounds lead renamed (a test comment cites the old name: stage 1b) |
| S1F | 2026-11-23 | 15 drift notes resolved, about 25 rows repointed | 7 code items to the close-out | 0 | Sonnet | 15 of 15 Resolved | review: stale text in ring-extension, examples, version-1, base-library; a garbled bullet; fixed |
| S1E | 2026-11-23 | about 240 citations in seq, core repr, io: about 115 repointed, 125 dropped | — | 0 | Sonnet | the check passes; PENDING emptied by the orchestrator (TODO.md, sizes.mjs fixed) | review: EC-169 dropped (restored), the check's exemption too wide (narrowed, self-check), a missing heading cited, two io comments; fixed |
