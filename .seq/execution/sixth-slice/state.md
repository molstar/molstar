# State: the sixth slice (scaffolds)

## Decided

- Stage 0 (the author, 2026-11-23): `.scaffolds/`; the execution folders kept in the tree; meso\*'s modules kept in
  Seq\*'s instance for now; `.v6/` moves later; a closed scaffold kept, editable, never cited by code (plan, "Decided
  at stage 0").
- meso\* removed from `packages/core` (the author, 2026-11-23: "this is something that exists at a different "level"
  and should not leak into the common libraries"): `351d9e9fc`.
- Started (the author, 2026-11-23: "Yes, start it").

## In progress

- 2026-11-23: stage 1 started. `.seq/_checks/designs-in-specs.ts` written (46 designs pending). Running: **S1A**
  (strong, `briefs/s1a-specs.md`) and **S1B** (strong, `briefs/s1b-specs.md`, merged `81a82d453`); S1A merged (`936242251`); **S1D** (Sonnet, `briefs/s1d-specs.md`) started beside it; **S1C** (strong, `briefs/s1c-specs.md`) started after S1A's merge.

## Findings by lane

| Lane | Merged | Spec rows | Design-only | Drift | Model | Items ticked | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| S1B | 2026-11-23 (`81a82d453`) | about 245 rows across 13 designs | 41 | 2 (viewer-models-behaviors, header-color-slot) | strong | 13 of 13 designs out of PENDING | review: the clipboard rule wrong (a behavior's param, fixed), the clock's port, a drift list, the tooltips row split, 3 rows missing (added); the examples package has no spec (S1D writes it); a HitTarget rule owed to core (S1A) |
| S1A | 2026-11-23 (`936242251`) | 15 designs, about 310 rows (core, repr, viewer) | about 50 | 5 (msgpack-codec, unknown-kinds, language-description, locations-and-queries, molecules) | strong | 15 of 15 out of PENDING | review: drift notes against existing TODO entries, security's bytes in flight (TODO line added), ambiguous pointers, an HLA example in core's repr spec, a wrong path; fixed; Targets are closed added to core for S1B's row |
