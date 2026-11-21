# State: the fifth slice

## Decided

- Gate 0 (the author, 2026-11-20: "the recommendation make sense"): every question as recommended (plan, Gate 0).
- Stage 0's drafts and their reviews run overnight, before the fourth slice's close; the author agrees them in the
  morning (the author, 2026-11-20: "Let's schedule the designs steps for an overnight run (we can do the design drafts
  and their reviews before we close slice 4 in the morning after i wake it)").

- The drafts (the author, 2026-11-21): "generally reasonable and we can fold in the reviews": each draft's "Review,
  unresolved" and the text it asks of another draft folded in before promotion.
- **Three genome formats built, the rest designed only** (the author, 2026-11-21: "Yes, those three"): `runtimes/bed`
  (BED, bedGraph, BEDPE), `runtimes/gff` (GFF3, GTF as transcripts) and `runtimes/vcf` with BGZF and tabix in
  `@molstar/io`; bigBed, 2bit, bgzipped FASTA (`.gzi`), PAF, PIF and chain stay designed, no runtime. Why: a runtime is
  maintenance and fixtures, a host brings its own through `withAdapter`. The dataset's synteny, dot plot and liftover
  views read the alignment from the `seq` format (written from the alignment truth in code); "open by name" reads GFF3
  names only; fetches 4, 5, 6, 9 and 10 dropped.

## Overnight: stage 0 (run by the fourth slice's loop)

The fourth slice's loop runs this while the author sleeps (its `state.md`, Next, points here). The lane limit is the
fourth slice's `loop.md` line (two), shared: one fifth-slice draft (or its review) beside one fourth-slice lane, while
both have work (the author, 2026-11-20: "you can continue with slice 4 too, just don't close it until i review it in
the morning").

1. **Order**, one at a time beside the fourth slice (two when it has nothing ready): GA, RR (the foundations), GF,
   LX (it reads GA's and RR's), RN, GD (they read the others). Briefs in `briefs/`; each writes one file in `drafts/` in the main checkout (no worktree,
   no code, no commits).
2. **Each draft**: on its report, an independent review by a strong reviewer (read only; it gets the draft, the
   designs it amends and the brief, never the lane's report): the design's rules followed, the self-check, Gate 0's
   answers built on, precedents, the interaction matrix, Defaults, the Configuration table, edge cases by id, conflicts
   with the other drafts; findings ranked. Then one revision round by the drafting lane (SendMessage). The orchestrator
   commits the revised draft with the review's unresolved findings listed at its end, under "Review, unresolved".
3. **Stop**: a draft that needs a decision the plan's Gate 0 does not give lists it as a gap with a recommendation; it
   never stops the night. A refused permission or a fetch (none is needed for a draft) stops that lane and waits.
4. **Morning**: a summary per draft (the decision, the gaps, the review's unresolved findings) for the author; the
   drafts are agreed and promoted to `.seq/design/` at stage 0's gate.

## Progress

| Draft | Lane | Drafted | Reviewed | Revised | Committed |
| --- | --- | --- | --- | --- | --- |
| genome-coordinates | GA | 2026-11-20 | 2026-11-20 (3 blocking, 8 to fix) | 2026-11-20 | yes |
| region-reads | RR | 2026-11-21 | 2026-11-21 (2 blocking, 8 to fix) | 2026-11-21 | yes |
| genome-formats | GF | 2026-11-21 | 2026-11-21 (3 blocking, 8 to fix) | 2026-11-21 | yes |
| across-sections | LX | 2026-11-21 | 2026-11-21 (2 blocking, 8 to fix) | 2026-11-21 | yes |
| ring-extension | RN | 2026-11-21 | 2026-11-21 (2 blocking, 7 to fix) | 2026-11-21 | yes |
| genomes-dataset | GD | 2026-11-21 | | n/a | yes |
