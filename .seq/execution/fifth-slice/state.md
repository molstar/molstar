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

- **The whole-designs review's twelve recommendations taken** (the author, 2026-11-21: "yes with the
  recommendations"; `reviews/whole-designs.md`): GF's dropped runtimes to a "designed only" appendix and their
  additions dropped; no `gzi` in io; the `seq` runtime serves sided systems and pairs from either end (LA after RG and
  FR); BEDPE indexed by both ends after a whole read, held rows in a byte budget; no doubled BEDPE (pairix noted for
  Hi-C); a per-view byte bound for tabix levels; a host's bigWig as the whole-genome remedy for 1D layers; the Manhattan
  family, `split` and EC-224's join to the close-out slice; MAF cut, EC-221 to the matrix, EC-226 to themes; fetch 3
  cut; RN's gene detail on GFF and `coarsenBins`; GD reviewed before promotion. Folded into the drafts with their
  reviews.

- **The fold-in's three** (the author, 2026-11-21: "yes"): each draft's "Gaps for the author" becomes "Decided" at
  promotion; the fifth slice's plan reconciled with the drafts when written (no BD, SQ or PF lanes; axis.md's tension
  "narrowed"); GTF names indexed beside GFF3's (`runtimes/gff`).

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

## Ready to start (2026-11-21, the fourth slice's finalize, step 5)

- **Promoted** (uncommitted at writing): `design/genome-coordinates.md` (GA), `region-reads.md` (RR),
  `genome-formats.md` (GF), `across-sections.md` (LX), `ring-extension.md` (RN), each with "Decided by the author"
  in place of its gaps and its references moved; the Genomes dataset (GD) as `examples.md`'s Genomes section, as the
  other datasets. Edited from their Consequences: axis (the tension narrowed), extension-openings, files, sources,
  examples, layouts, navigators, grids-and-rows, protocol, base-library; the README's table. The other designs the
  Consequences name are amended by the lanes that build them (plan, The designs). `drafts/README.md` marks the drafts
  promoted; the folder is kept as the record.
- **The plan written** ([plan.md](plan.md)): FX1 first (B6 and B7's checks, the services' structure section, the flake
  watch list, a bigWig story, `fixes.md`); stage 1 GC1 → SC, GC2 → GC3 (gate with F3, F9, F11, F14, F25); stage 2 BZ
  → BE, GG, VC; stage 3 RG → FR → LA; stage 4 GX, RX; the novice pass, M5 (with M4's carried measurements),
  finalize; FX2 to FX6 at the fix points.
- **`loop.md` in place** (the fourth slice's with its fold-back applied, two lanes at once); `fixes.md` the fix phase's
  queue.
- **RR's language part declared by GC1** (the author, 2026-11-21: "Yes, in GC1"): `pairs.source`, `grid.source`,
  `pairs.field`, `grid.field`, the `matrix` kind and layer type join GC1's one minor bump, as GF's and LX's params do;
  RG builds what reads them.

## Next

Start the loop: `/loop run .seq/execution/fifth-slice/loop.md`. Its first step writes `briefs/fx-fix-lane.md` (from the
fourth slice's) and `briefs/gc1-genome-language.md`, and starts FX1 beside BZ or GC1 (two lanes at once).
