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

- **BZ's two core changes past the surface accepted** (the author, 2026-11-21: "we don't need to touch prism; we will be
  rewriting that. the other accept"): core's internal `_inflateMember` export; `ungzip` throws on a truncated gzip
  instead of returning partial bytes (Mol*'s plugin fetch, core's assets and data source see it; Prism is left as it is,
  to be rewritten), recorded in `.v6/plans/breaking-v6-changes.md`.

- **The size bounds raised for language 0.6.0** (the author, 2026-11-21: "Raise the bounds"): reading a spec 240 kB
  (was 230), the minimal viewer 960 kB (was 940); headless 830 and default 1,330 unchanged. GC1's declarations alone
  took the minimal viewer to 944.6 and reading to 228.7 kB. Decided as recommended with it: the include floor named
  `min_length` (snake_case like every spec param); a sequence under an assembly inherits the assembly's alphabet.

- **Genome resolve loads lazily** (the author, 2026-11-21: "Lazy genome resolve"): assembly resolve, the names rule and
  the locus and block forms are a chunk loaded only when a spec holds an assembly (as navigators in `renderSvg`), so a
  viewer that never opens a genome pays nothing; the bounds stay (reading 240, headless 830, default 1,330, minimal 960
  kB). GC1's resolve stopped at headless 831.1 and minimal 960.1 kB before it.

- **Sizes for GC2's genome axis** (the author, 2026-11-21: "raise the bound, even 10% is still fine"): headless
  870 kB (was 830), the default viewer 1,400 (was 1,330), the minimal viewer 1,000 (was 960), about 5 % each, room for
  the slice; reading stays at 240. GC2 measured 833.6, 1,332.3 and 964.4 kB.
- **One lane at a time** (superseded 2026-11-22 by three) (the author, 2026-11-21: "slow down to just 1 active lane after current batch finishes, we
  are nearing usage limit"): after GC3 and GG finish, one agent at a time (a lane, then its reviewer, then its fix
  round), never two; VC waits for both.
- **Stage 1's gate passed** (the author, 2026-11-21: every story walked, the findings in `fixes.md`; the novice
  questions "Yes, all five"): F3 "Within 5 Å" on a whole chain selects the other chains' residues within 5 Å, its card
  naming it; F9 the colouring track's icon pressed, one at a time, its tooltip naming it; F11 an empty section names
  its issue in its own row (the empty-well rule); F14 the "n filtered" card names the filter and says a click clears;
  F25 a navigator window at the fit drawn as one quiet frame (no grips, no target), the story opening zoomed.
- **Three lanes at once** (the author, 2026-11-22: "we can run up to 3 lanes in parallel now"), replacing "one lane at
  a time"; `loop.md`'s line changed by the author's word.
- **The builder stays in core; a Features story showcases it** (the author, 2026-11-22: "keep it in core for now, but
  add a storybook example to features section that showcases it"): lane BS; no `seqx-builder` package for now.
- **The GenBank download, and the two drifts as recommended** (the author, 2026-11-22: "Yes download and go as
  recommended with the drifts"): GB fetches L09137.2's flat file once; the SVG export draws what resolves on a spec
  with an error (viewer-composition.md); a mixed run takes its majority category's fill (genome-formats.md, Defaults
  2). Both fixes to FX6.
- **Sizes for RG's region reads** (the author, 2026-11-22: "Ok increase the limit", after a breakdown was offered
  first): reading 260 kB (was 240), headless 910 (870), the default viewer 1,450 (1,400), the minimal viewer 1,040
  (1,000); RG measured 254.9, 897.2, 1,428 and 1,026. RG still removes an import reading should not have, if its
  breakdown finds one.
- **Sizes raised again** (the author, 2026-11-22: "Raise again", offered a trim after the review first): headless
  920 kB (was 910), the default viewer 1,460 (1,450); RG after the merged lanes measured 912.1 and 1,450.1. The review's
  `tables-region.ts` cut stays a minimal-code item for RG, no longer a size gate.
- **A GenBank lane in stage 2** (the author, 2026-11-21: "Add a GenBank lane", after "Try pUC19.gb" in Services ›
  Your own file showed only a refusal): GB after VC, strong, a short design section first (plan, Stage 2); GenBank
  leaves Out of scope.
- **A Nightingale showcase as a storybook dev dependency** (the author, 2026-11-21: "yes, we can do a dev dependency
  version. maybe add a phase for this to slice 5?"): lane NS in stage 4; the npm install is approved with it; the
  fourth slice's "Nightingale dropped" stands for the packages (no dependency, no lane there).
- **Sonnet for well-scoped lanes** (the author, 2026-11-21: "maybe use this more in general for well-scoped
  subtasks?"): fix lanes, storybook-only lanes, pure splits and fast reviews run on a Sonnet subagent; lanes building
  a design's core and their strong reviews keep the strong model; the orchestrator's look and checks gate both.

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

Resume with `/loop run .seq/execution/fifth-slice/loop.md` (a new session; the 2026-11-21 one hit its usage limit). In
order:

1. **VC** (stage 2's last lane): write `briefs/vc-vcf.md` from genome-formats.md ("VCF", Bounds, Consequences) on BE's
   and GG's shape (`be-bed.md`, `gg-gff.md`); strong model (untrusted input), its strong review, one fix round.
   ClinVar's recording (#2: the `.tbi` whole and TP53's pages) is approved (the author: "Fetch all four").
2. **GB** (stage 2, after VC): write `briefs/gb-genbank.md` from the plan's GB entry; the genome-formats.md
   "GenBank" section agreed in the brief first; ask before downloading L09137.2's flat file. Strong model.
3. **Stage 2's GATE** (plan): `pnpm build`, `pnpm check:workspace`, each runtime's size row.
4. ~~Stage 1's GATE~~: passed 2026-11-21 (Decided).
5. **RG** (stage 3, strong, alone in the store), then FR, LA; GX and RX in stage 4.
6. **Fix lanes** on Sonnet (Decided; the two class items in `fixes.md`, the hover contract and the empty well, on the strong model) as `fixes.md` gathers items; open now: the hover card (placement, wrapping, the density row's highlight, bare row labels), the navigator's readout covering the control, track sizing with a cap and its own scroll (a design first), the navigator's window and zoom strip (a look design first), the square's sweep without feedback, the palette control (shown where it does nothing, twice, its arrow), the overview cut at the right on a wide page, the tree's leftmost line the 3,000-contig speed (M5).

Agents: one at a time while usage is tight (Decided); the orchestrator's look and checks gate every merge.
`foldback.md` collects the slice's lessons for the close.

## In progress

- 2026-11-21: the loop started (`/loop run .seq/execution/fifth-slice/loop.md`). Running (B3): FX1 (fast; started
  with this commit) and BZ (strong; same). Both merged (BZ `167d24491`, FX1 `2afb32e91`). GC1 merged (`e3b3a98d0`). SC merged (`2fdc44e76`), GC2 merged (`db2364d09`). BE merged (`98028d3f6`), GC3 merged (`61e71891b`). GG merged (`9b8cf675b`). FX2 merged (`e678d2dbc`), NS merged (`190dfc135`), FX3 merged (`ded361691`), FX4 merged (`5bfebf366`). No lane running (2026-11-21, the session ended at its usage limit); next VC, then stage 1's gate (the author's), RG; open fixes in fixes.md.
- 2026-11-22: the loop resumed in a new session. BS merged (`5bf984d5b`), VC merged (`364acfb5e`), FX5 merged (2026-11-22), FX6 merged (2026-11-22), GB merged (2026-11-22). Running (B3, three at once from 2026-11-22): **RG** (strong; same commit, `briefs/rg-region-reads.md`; alone in the store) (2026-11-22 ~09:40: reported partway; ~11:30 step 4 done: grid on store tiles, the matrix layer, ring images as runs, hover texts and the card; no accidental import in reading (its growth is the region forms' check through `tables.ts`); headless 909.9 of 910 kB; two allocation regressions (the 40-state ring 16 B, the p53 pan 58 kB vs 54); ~12:00: the base merged in (no conflicts), both allocation regressions fixed (the pan 53.5/55.0 kB, the ring 0 B); headless 912.1 of 910 kB and the default viewer 1,450.1 of 1,450 after the merge: asked of the author; RG stopped after step 2. Not merged). Running from ~12:05: **RS** (Sonnet, `briefs/rs-region-stories.md`, on `seq/rg`: RG's step 5); RG's strong review done (~12:15: 3 bugs, 6 missing design, weak properties, 5 EC titles, the `tables-region.ts` duplication as the size cut); **RG's fix round** from ~12:20 (the same agent, core only); RS reported ~12:35 (`seq/rs`: generators, nine stories, the Design page, the host-adapter test; six core gaps as `it.fails` tests, given to RG's fix round: sources over an assembly, BEDPE through the store, a last partial tile, the matrix's zoom-in, `checkAdapter`'s `y`, held pairs' bins); RS adapting its stories to RG's fixes (~13:05); **FX7** (Sonnet, `briefs/fx7-small-fixes.md`, six items clear of RG's files) and **DL** (strong, `briefs/dl-look-designs.md`, two design drafts, no code) started ~13:10 (the author: "we have 1-2 free lanes");  Downloads stay asked before they are made (plan, How this slice runs).

## Findings by lane

One row a merged lane (loop, step 5, record): review and look findings by kind, and who looked (B2).

| Lane | Merged | Bug | Missing design | Look | Who looked | Notes |
| --- | --- | --- | --- | --- | --- | --- |
| BZ | 2026-11-21 | 6 (review: plain members allocating the rest of the file each; a quadratic header scan; chunks before `min_off`; a one-block chunk costing 0; cached blocks aliased; internal folders published) | 2 (`_inflateMember` past the surface, made private; the generic index's end rule credited to `tbx_parse1`) | 0 | nothing drawn | the author accepted the truncated-gzip throw; Prism untouched; no htslib fixture (not installed) |
| FX1 | 2026-11-21 | 7 (my look: letters under the labels' cap height (a flat 0.9); review: B7 unable to fire; the readout's ellipsis cutting the number; author numbering across insertion codes; the bigWig twin reading whole; the author: the legend hover sticking) | 2 (B6 as a list to check, not to exempt; unobserved ends as gaps) | 0 | the orchestrator: coding DNA (light; SVG letters 11.75 px against 12 px labels; legend fast passes), Globins rows (dark; readout) | the two checks in the loop's Checks; a 20 KB bigWig too small to measure pans (M5 needs a larger) |
| GC1 | 2026-11-21 | 15 (review: a failed chunk leaving a spec pending; unread params silent; ReDoS in include patterns; an excluded name unreported; an override's letters without an alphabet; duplicate names overwriting; an alias equal to another name dropped; the hint; HLA names misread; `(block :name)` outside assemblies; a vacuous property; three EC ids; blocks built twice; my checks: its own perf regressions (48 and 16 B a frame), the 100,000-sequence resolve at 53 ms against 50) | 1 (GC1's claim that its perf failures were the base's) | 1 (a bare number as a block name refused with "quote it": fixes.md) | nothing drawn | bounds raised (reading 240, minimal 960) and genome resolve lazy, by the author; reading at 239.2 kB, 0.8 kB left |
| SC | 2026-11-21 | 10 (my look: residues never drawn (a held unit's letters never reached the store); the legend differing before and after a drop; the hover card's raw address; review: the describer and copy over the static model; the perf case never hovering; a 100k pre-warm hiding `ViewerRails`' deopt loop; the overlay fix untested; a mistyped override dropped silently; the frame's and view's overlays keyed apart) | 0 | 0 | the orchestrator: the host's complex (light: hover, zoom to letters; dark: drop and restore) | a heavy theme ratio failed under load (GC2 building), passed alone on the branch and the base |
| GC2 | 2026-11-21 | 15 (my look: the hover card's concatenated offset; the whole ruler lit at the fit; labels cut to "c…" and the direction mark lost; no ruler ticks in the 25-block row; group labels among block labels; review: `data-assembly` off by an empty group; uncapped row group gaps; runs depending on which spec opened first; a lost companion load; `placesOf` linear; O(n²) block themes; label strings per layout; thin tests; surface past the design) + 1 found by its own property (a snapped block in a run read past its end) | 0 | 4 (my second look: fixes.md) | the orchestrator: row, B in A's order, two on a ring, 3,000 contigs (light and dark) | the author raised headless 870, default 1,400, minimal 1,000 kB; the lazy split undone with that room |
| BE | 2026-11-21 | 16 (my look: no strand on BED12 genes; a gene hovered and clicked as a position; no band labels at the fit; review: chrom.sizes quadratic (80k lines 16.9 s); BED12 blocks looser than the design and unbounded pieces; a truncated recorded page read as zeros and errors cached as no recording; the property's gaps (zero-length items, `head` past `track` lines); the BEDPE bound's name; every example's status line and a pages fetch per file; racy story waits; titles claiming more; BEDPE eager; the script's `../`) | 0 | 1 (a double-click on a long sequence jumps to bases: fixes.md) | the orchestrator: cytobands (light, dark: labels at the fit), BED12 genes (strands, hover, click), the misread file | the ranges layer changed (items at the fit), outside the brief's runtimes, after my look asked for labels; BED row 23.2 kB first load, 20.5 kB lazy |
| GC3 | 2026-11-21 | 24 (my look: the navigator's click doing nothing; no bands or zoom strip, misaligned; the square unlabelled and deaf, its readout ": 1–0"; the zoom strip's label clipped; review: history mirrored after a reverse; the square's wash past `maxCells` per frame; a stale stop block; a vacuous selection property; an O(blocks) history signature; duplicate overrides; the cursor jumping copies; the locus fallback's abort and hidden hits; the zoom strip's 16 B; leftovers; the author: header icons blinking along the ruler, the zoom strip cut, switching keeping a width, the rail dead at the fit, small chromosomes unlabelled) | 0 | 2 (bands wait for RG; menu entries static: TODO) | the orchestrator: a chromosome at a time (clicks, labels, hover), the square, a ring that switches, the row at the fit (dark) | the author's rule for switching; a heavy pan bound failed under load (GG running), passed alone on the branch and the base |
| GG | 2026-11-21 | 13 (my look: labels drawn as boxes past a model and over exons; a stray strip in Genes (`features` dropped by resolve); find-by-name not showing a find; review: a deep `Parent` chain and a 200,000-exon transcript overflowing the stack; `features: rows` quadratic; an encoded comma splitting `Parent`; a gene losing its mRNA; a GTF gene hit's label; long items making tiles O(n)) | 0 | 4 (the author and my second look: the empty actions button, the hover card's anchor, the Transcripts track's spacing, the SVG not repainted: fixes.md) | the orchestrator: gff3 tp-53 and find-by-name (light, dark) | a ranges-layer fix (lanes past the band) outside the brief, small, with its test; GFF row 25.8 kB |
| FX2 | 2026-11-21 | 6 fixed of 9 (the author's: a hidden track unnamed, clearing a sort, the SVG on a cleared selection, the hover card's anchor, lanes from a tall band's top, a double-click's step) | 0 | 3 open (the empty actions button not reproduced and closed by the author; block rulers at the fit and a ring's own hover mark open) | the orchestrator: gff3 tp-53 (card by its lane, lanes from the top, the Genes menu); no reviewer (the last agent task) | my check found the heavy hover walk drawing nothing after the lanes moved up: the walk now probes a track's top first (`59a20e919`) |
| NS | 2026-11-21 | 0 (a Sonnet subagent, the author's ask to see a weaker model on the library: three stories working first time, light and dark; the link one way as asked) | 0 | 2 (Nightingale's panel not column-aligned with Seq*'s; the NG consumer's SVG not shown: it is test-only, so Seq*'s own export stands as the reference) | the orchestrator: sequence-and-features (light), the-link (dark: DNA-binding selected, highlighted in Nightingale) | Nightingale 5.11.0 (three packages, MIT) as storybook dev dependencies; no biome rule added: the package boundary check already keeps it out of the other packages |
| FX3 | 2026-11-21 | 3 fixed (the rail beside a declared navigator and its place, designed first; a ring's own hover mark; an override naming nothing) + 1 not reproduced (block rulers at the fit) | 0 | 0 | the orchestrator: Globins rows family (dark: the rail beside the overview), three-thousand-contigs (the ring's hover box) | a Sonnet subagent; 5 tests changed for the new default; `railPlace: navigator` rows at the canvas's first or last row (a navigator is a host part); its storybooks kept waking it: stopped after its report |
| FX4 | 2026-11-21 | 2 built (the rail's place as a list; `pinned`) | 0 | 0 | the orchestrator: FamilyRailPinned (two rails; a drag on the pinned grip zooms and the top rail follows) | a Sonnet subagent; stopped after its report |
| BS | 2026-11-22 | 2 (review: the numbering's handle dropped in Handles; a tautological source test) | 1 (the SVG export blank on a spec with an error while the canvas draws: drift/svg-on-error.md) | 3 (my look: a stray "Builder Story" export; a chromosome called a locus; review: a hex in a comment) | the orchestrator: First view (light), Handles (light: hover card "a, hits, hello 3–8", click selects 3–8 on canvas and SVG, pointer off), Layers and themes, Check (dark), A genome view | Sonnet lane and Sonnet review; its guard found three element classes listed as stories (fixes.md); `vite/client` types added to the storybook tsconfig (hot, accepted) |
| VC | 2026-11-22 | 10 (review: three quadratic scans on untrusted input (fixedEnd 77 s on 8 MB, INFO `=`, Number=A/R per allele); merged runs losing their fill; a fix commit without its test; `#2` ids unstable under END=; QUAL `''` and `0x10` read as numbers; my look: marks grey not by category; review nits: IndexHint duplicated, views.ts at 978) | 2 (the zoom-in text not built: given to RG; mixed-category runs neutral: drift/vcf-mixed-runs.md) | 2 (the stories' text claiming a zoom-in text and an allele card: corrected; the allele's card: fixes.md) | the orchestrator: TP53 ClinVar (light; hover a 7 kb deletion: card and legend; double-click focus), long deletion (dark), chr17 at the fit; canvas pixels counted (the pane reported hidden) | the ClinVar licence row says the terms page was not fetched: to confirm before stage 2's gate; VCF row 25.2 kB (estimate 15, bound 30); heavy lane not rerun after the fix round (load 11), at the gate |
| FX5 | 2026-11-22 | 5 fixed of 6 (the readout over the control; the square's sweep, snapped to cells after review; the palette control's four bullets; the tree's root line; story exports) + 1 review (a key and a submenu losing the arrow in `packages/ui`'s menu) | 0 | 1 open (the overview cut at the right: not reproduced, hardened, open) | the orchestrator: BED cytobands (mid-drag readout by DOM boxes clear of window and pointer; the pane reported hidden), Themes (the Variants menu, Color submenu, a palette round trip) | Sonnet lane and review; `packages/ui`'s context menu changed (Prism has its own copy, unaffected); new surface `SquareNavigator.snapSweep`, `dragTo(at, from?)` |
| FX6 | 2026-11-22 | 2 fixed (the drifts as decided) + 1 review (categories in `Hit.state`, which rows use: `Hit.categories`) | 0 | 0 | the orchestrator: ClinVar TP53 (dark; canvas pixels by category), Builder › Check (the SVG drawn, the error in `<desc>`, by DOM) | Sonnet lane and review; export tests added (clean spec, escaping, the bound blank, ten errors); the bounds exception written into viewer-composition.md |
| GB | 2026-11-22 | 6 (review: an origin join past 2^31 overflowing Int32; the nesting bound off by one; an unclosed quote swallowing a feature silently; a repeated name indexed joined; qualifiers wrapped mid-word; my look: "1 sequences") | 2 (`order` read as a wrap: decided no; the open quote's rule: decided) | 1 (the twin's GenBank file without the sites: written as `protein_bind`) | the orchestrator: GenBank plasmid (light: hover bla across the origin, one item; click), Your own file (dark: Try pUC19.gb loads) | the lane refused the relayed download approval (right): the orchestrator fetched L09137.2 itself on the author's word; a click across a cyclic origin selects one piece (fixes.md, shared with GFF3); heavy lane not run in full (load), at the gate |
