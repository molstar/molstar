# State: the fifth slice

## Decided

- The values layer's look: variants first (the author, 2026-11-22: "is there a way to make the barcharts look nicer
  too?", then "Yes, variants first"); VL (`briefs/vl-values-look.md`) after RL merges.
- The rail look picked (the author, 2026-11-22): "let's go with D1 by default, but customizable to D2 as well"; "the
  hover card should not show on the rail". A, B, C and today's deleted. To RL step 2.
- Header actions show on header hover only (the author, 2026-11-22: "Yes, header only"): not on data or rail hover;
  also with the keyboard cursor on the header or its menu open. To RL.
- The default viewer raised to 1,480 kB (the author, 2026-11-22: "raise default to 1,480"), LA at 1,470.4.
- The rail's look: today's, polished (the author, 2026-11-22, after A–C: "can we try option D as 'todays' but with
  better color and knobs? in general, i like the simple look of today's one .. it just needs some polish"): RL adds D
  (the line and a heavier window in a mid ink, round knobs where the box zooms, quiet at the fit, inset ends).
- The rail's look: variants first (the author, 2026-11-22: "Yes, variants first"); RL (`briefs/rl-rail-look.md`) after
  FX10 merges: capsule, capsule in a groove, quiet until hovered, beside today's.
- Minimal raised to 1,060 kB (the author, 2026-11-22: "raise minimal to 1,060 too"), LA at 1,039.95 of 1,040.
- Sizes raised a third time (the author, 2026-11-22: "You can raise the limits"): headless 930, default 1,470 kB; the
  chromosome pan's fit bound 55 KB (FX9 measured 54.1).
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
- **The navigator's readout in one fixed place outside** (the author, 2026-11-22: "it's better if the label stays
  outside in top right instead of jumping around all the time … the same thing for the 'vertical maps' … 'above' the
  map"): top-right above a horizontal navigator, above a vertical one; into DL's navigator draft, then a fix lane.
- **DL's two drafts promoted** (the author, 2026-11-22: "Yes, as recommended", all twelve gaps): the window's look and
  readout into navigators.md, track sizing into layouts.md, each with its Design page line; the drafts kept as the
  record. Fix lanes: the navigator's (now, clear of RG), track sizing after RG (tracks*.ts are RG's).
- **Inline data over an assembly names its sequences** (the author, 2026-11-22: "Ok we can use this naming"):
  `a_name`/`b_name` for pairs, `x_name`/`y_name` for grids (`name` for a matrix's items), parallel to maps'
  `from_name`/`to_name`; a language param (a minor bump), built with the matrix's `source` typing in one language
  lane (GA3) after GA2 and FR merge; until then the info `seq.resolve.assembly.inline`.
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
- 2026-11-22: the loop resumed in a new session. BS merged (`5bf984d5b`), VC merged (`364acfb5e`), FX5 merged (2026-11-22), FX6 merged (2026-11-22), GB merged (2026-11-22), FX7 merged (`29546ad25`), RG with RS merged (`590151907`), FX8 merged (2026-11-22), GA2 merged (2026-11-22), FR merged (`0af3f0809`), GA3 merged (`64e2f87cc`). Running: **FX12** (Sonnet, `briefs/fx12-empty-wells.md`, the free slot, the author: "use the free slot for a fix lane"); **LA** (strong; steps 1, 2 and 2's rest done on `seq/la` (resolve, layout and render, per-side reads drawn once, blocks merged under a pixel, the hatch; the bridge chunk 20.8 kB, its own row raised by LA to 24 kB: the author's to confirm); its agent out of context ~16:35, a fresh strong agent on `seq/la` for steps 3 and 4), **RL** (Sonnet, `briefs/rl-rail-look.md`, three rail looks to pick from), FX10 merged (`38f147b49`), FX11 merged.

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
| RG | 2026-11-22 | 13 (review: items-mode pairs not drawn after a pan inside tiles; per-row strings and Sets each layout; a refused pyramid reading before refusing; margin wants counted toward the layer bound; my look and the author's: the zoom-in state never drawn for one want past its bound; RS's review: a grid's colours cached under a stale domain (the canvas saturated, the SVG not); the store settling while summaries build; RS's gaps: a last partial tile drawing nothing, the matrix's zoom-in, `checkAdapter` without `y`; two allocation regressions (the p53 pan 58 kB, the ring 16 B)) | 8 (coarser before zoom in unbuilt; no field sums; matrix shares; mergePairs in one task; a dropped level answered as items; the bin element unspecified; data over an assembly unread (GC2's, to GA2); BEDPE not through the store (to FR)) | 0 | the orchestrator: Pairs zoomed out (both rows), Hi-C by region (dark), Matrix samples (about 7 s at load 19), the zoom-in state (with RS); an earlier look of mine saw empty stories from a background tab (corrected) | sizes raised twice by the author (920, 1,460); `tables-region.ts` on one shared schema (−5 kB); held pairs' bins to GA2; `store.ts` 854 and `viewer.ts` 818 lines (a question); weak properties strengthened, 5 EC titles added |
| RS | 2026-11-22 (with RG) | 2 (the twins not composing as the stories; a closed colour domain hiding a core bug: removed) | 0 | 1 (synteny's 2,000 arcs a blob: about 80 blocks) | as RG's | Sonnet lane and review; generators deterministic on PCG; the genome ring and Hi-C over the genome on one system until GA2 |
| FX7 | 2026-11-22 | 5 fixed + 5 review (a per-hit row scan matching keys; zoom-to-hit across the origin; `leaderRange`'s complement; a wall-clock flake fix; float ties) | 0 | 0 | the orchestrator: ClinVar's 7 kb deletion pink on canvas and SVG (the tab fronted); the plasmid click not landed at this page size (tests stand) | Sonnet lane and review; the tie test passes on its parent (unproven); the tree's way back stopped on RG's files (open); navigator lines to FX8 |
| FX8 | 2026-11-22 | 2 fixed (the window's look; the readout's place and its lagging text) + 5 review (config as constants; short forms untested; the fit sweep untested with real pointers; closures per frame; `navigator.ts` at 1,500 lines, split to 679 and 853) | 0 | 0 | the orchestrator: BED cytobands (the readout above at rest; a stepwise sweep at the fit zoomed to 12.5–28.3 Mb, the readout equal; a one-shot drag too fast to sweep); my tab twice behind a lane's tab (checked `tabs_context`) | Sonnet lane and review; the lane's merge of the base refused (right: lanes never merge), merged by the orchestrator; a drag across 3,000 contigs allocates about 110 kB a frame (the drag rule's item); its 2 ms bound to the heavy lane |
| GA2 | 2026-11-22 | 6 (review: a symmetric grid's mirror asked of the wrong pair; per-frame gathering in pairs over an assembly (a summarized coverage getter allocating 40 B a call); the fit cache missing per block and freeing its signature; unknown names silent over two assemblies; my look: the genome ring a solid disc; its text stale) | 3 (inline data over an assembly dropped silently: said, columns decided for GA3; held pairs over blocks by device column; the whole-assembly level not cut into block pairs) | 1 (the ring still dense with uniform random links: fixes.md) | the orchestrator: cytobands over the genome, Hi-C over the genome, the genome ring (light and dark), the tab fronted | heavy lane unconfirmed (three wall-clock bounds under load 4–9), at the gate; `view.ts` 1,494 → 1,287 lines; headless 917.2 of 920 kB |
| FR | 2026-11-22 | 9 (review: GFF reading past `readCost` (a gene's hull); a whole VCF held in memory; the genotype line buffers uncounted; tiles not `level × 1024` at 4096; levels asked past the first that does not fit; io's BGZF in every runtime's first load; `numberAt` through strings; my look: the Genotypes legend a 0-to-1 ramp, not classes; an unneeded `as never`) | 4 (ClinVar's coarse levels need chr17 recorded past TP53 (a fetch: the author's); BED/GFF coarse stories need a BGZF writer in examples (hot `package.json`); GTF without transcripts not saying `extent` by region; a GFF tile's first `readCost` under-counts its flanks, unsaid) | 0 | the orchestrator: VCF Genotypes (five class swatches, canvas = SVG), Pairs zoomed out (held and BEDPE one picture), the tab fronted | heavy: the 600k pyramid bound failed at load 13 (FX9 running), at the gate; runtimes BED 20.0, GFF 23.8, VCF 22.3, GenBank 22.4 kB; headless 917.6 of 920 |
| GA3 | 2026-11-22 | 8 (review: inline pairs over an assembly dropping `id`, `label` and ride-along columns; a grid serving every sequence pair unbounded (3,400 scaffolds out of memory); ends past a length unchecked; unknown names said while the assembly is pending; rows with no name and matrix items past a length dropped silently; two assemblies' ranks interleaved; bins prebuilt for every level at resolve (now lazy by the halving rule); `binPairs` eager through the data barrel) | 3 (a level-1 source's sparse links over a genome "zoom in" at the fit (`it.fails`, the author's); a grid's inline values per sequence pair one dense x-major grid (a binned Hi-C cannot be written inline); region pairs coloured by kind unlike held pairs) | 0 | the orchestrator: links written in the spec on a chord ring, a matrix written in the spec over two chromosomes (canvas = SVG); my "22"/"M" label finding was the narrow blocks' shared-`chr` drop, by design | language 0.7.0 (`a_name`/`b_name` on pairs over a sequence, ride-along before, now refused: in the version history); headless 919.9 of 920, default 1,459.6 of 1,460 |
| FX9 | 2026-11-22 | 7 (review: wheel deltas lost within a frame; the tree entry shown and mislabelled in a wide frame, `lacks` true with the tree kept; EC-249's lane cropped and a negative rail; the stored top not clamped after lanes shrink; allocation on the wheel path; my look: the SVG export clipping 40 lanes to 334 px) | 12 (keys, rail drag, autoscroll, EC-247 ticks, focus, the announcement, the "+n" headers, `track.overflow` a param, bounds in `SpecBounds`, `bandTop` as a key, stacked blocks per block, settle batching: TODO) | 0 | the orchestrator: GFF3 TP53 light and dark (capped at 240 px, the rail, every lane in the SVG; one wheel tick scrolls only the lanes: my first sighting was five ticks past the end) | sizes raised by the author (930, 1,470; the pan's fit 55 KB); headless 923.3, default 1,464.2, minimal 1,039.3 of 1,040 |
| FX10 | 2026-11-22 | 6 (the author's: the lanes rail wheel-only, a follower's window shrunk by a body drag; review: dead `lanes.*` kinds, Page keys and wheel on the rail missing the band, a stale press when lanes change mid-drag, Escape pinning an unset top, a floored window clamped by its true width (its test passes on the parent)) | 4 (the text rail, a tracks section's scroll past `section.maxPx`, Home/End, the table-driven test: SX) | 0 | the orchestrator: dragging the GFF3 TP53 lanes rail scrolls to the end | headless 923.5, default 1,465.8, minimal 1,040.0 (of 930, 1,470, 1,060) |
| FX11 | 2026-11-22 | 8 (the author's: labels below a loading band, a scrolled band's labels unclipped; review: the reach step allocating every frame (p53's pan 55.5 KB past its bound), a latch never reset, duplicated packers, a second read budget and a silent refusal, tests passing on the parent) | 0 | 0 | the orchestrator: GFF3 TP53 loading (one jump, labels clipped), scrolled (no label past either edge) | the 3,000-contig leak check flaked under load, passes alone; headless 925.8, default 1,468.1, minimal 1,042.3 |
