# fourth-slice: integration in Seq\* (next), version 1's third milestone (2026-11-21)

## Asked

- Run the fourth slice's loop (the work id `fourth-slice`): milestone 3, integration. Molecules, maps and copies,
  queries by role, structure syncing on a fake 3D viewer, track headers, seq data ports, range reads and bigWig,
  archives, services, a test consumer in place of Nightingale, extension openings, built in `packages/seq` through
  lanes, then measured, reconciled, folded back and closed.
- Standing instructions, as they came:
  - **Lane limits**: "2-3 lanes, all either heavy or light depending on context", then "also we can run 2 lanes in
    parallel again now"; the order of HP1 and ST1 "order doesn't really matter, decisions as recommended".
  - **Public API**: "you can change public api of packages/seq, core/repr, packages/ui if needed. this has no other
    consumers beyond the work here yet"; each change listed in `gaps.md`, the rest of core and io still stops a lane.
  - **Overnight**: "you can continue with slice 4 too, just don't close it until i review it in the morning"; decisions
    a gate would ask made as recommended into `gaps.md`, beside the fifth slice's stage 0 drafts.
  - **Fixes**: small findings and the author's check of the third slice collect in `fixes.md` for fix lanes (FX1–FX7).
  - **To the close**: "you can run slice 4 loop all the way to the close. I will review the foldback.md file before we
    apply it".

## Decided

Kept in the designs as well. "As recommended" marks one accepted from the loop's recommendation; a quote is the
author's own words.

- **Gate 0**: "Accept all". Stage 0's designs agreed and promoted: structure syncing (a library, not a port; DR), seq
  data ports (SP), extension openings (XD); track header parts.
- **Molecules.** Core's standard components (ML1) "Accept"; a unit's positions from 0 ("From 0, as built"); the atom
  form monadic (ML2); `block.form` per unit (MR); a pair's end on a split copy no unit names "never on a guessed copy"
  (CU); T-DM1's placeholders "keep UNL placeholders" (no CCD component exists).
- **Range reads.** The surface in `drafts/bw-range-reads.md` "Accept", 64 KiB pages, "fail on change".
- **Drifts asked at the next gate**: MX's "ok for the drifts" (a section's `via` over every path; a component query
  pending); AC's and ST1's "Yes, as recommended" (`<alignment>/columns`; a map onto a described unit an overlay edge).
- **Stage 1's gate**: "Pass, FX3 fixes times"; submenus from a select param ("Commands own choices"); T-DM1's ids and
  the pairs layer's numeric columns "Both to FX3".
- **Sizes as bounds, moved with the feature behind them**: reading a spec 230 kB, headless 820 then 830 ("Move to 830 /
  930 kB"), the default viewer "Move to 1290 kB" then "Move to 1330 kB", the minimal viewer 940 kB (gaps.md).
- **`labels.minChars`**: an ellipsis inside before a label is dropped, at least 4 characters, configurable.
- **Module size**: SZ splits the four files over 1500 lines ("Yes, add all three", carried from the third slice).
- **The overnight gaps**: "for gaps go as recommended. download what you need for fixtures etc". Every `gaps.md`
  recommendation taken: a hover through the structure link at most 2 kB a lit unit; `header.slots`' rationale amended
  (one slot in the measured column); HD's drift (sections over a source's complex) to the next slice's view-model work;
  the viewer folder's budget unchanged; F1 and F22 to FX5, F5 kept (the host owns its camera); fetches allowed (SV's
  services, CDK2's PAE, each with its licence).
- **Public API changes accepted** (gaps.md, each "Keep? Recommended yes"): ST2's `JumpKind` `zoom` and `move`; ST3's
  `onHover` with `at`; HD's generic `withAdapter`; EX's extension surface (`SeqExtension`, `defineExtension`,
  `<id>:<name>` types, `extend()`, `KindSpec.extensionParam`/`fallbackParam` in core's repr); LE's `LegendSelection`
  and the compact form; CN's `navigator-picture` entry and `NavigatorModel.xRange`/`positionAt`; FX5's
  `LayerType.mark`/`ownHover`, `greens`; PC's `ranges.connector`; FX6's `submenuRestMs`; HT's `hiddenTracks`,
  `mountHiddenTracks`, `track.show {all}` and `@molstar/ui`'s `ContextMenuItem.items`. Every item card names its track
  (FX4, as recommended); a reveal goes out through the link (ST2, as recommended).
- **Hidden tracks are the application's**: "hiding a track moves the UI around. we need a better way" … "give this to
  control of the app layer and that is then responsible to showing the hidden tracks and ability to restore them"
  (track-header-parts amended; lane HT).
- **The use gates passed**: "i have looked at things here and there in the storybook, so for now user gates passed --
  we will need to do followups anyway as small things appear".
- **The fold-back merged**: "All 11" (below); the edge-case outcomes pasted into `design/edge-cases.md` (kept 39,
  changed 2: EC-66, 235; not met 5: EC-89, 171, 186, 201, 204).
- **For the fifth slice, decided during this one** (its `state.md`): Gate 0 "the recommendation make sense"; the six
  drafts "generally reasonable and we can fold in the reviews"; **three genome formats built** ("Yes, those three":
  BED/bedGraph/BEDPE, GFF3/GTF, VCF with BGZF and tabix in `@molstar/io`; bigBed, 2bit, `.gzi`, PAF, PIF and chain
  designed only); the whole-designs review's **twelve recommendations** ("yes with the recommendations"); **the
  fold-in's three** ("yes": gaps become Decided at promotion, the plan reconciled with the drafts, GTF names indexed).

## Delivered

35 lanes (36 merges, SZ in two parts; about 410 commits from AR's merge `7de920adf`), each through a brief, an
independent review in its own worktree, a fix round and, where it draws, a look before the merge.

- **Stage 1, molecules, maps, copies, queries**: AR `7de920adf`, NG `cb4c068b6`, ML1 `7a6e3e4a0`, DP `b59c12839`, BW
  `7de7f14ca`, ML2 `59fca56ee`, MC `9a329fd6c`, MK `b3b8912be`, QE `7ba97a5fd`, AC `f46985750`, MR `926ca4c7a`, CU
  `e627da916`, MX `97df4914d`.
- **Module split**: SZ part one `cb35a022b` (tracks.ts 2871 → 865, rows.ts, ranges.ts), part two `eee4718b4`
  (viewer.ts, models/view.ts); the file-size baseline empty.
- **Stage 2, structure syncing and track headers**: HP1 `42393a740`, ST1 `5dac66471`, HP2 `034ff29f0`, ST2
  `a5c8e1e0a`, ST3 `3e0144203`, HD `79a92c8fb`, HP3 `94722ac37`.
- **Stage 3, data and references**: SV `17b52cab9` (UniProt, PDBe SIFTS, AlphaFold DB from 28 kB of licensed
  fixtures); AR, NG and BW above.
- **Stage 4, extension openings**: EX `cf8331176`, LE `dc7bfa7b0`, CN `7a6abc713`.
- **Use-driven**: NP `34a944a17` (the novice pass's 29 findings), PC `e3c94d466` (a protein over its coding DNA, the
  connector), HT `5bf619ac6` (hidden tracks, flyout submenus).
- **Fix points**: FX1 `51d4fc141` (the third slice's 100 ms zoom frames: rows levels built across frames), FX2
  `05fe74453`, FX3 `00b7a5d24`, FX4 `29ea310a4`, FX5 `fd3b13b01`, FX6 `c1fc7530a` (the header menus own their
  events), FX7 `23480e1fb` (32 B a step boxed in the header actions). `fixes.md` empty at finalize, its rest in the TODO.
- **Language 0.5.0**: the slice's one bump (molecules, extensions, legends' and navigators' params, the connector).
- **Gates**: stage 1 (heavy 122; sizes 221.3 / 824.8 / 1289.9 / 922.6 kB), stage 2 overnight (heavy 123 at load 3.7;
  222.4 / 826.5 / 1289.3 / 925.3), stage 4 overnight (heavy 123 at load 2.6; the chromosome pan ratio 1.97× against
  1.9× once at load 4.9, green after; 223.1 / 806.0 / 1296.4 / 929.5), stage 3 after FX7 (medium 241, heavy 123 at
  load 3.9, both tiers, `pnpm build`, `check:workspace`).
- **Checks at the close**: quick 4,340 tests, medium 241, heavy 123; the allocation tiers (core 146 cases, Storybook 9)
  at 0 bytes a warm frame; sizes **223.4 / 810.3 / 1304.2 / 936.2 kB** against bounds of 230 / 830 / 1330 / 940
  (reading a spec, headless rendering, the default and the minimal viewer). One quick test flaky at the final gate:
  services-app's "draws what it has, waits for SIFTS" failed 2 of 3 full runs, green alone; to the fix lane **FX8**,
  running at the close.
- **M4** (Chrome at dpr 2, load 1.3 to 1.4, 300 frames a case, main-thread CPU with `requestAnimationFrame` flushed by
  the harness): the 2OCJ hover lighting 4 chains 0.7 / 1.2 ms median and p95 (cycling three rows: p95 2.6, max 7.5,
  over the 2 ms hover target, under a frame); 200 copies merged and split, hover, pan, zoom and row scroll all under
  1.2 ms p95 (max 5.7); gene by region and the chromosome's levels under 0.5 ms p95. **Not measured**: a bigWig pan (no
  story reads one), raster and compositor time (no real-vsync run), the split patch's click to drawn (one sample of
  33.6 ms with timer slack, against 16 ms; re-measure); 200 copies are CDK2's, not 2OCJ's.

## Learned

The fold-back's eleven (B1 to B11; `foldback.md`, merged "All 11"; applied in the fifth slice's `loop.md`, B1's
overlay rule in the viewer `_spec`, B6's and B7's checks built by the fifth slice's first fix lane):

1. A look is a gesture, not a screenshot: real pointer paths for each matrix cell; a DOM part over the canvas proves
   the view receives none of its events (`isViewEvent`).
2. The orchestrator looks itself at every merge, named in the findings row; no look, no merge.
3. The lane limit counts every running agent, listed with start times in `state.md`.
4. A pane's screenshot can be one step stale: two agreeing pictures, or one plus the DOM or pixels.
5. Timed checks run on a quiet machine (no lane, `uptime` under 3); a once-under-load failure rerun alone and watched;
   the allocation tier at the check of a lane on a frame's path.
6. Defaults sections by a check (a design naming a layer, layout or part; a baseline that only shrinks).
7. `features.md` stale rows by a check (28 found at this close, as the third slice's B5).
8. The review checks each `EC-<n>` the brief names against the lane's test titles (EC-235 built, unnamed).
9. The strong reviewer for every lane but a pure split.
10. The novice pass after the last lane that draws, its script with pointer paths onto menus.
11. Finalize starts with the unticked plan items listed and each run, cut or carried.

By kind (summed from the fold-back's table): **151 bugs, 51 missing designs, 34 looks**. Lane reviews found 126 of
the bugs in 33 of 35 lanes (the two SZ parts none); the orchestrator's looks 16 of the 18 lane-level looks; **seven of
the author's eight bugs were a menu's interactions**, the eighth a hover, all after screenshots had passed them.

## Left open

- **FX8**: the flaky services-app test, running at the close.
- **M4's carried measurements** (fifth slice): a bigWig pan in the browser (needs a story), real-vsync frames with the
  pane in view, the split patch's click to drawn.
- **The fold-back's checks** (B6, B7): the fifth slice's first fix lane.
- **Placed in `packages/seq/TODO.md`**: sections over a source's complex (HD's drift, the next view-model work); the
  structure section's fit in the services story; the quick lane's folders (viewer about 2.8 s, fewer jsdom files);
  looks not taken in a browser (the alignment's last column, the toolbar's dimming, the Antibody's rail); a rows gap
  mark; the host command's change kind; the pairs layer's numeric columns; the `pairwise` and `alignment` map kinds;
  tests watched at the gates (Prism's `commands`, `collapse`, the services story).
- **The novice pass's design questions** F3, F9, F11, F14, F25: the fifth slice's gate.
- **`render/navigator.ts` at 1473 lines**: a split before the fifth slice's navigators grow it; 18 files over 800.
- **Edge cases not met**: EC-89, 171, 186, 201, 204; eight kept without their number in a test.
- **`sideEffects`**: still in the TODO.

## Self-check

The principles changed this close through the fold-back: the look as gestures, the orchestrator's look, the counted
lane limit, the quiet machine (amending "The heavy lane runs alone"), Defaults and `features.md` by checks, EC ids
checked in review. They are consistent with what was built: each answers a failure this slice measured (eight author
bugs past screenshots, three lanes past a two-lane limit, every timing miss under load, no Defaults heading in any
design, 28 stale rows), and each lands as a loop step or a check. The fold-back's own negation holds: the menu bugs
were found by use with intent, which a scripted gesture list only approximates; the fifth slice's use will again find
faults off the matrix in whatever it adds. One tension carried: the novice pass's browser tool cannot send a pointer's
path, so B10 needs a tool that can.
