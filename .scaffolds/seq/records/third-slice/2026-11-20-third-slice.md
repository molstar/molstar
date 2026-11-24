# third-slice: figures in Seq\* (next), version 1's second milestone (2026-11-20)

## Asked

- Run the third slice's loop (the work id `third-slice`): milestone 2, figures. The axis, themes and legends, the
  layouts (text, blocks, rings), pairs, grids and rows, the navigators and stories, the datasets and their views, output,
  accessibility and use, built in `packages/seq` through lanes, then measured, reconciled, folded back and closed.
- Standing instructions, as they came:
  - **Lane limits**, changed five times: "let's keep 1heavy+2 light or heavy active", "slow down to just 1+1 lane",
    "slow down to just 1 lane after current batch is done", "let's pause after current tasks are done" (usage),
    "default is 1+1-2 lanes now"; for the fourth slice "2-3 lanes, all either heavy or light depending
    on context".
  - **Use as lanes merged** in place of a fixed early gate: the author used the Storybook as lanes merged and reported
    what they saw; the loop turned each finding into a decision and a lane. The novice pass came before it.
  - **Defaults**: "really good defaults for everything"; "configurable? reserve by default"; a ruler "not numbered" as
    the default. A default is a decision with a picture, not an accident.
  - **Nightingale** dropped to a test consumer of the resolved view (nothing of it is a dependency).
  - **The fourth slice planned early**, while this one ran ("you can also run another 1-2 lanes", the drafts of
    structure syncing and header parts as a design-only exception).
  - **Sizes**: "reduce render first, then move the size targets"; "give the bounds a little wiggle room".
  - **Worktrees**: "part of the slice finalize should be cleanup of the worktrees".
  - **The fix phase**: findings that need no design collect in the fourth slice's `fixes.md` and a fast lane takes them
    at the fix points.
  - **Module size**: "do we have guards in our principles that check a single file size with some reasonable limits?
    so we don't end up with 10k line superfiles?"

## Decided

Kept in the designs as well (decisions are copied, not moved). "As recommended" marks one accepted from the loop's
recommendation; a quote is the author's own words.

- **Layouts and labels.**
  - Plain defaults (`pieces` whole, `direction` the writer's); a collapsed intron zero width on every layout;
    a one-position block is a chip in every layout (as recommended).
  - One label rule for the ring and the line: "Inside if it fits, else beside"; ink "contrast inside, text beside".
  - A row column's cut value whole only under the pointer ("the full labels should only display when you actually hover
    on the species column"); "Leave it out of exports" for the rail's row.
  - A label that does not fit keeps at least four characters before its ellipsis, configurable
    (`labels.minChars`; fourth slice).
- **Heights.** A track's height fits its content in steps, tracks of one kind in a section sharing one: "configurable in
  the spec" and "overrideable via a port too" (`heights`: content, uniform, fixed, `ext:`; the `trackSizing` port). The
  defaults audit's numbers and all four smaller proposals as drafted.
- **Legend.** `channels: reserve | collapse`, reserve by default ("configurable? reserve by default"); the clade
  visuals "Light wedge"; "2 columns hidden" named in the chip, which stays (as recommended).
- **Navigators and rails.** The ruler rail is the default linear navigator; the rows rail idles on `whole` alone, as
  its design says (corrected); a ruler's `units: none | bases`; a plot's box fixed, a double-click "Zoom a step about
  it", "Never scroll the page".
- **Rows and grids.** "Small blocks exact" (under 64 rows counted exactly after the PCG fix); "Caret, not first copy"
  for a pair's end on a split copy, then CE parked ("should we wait until we can address copies later?").
- **Numbering.** `unnumbered: mark | fallback`, from "customizable, ... 'not numbered' as default".
- **Linked pairs.** "Behavior, on by default": pairs with an end in a hovered or selected range are lit; pairs
  answer `items`, a highlight channel's `elements` (the public API, accepted).
- **PCG32.** "Fix PCG in core": a real PCG32 in `@molstar/core`, Mol\*'s distinct colors regenerated once.
- **The minimal viewer** injectable and opt in ("Injectable, opt-in"); the default composition's size left as a guide
  then ("i guess that's fine").
- **Sizes as bounds**, after RR: 195 / 780 / 1250 / 880 kB (reading a spec, headless rendering, default viewer,
  minimal viewer), moved only by the author, as a test budget is ("give the bounds a little wiggle room");
  `sideEffects` not now (TODO).
- **Module size**: "Yes, add all three" (§module-size, `file-size.ts` with a baseline that only shrinks, the split
  lane SZ in the fourth slice).
- **The five close answers**, all fixed or accepted ("Fix all three" for `reveal.opensClades`, `navigator.idle`,
  `rows.reduce`; "Accept both" for EC-63 and EC-138's amendments).
- **Designs agreed**: the molecules design and the seq data ports ("Other recommendations from the design ok"), files
  (as recommended), structure syncing and header parts for the fourth slice ("Yes, after folding"); a Rondo-like genome
  ring as a fifth-slice extension ("Rondo-like, extension"); the gaps (palette split, `classes-safe`, the dark cell
  edge, `view.motion`) "Yes, as recommended".
- **The fold-back merged**: "add a note about the sizes to the slice 3 foldback too. and then you can apply the
  foldback." Applied in `c435eb293` and `116cebde1`, §module-size among them.

## Delivered

About sixty lanes (73 merges, about 650 commits from C0's merge `6d4b96898`), each through a brief, an independent
review, a fix round and, where it draws, a look in light and dark before the merge.

- **Axis, themes, small layers** (stages 0 to 1): C0 `6d4b96898`, AX `1e7e87870`, SL `2457f3ab9`, TH `1772436b9`.
- **Layouts** (stage 2): TX `c6282f553`, TA `24db27a58`, TV `be56ae0be`, TE `78f20c5fd`, RB `4e9a5ded6`, G2 `e2aac7989`.
- **Rows, pairs, grids, navigators, stories** (stages 3 to 4): RW `98d1af9e2`, SY `6b36ad68d`, CH `a5252eabd`, PG
  `e99c48af6`, NV `96a649d1f`, GD `b69771348`, QT `61f351a30`.
- **Files, interaction, datasets, output** (stages 5 to 6): SO `21ac00d01`, FS-A1 `73ed34d03`, FS-A2 `1344c62ee`, the
  hover and cursor lanes (HB, DX, CR, HL, IZ, RK), the datasets (DP `160e49a2e`, DT `57244533e`, DG `33dafa378`, DC
  `068ae8f2b`, DA `c7ffb2b32`, CS `ce6566143`, MF `a455d18a7`), output and accessibility (OP `a1eb09c83`, AY
  `fc77e3c74`), the novice pass's fixes (NP `c813defa1`), HC `16112a9b7`, the edge-case tests (AT `f88979d0c`, ET
  `50956fce0`).
- **Use-driven lanes**: CF `ecb040933`, MV `9bb1b35d5`, UN `2439e9ba8`, PG-PCG `a4da637c3`, RX `c55db3098`, LP
  `099493fb3`, DB `6e290d105`, PZ `323060146`, PF `413204753`, DS `acda53035`, LX `d4fe8931a`.
- **Measure and sizes**: EA `f289450ce` (the emphasis frame's 800 B), RR `1dcb175b6` (the minimal viewer 853.3 to
  801.7 kB), the bounds in `c5764af31`.
- **Parked**: CE (a pair's end on a split copy) is **not merged**: branch `seq/ce-copy-ends`, kept with its review in
  `followups.md` until units exist in the fourth slice. Nothing else unmerged.
- **Core additions**: PCG32 (PG-PCG), the palette split, the records runtime and MessagePack (FS-A1).
- **Checks at the close**: quick 2465 tests, medium 110, heavy 109; the allocation tiers (core 135 cases, Storybook 5)
  at 0 bytes a warm frame; typecheck, lint, format, `pnpm build`, `check:workspace` and the `.seq` checks green.
- **M3** (measured on the merged tree, Node and Chrome at dpr 2): every time target passes: pan 0.1 / 0.2 ms median and
  p95 on 10,000 x 300 rows, row scroll 2.1 / 2.9, zoom 1.2 / 2.5 (6.2 worst), a ring's pan and zoom under 2.5 ms,
  text reflow 1.6 / 3.3, a story's 647 frames 0.6 / 2.1, SVG export 3.3 ms, PNG 4 to 6 ms warm. Two misses: the
  viewer's warm emphasis frame at 800 B (fixed by EA, 0 B after) and **two frames of about 100 ms in every 120-step
  zoom** of the large alignment (to the fourth slice). Sizes at M3 failed every row (175.4 / 712.0 / 1142.4 / 853.1
  kB); RR reduced `render`, and the bounds were then set from what was measured, with room.

## Learned

The fold-back's ten (B1 to B10; the work folder's `foldback.md`, merged):

1. Defaults are design, decided with a picture before the build: a Defaults section a design gate reviews.
2. Configuration tables are checked against `DefaultRenderConfig` and `DefaultViewerConfig` by a check, not kept by
   hand (SR found 17 stale rows and six items the code never had).
3. A drift note is a gap: asked at the next gate, never parked for the close.
4. An edge case a lane builds is tested by name in that lane; `ec-ids` lists the untested.
5. The Features page and `features.md` move with the merges, in the lane's own change.
6. The findings table is kept for every lane, not only the early stages.
7. The quick lane's time is the lane's to keep: it moves its own slow tests before its review.
8. Parallel lanes and their reviewers get worktrees of their own; the lane limit is one line in the loop file.
9. A seeded generator carries a reference-vector test against its published algorithm (core's PCG was biased).
10. Files and bundles grow one feature at a time: §module-size, `file-size.ts`, the split lane, sizes as bounds.

By kind, over the gates and the use (counted in the fold-back): **98 bugs, 27 missing designs, 33 looks**. Lane reviews
found 82 of the bugs; the novice pass (21 findings) and the author's use found 16 of the missing designs and 14 of the
looks, nearly all of them defaults: how a figure looks before anyone sets anything.

## Left open

- **The author's look at what changed**, before the fourth slice builds on it: the use gate counted what the author
  used as lanes merged, so a last look is at the figures changed since; its findings go to the fourth slice's fix stage
  (`fixes.md`), not to a lane here.
- **`fixes.md`'s items carried** (fourth slice): ring letters off their cells, the Reverse track's empty label lane,
  dense leaders when a ring group opens, the Antibody's IMGT label as a second line, a moved label in a lower lane
  dropped, Alignment Rows And Columns collapsing when the frame resizes.
- **`sideEffects`** on the packages: not now; in `packages/seq/TODO.md` with the measured gain.
- **The 100 ms zoom frames** on 10,000 rows: build the new level's tiles across frames, keep the old level drawn.
- **Moved to the fourth slice**: the seq data ports and structure syncing on a fake 3D viewer, track header parts,
  legends' library surface, custom navigators, units and CE's copies, molecules (the `molecule` kind), the file split
  lane SZ over the four files above 1500 lines, `labels.minChars`, the test consumer replacing Nightingale (NG), the
  novice pass on a host's integration task.

## Self-check

The principles changed this close: the fold-back's rules (a Defaults section, checked tables, drift notes as gaps,
tests by name, the Features page with the merges, the findings table per lane, the lane's quick-lane time, reviewers'
worktrees, a seeded generator's reference test) and §module-size with sizes as bounds. They are consistent with what was
built: every rule answers a failure this slice produced and measured (most of the missing designs and looks the author and the novice pass found were
defaults; 17 stale table rows; a 2882-line `tracks.ts`; bundles at twice their targets), and each has a mechanical
check or a step in the loop, not a wish. The fold-back's own negation holds: this was the slice of figures, so the
Defaults section may be ceremony in the fourth slice; B1 is scoped to designs that add a layer, a layout or a part,
and the novice pass's task becomes a host's. One tension carried: the file-size check bounds lines, not complexity.
