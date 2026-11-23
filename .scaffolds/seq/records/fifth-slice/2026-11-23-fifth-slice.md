# fifth-slice: Genomes in Seq\* (next), version 1's fourth milestone (2026-11-23)

## Asked

- Run the fifth slice's loop (work id `fifth-slice`): the Genomes milestone, built in `packages/seq` through lanes, then
  measured, reconciled, folded back and closed. Six designs: genome coordinates, region reads, genome formats, across
  sections, the genome ring as an extension, tooltips; the Genomes dataset as their example.
- Standing instructions, as they came:
  - **Lanes**: one at a time ("slow down to just 1 active lane"), then "we can run up to 3 lanes in parallel now";
    "maybe use this more in general for well-scoped subtasks?" (Sonnet for scoped lanes).
  - **Sizes**: soft bounds, trimmed or lazy-loaded first, then raised by the orchestrator without asking; performance
    and allocation are hard gates.
  - **Fixes**: "lets have a 'fix phase' for slice 5 too": small findings in `fixes.md` for fix lanes (FX1–FX24);
    "we should not be fixing this one-by-one": a class by a rule and a sweep.
  - **The ring**: "a real window, as recommended" (the square); "make it RX3 after RX2" (kinds of links); "let's pause
    the loop after RX2 is done", then "can we resume the loop?".
  - **Scope**: the Prism swap and the old package's removal deferred ("i was planning to rewrite whole prism anyway");
    `.seq/` reconciled into a common place for scaffolds as the sixth slice; the storybook's Docs taken early ("i think
    we can write the slice 6 docs now ?").

## Decided

Kept in the designs as well (each in `state.md`, Decided, with its quote). "As recommended" marks a recommendation
accepted.

- **Gate 0** "the recommendation make sense"; the drafts "generally reasonable and we can fold in the reviews"; the
  whole-designs review's twelve recommendations ("yes with the recommendations"); three genome formats built ("Yes,
  those three"), the rest designed only; a GenBank lane ("Add a GenBank lane").
- **Stage gates**: stage 1 passed after the author walked every story ("Yes, all five" to the novice decisions);
  stages 2 and 3 "as recommended" (the ClinVar slice committed with attribution; genotype codes on tokens of their own;
  a scaled leak control in the shared harness; `formatGenomic` parsing once); stage 4 passed by the checks, its use
  gate **carried to the next plan** ("Carry to next plan").
- **Looks decided by the author**: the values layer "C with the bars"; the rail "D1 by default, but customizable to
  D2"; header actions on header hover only; tracks share the tallest height "only when close"; the dot plot's gaps
  tinted; the navigator's readout in one fixed place outside; letters centred on the whole cell; the square "a real
  window"; ring widths "Per layer, said".
- **The engine, not the game**: "the base library is like a "game engine", not the "game itself"" (the design README,
  ring-extension's Intent, the Docs' introduction).
- **Sizes raised**, each with its feature: reading 265, headless 980, default 1,545, minimal 1,090 kB, the bridge chunk
  28, the ring 45 (soft-size rule; the author's own raises quoted in `state.md`).
- **The fold-back merged**: "Apply all ten", B10 as the author refined it ("how about instead we mocked the
  interactions too during the design stage?").

## Delivered

56 lanes, about 1,050 commits from the plan (`0669f8fab`), each through a brief, an independent review, a fix round and,
where it draws, the orchestrator's look in a visible pane before the merge.

- **Stage 1, coordinates**: GC1 `e3b3a98d0` (language 0.6.0), SC `2fdc44e76`, GC2 `db2364d09`, GC3 `61e71891b`.
- **Stage 2, formats**: BZ `167d24491` (BGZF and tabix in `@molstar/io`), BE `98028d3f6`, GG `9b8cf675b`, VC
  `364acfb5e`, GB `36d87ca30`, FR `0af3f0809`.
- **Stage 3, region reads and across sections**: RG with RS `590151907`, GA2 `dc911cfda`, GA3 `64e2f87cc`, LA
  `f67063de4`, LA2 `597d42deb`, LA3 `fe30f34c6`, LA4 `344ed99f1`.
- **The viewer's looks and parts**: BS `5bf984d5b`, RL `12f58e1c7`, VL `385c98c60`, SX `7e6d0649a`, VP `ec8916c85`, PF
  `ebc7b0ad6`, HC `50286b413`, TT `5b8a02f89` (tooltips, 0.9.0).
- **Stage 4, the dataset and the ring**: GX `a1b4f4e18` (the Genomes dataset), RX `ed692a648` (the ring, 0.10.0), RX2
  `38c07b835` (the square a window, menus), RX3 `6507a81fb` (kinds of links).
- **Showcases and docs**: NS `190dfc135` (Nightingale), DX `c9d120931` (the storybook's Docs).
- **Fix lanes**: FX1 `2afb32e91`, FX2 `e678d2dbc`, FX3 `ded361691`, FX4 `5bfebf366`, FX5 `95c7104a4`, FX6 `fa175dce8`,
  FX7 `29546ad25`, FX8 `f04196e81`, FX9 `1b651968d` (the hover contract and its sweep), FX10 `38f147b49`, FX11
  `4d664bcea`, FX12 `f6f8823be` (no track an empty well), FX13 `668bd159d`, FX14 `d219f79b7`, FX15 `3f2d2bc73`, FX16
  `8c07380aa`, FX17 `1857af215`, FX18 `684407ba6`, FX19 `2b4f75eb8`, FX20 `a7ac59914`, FX21 `471014f35`, FX22
  `b508aad27`, FX23 `9c4b1f910`, FX24 `b386475be` (the novice pass's bugs).
- **Language 0.10.0**: 0.6.0 assemblies, 0.7.0 inline data by name, 0.8.0 points, 0.9.0 tooltips, 0.10.0 a ring's
  `center`; no migration (pre-release).
- **Novice passes** (two, Sonnet): 36 findings, 5 bugs (2 fixed by FX24, 3 not bugs), 14 design questions, 15 looks
  (`novice-findings.md`).
- **M5** (`measurements.md`): every frame target met at dpr 2 but the 200-copy split patch (median 16.2 ms against 16)
  and one Signal zoom frame (25.7 ms); real vsync without drops; the bigWig pan measured (M4's gap).
- **Checks at the close** (load 1.8 at the start): quick 5,088 tests (`packages/seq`, `cli/seqx-render`, core's `util`),
  medium 676, `tsc -b`, the test typechecks, lint, format, the `.seq` checks, `pnpm build`, `check:workspace`, the
  allocation tiers (core 245, Storybook 22). Sizes **260.2 / 974.3 / 1,539.7 / 1,086.9 kB** against 265 / 980 / 1,545 /
  1,090 (reading a spec, headless, the default and the minimal viewer); the ring 43.9 of 45. Heavy 173 of 176: three
  tight bounds failed under the lane's own contention (the interval lookup, the pairs pyramid's slices, the field
  theme's 1.1× ratio) and passed alone at load 2.9 (B7).

## Learned

The fold-back's ten (B1–B10; `foldback.md`, "Apply all ten"; in `principles/README.md`, "Folded back from the fifth
slice", the design README, the viewer `_spec` and `loop.md`'s close section):

1. A lane's model follows its scope (21 Sonnet lanes, none reverted).
2. A repeated finding is a class: its rule and a sweep (the hover contract, the card over blocks).
3. A brief's items are checked as landed, not as reported (`section.center` and `PairsReader` were reported, not built).
4. The orchestrator's files stay out of lanes' reach; it stages by path (a lane overwrote its checks script; its own
   `commit -a` swept a temporary file, reverted).
5. The engine, not the game: application-level work labelled.
6. A novice finding is confirmed before it is placed (3 of 5 bugs were not bugs).
7. Tight heavy bounds run isolated (two failed under the lane's own contention at load 1.2).
8. A plan's size line is a budget, not "not grown" (the minimal viewer +150.7 kB).
9. A test over lazy parts awaits them (the hover card's flake).
10. The picture is agreed before it is built: an interactive mock at the design stage, thrown away at the merge.

By kind: lane reviews and looks **427 bugs, 69 missing designs, 42 looks** in 56 lanes (53 with bugs). The author's
use found what nothing else did: almost all of it interaction on the ring and its square, imaginations that differed
from the designs' text (B10).

## Left open

- **The use gate** and the novice passes' **14 design questions**: carried to the next plan (the author's word for the
  gate; the questions with it).
- **The minimal viewer's +150.7 kB**: a question for the next plan (accept as the budget, or trim).
- **Carried from `fixes.md`**: the drag rule (a drag's work bounded by the screen; the PAE drag, the 3,000 contigs); the
  stage 1 novice decisions F3, F9, F11, F14.
- **The fold-back's code changes**: `features-planned` over a brief's items (B3); the heavy bounds isolated in
  `vitest.config.ts` (B7: the interval lookup, the pairs pyramid, the field theme's ratio).
- **Feature rows**: five moves out of `deferred` (two assemblies on one axis, layers across sections, variants by
  sample, gene hierarchies, the name index) are the author's.
- **In `packages/seq/TODO.md`**: M5's leftovers (the split patch, the Signal zoom frame, Cmd+wheel on the ring), the
  novice passes' looks, the ring's wheel direction, a canvas–SVG twin sweep, the hover card's place, coarse row levels,
  a bare number as a block name, four files over 800 lines.
- **Deferred, separate work**: Prism on the new Seq\* (`sixth-slice/drafts/prism-on-new-seq.md`), then removing
  `packages/sequences` and renaming `seqx` to `seq`.
- **Next**: the sixth slice (scaffolds, draft), then version 1's close-out.

## Self-check

The principles changed at this close through the fold-back. Each rule answers a failure the slice measured: unlanded
items, lanes reaching the orchestrator's files, tool artefacts taken for bugs, bounds failing under contention, a size
line nobody could meet. Its strongest change, B10, answers the one kind no rule had caught: the author's picture of a
view differing from the design's words. The negation holds: a mock agreed at the gate still pictures only what the
design names, and real data at scale will show what a crude mock cannot. The use gate stays, carried, to catch that.
