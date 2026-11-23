# Finalize's reconcile, prepared (SR, 2026-11-19)

Lane SR's read-only list (`briefs/sr-reconcile-prep.md`), at `packages-seq` before RK and CS merged. Finalize step 1
applies it; "after RK" and "after CS" rows are checked again then. Its design questions are in
`drift/finalize-config.md`.

| List | Now | After RK / CS |
|---|---|---|
| 1 Configuration tables | 17 | 2 after CS |
| 2 features.md | 41 rows | 3 after RK, 1 after CS |
| 3 status.md | 7 | 1 |
| 4 edge-cases.md | 6 (one a missing section, about 50 ECs) | 2 |
| 5 drift | 1 note, 3 edits | 0 |

## 1. Configuration tables

Viewer (`packages/seq/viewer/_spec/README.md`, Configuration) against `ViewerConfig` (`viewer/src/config.ts`):

1. `ui.canvas` and `ui.history` have no rows (`viewer-spec.ts:80-109`; the prose covers them). Add.
2. `bounds.themeDepth` (4, `seq.lang.theme.depth`; `core/src/lang/bounds.ts:30-41`) has no row. Add.
3. `RenderSvgOptions` (`core/src/render/export.ts:264-330`: `legend`, `cvdSafe`, `folds`, `navigators`, `navigator`,
   `window`, `viewport`, `viewports`, `idPrefix`, …) has no table though the designs cite "option of `renderSvg`".
   Add a Call options table to the core `_spec`.

Core (`packages/seq/core/_spec/README.md`, `RenderConfig`) against `core/src/render/config.ts`:

4. `heights` lacks `rows` 240 and `logo` 48 (config.ts:46-47, 163-164).
5. `text.heights` lacks `rows` 18 and `logo` 24 (config.ts:251-252).
6. No rows for `rows.labelMinPx` 8, `rows.pinInPlace` true, `rows.columns.{floorPx 32, gapPx 6, dataMinPx 200,
   titlePx 16}` (config.ts:117-125, 221).
7. No row for `imageBytes` (64 MB, config.ts:116).
8. `navigator.*`, `pairs.*`, `grid.*`, `text.*`, `rings.*`, `emphasis.highlightTint`, `trackPadding`, `trackShape`:
   match. `ResolveConfig` (core `_spec`:815): matches.
9. `examples.md` (Configuration's preamble) promises the options "in the package's Configuration table", but
   `packages/seq/examples` has no `_spec`: create one with the generator and record-script options, or reword
   (`drift/finalize-config.md` 6). `chromosome.*`: after CS.

Designs against the code:

10. grids-and-rows.md:553 `rows.clades.cutLine` not built (no item, no drag). Drift 2.
11. grids-and-rows.md:554 and layouts.md `reveal.opensClades`: not an item; the code always opens
    (`viewer/src/behaviors/rows.ts:271`). Drift 2.
12. grids-and-rows.md:549 `rows.reduce`: not an item; the reducer is fixed (`core/src/render/rows.ts:1336`). Drift 4.
13. navigators.md `navigator.idle`: not an item; idle is quiet only. Drift 3.
14. grids-and-rows.md:549 spells `rows.pin.inPlace`; the code `rows.pinInPlace`. Drift 1.
15. layouts.md's `text.*` rows lack `text.scale` 12, `text.numberGap` 8, `text.heights` (config.ts:243-252); its
    `figure` preset is not built. Drift 5.
16. themes.md: ESPript `:gaps`, `:similarity`, `:mode` and `themeDepth` built; fine.

After RK: `ViewerConfig.cursor` and `view.keys` (Page, Cmd+arrows, Enter, Space). After CS: the `chromosome.*` row's
home.

## 2. features.md (line numbers as of `4f2da9db5`)

- 46 Numbering schemes says "built:", not a status word: `done` for labels, the offset and segments as remainder.
- 45 Indices and labels: the Antibody view exists (DA): drop the remainder.
- 47 Cyclic systems: `done` (DP).
- 32 Language as data, 35 Validation: the Antibody remainder: verify in `examples/src/antibody` and drop.
- 37 Encodings: MessagePack round trip of every view's fixtures: verify.
- 51 Direction on every run: "reversed pairs rows not" stays.
- 56 Pieces and layout maps, 57 Letters through a map, 58 A section's direction, 59 Collapsed introns: `done` (AX, SL).
- 66 Ligands: `building: chips and arcs (DC, HB)`; verify the remainder against `examples/src/complex/ligands.ts`.
- 67 Alignments, 68 Trees: `done` once EC-38 to EC-44 have outcomes.
- 69 Annotations: "building: ranges, values, pairs, grid, rows; matrix deferred".
- 72 Formats in the core: add aligned FASTA, Stockholm, records (tsv, csv, json, msgpack, gff3), alignment; A3M and
  Clustal remain.
- 75 Structure files: stays planned.
- 77 Example datasets: building, p53's chromosome not (then `done` after CS).
- 92 Snapshots: SY built the player and a step's view; verify the p53 story.
- 104 Self-describing sources: "assembly not" only.
- 116 Format runtimes: list `records`, `alignment`, `aligned-fasta`, `stockholm`; the conformance suite: verify.
- 125 Letters: `done` once EC-116/117/120–122 have outcomes.
- 126 Ranges: line styles and arrowheads built; verify split and collapse.
- 127 Values: `done` (line and area, EC-57).
- 129, 131, 133–140 (variants, pairs, grid, logo, secondary, ruler, frames, rows, clades, lazy rows): `done` once their
  outcomes exist; 140's cursor part after RK.
- 142–148 (themes, schemes, shared domains, legends, color-blind substitution, column measures, sizes): `done` (TH, GD).
- 157 Tracks layout: groups not found: verify.
- 158–165, 168 (text, alignment figures, blocks, rings, navigators, ruler rail, layout interface, copies): `done`; 163
  per NV's remainder.
- 176 Selection: only codon snapping (EC-149) remains.
- 177 Keyboard cursor: after RK.
- 184 Copy the selection's letters: `done` (`viewer/src/commands.ts:45`).
- 193 SVG and PNG: PNG built (OP); 194 Render CLI: `done` (`cli/seqx-render`).
- 196 UX foundations: motion built.
- 203 Ports: add the animation port and the `hoverCard` behavior's port.
- 206, 211: still true (one registry kind) until MV.

## 3. status.md

1. "Built" stops at the second slice: add the third slice's lanes and merges (from `state.md`).
2. "Next" is the second slice's: RK, CS, MV, stage 7, finalize; the fourth slice; Genomes the fifth.
3. Language version says 0.3.0; the code is 0.4.0 (`core/src/lang/language.ts:19`).
4. Open item 1 and 4 stale (navigators and stories built); the third slice's open questions missing.
5. Measured sizes are the second slice's (reset at finalize).
6. Provisional 3: the source protocol agreed with files.md.
7. `packages/seq/TODO.md` Next (lines 7-15): stale the same way.

## 4. edge-cases.md

1. No "Outcomes of the third slice": write it (kept / changed / not met, each with its test). ECs cited in
   `packages/seq` with no outcome row: 11, 20, 33, 35, 38–41, 43, 44, 50, 51, 57, 59, 62–65, 67, 87, 92, 94–97,
   99–103, 106–109, 112, 113, 115–118, 128, 129, 132, 137, 138, 145, 195, 204, 206, 209, 232, 233. Tests naming an EC:
   EC-33 `render/_test/axis-perf`, EC-38/39/40/43 `data/_test/alignment-readers`, EC-57 `render/_test/values-line`,
   EC-67 `render/_test/secondary`, EC-87 `data/_test/records`, EC-100 `globins-stories`, EC-112/113
   `viewer/_test/view`, EC-128/232 `render/_test/ranges-lines`, EC-170 `data/_test/{store,keys}`; the rest need their
   test found by behavior or "not met".
2. The second slice's EC-90 row and its section's intro name `packages/io/src/repr/_test`, which does not exist: the
   codec's tests are `packages/core/src/repr/codecs/msgpack/_test`.
3. EC-85 (a step's view) now met (`core/src/story/state.ts:28,100`).
4. EC-126 "(in part)": complete with `cli/seqx-render/src/_test`.
5. EC-178, EC-149 stay open (TODO.md). EC-135 built (`theme-language.ts:138`); EC-130 check `secondary.test.ts`
   names it.
6. After RK: EC-154, EC-155. After CS: EC-169 to EC-175.

## 5. drift/themes-spellings.md

1. themes.md:152-162 still writes `=`, `!=`, `<`, `<=`, `>`, `>=` and `(>= (field plddt) 90)`: the words (`eq` … `ge`)
   as built and decided, symbol aliases later.
2. themes.md:212, 451, 562 write `non-positive`: `non_positive`.
3. A data column named `kind`: the core `_spec` says it (48, 58); the builder's docs do not. Add one sentence there, or
   accept the `_spec` as the home. Then the note is deleted.
