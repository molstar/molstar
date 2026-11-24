# The close-out's novice pass (2026-11-24, `packages-seq` at `99a61758d`)

Brief: [briefs/novice-pass.md](briefs/novice-pass.md). Three Sonnet agents, one group each, run one at a time in a
visible pane (a first start stopped on a hidden pane, as the brief says): **A** legends, layers and track groups (after
FX32), **B** finding things and reading values, **C** the genome and the ring. Every task finished but four: hide one
colour and a merged row's layers at once (A13, A3), a residue by number (B2), every WRAP53 transcript on the canvas
(B7).

**Reproduced by the orchestrator** in the pane before placing: C1, C3 (not reproduced: a click opens the third level),
C2 (not reproduced: Back stays disabled after ring clicks), C4 (by design: a piece shows the pair), C5, B4, B7; B1 read
against `search.md` (a step adds history only when it pans: by design); A6 is `packages/seq/TODO.md`'s FX29 followup.

**Counts**: 43 findings; 9 bugs (2 not reproduced, 1 by design), 12 missing design, 22 looks.

## Placed

| # | What | Kind | Placed |
|---|---|---|---|
| C5 | A partner or a piece opened from the stem menu is no history step: Back stays disabled (`ring-extension.md`, The stem menu: "shows the pair, one history step") | bug | FX34 (`293155bb1`) |
| B4, A8 | The stories' SVG twin misses state changes: search's Clear leaves its tints; a group folded on the canvas stays open in the twin after an undo | bug | FX34 (`293155bb1`) |
| A6 | A folded group's summary layer titled by its raw id ("ranges#summary", "secondary#summary") in its legend row and eye, right after a fold its own name | bug | FX34 (`293155bb1`, from the TODO) |
| A4 | Long Legend's track named "letters" titles its row "letters letters" | bug | FX34 (`293155bb1`) |
| B7, B10 | GFF3 Find By Name: the canvas's Transcripts track capped with an inner scroll while the SVG draws every lane (By Data's case before FX30); GTF's result line cut at the right | look | FX34 (the stories at `height: 'content'`; the line wraps) |
| B12 | Manhattan's card "chr 11,984, index 11983": an index, one off the position, jargon to a reader | look | FX34 (the card's position only) |
| C1, C9 | Nothing on the ring opens a chromosome by click or double-click; the stories never name the gestures (Cmd/Ctrl+wheel, the label's "Show only", the square) nor the partner keys | missing design | D1 (agreed) → FX35 |
| A13 | A legend swatch's click highlights its colour; a reader expects it to hide that colour | missing design | D2 (agreed) → FX35 |
| B8 | No way from the track's menu to lift its cap (every lane shown): `height: 'content'` is the spec's only | missing design | D3 (agreed) → FX35 |
| B11 | Manhattan unreadable: pale small dots, no y axis, the lines' values (7.3, 5) only in the description | missing design | D4 (agreed) → FX35 |
| C10 | A partner step barely shows on the ring (a tick, the card far off) and is said in the stem level's bins | missing design | TODO |
| B2, B3 | Search: a residue number gives no word; IUPAC's nucleotide letters applied to a protein without a hint | missing design | TODO |
| A3 | The Ligands story's seven ligand layers have no legend (no per-layer or row hide) | missing design | TODO (the story) |
| A7 | Folding a group by its header is no history step; only the header again undoes it | missing design | TODO |
| B1 | Search steps that pan nothing add no Back step (by design, `search.md`); a reader expects Back after stepping | missing design | TODO (a hint) |
| C3, C2 | Not reproduced | — | none |
| C4 | A piece opens the whole pair (by design: "Each row and piece shows the pair") | — | none |
| the rest | Looks: C6 a third level over the ring's labels; C7 the menus' fade; C8 the story's menu far from its stem; C11 the keys panel's close scrolled away; C12 a cell's hover and header in opposite orders; C13 the page scrolls after a view change; C14 translocations in the selection's brown; C15 a ribbon's card on a second move, repeated strand rows; C16 genome names absent from headers, a tooltip broken in a number; C17 dim links in dark; B5 a card over the track beneath; B6 an origin-spanning hit hard to see; B9 Back's view late, Forward lit early; B13 "–" for minus strand, a chromosome run into its number; A1 a hidden letters layer's row blank on Legend Placements, its label silent (As A Complex's label says "1 layer hidden"); A2 legend rows reorder under the pointer, "1 layer hidden" not marked clickable; A5 "+28" expands all three legends and scrolls the panel; A9 a folded group's strip overlaps its members' marks; A10 the focus residue not named, Back's animation zooming in first; A11 the hidden-tracks tip over the first row; A12 "Collapse or expand hidden stretches" is jargon | look | `packages/seq/TODO.md` (one line) |

## Gaps for the author (decided 2026-11-24: "Agree all four")

- **D1, opening a chromosome on the ring**: recommend a double-click on a chromosome's label or arc opens it (the cut,
  one undo step, as the label menu's "Show only"), the genome's own gesture (genome-coordinates.md: "double-click on its
  label: fit to it"); and each ring story's text names the gestures (Cmd/Ctrl+wheel, "Show only", the square, the
  partner keys).
- **D2, a swatch's click**: recommend keeping the highlight (CL1's answer A, agreed) and saying it: the swatch's tip
  "Click to highlight · the eye hides the layer".
- **D3, every lane from the menu**: recommend a track menu entry "Show every lane" / "Cap the lanes" setting the
  track's `height` to `'content'` or back (one patch, one undo step), on a track whose lanes scroll inside.
- **D4, the Manhattan's numbers**: recommend a reference line's label carries its value ("genome-wide · 7.3") when it
  has a label; the values layer's y axis is a design of its own (TODO), and the story's dots take the layer's `size`.
