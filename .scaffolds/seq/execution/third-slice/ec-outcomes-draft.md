# The third slice's edge-case outcomes, drafted (EO, 2026-11-19)

Lane EO's read-only draft (`briefs/eo-ec-outcomes.md`), for finalize to paste into `.scaffolds/seq/design/edge-cases.md` under
"## Outcomes of the third slice", after RK (EC-154, 155) and CS (EC-169–175) fill their rows and ET adds tests for the
built-but-untested ones. EO ran 34 of the cited files (497 passed); the heavy files, `small-review`, the CLI's and two
storybook tests it cites were read, not run. Paths relative to `packages/seq/` (EC-126's to the repository root).
Counts: kept 45 (12 with a stated gap), changed 2, not met 9. Its design questions: `drift/ec-contradictions.md`.

| EC | Outcome | Test |
|---|---|---|
| 11 | kept | `core/src/resolve/_test/encodings.test.ts` "a segment crossing the origin is split when built (EC-11)"; `model.test.ts` "an offset map onto a cyclic target wraps, split at the origin (EC-11)" |
| 20 | kept | `core/src/resolve/_test/encodings.test.ts` "maps: one-to-many backward (EC-20)" |
| 33 | kept (the 16 ms bound heavy only) | `core/src/render/_test/axis-perf.test.ts` "a pan frame lays out and draws under 16 ms" |
| 35 | kept | `core/src/query/_test/forms.test.ts` "row: an alignment's rows by key, positions of its rows system"; `core/src/render/_test/dynamic.test.ts` "a transcript keeps its item system while its domains show" |
| 38 | kept (A3M lowercase not built) | `core/src/data/_test/alignment-readers.test.ts` "Stockholm: … insert columns from RF", "aligned FASTA: lowercase-only columns are inserts" |
| 39 | kept | `core/src/data/_test/alignment-readers.test.ts` "Stockholm: … numbers from labels" |
| 40 | kept (Clustal not built) | `core/src/data/_test/alignment-readers.test.ts` "a row of another width is said and left out", "the most common width wins", "fuzz: malformed text is reported, never thrown" |
| 41 | kept | `examples/src/_test/globins.test.ts` "the lazy 10,000,000 rows: the rows in view read alone, drawn and hit"; `core/src/render/_test/source-rows.test.ts` "waits while the rows have not arrived" |
| 43 | kept | `core/src/data/_test/alignment-readers.test.ts` "a leaf matches its row by label, else by the label before its last `/` (EC-43)" |
| 44 | kept except 1,000 × 2,000 (no test) | `core/src/render/_test/rows-time.test.ts` (heavy); `examples/src/_test/globins.test.ts` "the 10,000 rows across a tile edge"; `core/src/render/_test/properties.test.ts` (3M positions); `core/src/render/_test/text-time.test.ts` (heavy, 35,000 residues) |
| 50 | kept | `core/src/data/_test/records.test.ts` "reads values sparse by position or a row each, `__proto__` and `constructor` ordinary names" |
| 51 | kept | `core/src/data/_test/records.test.ts` "range: past the letters left out or clipped, counted in one issue, hinting the other base" |
| 57 | kept | `core/src/render/_test/values-line.test.ts` "is a path a run of values …", "a color that changes does so halfway between two positions" |
| 59 | not met: two grids on one block pair merged later over earlier, the size cap per pair | none |
| 62 | kept (folding contacts by copy group untested) | `core/src/render/_test/pairs-grids.test.ts` "a cell's hover lights the arcs ending at either position"; `storybook/src/_stories/examples/_test/pairs-grids-stories.test.ts` "the contacts as arcs between the chains' blocks: a hover lights both ends" |
| 63 | not met: dots sliced at one cell (`render/pairs.ts` `sideBySide`); no glyphs, no "+" | none |
| 64 | kept (one position per heavy atom, hydrogens rejected untested) | `examples/src/_test/complex.test.ts` "copies merged ×2, the ATPs after the polymers …"; `core/src/render/_test/blocks-chips.test.ts` "a row: the chip hit, hovered and selected as its one position" |
| 65 | kept | `core/src/render/_test/blocks-layout.test.ts` "copies merged into one block labeled ×N …"; `core/src/render/_test/hover-blocks.test.ts` "CDK2 position … lights the CDK2 copies, no cyclin" |
| 67 | kept | `core/src/render/_test/secondary.test.ts` "200 identical copies draw what one does", "copies that differ in kind or extent are all drawn" |
| 87 | kept | `core/src/data/_test/records.test.ts` "one severity a code: every issue has its catalog's (EC-87)" |
| 92 | changed: text cells kept as written, not trimmed; numbers parse past spaces; a duplicate header refused | `core/src/data/_test/records.test.ts` "keeps text cells as written: no trimming …", "reads equal tables from TSV, CSV, JSON and MessagePack", "fuzz" |
| 94 | kept | `core/src/render/_test/text-layout.test.ts` "lines of whole groups of ten, a number a line, a band collapsed where no layer drew" |
| 95 | kept | `core/src/render/_test/text-layout.test.ts` "numbers at both ends: … inline group labels in the gap after their group" |
| 96 | kept | `core/src/render/_test/text-layout.test.ts` "overlapping ranges in a thin band are its density strip …" |
| 97 | kept | `core/src/render/_test/rows.test.ts` "in text (an alignment figure) every row is on every line …"; `rows-layer.test.ts` "the label column widened by the tree's width …"; `rows-columns.test.ts` "a text figure draws the same columns on every line" |
| 99 | not met: no test of ligand chips flowing with the words in a text section | none |
| 100 | kept (a collapsed clade framing its wedge untested) | `core/src/render/_test/text-layout.test.ts` "frames across tracks: … open where a line break cuts it …"; `storybook/…/globins-stories.test.ts` "family and folded" |
| 101 | kept (own units and the first alignment as default blocks untested) | `examples/src/_test/complex.test.ts` "copies merged ×2, the ATPs after the polymers" |
| 102 | kept | `core/src/render/_test/blocks-layout.test.ts` "a row of many blocks: the gaps never take more than a quarter of the width", "a stack: …" |
| 103 | not met: track groups and summaries not built (features.md) | none |
| 106 | kept | `core/src/render/_test/rows-layer.test.ts` "a field order …"; `rows.test.ts` "in tree order folds clades from `collapsed` …"; `rows-columns.test.ts` "an order by a column's field …"; `filter.test.ts` "leaves the rows out as the items are read" |
| 107 | kept | `core/src/render/_test/rows-columns.test.ts` "yield to a narrow view: right ones first …" |
| 108 | kept | `core/src/render/_test/rings-layout.test.ts` "places the section round its whole axis …", "labels outside, each hit as its item", "no two labels' hit boxes overlap" |
| 109 | kept (the scrollbar-flicker rule untested) | `core/src/render/_test/pairs-grids.test.ts` "places one strip, the plot … square cells at aspect 1"; `storybook/…/pairs-grids-stories.test.ts` "the plot's x is the letters' x …" |
| 112 | kept | `viewer/src/_test/view.test.ts` "a move before the section has a width is kept as an intent, then placed", "a section hidden for lack of width keeps its zoom …" |
| 113 | kept | `viewer/src/_test/view.test.ts` "scroll groups and follow (EC-113, EC-114)"; `follow.test.ts` "the protein follows the gene at scale 1/3 …" |
| 115 | not met: a stack's zoom about a point untested (the row's is) | none |
| 116 | not met: the letters, dots and bars steps exist, no test of "a row shows text for all its columns or none" | none |
| 117 | kept | `core/src/render/_test/axis.test.ts` "three-letter codes where the narrowest part in view has room …" |
| 118 | kept (Complex only) | `storybook/…/complex-stories.test.ts` "confidence: … cyclin A2's first 172 gray" |
| 128 | kept | `core/src/render/_test/ranges-lines.test.ts` "each style draws its mark …", "an invalid style is drawn solid and said once a layout" |
| 129 | kept (labels under an apex in reach untested) | `core/src/render/_test/pairs-grids.test.ts` "nested arcs nest", "an end past the window is taken `reach` past its edge …" |
| 132 | kept ("+n more" untested) | `core/src/theme/_test/theme.test.ts` "merges legends that say the same …" |
| 137 | kept | `core/src/render/_test/themes.test.ts` "the size channel shrinks a mark, never its lane or its hit box" |
| 138 | changed: an external hover projected onto every section; the ≥3:1 hover outline untested | `core/src/render/_test/emphasis.test.ts` "a hover from elsewhere … is projected onto every section" |
| 145 | kept | `core/src/render/_test/pairs-grids.test.ts` "a hover on a cell: a crosshair …"; `storybook/…/complex-stories.test.ts` "a rect dragged across the chain boundary is split at it" |
| 195 | not met: a streamed sequence's default numbering untested | none |
| 204 | kept | `core/src/theme/_test/theme.test.ts` "derives swatches … says hashed values once" |
| 206 | kept | `core/src/render/_test/pairs-grids.test.ts` "a stretch runs with its partner when both ends are on one strand …" |
| 209 | not met (deferred): the dot plot's y on another section kept; ribbons between sections not built | `core/src/render/_test/pairs-grids.test.ts` (the dot plot only) |
| 232 | kept | `core/src/render/_test/ranges-lines.test.ts` "`end` is on the mark holding the item's end …" |
| 233 | kept | `core/src/render/_test/axis.test.ts` "`never` keeps the stored letters on a reversed axis …"; `viewer/src/_test/layouts.test.ts` "copies a reversed DNA section's letters complemented …" |
| 85 | kept: now met | `core/src/story/_test/story.test.ts` "a view channel pans when the windows overlap …"; `viewer/src/_test/story.test.ts` "moves every leading section, each its own way" |
| 126 | kept: complete | `cli/seqx-render/src/_test/cli.test.ts` "exits 4 …", "exits 1 …", "prints its usage … (exit 2)", "a link that leads outside the root" |
| 130 | kept: fixed | `core/src/render/_test/secondary.test.ts` "secondary: the arrowhead (EC-130)" |
| 135 | not met: built (`core/src/lang/theme-language.ts:138`, `core/src/theme/schemes.ts:275`), untested | none |
| 154, 155 | pending (RK) | |
| 169–175 | pending (CS) | |
