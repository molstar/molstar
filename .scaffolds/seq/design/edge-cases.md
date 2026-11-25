# Edge cases: kept, not smoothed over

## Intent

The author, 2026-11-17: "i would like to keep the richness, the editing experiences etc. basically watch out for "edge
cases", so the 2nd implementation doesn't smooth over them".

A simpler implementation tends to become simpler by dropping what was hard: a codon split by an intron, a feature
across a plasmid's origin, a stale hover owner, a label that is not in index order. Those cases are much of what made
the previous implementation good. Simplicity here is fewer mechanisms, not fewer cases handled.

## Decision

- **The catalog** is [sources/seqstar/edge-cases.md](../sources/seqstar/edge-cases.md): 244 numbered items (`EC-n`)
  across coordinates and maps, data, format and patches, layout, rendering, interaction and editing, streaming and
  integration (EC-1 to EC-184, from the previous implementation), and the cases the JBrowse reply raised by drawing
  specs in another viewer and probing the format (EC-185 to EC-244).
- **Every item has one of three outcomes**, decided when its area is built, never by omission:
  - **kept**: a test or an example view shows it (the test name or the story cites the behavior, not the number);
  - **changed**: the new design handles it differently, and says how (confirmed by the author, 2026-11-17, "ok": EC-133, palette order by first occurrence, is
    changed by fetch-by-region; EC-55's doc conflict is settled as "drawn at the value, the domain bounds the height");
  - **dropped**: only by the author's decision, recorded in the design of its area with the reason.
- **Known flaws are fixed, not ported** (confirmed by the author, 2026-11-17): EC-130 (minus-strand secondary structure), EC-135 (ESPript and gaps), EC-151
  (cross-block text selection, deferred before), EC-181 (repeated ids across structures).
- **The reply's cases are cases, not wishes**: what it asked for that the designs take (EC-185 to EC-206, EC-235,
  EC-236) is built and tested like any other item; what is deferred (EC-209, EC-213, EC-223 and the smaller layer
  features) stays in the catalog with its trigger; EC-238 to EC-244 are tested as fallbacks a reader without them
  must get predictably.
- **An area is not done** until each of its items has an outcome. The area's plan lists them; its record says the
  outcome of each.
- **The catalog lives in the source's file, and only there**: this design refers to it, and every `EC-<n>` cited in
  `packages/seq` (a test, a comment, a `_spec`) resolves against it, checked by `.scaffolds/seq/_checks/ec-ids.ts` at the gate.
  Why: the ids are numbers only, so code cites them without a `.scaffolds/seq/` path, and an id looked up in the wrong file was
  once recorded as missing when it was not.
- **At the scaffold's close**, the kept items live on as the packages' tests and their `_spec` statements; the catalog
  goes with the scaffold, so the kept items it numbers need a home in a `_spec` first, or the citations go.

## Consequences

- The simplification shows up as fewer mechanisms carrying the same cases: one segment-table map type for every map
  kind (EC-9 to EC-32), one transient overlay for previews instead of four document layers (EC-163), one source protocol
  for streamed and inline data (EC-169 to EC-180).
- The examples are chosen to exercise the catalog: few datasets, many views ([examples.md](examples.md)).

## Outcomes of the first slice

Recorded for the items the first slice built and tested. A test cites the behavior and the item's id; the places are in
`packages/seq` (`core/src/{query,resolve,view,data,render}/_test`, `viewer/src/_test`). Every other item has no outcome
yet: its area is not built.

**Kept** (a test shows the behavior):

| Items | Behavior | Tests |
|---|---|---|
| 1–8 | `index` 0-based, `ordinal` 1-based and the default; label order (integer, no suffix, letters, `.n`); a numbering that skips zero; labels outside int32 refused, a repeated label names its first position; insertion codes; a missing end snapped inward with an info; a backwards range an error in index order; a unit's numbering before its molecule's | `query/_test/{numbering,forms}` |
| 9–12, 14 | a wrapped range is two runs, on every encoding; a wrapped offset image never runs past the end; a segment crossing the origin is split when built; a follower on a cyclic system takes the shortest arc; a cyclic label range is shorter than a lap | `resolve/_test/{encodings,model}`, `view/_test/view`, `query/_test/forms`, `lang/_test/lang` |
| 15–18, 20, 21 | codon maps on both strands with junctions inside a codon; codon status (IUPAC sets, start as M, Sec, Pyl, stops pair with nothing); a coding sequence that is not whole codons reported; maps one-to-many backward; a map that pairs nothing reported, a partial one not | `resolve/_test/{encodings,model}` |
| 19, 22–27 | at most one codon step; shortest path, equal compositions counted once; two maps joining the same systems differently are ambiguous; at most four counted steps and a bounded number of candidates; the copy rule is part of the search; a path ending at a ranged unit is clipped; a protein on a chromosome is as sparse as its coding sequence | `resolve/_test/paths` |
| 28–30 | a boundary lands the same on every encoding; closure lifts a nucleotide to its whole codon | `resolve/_test/{encodings,geometry,model}`, `viewer/_test/hover` |
| 81 | a change made against another revision is refused | `viewer/_test/spec` |
| 111–114 | a viewport keeps its center and zoom across updates and widths; scroll groups share the range on equal axes; followers never lead (changed in the second slice: below) | `viewer/_test/{view,follow,navigate}` |
| 125 | an export has its own want set and never touches the live view's; a timeout says unsettled | `data/_test/store`, `render/_test/export` |
| 141, 144 | one hover owner, a stale clear does nothing; a nucleotide lights its residue's whole codon, a split residue every piece | `viewer/_test/hover` |
| 143 | a resting pointer is hit again once per frame in which the layout was rebuilt (`pointer.rehit`); the picture follows a pan under it | `viewer/_test/viewer` |
| 169, 170, 172–174 | pending is never empty; the coarsest listed level no wider; requests ordered view, margin, prefetch; cancellation after a grace, at once when view work waits; failures retried with backoff then cooled down, malformed ones not retried | `data/_test/{store,keys}` |
| 176, 177, 179 | nodes by (source, key), a conflict keeps the held node, eviction quiet; an item system made on demand; nested sources stop at the depth limit | `data/_test/{store,regions}`, `resolve/_test/overlay`, `render/_test/dynamic` |
| 175 | streamed data never enters the spec: the store holds it and the overlay reads it, by construction | `render/_test/dynamic` |

**Changed** (the design handles it differently):

- EC-9, EC-10: in a spec a range with `start > end` wraps only on a cyclic system, and is an error elsewhere.
- EC-143: the catalog says once per frame when the view moves; built as once per frame in which the layout was rebuilt,
  which covers a pan, a zoom and a data arrival without a hit test per pointer event.

**Not met yet** (open, in `packages/seq/TODO.md`):

- EC-178: every template re-expands on each overlay version, not only the arrival's own (cheap at this size).
- Projected hover (linked-section outlines, EC-141's family): the hover model projects onto linked sections, the draw lists do not yet draw it. (Met in the second slice: below.)

## Outcomes of the second slice

Recorded for what the second slice built: editing (patches, history, previews, story steps), addressing (refs, links),
selection, highlights and the cursor, the security bounds, MessagePack, followers that lead, the view's history. Places
in `packages/seq` as above, plus `storybook/src/_stories/{features,examples}/**/_test`; the tree core's in
`@molstar/core` (`packages/core/src/repr/patch/_test`) and the MessagePack codec's in `@molstar/core` (`packages/core/src/repr/codecs/msgpack/_test`).

**Kept**:

| Items | Behavior | Tests |
|---|---|---|
| 76, 86 | nesting capped (256, queries 64), a 20,000-deep spec refused with its issue alone and no stack spent, in the whole check and a patch's; the first 10 places of a code, 100 issues a node, the rest counted | `lang/_test/bounds`, `viewer/_test/bounds` |
| 77 | a node built in code that is no record with a kind, or whose children are no list, is refused, never thrown on; a malformed op tree likewise | `lang/_test/bounds`; core `repr/patch/_test/patch` |
| 78, 80 | a patch's ops checked, applied all or nothing, the result validated, issues at the op's path; a removal leaving a ref dangling refused | core `repr/patch/_test/{patch,conformance}`, `lang/_test/patches` |
| 81 | a stale patch refused; rebased through `merge3` only when the host opts in | `viewer/_test/spec` |
| 82 | untouched subtrees shared; exact inverses; `diff` then apply, both ways, on seeded random trees | core `repr/patch/_test/{patch,fuzz}` |
| 141 (projected) | a hover and a selection on the protein drawn on the gene, lighter, in the canvas and the export | `viewer/_test/selection` |
| 148 | click sets, Shift extends the anchor's range, Cmd/Ctrl toggles; a double-click's second click selects nothing | `viewer/_test/selection` |
| 155 | key chords match their modifiers exactly | `viewer/_test/{bindings,history}` |
| 157 | highlight channels set and cleared independently; a link's hover on its own channel | `viewer/_test/addressing` |
| 159 | focus zooms to the address; reveal pans only what is not shown | `viewer/_test/addressing` |
| 162 | a bookmark is one patch, one undo step | `viewer/_test/addressing`, `features/addressing/_test/bookmarks` |
| 182 | a link hovers on its own channel, leaving the pointer's hover alone | `features/addressing/_test/links`, `viewer/_test/addressing` |
| 126 (in part) | the headless policy reads files under its root only and no network unless hosts are listed (the CLI and its exit codes wait for the render CLI) | `context/_test/headless-policy` |
| 127 (in part) | a selected item of zero length is a hairline at its place (the pinned caret waits for boundary layers) | `render/_test/emphasis` |

**Changed**:

- EC-84: a story step has no base: it applies over whatever the spec is; one that no longer applies is shown as not
  applying, never refused (stories.md); `viewer/_test/spec`, `lang/_test/snapshots`. The player is not built.
- EC-90: non-finite numbers and −0 are an extension (`0x10`), not nil, so they round-trip exactly; core `repr/codecs/msgpack/_test/msgpack`.
- EC-114: a follower leads: a move asked of it moves its leader through the follow (the author, 2026-11-18;
  `follow.followerMoves`, `ignore` the catalog's rule); `viewer/_test/follow`.
- EC-163: one transient preview (a patch over the root) and a story step over it replace the four document layers; a
  preview that no longer applies is dropped; `viewer/_test/spec`.
- EC-164: undo restores a recorded root (exact by construction), so no inverse can fail to apply; undo and redo
  through the host's history in `viewer/_test/{spec,history}`, the inverse laws in core `repr/patch/_test/fuzz`.

**Not met yet** (open, in `packages/seq/TODO.md`):

- EC-85 (a step's view): steps carry no view yet (stories.md: language fields). (Met in the third slice: below.)
- EC-149: drag snapping to whole codons over nucleotides (items snap). (Met in the close-out's CL9: `viewer/src/_test/snap-codons`.)
- EC-154: the cursor through the codon map only; block gaps and lazy rows wait for those layouts. (Lazy rows met in the third slice; block gaps in the close-out's CL9: `viewer/src/_test/cursor-blocks`.)
- EC-178 as above.

## Outcomes of the third slice

Recorded for what the third slice built: the axis (direction, letters through a map, collapsed introns), themes, the
text, blocks and rings layouts, navigators, pairs, grids and dot plots, the logo, the rows layer (columns, clades, lazy
rows), the files and their readers, the six datasets and p53's chromosome. Places in `packages/seq` as above
(`core/src/render/_test` is `render/_test`, and so on), plus `examples/src/_test` and
`storybook/src/_stories/examples/_test` (`…/` below); EC-126's in `cli/seqx-render/src/_test`. Drafted by lane EO
against the tests, completed by ET's (EC-63, 99, 115, 116, 138, 195), CF's (EC-135), RK's (EC-154, 155) and CS's
(EC-169 to 175).

**Kept**:

| Items | Outcome | Tests |
|---|---|---|
| 11 | kept | `resolve/_test/encodings` "a segment crossing the origin is split when built (EC-11)"; `resolve/_test/model` "an offset map onto a cyclic sequence wraps, split at the origin (EC-11)" |
| 20 | kept | `resolve/_test/encodings` "maps: one-to-many backward (EC-20)" |
| 33 | kept (the 16 ms bound in the heavy lane only) | `render/_test/axis-perf` "a pan frame lays out and draws under 16 ms" |
| 35 | kept | `query/_test/forms` "row: an alignment's rows by key, positions of its rows system"; `render/_test/dynamic` "a nested source is not requested until a template needs it; then its domains show through the codon map" |
| 38 | kept (A3M's lowercase not built) | `data/_test/alignment-readers` "Stockholm: … insert columns from RF", "aligned FASTA: lowercase-only columns are inserts" |
| 39 | kept | `data/_test/alignment-readers` "Stockholm: … numbers from labels" |
| 40 | kept (Clustal not built) | `data/_test/alignment-readers` "a row of another width is said and left out", "the most common width wins", "fuzz: malformed text is reported, never thrown" |
| 41 | kept | `examples/src/_test/globins` "the lazy 10,000,000 rows: the rows in view read alone, drawn and hit"; `render/_test/source-rows` "waits while the rows have not arrived" |
| 43 | kept | `data/_test/alignment-readers` "a leaf matches its row by label, else by the label before its last `/` (EC-43)" |
| 44 | kept, but 1,000 × 2,000 (no test) | `render/_test/rows-time` (heavy); `examples/src/_test/globins` "the 10,000 rows across a tile edge"; `render/_test/properties` (3M positions); `render/_test/text-time` (heavy, 35,000 residues) |
| 50 | kept | `data/_test/records` "reads values sparse by position or a row each, `__proto__` and `constructor` ordinary names" |
| 51 | kept | `data/_test/records` "range: past the letters left out or clipped, counted in one issue, hinting the other base" |
| 57 | kept | `render/_test/values-line` "is a path a run of values …", "a color that changes does so halfway between two positions" |
| 62 | kept (folding contacts by copy group untested) | `render/_test/pairs-grids` "a cell's hover lights the arcs ending at either position"; `…/pairs-grids-stories` "the contacts as arcs between the chains' blocks: a hover lights both ends" |
| 64 | kept (one position per heavy atom, hydrogens refused, untested) | `examples/src/_test/complex` "copies merged ×2, the ATPs after the polymers …"; `render/_test/blocks-chips` "a row: the chip hit, hovered and selected as its one position" |
| 65 | kept | `render/_test/blocks-layout` "copies merged into one block labeled ×N …"; `render/_test/hover-blocks` "CDK2 position … lights the CDK2 copies, no cyclin" |
| 67 | kept | `render/_test/secondary` "200 identical copies draw what one does", "copies that differ in kind or extent are all drawn" |
| 85 | met (not met in the second slice) | `story/_test/story` "a view channel pans when the windows overlap …"; `viewer/src/_test/story` "moves every leading section, each its own way" |
| 87 | kept | `data/_test/records` "one severity a code: every issue has its catalog's (EC-87)" |
| 94 | kept | `render/_test/text-layout` "lines of whole groups of ten, a number a line, a band collapsed where no layer drew" |
| 95 | kept | `render/_test/text-layout` "numbers at both ends: … inline group labels in the gap after their group" |
| 96 | kept | `render/_test/text-layout` "overlapping ranges in a thin band are its density strip …" |
| 97 | kept | `render/_test/rows` "in text (an alignment figure) every row is on every line …"; `render/_test/rows-layer` "the label column widened by the tree's width …"; `render/_test/rows-columns` "a text figure draws the same columns on every line" |
| 99 | kept | `render/_test/ec-blocks` "the chips share a line like words, and wrap onto the next line with them when the width runs out" |
| 100 | kept (a collapsed clade framing its wedge untested) | `render/_test/text-layout` "frames across tracks: … open where a line break cuts it …"; `…/globins-stories` "family and folded" |
| 101 | kept (own units and the first alignment as default blocks untested) | `examples/src/_test/complex` "copies merged ×2, the ATPs after the polymers" |
| 102 | kept | `render/_test/blocks-layout` "a row of many blocks: the gaps never take more than a quarter of the width", "a stack: …" |
| 106 | kept | `render/_test/rows-layer` "a field order …"; `render/_test/rows` "in tree order folds clades from `collapsed` …"; `render/_test/rows-columns` "an order by a column's field …"; `render/_test/filter` "leaves the rows out as the items are read" |
| 107 | kept | `render/_test/rows-columns` "yield to a narrow view: right ones first …" |
| 108 | kept | `render/_test/rings-layout` "places the section round its whole axis …", "labels outside, each hit as its item", "no two labels' hit boxes overlap" |
| 109 | kept (the scrollbar-flicker rule untested) | `render/_test/pairs-grids` "places one strip, the plot … square cells at aspect 1"; `…/pairs-grids-stories` "the plot's x is the letters' x …" |
| 112 | kept | `viewer/src/_test/view` "a move before the section has a width is kept as an intent, then placed", "a section hidden for lack of width keeps its zoom …" |
| 113 | kept | `viewer/src/_test/view` "scroll groups and follow (EC-113, EC-114)"; `viewer/src/_test/follow` "the protein follows the gene at scale 1/3 …" |
| 115 | kept | `render/_test/ec-blocks` "the block under the pointer stays under it through successive zooms in" |
| 116 | kept | `render/_test/ec-letters-density` "the steps come in order, and the letters step starts where the smallest letter fits its cell", "text for all the columns in view or none: never a lettered few, from the widest to the cells", "then dots a position, then one bar a run: …" |
| 117 | kept | `render/_test/axis` "three-letter codes where the narrowest part in view has room …" |
| 118 | kept (the Complex only) | `…/complex-stories` "confidence: … cyclin A2's first 172 gray" |
| 126 | kept, now whole (in part in the second slice) | `cli/seqx-render/src/_test/cli` "exits 4 …", "exits 1 …", "prints its usage … (exit 2)", "a link that leads outside the root" |
| 128 | kept | `render/_test/ranges-lines` "each style draws its mark …", "an invalid style is drawn solid and said once a layout" |
| 129 | kept (labels under an apex in reach untested) | `render/_test/pairs-grids` "nested arcs nest", "an end past the window is taken `reach` past its edge …" |
| 130 | kept | `render/_test/secondary` "secondary: the arrowhead (EC-130)" |
| 132 | kept ("+n more" untested) | `theme/_test/theme` "merges legends that say the same …" |
| 135 | kept | `render/_test/rows-layer` "ESPript's `gaps` reaches the column share (EC-135): …" (the scheme in `lang/theme-language.ts` and `theme/schemes.ts`) |
| 137 | kept | `render/_test/themes` "the size channel shrinks a mark, never its lane or its hit box" |
| 145 | kept | `render/_test/pairs-grids` "a hover on a cell: a crosshair …"; `…/complex-stories` "a rect dragged across the chain boundary is split at it" |
| 155 | kept (as in the second slice); the rows' page keys taken only while a cursor shows | `viewer/src/_test/cursor-rows` "with no cursor shown Page Down is the rows' scroll, Cmd/Ctrl+↓ the page's"; `viewer/src/_test/mount` "Space is the page's when it is no modifier and there is no cursor; …" |
| 169, 170, 172–175 | kept (as in the first slice); a chromosome read by region asks only for its view and margin | `data/_test/{store,keys}`; `examples/src/_test/chromosome` "opens fit to the chromosome and lays out to TP53's letters, asking only for the view and its margin" |
| 195 | kept | `resolve/_test/ec-streamed-numbering` "is stated by the model without data: …", "a label is read 1-based by default, …" |
| 204 | kept | `theme/_test/theme` "derives swatches … says hashed values once" |
| 206 | kept | `render/_test/pairs-grids` "a stretch runs with its partner when both ends are on one strand …" |
| 232 | kept | `render/_test/ranges-lines` "`end` is on the mark holding the item's end …" |
| 233 | kept | `render/_test/axis` "`never` keeps the stored letters on a reversed axis …"; `viewer/src/_test/layouts` "copies a reversed DNA section's letters complemented …" |

**Changed**:

- EC-92: a text cell is kept as written, never trimmed (the RFC's rule, [files.md](files.md): trimming loses meaning,
  EC-49); a number parses past the spaces around it; a header twice is refused. `data/_test/records` "keeps text cells
  as written: no trimming …", "reads equal tables from TSV, CSV, JSON and MessagePack", "fuzz".
- EC-138: a hover from elsewhere (a linked section, a structure) is projected onto every section as the viewer's own
  ([viewer-models.md](viewer-models.md), Projection: what linking is for), not left undrawn; `render/_test/emphasis` "a
  hover from elsewhere … is projected onto every section". Its outline holds 3:1 against the surfaces it is drawn on
  (the page, the track, the sequence's cells, the missing and pending fills), the residue and nucleotide class fills and
  a range's fill, in both themes: `render/_test/ec-hover-contrast` "seq-hover holds 3:1 against each surface and fill
  it rings". The EC asks 3:1 against the track (what the outline sits on), not every fill; a saturated data fill it
  fails (the categorical, ramp and secondary colours) gets a double outline, light and dark (decided by the author,
  2026-11-20: "Accept both"; the double outline not built yet, below).
- EC-63: kinds at one cell are drawn as slices, one a kind, each its own hit, not glyphs with "+" (decided by the
  author, 2026-11-20: "Accept both"; [grids-and-rows.md](grids-and-rows.md)); `render/_test/ec-pairs-kinds` "the cell
  is cut into a slice a kind, in order, each slice hit as its own pair". Glyphs are in `packages/seq/TODO.md`.
- EC-171: the zoom gate is the store's bound on a want's tiles (`maxTilesPerWant`): a want past it is not asked, said
  once (`seq.store.too_many_tiles`) and drawn pending, while a coarser level shows density; there is no "load anyway".
  `data/_test/store` "are bounded: lengths clip, a want over the tile cap is an issue and not asked";
  `render/_test/inline-tiles` "a view past a shorter source's end is drawn empty, not pending; past the tile bound it is
  pending".

**Not met yet** (open, in `packages/seq/TODO.md` or the next slices' plans):

- EC-59: two grids on one block pair merged, later over earlier, and a size cap per pair: not built.
- EC-103: track groups and summaries: not built.
- EC-138's double outline over a saturated data fill the hover outline fails: above.
- EC-149: drag snapping to whole codons over nucleotides (as in the second slice). Met in the close-out's CL9: a split
  codon snaps to its pieces, the intron left out (`viewer/src/_test/snap-codons`).
- EC-154: the keyboard cursor steps a million lazy rows by index, and hidden and folded rows by the nearest shown
  (`viewer/src/_test/cursor-rows` "steps by index over rows never read: …", "a fold puts the cursor on the summary row;
  a hidden row's on the nearest shown: …"), through reversed and collapsed axes (`viewer/src/_test/axis-commands`
  "through each reader: …"); it does not cross a gap between blocks of different systems
  (`viewer/src/_test/blocks` "the cursor does not step another system's positions along the section's units"). Met in
  the close-out's CL9: the cursor crosses to the next block's first position, along x and down a grid's y, skipping
  empty blocks (`viewer/src/_test/cursor-blocks`).
- EC-178: as in the first slice.
- EC-209: the dot plot's y on another section is built; ribbons between sections are the Genomes milestone's (layers
  across sections, deferred).
- A folded clade's summary row has no accessible name (`render/rows.ts` `rowItems` skips it), against
  [grids-and-rows.md](grids-and-rows.md) (Accessibility: a node announced with its state): in `packages/seq/TODO.md`.

## Outcomes of the fourth slice

Recorded at the fourth slice's close (2026-11-21) from the tests on `packages-seq`. Counts: kept 39 (8 not named by a
test; 4 with a stated gap: EC-51, 67, 166, 238–244), changed 2 (EC-66, 235), not met 5 (EC-89, 171, 186, 201, 204).
Paths relative to `packages/seq/`.

| EC | Outcome | Test |
|---|---|---|
| 1–8 | kept, on a structure's author numbering too: negative and zero numbers, insertion codes, a repeated label naming its first position, a unit's own numbering before its molecule's | `core/src/structure/_test/facts.test.ts` "negative and zero author numbers, insertion codes, a label given twice … (EC-1 to EC-8)", "the unit's own author numbering beats its molecule's … (EC-8)"; `core/src/query/_test/numbering.test.ts` "(EC-2)", "(EC-4)", "(EC-8)" |
| 3 | kept, by `offset` and `segments` numberings (MK, NB's leftovers) | `core/src/resolve/_test/maps-segments.test.ts` "an offset numbers a position plus `by`, answering a label without a scan (EC-3)" |
| 11 | kept on `segments` maps | `core/src/resolve/_test/maps-segments.test.ts` "onto a cyclic target a run wraps, split at the origin (EC-11; …)" |
| 13 | kept: a unit never wraps | `core/src/lang/_test/molecules.test.ts` "a unit never wraps (EC-13) …" |
| 20, 21 | kept on `segments` and `structure` maps (runs overlapping on `to` refused for `structure`) | `core/src/resolve/_test/maps-segments.test.ts` "(EC-20; …)", "(EC-21; …)"; `core/src/lang/_test/maps.test.ts` "a structure map's runs overlapping on `to` are refused … (EC-20)" |
| 22–27 | kept on every encoding with a ranked pair: `priority` never lengthens a path, a tie stays ambiguous with a hint naming `priority` and `via`, a ranked map down to each copy, a ranged unit clipped, a composed path as sparse as its maps | `core/src/resolve/_test/maps-choice.test.ts` "two maps of one pair: priority (EC-22 to EC-27, every encoding)" and its cases; "paths: properties at small seeds (EC-22, EC-24)" |
| 35 | kept: a molecule its residues and atoms systems, a unit its molecule's positions from 0 | `core/src/resolve/_test/molecules.test.ts` "… from 0 (EC-35)"; `core/src/resolve/_test/overlay-molecules.test.ts` "(EC-35)"; `core/src/query/_test/elements.test.ts` "(EC-35)" |
| 36 | kept, not named: a node's identity by key across a reattach | `core/src/data/_test/reattach.test.ts` "keeps a node's identity by key and takes its new content …" |
| 48 | kept: pairs, ranges and maps on molecules, units and an alignment's columns | `core/src/lang/_test/molecules.test.ts` "take molecules and units as targets, and a complex as a section (EC-48)"; `core/src/render/_test/pairs-copy-ends.test.ts` "a pair's ends on copies (EC-48, EC-65)" |
| 51 | kept, not named for the services: the third slice's records test holds the clip; SV's UniProt runtime reads through it | `core/src/data/_test/records.test.ts` (the third slice's row) |
| 62 | kept: a glycan's links found from either end, its Asn by link | `core/src/query/_test/elements.test.ts` "linked: the Asn and what the links join to it … (EC-62)"; `core/src/query/_test/elements-pdb.test.ts` "3SGJ: … (EC-62)" |
| 64 | kept: one position per heavy atom; a hydrogen left out and said, a host's and a spec's alike | `core/src/data/_test/seq-molecules.test.ts` "a host's molecule with a hydrogen is said … (EC-64)"; `core/src/io/_test/structure.test.ts` "leaves a molecule's hydrogens out, said (EC-64)"; `core/src/query/_test/elements-pdb.test.ts` "4HHB: four irons, a heme each … (EC-64, EC-65)" |
| 65 | kept: copies by unit, merged and split, from a file's entities and chains | `core/src/render/_test/pairs-copy-ends.test.ts`, `core/src/view/_test/complex-blocks.test.ts` "… copies merged (EC-65)", `core/src/runtimes/_test/mmcif-pdb.test.ts` "(EC-118, EC-65)" |
| 66 | changed: a pair's end on a split copy that no unit names is drawn on no copy, never a guessed one, said once (`seq.render.pairs.copy`; the author, CU); streamed data on a merged block is not tested | `core/src/render/_test/pairs-copy-ends.test.ts` "(1) contacts on the molecules naming no units, split: on no copy, never a guessed one, said once" (not named) |
| 67 | kept with molecules: 200 merged copies and a glycan chip at frame time | `core/src/render/_test/molecules-time.test.ts` "200 merged copies and a glycan chip: a pan and a zoom frame (EC-67)" (the bound heavy only); `core/src/resolve/_test/molecules.test.ts` "200 merged copies (EC-67)" |
| 89 | not met: a stale query text beside its expression is not reported (no query text is stored beside the expression) | none |
| 91 | kept: archives bounded by entries, size, ratio and names, read and written alike | `core/src/data/_test/archive.test.ts` "EC-91: refuses to write past the bounds it reads with", "writing an archive: what is written opens (EC-91)" and the bounds' cases |
| 95 | kept: a block's label inline on its first text line, the start number clear of it (LB's leftover) | `core/src/render/_test/molecules-draw.test.ts` "in text, a block's label inline on its first line … (EC-95)", "… clear of its inline label (EC-95)" |
| 99 | kept: a complex's chips flow like words in a text section (the third slice's not met) | `core/src/render/_test/molecules-draw.test.ts` "in text, the molecules' chips flow like words on one line, each labeled (EC-99)"; `core/src/render/_test/ec-blocks.test.ts` "EC-99: ligand chips in a text section" |
| 101 | kept: a block a unit, polymers first, then branched, ligand, ion; water left out | `core/src/view/_test/complex-blocks.test.ts` "a block a unit: … water left out (EC-101) …" |
| 110 | kept, not named: a custom picture drawn, placed and hit as a built-in's; a host-drawn navigator over `NavigatorModel` | `core/src/render/_test/navigator-picture.test.ts` "is drawn as its tracks, its window placed and hit as any picture's" |
| 118 | kept: observed runs, author numbering and insertion codes from `atom_site`, the first model; colors by run | `core/src/runtimes/_test/mmcif.test.ts` "a unit's observed runs and auth numbering from atom_site … (EC-118, …)"; `core/src/structure/_test/colors.test.ts` "(EC-65, EC-118)" |
| 132 | kept, not named: "+n" past `maxItems` in the compact form, merged legends | `viewer/src/_test/legend-options.test.ts` "past maxItems the compact group says +n …"; `core/src/theme/_test/theme.test.ts` "merges legends that say the same …" |
| 141 | kept for the structure link: a late hover from the other side never takes the pointer's, a late leave clears nothing | `storybook/src/_stories/features/structure-link/_test/structure-link-stories.test.ts` "latency: … (EC-141)", "a late hover of the fake, older than the pointer's own, is dropped … (EC-141)" |
| 146 | kept (refused by default, as designed): a header's hover lights none of the track's items | `core/src/render/_test/headers.test.ts` "EC-146: a header's hover lights none of the track's items …" |
| 155 | kept for header parts: roving focus, Esc back to the canvas, Tab never trapped | `viewer/src/_test/header-overlay.test.ts` "a part at each box … (EC-155)" |
| 156 | kept: the menu key once per ask | `viewer/src/_test/headers.test.ts` "the menu key once per ask … (EC-156)" |
| 166 | kept with a gap, not named: the card's items capped (`hover.cardItems`) and its track named; bridging the gap into the card is untested | `viewer/src/_test/hover-card.test.ts` |
| 167 | kept: the size gate lights whole units past `maxUnits`, a run on each below | `core/src/structure/_test/units.test.ts` "the size gate (EC-167)"; `core/src/structure/_test/hover-time.test.ts` "a selection on 4,096 copies … (EC-167)" |
| 169 | kept on the slice's new reads: a unit past 0 on a source's sequence, a unit-level property, a map onto a source's unit pending | `core/src/render/_test/unit-letters.test.ts` "a unit past 0 on a source's sequence (EC-169)"; `core/src/query/_test/extensions.test.ts` "a unit-level property pending (EC-169)" |
| 174 | kept for services: a failing service leaves the data drawn and says which; an unrecorded location said, never retried | `storybook/src/_stories/support/_test/services-app.test.ts` "a service that fails leaves the data drawn and says which … (EC-174 …)", "… NotRecorded, said so, never retried (EC-174) …" |
| 175 | kept, not named: a host's data reopened in place enters no patch and no history (seq data ports' rule) | `core/src/data/_test/adapters.test.ts` "a chain removed: the source reattached …" |
| 176 | kept: a node's identity by key across an update, its content new | `core/src/data/_test/reattach.test.ts` "(EC-176)"; `core/src/data/_test/adapters.test.ts` "a chain removed: … (EC-176)" |
| 181 | kept: a host's edit renames nothing | `core/src/io/_test/structure.test.ts` "keeps every other key and node when a chain is removed (… EC-181)" |
| 183 | kept: a second consumer draws what core draws from the resolved view and the store alone | `examples/src/_test/consumer/consumer.test.ts` "EC-183 the second consumer draws what core draws" and its layers |
| 185 | kept: a source's `ref` reaches the host's adapter, alone or beside a uri | `core/src/data/_test/adapters.test.ts` "… a `ref` alone reaches it (EC-185)"; `core/src/lang/_test/molecules.test.ts` "(EC-185)" |
| 186 | not met: bigWig is self-indexed; no format with a separate index was added this slice | none |
| 200 | kept for a dropped file: a name with a URL's characters is one location | `storybook/src/_stories/support/_test/dropped-file.test.ts` "a name with characters of a URL is one location (EC-200 …)" |
| 201 | not met for bigWig: its runtime reads by region and zoom level, no open-domain view over it | none |
| 207 | kept: several maps between one pair without every path ambiguous (by `priority` and `via`) | `core/src/resolve/_test/maps-choice.test.ts` (EC-22 to EC-27's cases); not named |
| 229 | kept: a query over data read by region stays pending while its tile has not arrived | `core/src/query/_test/extensions.test.ts` "pending while the tile has not arrived … (EC-229 …)" |
| 235 | changed: an extension's unknown `<id>:<name>` type is drawn as its declared fallback, said once per extension (`seq.lang.layer.fallback`, `seq.view.layer.extension`); a newer bare built-in type likewise, a warning, never the enum's error; not named by its number | `core/src/context/_test/extension.test.ts` "is drawn by a composition that has it, and as its fallback by one without it", "draws a newer bare type as its fallback, a warning and never the enum's error"; `storybook/…/extensibility-stories.test.ts` "draws the fallback without the extension …" |
| 237 | kept for track headers: moving and showing tracks through the menu | `viewer/src/_test/headers.test.ts` "seq:move-up and seq:move-down … (EC-237)", "seq:hidden-tracks: … (EC-237)" |
| 238–244 | kept for the test consumer's subset (sections on alignment columns and chains, placement through a path, composited layers); not named by number | `examples/src/_test/consumer/consumer.test.ts` (EC-183's cases) |
| 36, 51, 110, 132, 175, 207, 235, 238–244 | **to name**: the tests above hold them without their number; a test-only lane (as ET) adds `(EC-<n>)` to their titles | — |
| 171 | not met: no zoom gate with "load anyway" over bigWig | none |
| 204 | not met for the new legends: painted vs vocabulary is not said by `LegendSelection` | none |

## Outcomes of the close-out

Recorded at the close-out's reconcile (2026-11-24) for the items of built rows that no test names by number and that had
no outcome yet (`ec-ids.ts --worklist`). Each was matched by behavior against the tests on `packages-seq`; a clause no
test asserts is said in its row. Counts of the 56: kept 39 (11 of them in part: the clause named), changed 10, not met 7
(EC-34, 52, 98, 104, 105, 158, and 168 in part). No test was written to cite a number. Paths relative to
`packages/seq/`; the generic repr tests are in `packages/core/src/repr`.

| EC | Outcome | Test |
|---|---|---|
| 31 | kept: a hidden stretch and a reversed block are layout (axis breaks with hidden ranges, a direction per block), two pieces of one unit are blocks over one system | `core/src/view/_test/axis.test.ts` "has its breaks with their hidden ranges"; `core/src/view/_test/genome-blocks.test.ts` "keeps a reversed copy apart from a forward one" |
| 34 | not met: a block past 2^31−1 positions left out with a "use a range" message; only the numbering limit and the refusal of a contig past int32 at open exist | related: `core/src/runtimes/_test/vcf.test.ts` "a contig longer than 2^31 − 1 is refused by name" |
| 42 | kept by the design, not tested for a row bound twice: paths search bound rows only; the generic ambiguity rule is tested | `core/src/resolve/_test/paths.test.ts` "two maps joining the same systems differently are ambiguous (EC-24)" |
| 45 | kept in part: rows sharing a key are one item; "where rows overlap, the last covering row wins" has no test | `core/src/render/_test/source-wrap.test.ts` "is one item: both pieces answer to one key" |
| 46 | kept: a zero-length item meets the position it is at and draws as a hairline when selected | `core/src/data/_test/intervals.test.ts` "a zero-length interval meets the position it is at"; `core/src/render/_test/emphasis.test.ts` "a selected item of zero length is a hairline at its place" |
| 47 | kept: misfit rows left out, one issue a cause with its count and first places; a table with no rows is a valid empty source | `core/src/data/_test/records.test.ts` "row-level misfits: each left out, one issue a cause" |
| 49 | kept: `fields` maps roles by column or path, text cells keep their text (`01` stays `01`); an explicit swap is not asserted | `core/src/data/_test/records.test.ts` "keeps text cells as written", "maps roles and fields by column or path" |
| 52 | not met: no test for a comparison with a missing (NaN) value being false and its `not` true; NaN is kept as missing in values and aggregates | related: `core/src/lang/_test/lang.test.ts` "keeps NaN in values through JSON" |
| 53 | kept in part: numbers exact through JSON, NaN and −0 preserved; the Float32 shortest-decimal comparison has no test | `core/src/data/_test/seq-format.test.ts` "round-trips through JSON: typed columns, NaN, -0" |
| 54 | kept: a values node holds one array per target; a length that is not one per position is `seq.lang.values.length` | `core/src/lang/_test/lang.test.ts` "is reported, with its severity and a place" |
| 55 | kept by the design, not tested: values draw at their value, the domain bounds the height (`themes.md`) | related: `core/src/render/_test/values-reference.test.ts` "a closed domain’s outside line is not drawn" |
| 56 | changed: shared domains are a named `scale` node; an open domain is the union over the layers naming it at settle, a closed end does not move, a diverging scale is set by `mid` | `core/src/render/_test/themes.test.ts` "unions an open domain over the layers naming it at settle"; `core/src/theme/_test/scales.test.ts` "makes a numeric field’s default palette diverge" |
| 58 | kept: a grid is x-major, a wrong-size grid is dropped saying the count expected, typed columns kept by identity | `core/src/view/_test/grid.test.ts` "drops a grid of the wrong size, saying the count it expected" |
| 60 | kept in part: zoomed out a bin draws its majority at its share's alpha; the tie rule and a `maximum` density have no test | `core/src/render/_test/matrix.test.ts` "zoomed out, draws the bins as an image tile of bins × states" |
| 61 | changed: a matrix has one states system per entry, so its items share one stride | `core/src/resolve/_test/region-plans.test.ts` "give a matrix a states system, pending until its source's coverage lists them" |
| 68 | changed: a reference is an address, not the old ref syntax; an id is required on every child but the root | `core/src/lang/_test/lang.test.ts` "require an id on every child of the spec" |
| 69 | kept: the kind a ref param takes is checked statically; a ref to no node is reported on resolve | `core/src/lang/_test/axis.test.ts` "check the kinds their references name"; `core/src/resolve/_test/model.test.ts` "a reference to no node is reported on an unchecked spec" |
| 70 | kept: a ref cycle is a warning on each ref in it, the rest of the check runs | `core/src/lang/_test/refs.test.ts` "reports a ref naming nothing and a chain that comes back, as warnings" |
| 71 | kept: dependencies are the typed ref params only | `core/src/lang/_test/axis.test.ts` "are followed by a change: a node `along` or `pieces` names, removed, is reported for the section" |
| 72 | kept: a source no runtime reads is reported and stays unavailable; an extension a composition lacks is drawn as its fallback | `core/src/resolve/_test/model.test.ts` "reports a source no runtime reads, and keeps it unavailable"; `core/src/context/_test/extension.test.ts` "is drawn by a composition that has it, and as its fallback by one without it" |
| 73 | kept: a newer minor's unknown kind is a warning drawn by fallback, a newer major is refused alone | `core/src/lang/_test/lang.test.ts` "a newer minor’s unknown kind is a warning", "a newer major version is reported alone" |
| 79 | changed: no `first` or `last` anchors; a place is `{ after }` or `{ last }`, `set` and `unset` act on one key | `core/src/lang/_test/patches.test.ts` "hold the core's cases and Seq's own" |
| 83 | kept by the spec, not tested for a duplicated id: `diff` falls back to a `replace` of the parent | `core/src/repr/patch/_test/patch.test.ts` "falls back to a replace" |
| 88 | kept by the code, not tested: an unknown key and a bad value carry hints listing the keys and values; tests assert the code and path | `core/src/lang/_test/lang.test.ts` "give the description's codes first" (code only) |
| 98 | not met: `align: center` and `max_width` are implemented, no test asserts that hits and margins move with them | none |
| 104 | not met: split rows are deferred with Ranges' split (`track-groups.md`) | none |
| 105 | not met: deferred together with EC-104 | none |
| 120 | changed: no glyph atlas; letters are drawn as glyph runs through `fillText`, fonts and glyph strings cached | `viewer/src/canvas/_test/backend.test.ts` "keeps colors and glyph strings: drawing the same lists twice makes the same calls" |
| 121 | changed: no atlas pages to upload; the font is set only when it changes, nothing re-measured by zooming | `viewer/src/canvas/_test/backend.test.ts` "measures in the CSS font the backend draws with, setting the font only when it changes" |
| 122 | kept in part: a letter is centered in its whole cell, edge and wider-than-view cells included; per role, size and DPR is not asserted | `core/src/render/_test/letters-whole-cell.test.ts` "a cell wider than the view centres its letter off-screen" |
| 123 | kept in part: the hit is the layout's own, a pixel hits the cell drawn on it; a glyph's overhang into a gap is not asserted | `viewer/src/_test/viewer.test.ts` "a hit under the pointer is the layout's: a letters cell with its position" |
| 124 | kept: the SVG and PNG share the draw list and a bundled metrics table, the PNG writes valid CRCs | `core/src/render/_test/png.test.ts` "writes the signature, IHDR, one IDAT of stored rows and IEND, each with its CRC"; `core/src/render/_test/backends.test.ts` "measures with the bundled table" |
| 133 | changed: the palette order is the section's domain, known up front, never the first occurrence by region | `core/src/theme/_test/theme.test.ts` "compiles over one domain known up front share its colors" |
| 134 | kept: an unlisted value takes the list color at its FNV-1a hash over UTF-8; a dictionary scheme gives an unlisted letter the missing color | `core/src/theme/_test/theme.test.ts` "an unlisted value takes the color at its FNV-1a hash over UTF-8" |
| 136 | kept in part: dark variants derived, only fitting schemes offered, the reader's substitution writes nothing into the spec; "a pick survives patches" is not tested | `core/src/theme/_test/theme.test.ts` "a scheme's literal colors take a derived dark variant"; `viewer/src/_test/legends.test.ts` "offers the schemes that fit a layer" |
| 139 | kept in part: ink by the contrast rule, chip ink at 4.5:1 in both; the 0.179 threshold is not asserted | `core/src/theme/_test/ink.test.ts` "equals the contrast rule, first and memoized, for any pair of inks" |
| 140 | changed: no band copy; a pan redraws from the lists and equals a fresh renderer's; hide and show equal to a fresh mount is not tested | `core/src/render/_test/pairs-pan.test.ts` "(3) a pan's frame equals a fresh renderer's" |
| 142 | kept: a hover is an element in its own system, projected onto another section when asked, once per version | `viewer/src/_test/hover.test.ts` "is computed when asked, once per version; configurable" |
| 147 | kept in part: the row band is drawn in the hover pass from the emphasis input; that it changes no state is not asserted | `core/src/render/_test/matrix-emphasis.test.ts` "the hover's tint is the config's" |
| 150 | kept: a drag held past the box edge autoscrolls by the clock and stops at the release; growth with distance is not asserted | `storybook/src/_stories/features/text/_test/text-viewer.test.ts` "a drag-select held past the box's bottom scrolls it, by the clock, and extends" |
| 151 | kept in part: a text drag selects in reading order across lines (a known flaw fixed: cross-block selection is no longer deferred); "never pans" is not asserted | `storybook/src/_stories/features/text/_test/text-viewer.test.ts` "a drag selects in reading order across lines" |
| 152 | kept in part: a press always ends (release elsewhere, cancel, blur, lost capture); one pointer at a time and clamping are not asserted | `viewer/src/_test/drag-release.test.ts` "a second press after a lost release starts a new drag" |
| 153 | kept: input is a binding table by gesture, modifiers and hit target; an unbound drag does nothing | `viewer/src/_test/bindings.test.ts` "the default preset resolves each gesture by the hit’s target" |
| 158 | not met: no dwell or pointer guard on reveal-on-highlight; only the explicit reveal exists (panning only) | related: the search's "next reveals, panning only" |
| 168 | not met in part: `show` switches views and gives the other view's sections viewports; "prefer a view that holds the extents directly" is not tested | `viewer/src/_test/follow.test.ts` "show switches views" |
| 178 | kept: a changed overlay version changes the frame key only; the model and resolved view are reused | `core/src/render/_test/dynamic.test.ts` "the overlay’s version is part of the key; the model and the resolved view stay" |
| 180 | changed: a codon is cells on the coding axis, one cell per part of a split codon, no `\|` label; a minus strand reads from the highest exon down | `core/src/render/_test/axis.test.ts` "a residue a codon wide on the gene’s coding axis"; `core/src/resolve/_test/encodings.test.ts` "minus strand: the coding sequence starts at the highest exon and reads downward" |
| 192 | kept: `limits.maxBytes` refuses a file before it is read, `seq.source.too_large` | `core/src/runtimes/_test/genbank-bounds.test.ts` "its bytes (`limits.maxBytes`), before a line is read" |
| 193 | kept: lane heights, fonts and row pixels are the render config's; a track asks for pixels by its own `height` | `core/src/render/_test/track-heights.test.ts` "a track’s own height overrides" |
| 194 | kept: named schemes and token palettes have values in both light and dark tables | `core/src/theme/_test/conformance.test.ts` "colors in each table, ink at 4.5:1 on each, a provenance" |
| 197 | kept by the format design, not tested as a stated 0-based code: a category column is a dictionary and int32 codes | related: `core/src/theme/_test/theme.test.ts` "a dictionary-coded value’s color is its value’s" |
| 198 | kept by the code, not tested: a ranges item's label is its `label`, else its id | none |
| 199 | kept: each type states its fallback; a ranges layer takes pairs as marks, a values layer does not | `core/src/lang/_test/molecules.test.ts` "a ranges layer takes pairs (each end a mark at its site); a values layer does not" |
| 226 | changed: no sorted-input flag; lanes pack first-fit by start over the whole data whatever the file order (lane order itself not asserted) | `core/src/render/_test/ranges-reach.test.ts` "is as tall as the whole region’s lanes at the first settled frame" |
| 234 | kept in part: leading ruler tracks are a fixed head; an absent or lower ruler is not tested | `core/src/render/_test/section-scroll.test.ts` "scrolls the tracks below the ruler by `top`, the ruler fixed" |
| 236 | kept: `diff` writes an unset param back as a `set` and the patch round-trips; snapshots and metadata are not asserted separately | `core/src/lang/_test/molecules.test.ts` "unset on its node, and set again" |

## Notes for the close

- `ec-ids.ts --worklist` still lists many earlier rows' items (most from before the coverage rule, in its baseline);
  this slice moved no row to `done`, so the check stays green; the rows it moved to building with new cases (Archives,
  Fallbacks, Host data, A second consumer) have their items in the table above.
- EC-63, 115, 116, 135, 195 were the third slice's not-met items: tested before this slice (plan, Gate 0, 10); EC-59 and
  EC-209 to the fifth slice; EC-99 kept here (MR).

## In the specs

| Decision | Where |
|---|---|
| The catalog, its numbering, where it lives, and `ec-ids.ts` checking every cited id | design-only: the scaffold's catalog; the ids a spec cites are checked by the scaffold |
| Three outcomes (kept, changed, dropped), decided when an area is built | design-only: a process; each outcome is a test or a spec rule of its area |
| Known flaws fixed, not ported (EC-130, 135, 151, 181) | design-only: the catalog's outcomes; the fixes are in the areas' specs |
| The reply's cases are cases (EC-185 to EC-244) | design-only: the catalog's outcomes; the areas' specs hold the behaviors |
| An area is not done until each item has an outcome | design-only: a process (status and the plan's gates) |
| At the scaffold's close the kept items need a home in a `_spec` | design-only: the close's work (the coverage worklist), not a package rule |
