# The fourth slice's edge-case outcomes, drafted (finalize, 2026-11-21)

Drafted at finalize from the tests on `packages-seq` (`5bf619ac6`, HT merged), for the close to paste into
`.scaffolds/seq/design/edge-cases.md` under "## Outcomes of the fourth slice". The items are those the slice's features cover
(features.md's Cases of the rows it moved, the plan's named items): molecules and units, structure maps and the link,
seq data ports, archives, the test consumer, extension openings, track headers, services. Tests found by the name they
cite (`EC-<n>` in the title or a comment of a `_test/` file); where no test names the item, the test that holds the
behavior is cited and the row says "not named". Read, not run here (the gate runs them). Paths relative to
`packages/seq/`. Items the earlier slices recorded and this slice left as they were are not repeated (EC-15–18, 28–30,
63, 170–173).

Counts, by row: **kept 39** (8 not named by a test; 4 with a stated gap: EC-51, 67, 166, 238–244), **changed 2**
(EC-66, 235), **not met 5** (EC-89, 171, 186, 201, 204).

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

## Notes for the close

- `ec-ids.ts --worklist` still lists many earlier rows' items (most from before the coverage rule, in its baseline);
  this slice moved no row to `done`, so the check stays green; the rows it moved to building with new cases (Archives,
  Fallbacks, Host data, A second consumer) have their items in the table above.
- EC-63, 115, 116, 135, 195 were the third slice's not-met items: tested before this slice (plan, Gate 0, 10); EC-59 and
  EC-209 to the fifth slice; EC-99 kept here (MR).
