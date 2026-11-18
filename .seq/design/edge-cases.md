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
  `packages/seq` (a test, a comment, a `_spec`) resolves against it, checked by `.seq/_checks/ec-ids.ts` at the gate.
  Why: the ids are numbers only, so code cites them without a `.seq/` path, and an id looked up in the wrong file was
  once recorded as missing when it was not.
- **Before `.seq/` is deleted**, the kept items live on as the packages' tests and their `_spec` statements; the catalog
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
| 111–114 | a viewport keeps its center and zoom across updates and widths; scroll groups share the range on equal axes; followers never lead | `viewer/_test/{view,follow,navigate}` |
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
- Projected hover (linked-section outlines, EC-141's family): the hover model projects onto linked sections, the draw lists do not yet draw it.
