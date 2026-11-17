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
  - **changed**: the new design handles it differently, and says how (proposed so far, to confirm: EC-133, palette order by first occurrence, is
    changed by fetch-by-region; EC-55's doc conflict is settled as "drawn at the value, the domain bounds the height");
  - **dropped**: only by the author's decision, recorded in the design of its area with the reason.
- **Known flaws are fixed, not ported** (recommended; the author confirms): EC-130 (minus-strand secondary structure), EC-135 (ESPript and gaps), EC-151
  (cross-block text selection, deferred before), EC-181 (repeated ids across structures).
- **The reply's cases are cases, not wishes**: what it asked for that the designs take (EC-185 to EC-206, EC-235,
  EC-236) is built and tested like any other item; what is deferred (EC-209, EC-213, EC-223 and the smaller layer
  features) stays in the catalog with its trigger; EC-238 to EC-244 are tested as fallbacks a reader without them
  must get predictably.
- **An area is not done** until each of its items has an outcome. The area's plan lists them; its record says the
  outcome of each.
- **Before `.seq/` is deleted**, the kept items live on as the packages' tests and their `_spec` statements; the catalog
  goes with the scaffold.

## Consequences

- The simplification shows up as fewer mechanisms carrying the same cases: one segment-table map type for every map
  kind (EC-9 to EC-32), one transient overlay for previews instead of four document layers (EC-163), one source protocol
  for streamed and inline data (EC-169 to EC-180).
- The examples are chosen to exercise the catalog: few datasets, many views ([examples.md](examples.md)).
