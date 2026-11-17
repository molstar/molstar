# Drift: query decisions the orchestrator made from lane C's review (to confirm with the author)

- `(chain :id A)` returns Nodes (a data-made system addressed by segments), so the structure path
  `/data/1tup/(chain :id A)/(at …)` works; locations-and-queries.md listed `chain` under Elements. Edit the design.
- `.n` labels are written quoted (`(at :label "100.1")`): the address atom rule reads `100.1` as a number.
- A boundary at p overlaps [s, e) when s ≤ p < e; a boundary at the system's end overlaps a range ending there (EC-46).
- Location equality is structural (rects are not merged); documented, not normalized.
