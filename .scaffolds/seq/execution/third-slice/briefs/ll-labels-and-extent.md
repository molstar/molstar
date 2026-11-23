# Lane LL: item labels inside or beside, and an origin-crossing extent (fast)

The author (2026-11-19, the plasmid's ring and linear views): the labels read differently in the two views. Decided:
"Inside if it fits, else beside" (`.seq/design/layouts.md`, "Labels inside if they fit, else beside"; the Design page
`layouts.mdx`, Rings). Read that section and `ctx.labels` under Layer types first.

## Do

1. One label rule in both layouts. A labeled item's text inside its mark when it fits (the mark's length less a
   padding, at the track's text size); else beside: a ring places it outside with a leader (as now); a linear track
   places it in a label lane just above the track at the item's start, a short leader when it moved. The same sink
   (`ctx.labels`) and placement time (when laid out, never on a window move or hover), the same drop priority (the
   layer's order, then the longer item), the same text style inside and beside. A placed label's hit is its item's.
   Linear placement is per strip and must stay inside the pan budget (no relayout on a pan in the same level: place
   once a level or a layout, as the ring does; say what you chose).
2. The extent bug: the hover card of bla (AmpR) on pUC19c says "Extent pUC19c 1–2,686", which is wrong: bla crosses
   the origin. Find where (`viewer/src/describe.ts` item card, `core/src/render/accessibility.ts` `spanText`/`itemText`,
   or the item read that gives it start and end) and say an origin-crossing item as its two runs or as start–end
   through the origin (e.g. "2,486–2,686, 1–175"; pick what EC-9 and the selection's words say and match them).
   Test: the card and the live region's words for an origin-crossing item.

## Tests

- Linear: a long feature's label inside its arrow; a short one's in the lane above, with a leader when moved; two
  overlapping short ones spread or one dropped by priority; the hit of each label is its item.
- Ring: unchanged behavior for labels that do not fit; a label that fits its arc drawn inside (if the arc is long
  enough at the ring's radius).
- The plasmid story's test: the same items labeled in the ring and the linear detail with the same text style.
- Each fix's test fails on its parent.

## Look

Features › Plasmid (ring and linear detail) in light and dark at dpr 2: labels readable, no label over another or over
a mark it is not inside.

## Rules

The lane rules of `.seq/execution/third-slice/loop.md` (hot files never staged, headers, changelog file to report, not
add). Report: checks with exit codes and counts, sizes, the pan perf (the plasmid's linear detail), what you chose for
linear placement time.
