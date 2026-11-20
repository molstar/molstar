# Fixes: the author's findings and the looks' notes, collected

Findings from use that need no new design collect here, and a fix lane (FX1, FX2, …, fast unless the cause is in a
model) takes the open ones at each fix point in [plan.md](plan.md): before stage 1, after each stage's gate, and before
M4. A finding that needs a design goes to a stage 0 style draft instead; one that blocks the author's use jumps the
queue. Each line: what is seen, where (the story), where it came from, its kind (bug, look, missing design).

## Open

- **Ring letters off their cells** (look): on a ring the letters drift from the cells they name, unevenly spaced and
  crowding in places (UUUAG near 10, GGGAG near 20), not centred on their cell's angle and the track's mid radius. Where:
  tRNA › Views › Chords (and every ring with letters). The author, 2026-11-20. Check the placement against the cell's
  centre angle and the glyph's centre (not its baseline start), in the canvas and the SVG.
- **The Reverse track reserves an empty label lane** (look): Plasmid › Ring › Detail with a strip. Third slice's DS look.
- **Dense leaders when a ring group opens** (look): Plasmid, a group of more than 6 labels hovered. Third slice's DS look.
- **The Antibody's IMGT label reads as a second line** (look): Antibody, the ruler's numbering label. Third slice's look.
- **A moved label of an item in a lower lane is dropped** (bug): leaders exist for lane 0 only (LX). Third slice's LX
  review.

## Fixed

(none yet; each line moves here with its lane and commit)
