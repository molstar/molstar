# Lane SX: a section scrolls past its cap, with the one rail

**Base**: `packages-seq` at the orchestrator's head after FX10, FX11 and RL's pick merge (`git reset --hard
packages-seq` first). Branch `seq/sx`. Model: strong (a scroller abstraction over three geometries, the section
layout, the wheel chain, export).

**Why**: the author, 2026-11-22: "not just text sections. imagine having multiple tracks with standard features, these
should be scrollable too. and 600px should definitely be configurable"; "we should use the same rail everywhere".
layouts.md: "A slim rail" (one rail for every box that scrolls) and "A section scrolls past its cap, whatever its
layout". FX10 made the lanes rail the rows rail; it stopped at the rest: a text section scrolls by a range anchor
(`viewport.start`), not by px or a band top, so the rail's gesture needs a scroller abstraction, and a tracks section
has no scroll at all.

**Task**:

1. **A scroller**: the one rail's gesture (`rowsRailGesture` and its hit) drives a small interface (its extent, its
   window, set, page, cancel), implemented by a rows band, a lanes band, a text section (lines, from its geometry) and a
   tracks section (px, the section viewport's `top`). One drawing, one hit, one gesture; grips only where the box
   zooms.
2. **A tracks section past its cap**: its box the section's `height` (spec param: check its reader covers tracks; widen
   doc and reader, no language bump; a new param stops the lane), else its content up to `section.maxPx` (render
   config, default 600; `text.maxHeight` folded in, one key, an alias if anything outside seq reads it); the ruler and
   navigators fixed above; the wheel chain (a lanes or rows band first, then the section, then the page); an export
   showing every track.
3. **Keys and the floor**: Page Up/Down and Home/End with the cursor in any scroller (resolve the binding conflict
   FX10 met by one command over the scroller under the cursor); the window's floor on keys and wheel as on drags.

**Tests**: one table-driven test across the four scrollers answering the same gesture the same way; a section of 30
gene tracks at 600, at a host's 300 and at a spec `height` of 200, scrolling by wheel, rail and keys, the ruler fixed;
the wheel chain; an export with every track; a warm scrolled frame 0 B (both perf tiers). Each fails on its parent.

**Where**: a worktree on `seq/sx`. Headers (`@author David Sehnal <david.sehnal@gmail.com>`, 2026); hot files never
staged: list them. Stories: a `Features/Layouts` story of a tracks section past its cap (light and dark, canvas and
SVG, headless twin). A file past 800 lines is a question. Sizes (headless 930, default 1,470, minimal 1,060 kB): a
bound passed stops the lane with a report.

**Checks**: as `fx10-lanes-rail.md`. **Report**: under 350 words.
