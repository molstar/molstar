# Lane NP: the novice pass's bugs, stories and looks

**Goal**: the findings [novice-findings.md](../novice-findings.md) places on NP (the numbers are its rows). Reproduce
each first in Storybook (your own server on port 6017 from your worktree); a finding that does not reproduce is
reported as such, not "fixed".

**Bugs**
- 12, 13: a story's SVG twin draws a step behind the canvas after an edit or an undo (Complex › Views › Ligands' hide
  and undo; Globins › Rows › Family's fold). Find why the twin renders from the previous state and fix it where the
  twins are built (the stories' shared code), with a test.
- 16: Enter in Features › Query's query box applies it, as the Select button does.
- 8: Globins › Rows › Lazy: after a click into the rows, Page Down moves a window of rows (grids-and-rows, The keyboard
  cursor). If the cursor needs focus the click does not give, fix the focus; a test of the key after a pointer click.
- 19: Features › Highlights › Overlapping channels: the Host toggle takes one click.
- 15: Antibody's Kabat ruler draws "88 82B" overlapping at the insertions: the ruler's tick labels never overlap
  (thin them by the label's width, insertion labels included).

**Stories** (missing design in the page, not the core)
- 4: p53 › Gene by region opens on a region, so its next and previous move; a disabled button says why (its title).
- 7: Globins › Hemoglobin gains the split and merge toggle the Editing story has (or links it).
- 14: Antibody › Design numberings gains the rulers' switch Schemes has.
- 18: Features › Navigators › Strip opens zoomed in, so its window can be dragged.

**Looks**
- 5: Globins › Family's logo tall enough to read at the fit (a track height in the spec).
- 6: the tree's clade hit target at least 8 px (the junction's), in the core's tree layer, with a test.
- 9: the plasmid ring names its strands (the two lanes' labels or the legend's words), and the selection tint does not
  cover the ruler's numbers in the hole (drawn under them).
- 11: the tRNA dot plot's light theme: the low covariation as quiet as in dark (the ramp's low end against the
  surface), so the pairs stand out.
- 19: the host highlight's outline at 3:1 against the light surface (the channel's color or its token).

**Not in scope**: the hover card and keys hints (the author's call), the site labels' fan (9) and long legends (17):
add both to `packages/seq/TODO.md`.

**Model**: fast (sonnet). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B seq/np-novice-fixes
<that commit>`. Nothing else runs in `packages/seq` beside you.

**Checks**: `pnpm test:unit packages/seq/` (and each folder's Duration: core, viewer and storybook sit near the 2 s
guide; add no slow quick test); `pnpm test:medium packages/seq/`; `npx tsc -b packages/seq/core packages/seq/viewer
packages/seq/examples`; `pnpm -C packages/seq/<core|viewer|examples> typecheck:test`; `pnpm -C packages/seq/storybook
typecheck`; `pnpm -C packages/seq/core test:perf`; `pnpm lint`; `pnpm format:check`. Look at each fix in light and dark.

**Headers**, **hot files** and **stop when**: as every lane's (see [qt-quick-lane.md](qt-quick-lane.md)).

**Report**: commits; per finding: reproduced or not, the fix, its test; checks; folder durations; what you saw.
