# Lane DL: two short designs the fix lanes wait for

**Model**: strong. **No code, no worktree, no commits**: write two drafts in the main checkout's
`.seq/execution/fifth-slice/drafts/` (the stage 0 drafts' rule: one file each, the orchestrator commits them). Read
first: `.seq/design/README.md` (how a design is written: Intent, Precedents, Decision, Defaults with a picture, the
interaction matrix's cells, Configuration, Tensions, Self-check, Consequences), `ux-foundations.md`, `navigators.md`,
`layouts.md`, `grids-and-rows.md` (the rows' scroll and rail), and the two items in `fixes.md` (Open) named below with
the author's words. Look at the code that draws them today (`packages/seq/core/src/render/navigator*.ts`,
`render/ranges*.ts`, the viewer's navigator mount) to say what changes; quote file:line.

1. **`drafts/navigator-window-look.md`**, an amendment to navigators.md, from "The navigator's window and zoom strip
   look inconsistent and rough": one rule a reader can see for when the window is a box with grips and when the zoom
   strip appears (or why always one), one stroke weight for the window, the funnel and the strip, their fill and
   border from the same tokens (the loudest mark at rest), the funnel thin and muted or dropped; ASCII pictures of the
   box, the strip, a window at the fit (F25: a quiet frame) and a window over the cyclic origin; the cells of the
   interaction matrix it changes; Configuration rows (`navigator.zoomStrip.minPx` and anything new) with whose and why.
2. **`drafts/track-sizing.md`**, an amendment to layouts.md / ux-foundations' heights, from "Sizing a track of many
   lanes: content height, a cap, then its own scroll": a ranges track sized to its lanes as known (growing as a
   source's rows arrive, never a blank band), up to a cap (`track.maxPx` or a share of the view: pick one and say why),
   past the cap its lanes scrolling inside the track (wheel, keys, a slim rail like the rows', the header and labels
   fixed), the first row labelled; how this meets the empty-well rule (`fixes.md`), export (an SVG of a capped track:
   all lanes or the cap?), the rows layout, and a ring (no inner scroll on a ring: say what happens); precedents
   (IGV's expanded/squished, JBrowse's max height, UCSC's dense/pack); edge cases by EC id where the catalog has them
   (`.seq/design/edge-cases.md` or wherever the catalog lives; new ones proposed as `EC-new-n`).

Each draft ends with "Gaps for the author" (yes/no questions, each with a recommendation) and its self-check. Under
about 250 lines each. Report under 200 words: the two files, the decisions each proposes, the gaps.
