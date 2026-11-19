# Defaults audit (lane DA)

Base c55db3098 (branch `seq/da-defaults-audit`). Storybook on port 6028, headless Chrome 1200x1000 at dpr 2, light and
dark, every story under `Seq* (next)` Examples and Features: 98 stories, 196 shots, none failed to render. Shots:
`/private/tmp/claude-501/-Users-davinci001--claude-spares-spare-BGbiLE/3fb68e09-fc5b-4c2f-a070-ae4f27a85f49/scratchpad/da-shots/<story-id>.<light|dark>.png`
(ids from `/index.json`; below, the `seq-next-` prefix is dropped).

Coverage honesty: I looked at nine shots in full (examples-antibody-views--schemes, examples-p53-tracks--tracks,
examples-plasmid-ring--ring, examples-complex-views--confidence (dark), examples-globins-rows--family,
examples-trna-views--arcs, features-layers--variants (dark), examples-p53-gene-by-region, features-text--text) and read
the code defaults. "Also" lists stories that share the default by construction, not confirmed by eye. Findings are
ranked by first-impression gain. Config: `packages/seq/core/src/render/config.ts` (heights 262-274, trackGap 261,
sectionGap 262, laneHeight 278, labelFont 284, itemFont 285, text.bandGap 338).

## 1. A ranges track is a fixed 36 px band whatever it holds (decided: layouts.md "fits its content, in steps")

- Stories: examples-antibody-views--schemes (4 CDR bands), examples-p53-tracks (Domains), examples-p53-gene-by-region
  (eight bands, "Transcript W1" entirely empty), examples-trna-views--arcs ("Modified"), examples-plasmid-ring (Forward
  and Reverse strand), examples-complex-views--confidence ("In 1FIN"). Also: the other antibody views, features-layers--*.
- What: one lane of marks sits at the top of a tall grey well, two thirds of it empty. Worst in antibody and
  gene-by-region. In "In 1FIN" the bar hugs the top edge.
- Cause: `heights.ranges: 36` (config.ts:263), `variants: 32`; type table lookup `tracks.ts:1699`; ranges place marks
  from the lane grid at the band top (`ranges.ts:667-673`, `lay.mlane[j] * laneHeight`), not centred in the band.
- Proposal (numbers for the decided rule): lane 14 px (mark 12, gap 2; today `laneHeight` 12 gives a 10 px mark), padding
  4 px top and bottom, step one lane, minimum one lane plus padding (22 px), cap 5 lanes (78 px); so 1 lane 22, 2 lanes
  36, 3 lanes 50, 4 lanes 64, 5 lanes 78 (today's 36 is the 2-lane case). Marks centred in the band. Same-kind tracks of a
  section share the tallest (two CDR tracks, all transcript tracks equal). An empty track takes the minimum, not zero (rows
  stay aligned; hiding all-empty tracks is a spec choice). Values (32), letters (20), ruler (22), pairs (64) stay fixed:
  they have no lanes.
- Why these numbers: the 22 px minimum equals the ruler band, so a one-lane track and a ruler look like siblings; 4 px
  padding clears the squircle corners (`squircleMin` 6); the 5-lane cap keeps a dense track under 80 px.

## 2. Item labels clip mid-word ("Transactivatio", "CDR-L", a blank blue mark) (decided: layouts.md "Labels inside if they fit, else beside")

- Stories: examples-antibody-views--schemes (CDR-L1 "(K", an unlabelled Heavy CDR1 Kabat mark next to its labelled IMGT
  twin), examples-p53-tracks and features-layers--variants (Domains: "Transactivatio", "Proli", "Nuclear", "C-ter",
  "Tetram"). Also: every ranges story with labels.
- What: the label is drawn over the mark and cut at the mark's edge; half a word reads as a bug.
- Cause: `ranges.ts:547-560` ("labels inside marks wide enough, faded past their width") and the `laneHeight >= 10` gate;
  `itemFont` 10px.
- Proposal: the decided rule as the default: whole text inside if it fits (mark length less padding), else beside in a
  label lane above, with a short leader if moved; no fade-clip. Interim until the lane exists: draw whole or not at all
  (hover card names it). The content height (finding 1) counts the label lane (one extra lane when any label goes
  beside), so the two are one change.

## 3. Track labels in the left column are hard-clipped at a fixed 112 px

- Stories: examples-p53-gene-by-region ("Domains T1-domain", "Domains T2-domain", Canvas and SVG), features-layers--variants
  ("On the gene <- reve"). examples-complex-views--confidence shows the column grown to fit "Secondary (UniProt)" in the
  Canvas, so authors set it by hand.
- Cause: `labelWidth: 112`, `labelFont: 12px sans` (config.ts:259, 284); the clip has no ellipsis (the globins species
  column does use one, "Glyc...", so the idiom exists).
- Proposal: ellipsis on overflow; default `labelWidth` measured as the widest track label plus `labelGap`, floor 80, cap
  160 px. A fixed 112 suits "Sequence" and clips real names. No wrap: it breaks one-row-height.

## 4. Sections run into each other

- Stories: examples-antibody-views--schemes (Light chain's last band to Heavy's ruler is 8 px, barely a break),
  examples-complex-views--confidence (second chain starts 8 px under the first), features-text--text (the stacked track
  names in the label column are 2 px apart and touch).
- Cause: `sectionGap: 8`, `trackGap: 4`, `text.bandGap: 2` (config.ts:261-262, 338).
- Proposal: `sectionGap` 8 -> 20 (a section is a unit; grouping by proximity), `trackGap` stays 4, text `bandGap` 2 -> 3.
  Finding 1 shortens sections, so 20 px costs little.

## 5. No default colour per layer: unstyled ranges are one steel blue

- Stories: examples-p53-gene-by-region (eight bands, one grey-blue #7b8fa6-like), examples-plasmid-ring ("pUC19c" and
  letters fill, a grey well), features-files--* tables. Reads as "placeholder".
- Proposal: a layer with no colour gets the next colour of the theme's categorical palette by layer order (per layer, not
  per item), so two layers in a track differ and a gallery looks designed. Cheapest large "pretty by default" gain; ask: it
  changes the rendered output of every colourless spec (a host's `Library, not a viewer` concern is only about analyses,
  not looks).

## 6. Legend: group titles look like chips

- Stories: examples-antibody-views--schemes ("Sequence" + 6 chips and "CDRs" + 3 on one line, no separator),
  examples-trna-views--arcs (two lines, three titles), examples-plasmid-ring ("Forward strand ... Reverse strand" inline).
- Proposal: one group per line; title muted, weight 600, 11 px; 12 px between groups and 12 px under the figure. Why:
  grouping must be visible without reading; today the title is the same size and colour as a chip label.

## 7. Contrast of muted text and pale scheme tints (light)

- Stories: all light shots. Ruler tick labels and "not numbered" (muted grey on the beige well) look about 3:1 by eye (not
  measured); the scheme pastels (hydrophobic vs aromatic vs polar) are close in the legend chips.
- Proposal: measure with the theme tokens first; if under 4.5:1, tick labels at the text colour 75%; scheme tint for the
  zoomed-out density strip up one step. Dark theme: no problem seen (confidence, variants dark are clear; the sequence band
  is dark and desaturated, the data colours lead, as wanted).

## 8. The viewport bar under the figure is the loudest mark

- Stories: those with a ruler (antibody, confidence, p53 tracks, plasmid ring, gene-by-region): a 2 px full-width coloured
  line with end caps (brown light, blue dark), drawn under the last track (p53 tracks: below Conservation, away from the
  ruler). It competes with data; unexplained.
- Cause: not traced (the navigator/window bar styling, navigator section of config.ts). Proposal: 1 px, muted at rest,
  accent only while hovered or dragged. Not decided anywhere I found: ask the author.

## 9. Ring: dense leader fan and Canvas/SVG mismatch

- Stories: examples-plasmid-ring--ring (+ detail-with-a-strip, rotated-linear-view, features-across-the-origin).
- What: 13 restriction-site labels (NdeI ... HindIII) fan out from a ~10 degree arc, leaders crossing the ring; Canvas ring
  is drawn with a tan outer window circle that the SVG does not have, and the SVG ring and the strip below it sit 40-50 px
  higher than the Canvas ones in the same story.
- Cause: not traced (the label sink, `rings.leaderBend 12`, `labelGap 8`; the window ring is the "followed" window).
- Proposal: when more than ~6 labels fall within 12 degrees, draw one grouped label ("13 sites", expand on hover or
  zoom); the window circle only when the window is smaller than the ring; the SVG host reserves the same top row as the
  Canvas.

## 10. Text layout

- Stories: features-text--text (+ centered, under-tracks), complex-views--text.
- What: the three track names ("p53 / Secondary structure / Domains") stacked at the first line only, tight; number labels
  (10, 20, ...) sit close above the letters.
- Proposal: `lineGap` 8 -> 10 (config.ts:337); drop the stacked names (the legend carries them) or put each at its band.

## 11. Globins rows

- Stories: examples-globins-rows--family (also folded, figure, scale, lazy).
- What: the brown "Species, Identity hidden" header text is state shown as a column title; the logo's bottom stack
  of letters overlaps into an unreadable pile; the SVG species column is "Spec..." / "Glyc..." (6 chars).
- Proposal: hidden columns show no header text (a small chip with a tooltip); the logo draws letters only above a minimum
  glyph height (about 9 px) and a flat scheme-coloured block below; species width measured to a cap (as finding 3).

## 12. Story scaffolding, not library defaults

Pairs of Canvas/SVG have 100+ px of empty space between description and figure and below the legend; the SVG twin starts
lower than the Canvas (navigator bar space on Canvas only). Story-side, listed because a reader judges the first view.

## Fine as is

Letters (20 px, 9-14 px steps), ruler band (22 px), values strip (32), arcs (64), keys and undo chips, code panel, the
AlphaFold bin colours, ring hole title. Light and dark invert correctly in every sampled story.

## Order of work

1. Content heights with the label lane (1 and 2): one change, touches ~40 stories.
2. Measured label column (3), section gap (4), per-layer default colour (5).
3. Legend grouping (6), contrast check (7), viewport bar (8), ring cluster (9), text and globins spacing (10, 11).

Stories visited: all 98 `seq-next-examples-*` and `seq-next-features-*` stories in light and dark; failed to render:
none. Nine opened and read in full (listed at the top); the others were rendered and saved, not individually inspected.
