# Values marks: reference lines on a values layer (points and the log scale are built)

<!-- Agreed by the author, 2026-11-24 ("Agree all"): close-out CL8; features 130 (reference lines). A layer param
     within the `0.11.0` bump track-groups.md makes. Sections marked **Provisional** say what would settle them. -->

## Intent

A reader of a genome-wide association plot sees the SNPs as dots at their height (-log10 p) and a line at the
genome-wide significance threshold (7.3, p = 5e-8): the dots over it are the finding. A reader of a signal track sees
a zero line or the mean it is read against. A writer says where the line is, once, in the spec, and it stays right at
every zoom, on a linear or a log scale, in light and dark.

What row 130 says, and where each part stands:

| Part (EC) | State | Closed by |
|---|---|---|
| Points at a height (216) | **built** (the fifth slice's VP: `encoding: 'points'`, `render/values-marks.ts`, language `0.8.0`): a dot a value, a thin min-to-max mark a device column below a pixel a value; a coarse level keeps min and max (core spec, Data, summaries), so a zoomed-out bin never averages one significant SNP away | nothing |
| Log scale (215) | **built** (`scale.type: 'log'` on y and on colors, `non_positive`, the log legend; EC-215 in themes.md, "Palettes and scales") | nothing |
| Reference line (217) | **built** (CL8): `reference` draws a line at a chosen value; the baseline hairline stays | this design |
| A scale's title and units (217) | **built** as `scale.label` and `scale.unit` for the legend; the values layer draws no y axis, and this design adds none | out of version 1 (Tensions) |
| Split by sample (214), colour by chromosome (218), LD (219) | other rows | not here |

Library work: one layer param, the values layer's draw (a line and an optional label), the SVG export, the check.
Application-level: nothing; the threshold's value is the writer's.

## Precedents

- UCSC browser (a wiggle track's "Draw y=0 line", `yLineOnOff` and `yLineMark`: one horizontal line at a value, drawn
  over the data): taken, a line at a value, off by default, a value not a keyword; not taken, one line only, and the
  line having no name.
- IGV (a track's data range with a baseline value, log scale; its bar charts rise from the baseline; a track has no
  threshold line): taken, a line is read through the track's own scale (log included); not taken, the data range as a
  dialog (a `scale` node and `domain` do it already, EC-201).
- JBrowse's wiggle tracks (`bicolor` pivot, `scaleType: log`, autoscale; no reference line of their own, a pivot
  colours either side): taken, a log scale shared with the track's other layers; not taken, a pivot that recolours the
  data (a theme's rule does that: `(rules (> value 7.3) …)` colours the dots over the threshold).
- LocusZoom (Manhattan plots: a horizontal "genome-wide significance" line at 7.301 and a "suggestive" one at 5, dashed,
  labelled in its legend; the points above are not recoloured): taken, two lines on one plot, dashed, a label; not
  taken, the label in a legend (here at the line, or the accessible text).

## Decision

### 1. A `reference` param on the values layer (a language change, one minor bump)

A `values` layer takes `reference`: a list of up to four lines, each a record `{ at, label?, color?, style? }`.

- `at` (number, required): the value, in the y scale's values (a log scale reads it as a value, 5e-8 on a log scale
  stays 5e-8).
- `label` (string): its name for people ("genome-wide", "mean"), drawn at the line and read by the accessible text.
- `color` (a color, a token or `light-dark(a, b)`; as a layer's `color`, never a theme form: a line is one colour).
  Absent: `msui-text` (the author, 2026-11-24: 3:1 against the ground and `seq-value`).
- `style` (`solid`, `dashed`, `dotted`, the ranges' `line` words; absent: `dashed`).

The line is a hairline across the whole well (1 CSS px, never thinner than a device pixel), drawn **over the bars, the
profile and the dots and under the hover**; it is never hit and never hovered (a line is not an item; its value is in
the label and the accessible text). The baseline hairline stays as it is; a line at the baseline's value draws over it
and wins.

### 2. Through the layer's scale, like a value

A line's y is `ys` of its value over the layer's domain, with the headroom: a line at the domain's top sits where a
value there sits, so a dot at the threshold and its line coincide. On a log scale the line at `at` is where a value
`at` is; `at <= 0` there is not on the scale (the layer's non-positive rule: no line).

### 3. A line outside the domain, and an open domain (**Provisional**)

- A **closed** domain (the layer's `domain`, or the scale node's both ends): a line outside it is not drawn, said once
  per layer (`seq.render.reference.outside`, info, the line's index and value). Why: clamping to the well's edge
  would claim a position it does not have.
- An **open** end: the line's value is part of the data the end is computed from (the union of EC-56 gains the lines),
  so a Manhattan track whose tallest SNP is 6.1 still shows its 7.3 line, and a plot with no significant SNP is
  seen to have none. Why: the threshold's absence from a plot is the finding. Cost: lines move a shared scale's
  domain, so two tracks on one scale agree.
  **What settles it**: the author's reading of the Genomes story with a line at 7.3 over data to 6 (a flat plot with
  a visible line, or the line missing).

### 4. The label (**Provisional**)

A label draws at the well's right end, above the line, 10 px in `msui-text-muted`, over a patch of the well's ground,
only when the well is at least 40 CSS px tall and the line is not within one label height of the well's top (else
below the line); two labels that would overlap, the later one dropped (its line stays). The text is "label · value"
("genome-wide · 7.3", the value as the legends write one); a line without a label draws the value alone (the author,
2026-11-24, close-out B11: the plot has no y axis, the value is read at the line). The accessible text is unchanged. **What settles it**: the Genomes story at the default well height (the label legible and out of the dots' way),
or the author preferring the label in the accessible text and the hover of the whole track only.

### 5. Marks, export, accessibility

- A solid line is a rect (as the baseline's, no allocation warm); dashed and dotted are a path with a dash (the draw
  lists have one), so the SVG export draws a dashed line and a canvas draws the same.
- The accessible text of the layer gains "reference line at 7.3 (genome-wide)" per drawn line, after its values.
- A warm frame of each encoding with two lines allocates nothing (a case in the perf tier beside the encodings').

### 6. Not here

A fill between two values (a band), a line from a field (a per-row threshold), a line named from another layer's data,
a baseline other than zero or the domain's low, a y axis with ticks, a log axis's decades. Version 1's need is a
threshold, a zero and a mean.

## Defaults

With no `reference` the layer is as built: no line. With a line, no other setting:

```
 light                                         dark
 ┌──────────────────────────────────────┐     ┌──────────────────────────────────────┐
 │                       ▪                │     │                       ▪                │
 │ - - - - - - - - - - - - - - - genome-wide    │ - - - - - - - - - - - - - - - genome-wide
 │        ▪        ▪   ▪       ▪        │     │        ▪        ▪   ▪       ▪        │
 │ ▪ ▪▪ ▪▪▪ ▪▪ ▪▪▪ ▪ ▪▪ ▪▪ ▪▪▪ ▪▪ ▪ ▪▪ ▪ │     │ ▪ ▪▪ ▪▪▪ ▪▪ ▪▪▪ ▪ ▪▪ ▪▪ ▪▪▪ ▪▪ ▪ ▪▪ ▪ │
 └──────────────────────────────────────┘     └──────────────────────────────────────┘
   baseline hairline at 0 (seq-sequence)        line: msui-text, dashed 4 on 4 off, 1 px
```

- **Light and dark**: tokens only (`msui-text` for the line, `msui-text-muted` for the label, the well's ground behind the
  label); the line reads against the well's ground at 3:1 in both themes and against `seq-value` at 3:1 in light (a test
  measures them). In dark `seq-value` is a light blue: no ink is 3:1 against it and the near-black ground (1.75:1 for
  `msui-text`): carried. A writer's own colour is theirs to check.
- **A line outside the domain**: closed domain `[0, 6]`, line at 7.3: not drawn, one info. Open domain: the domain
  becomes `[0, 7.3]` plus headroom and the line is drawn (Decision 3).
- **Two lines**: 7.3 ("genome-wide", dashed) and 5 ("suggestive", dotted, `color` set): two hairlines, two labels at
  the right end; when the wells' heights bring them within a label's height the second label is dropped, the line stays.
- **A log scale**: values 1 to 10^4, line at 100: halfway up the headroomed span, between the decades 10^1 and 10^3
  at their even steps; a line at 0 or below: none, the layer's log info already says why (`seq.render.scale.log`).
- **A line at the baseline's value** (0): over the baseline hairline, dashed in the muted ink over the solid one; the
  writer who wants only the baseline sets none.
- **Narrow and crowded**: the label cut to the well's right third with an ellipsis; a well under 40 px shows lines and
  no labels. **Zoomed out to a genome**: unchanged (a line is in values, not positions). **Empty** (no values in
  view): the line draws (the threshold is not data), the well's pending and failed spans as built.

**Mock**: none (no gesture). **Precedent**: LocusZoom's genome-wide line. **Reference picture**: the Genomes dataset's
"Manhattan" story (the lane builds it from the fifth slice's GWAS generator step, deferred with EC-216), light and dark,
canvas and SVG.

## Configuration

| Item | Default | Why | Alternatives |
|---|---|---|---|
| a values layer's `reference` (spec) | absent: no line | a figure draws what its writer asks; a zero line is the baseline already | a scale node's `lines` (Rejected) |
| a line's `at` (spec) | required | a threshold is a number the writer knows | a field's value (not built) |
| a line's `label` (spec) | absent: none | a label is a claim | the track's title (the writer's) |
| a line's `color` (spec) | `msui-text` | the series colours (by chromosome, by value) are the data: the line must not borrow one | any color or token |
| a line's `style` (spec) | `dashed` | a threshold is conventionally dashed and apart from the baseline's solid | `solid`, `dotted` |
| lines per layer (language bound) | 4 | a threshold, a suggestive one, a zero, a mean | any (not built) |
| `values.reference.widthPx` (render config) | 1 | the baseline's and the hairlines' weight | any (not built) |
| `values.reference.dashPx`, `.dotPx` (render config) | 4, 1.5 | a dashed and a dotted line, on then as long off (the draw lists' dash is one number: 4 on 3 off was not built) | any (not built) |
| `values.reference.labelPx` (render config) | 10 | the quiet weight of the axis labels | any (not built) |
| `values.reference.labelMinWellPx` (render config) | 40 | under it, text hides the data | any (not built) |

## Rejected alternatives

- **A new node kind** (`reference` or `guide` under the section, drawn across layers): a line belongs to one scale
  and one y; a node would need an address, a track and a slot, a language change larger than a param for no reader's
  gain, and a threshold drawn by two layers on one scale is two params.
- **Lines on the `scale` node** (`scale.lines`: every layer naming it draws them): attractive (the threshold written
  once for a Manhattan track and its bars sharing the scale), but a scale also colours legends and sizes, a bars layer
  and a points layer would both draw it, and the legend would need to list it. **Provisional**: the lane keeps the
  layer param; what would change it is a story where two layers on one scale both need the line.
- **A layer on the track** (`values` over a constant series): the line would be a data layer, unreachable by the
  scale's open end (Decision 3), cut by the data's gaps and bins, and a legend entry.
- **A threshold as a theme rule only** (colour above 7.3): shows the dots over it, not the line; the rule stays
  available beside it.
- **A single `reference: number`**: two lines is LocusZoom's default and a zero plus a mean is as common.
- **A theme form for the line's colour**: a line is one colour (Configuration).
- **Splitting points, log and lines into three designs**: points and log are built; only the line is new.

## Consequences

**The language** (a change: its own lane, the author's word before it starts):

- `lang/language.ts`: the layer's description gains `reference` (a `list` of `record`, `max-length` 4, fields `at`
  (number, required), `label` (string), `color` (string), `style` (enum `solid`, `dashed`, `dotted`)), its doc "For
  `values`: horizontal lines at values of the layer's y (a threshold, a zero, a mean), drawn over the marks."
- `SeqLanguageVersion` takes the next minor. Both landed in the same unreleased `0.11.0`, the track-groups bump first
  and `reference` after it. Its history entry, the builder's types (`values`' options),
  the catalogs and `config-keys` tests.
- `lang/check.ts`: a `reference` on another layer type is the existing "param for another type" warning; `at` a
  non-number, a colour that is not a color or token, a style outside the enum, a fifth line: the description's own
  checks. New codes: `seq.render.reference.outside` (info; the line's index, value and the domain). A test checks the
  catalog.
- Nothing else in the language: no node, no kind, no slot, no theme form. A file at the earlier version never has it.

**The render**:

- `render/values-reference.ts` (new: `values.ts` was at 775 of 800) holds `drawReference(lay, out)`, called after the data
  pass by each encoding, the labels' pass and the resolved lines (`ValuesLay.ref`: value, y once per layout, color, dash,
  label, in preallocated arrays); `values-domain.ts` takes the open domain, which includes the lines for an open end
  (Decision 3).
- The SVG export and the PNG: through the draw lists, nothing of their own; the accessible text (`accessibility.ts`).
- The legend and tooltip: unchanged.

**Tests** (core, `render/_test/values-look.test.ts` and a new `values-reference.test.ts`, quick lane; one medium for
the story): a line at a value sits where a value at it sits (bars, line, points, log); outside a closed domain: no
line, one info; an open domain includes it; at or under zero on a log scale: none; two lines; the label dropped
under 40 px and when overlapping; dashed in SVG; light and dark contrast; warm allocation 0 B (the perf tier's case
beside the encodings'); the check's new rows; the catalog test.

**Specs and examples**: the core spec's Render › Layer types › `values` and its Configuration table; the language
section's layer param list and version history; `themes.md`'s Values row and `layouts.md`'s "values look" cite it.
`examples/src/genomes` gains the "Manhattan" story (points at a height, a line at 7.3 and one at 5, over the GWAS
generator step; coloured by chromosome is row 218's, not here). `features.md` row 130 and its parts-not-built entry
close; `examples.md`'s deferred list loses "a reference line".

**Edge cases**: EC-217's reference-line part closes (kept, a test named); EC-217's title and units stay with the
legend (`scale.label`, `scale.unit`), the axis part is out of version 1. EC-216 and EC-215 are kept (built; the lane
adds a test that a coarse level's bin keeps its one significant SNP if no test names it today: **Provisional**, settled
by a grep of `values-look` and the summary tests for a "max over a bin" case).

**Matrix**: no gesture, no cell.

## Tensions

- **A line moves an open domain** (Decision 3) against "values are drawn at their value; the domain bounds, never
  the value" (EC-55): the lines are not values, but they shape the bound. Settled for the threshold's visibility;
  carried as Provisional.
- **A label is text in a layer that has none**: the values layer has drawn no text; a label brings a font, a clash
  rule and a width. Settled narrowly (one size, a minimum well, dropped on overlap); the fallback is the accessible
  text alone.
- **No y axis**: a reader cannot read 7.3 off a plot without ticks except by the label. Carried: a values layer's axis
  is a larger design (a scale's title, ticks, units: EC-217's rest) and not version 1's.
- **-log10 p on a log scale reads backwards**: p-values on a log scale put small p low; the Manhattan convention is a
  column of -log10 p on a linear scale, which the writer's data supplies (a transform of the column is the data
  step's, not the layer's). Carried; a reversed y is not built.
- **A language change this late**: one optional param on an existing layer, a minor bump, no new kind; the smallest
  that fits (Rejected: node, scale lines).

## Self-check

- **The claim**: a threshold line is a property of one values layer's y, so it is a param of that layer.
- **Negated**: the line is a property of the figure across tracks (one threshold over several tracks) or of the data
  (an annotation), not of one layer's y. The negated world explains a section-level node (rejected: each track has
  its own y and domain, and a Manhattan has one track), or a values layer of constants (rejected: data semantics,
  gaps, bins). The layer param explains the one-track need and the shared-scale case by repeating the value (carried
  as the `scale.lines` alternative, Provisional).
- **Its inverse**: a writer who wants no baseline: not built, not asked; the baseline hairline is the layer's.
- **Built**: points (VP), log (EC-215), baseline. the reference line (CL8). **Carried**: Tensions above.

## In the specs

| Decision | Where |
|---|---|
| the `reference` param (the list, its fields, bound 4, absent: none; the layer params' list is in Language's prose) | core › Render › Layer types |
| the `reference` param's version (`0.11.0`) | core › Versions |
| the line is drawn over the data and under the hover, a hairline, never hit, an outside line not drawn | core › Render › Layer types |
| the line reads through the layer's scale (log included, `at <= 0` none); an open end includes the lines (Provisional) | core › Render › Layer types |
| an open end includes the lines (Provisional) | core › Themes |
| the label (right end, 10 px, 40 px well, dropped on overlap) (Provisional) | core › Render › Layer types |
| `seq.render.reference.outside` | core › Render › Issues |
| the render config keys `values.reference.*` | core › Render › Configuration |
| the accessible text | core › Render › Export |
| the Manhattan story | design-only: an example, listed in `examples.md` |
