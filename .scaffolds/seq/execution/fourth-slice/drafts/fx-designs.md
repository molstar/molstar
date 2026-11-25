# Fix-item designs (FX drafts)

Seven items from `fixes.md` that asked for a design. Each: what is seen and why (code), options, a recommendation, the
Configuration row if one is added (a second use, else a plain default), and whether it is a fix after all. Causes marked
"read, not run" come from reading the code; the lane that builds it confirms them with a failing test first.

## 1. The rows rail and the vertical overview together

**Seen.** A declared vertical overview replaces the rows rail (navigators.md, "The rows rail": "a declared overview
replaces it on the rows it binds"). The author wanted both ("also have the thin slider"). The sentence is the cause: the
rule that lets a declared navigator replace the slim ruler rail was copied onto the rows rail.

**Why the two differ.** The ruler rail and a strip are two pictures of one axis in the same place; replacing is right.
The rows rail is 12 px at the band's right edge, the overview a 96 px column: the rail is the precise, always-there
control (grips, `[` `]`), the overview the picture. They agree on the window already ("several navigators on one
section must agree on the window", navigators.md self-check).

**Options.**

- A. A param on the overview (`rail: true`). Cost: one more param on a node whose params are the form's presets; every
  host that wants both repeats it.
- B. A declared rail form: the rows rail becomes `form: 'rail'` with `axes: rows`, declarable like any navigator; the
  overview and it side by side. Cost: the default rail would then also be a declared node, two ways to say one thing.
- C. A viewer-wide rule: the rows rail is drawn beside a declared overview unless the host turns it off. Cost: a config
  item; benefit: no spec change, existing specs gain the rail, and a narrow frame (item 2) has its fallback in place.

**Recommend C.** The spec says what the figure shows (the overview); the rail is a control, and controls are the viewer's
(as `navigator.form` is). Keep the replacement rule for the positions axis only. Both draw against one window, so a drag
on either moves the other in the same frame.

| Item | Where | Default | Why | Alternatives |
|---|---|---|---|---|
| `navigator.rowsRail` | `ViewerConfig` | `always` | the author's choice; the overview is a picture, the rail a precise control | `unlessOverview` (the old rule), `none` |

Also: the overview's picture goes to the left of the rail and the rail stays at the band's right edge beside the
numbers; the rail is idle by its own rule when every row shows. Not a fix: it changes navigators.md. The fix lane's
stories (both forms, then together) are right as written.

## 2. A narrow frame

**Seen.** Vertical overview fixed at 96 px (`NavigatorSizes.rows`, `render/navigator.ts`), the label column 80 to 160 px
(layouts.md), tree and number columns; at 480 px the plot has a few pixels. The plot box is fixed once the columns are
subtracted (grids-and-rows.md), so the last thing laid out is what is left.

**What gives way, in order** (each is a rule with a floor, the ruler's rule: never a column squeezed to nothing):

1. The overview narrows from 96 toward 48 px (the picture is a summary; 48 px keeps the window and handles hittable).
2. Below a plot width (here 160 px) the overview is hidden and the rows rail stays (item 1 makes the rail always there:
   the way to move is never lost, the same rule as a hidden ruler falling back to the slim rail).
3. The label column shrinks to its 80 px floor; the tree column collapses to its toggle glyph.
4. If the plot is still under the floor, the plot stays and the frame scrolls horizontally (never a hairline plot).

The order is by what is cheapest to lose: a miniature before names, names before data.

**Options.** (a) the order above; (b) only hide the overview at a fixed frame width (breaks in a narrow host panel with a
wide window); (c) the overview a share of the frame (3 columns of different widths per frame: unpredictable).
**Recommend (a), keyed on the plot's width, not the frame's.**

| Item | Where | Default | Why | Alternatives |
|---|---|---|---|---|
| `navigator.overview.minPlotPx` | `RenderConfig` preset | 160 | below it the data is the thing to keep | any px; 0 never hides |
| `navigator.overview.minPx` | `RenderConfig` preset | 48 | the narrowest picture that still shows the window | any px |

Decision change (navigators.md sizes, grids-and-rows.md column order); the overview's own size (96) stays the default.

## 3. The Reverse track reserves an empty label lane

**Seen (read, not run).** `laneGeometry` (`render/ranges.ts`) adds the label lane when `content.labels`, which
`labelsBeside` sets from the whole data at the fit view (sizing.ts `LayerContent.labels`). The content sizing gives
tracks of one kind the tallest height (`ContentSizing`), so Forward (labels beside) and Reverse (none, or not placed)
share a height; and the strip wraps the linear detail, so a label that needed the lane on one wrapped line leaves the
lane on every line. This is the fit-time rule working as written: the lane must not appear with a pan (layouts.md,
"A track's height by default fits its content").

**Options.**

- A. Keep the rule, accept the empty lane (the cost: a band with a 14 px gap, "feels empty", the author's own earlier
  complaint about the 72 px band).
- B. Reserve the lane per track by its own labels, not the shared kind height: sharing is for heights (lanes), not for
  where the block sits; Reverse's block is centred in its height with no label lane, as `laneGeometry` already does when
  `content.labels` is false. Cost: none for the pan rule (still decided by the whole data); the two tracks' marks no
  longer sit at one y. Whether that happens today needs the failing test.
- C. Reserve per line in a strip (cost: heights vary with the line, the very thing the rule forbids).

**Recommend B**, after the lane confirms which of two causes it is (shared kind height reading Forward's lane, or the
lane set from labels the strip's lines place elsewhere). Either way a fix, not a design change: layouts.md already says
"a label lane counts as a lane when a label goes beside its mark", of that track. If the lane confirms the shared height
carries the lane, one sentence goes into layouts.md ("tracks share the tallest lanes; each its own label lane"). No
configuration row: no second use.

## 4. Dense leaders when a ring group opens

**Seen (read, not run).** The fan is the full label set placed as if no group existed (`layoutLabels`, `layout/rings.ts`):
more than six labels within 12° (`rings.cluster.over`, `.degrees`) were found dense because they are close, so placed
singly they are still close: each takes a line apart (`placeLabels`, spread by `rings.labelGap`) and the leaders fan over
the ruler with `rings.leaderBend`. The group opens to every member, so a group of 13 gives 13 leaders in the sector of
12°: dense by construction; the precompute is correct, the picture is what a full fan is.

**Options.**

- A. Fan wider: spread the opened members over more arc than their own (a fan pivots around the group's middle angle,
  neighbours of the group keep their place since the fan is precomputed with them). Cost: a member's label stands away
  from its mark; leaders longer; risk of overlapping a neighbouring group's labels.
- B. Open as a list in the hover card, not a fan (the card already lists the members): the fan only for groups up to a
  second bound (`cluster.fanMax`, 12), above it the group stays one label and the card is the opening. Cost: the author's
  "relayout" is only partly done for large groups.
- C. Two columns: a large fan's labels in two staggered radii (alternating `labelGap` rings), leaders short and
  non-crossing. Cost: new layout code; a ring with a second label radius.

**Recommend B with A's wider pivot up to the bound**: a fan stays only while it can be read (a line apart over an arc
that does not reach a neighbour), else the card. The fan's rule "a leader never crosses another label or mark" (layouts.md)
already drops labels it cannot place; B makes the drop explicit instead of a ragged fan. Not a fix: a decision (what a
group above the bound does).

| Item | Where | Default | Why | Alternatives |
|---|---|---|---|---|
| `rings.cluster.fanMax` | `RenderConfig` preset | 12 | past a dozen lines a fan is a hairball; the card lists them | any count; `Infinity` the present behaviour |

## 5. The ruler's choice among insertion labels

**Seen.** `fit` (`render/ruler.ts`) draws majors first, then minors, and a minor label only `if (w > lay.axis.scale)
continue;`: only where it fits its own position's width. At the closest zoom `112.1` at `10px mono` is about 30 px, wider
than a position, and so is never shown even where its neighbours' labels (`112.2`, `112.3`) are also suppressed, leaving
the comb of ticks with no number. The rule was to stop an insertion crowding a round value. Short CDR labels (H1 by
Kabat) at the fit view are the same family: a label that fits its mark only sometimes (layouts.md, "Labels inside if
they fit, else beside" and the ellipsis rule already in `fixes.md`), not the ruler's.

**Options.**

- A. Drop the scale rule; keep only the clearing test against what is drawn (majors first, then minors in position
  order). Cost: at a zoom where three insertion labels are a residue apart, the first wins and the other two are
  skipped; which one is by order, not by use.
- B. A: and order the minors by priority: the one nearest a drawn major first, then the shortest label, then position.
  The reader gets a number between every pair of majors, and `112.1` (the first insertion after 112) is labelled before
  `112.3`. Cost: a small sort of the minors at layout (no allocation: the scratch array exists in `taken`).
- C. Stagger: draw overlapping insertion labels on a second row under the first. Cost: the ruler's height changes with the
  zoom (violating the fixed-height rule).

**Recommend B.** One rule for the ruler: every label is placed if it clears what is drawn; the order decides who wins.
Because the labels are decided at layout (never on a hover), no jitter; the hit card names the rest. A fix to the rule in
axis.md/ruler section of the core spec ("an insertion's label only where it fits its own position" is replaced), not a
new feature. No configuration row (the gap `ruler.labelGap` exists).

## 6. Where an opened source's issues are reported; `seq.resolve.length.unknown` before open

**Seen.** Two cases, one gap. (a) A section over a source sequence needs a sequence named after the source's system; if
none is, the length is unknown, nothing lays out, and nothing says why: `seq.resolve.length.unknown` is raised only for
a wrapped row or an exon on a cyclic sequence (`resolve/model.ts` lines near 746 and 842), never for "no length at all".
(b) The same code is raised before a source opens (a source's coverage is not yet known) and then clears: it flashes in
the viewer's issues while a source opens (viewer spec, "Data": "an opened source's own `issues` are reported there").

**Where reports live.** Static model issues come out of core's resolve (`model.issues`); the viewer's store holds
sources' issues by address (`models/store.ts`, `issues`) and `viewer.layoutIssues` is the frame's. The viewer is the only
holder of "pending" (data-store.md, entry states `pending`, `loaded`, `empty`, `failed`).

**Options.**

- A. Core reports a missing length always, as `seq.resolve.length.unknown`, and the viewer filters it while the source is
  `pending`. Cost: a viewer rule that knows one code.
- B. Core reports the length only once the resolver is told the sources are settled (a flag of the resolve call: the
  static model has `lengthsKnown`); before, it says nothing; after, the code. Cost: an input to a pure function; the
  model changes twice per open.
- C. Two codes: `seq.resolve.length.unknown` (a cyclic wrap waits for the length) and `seq.resolve.system.unnamed` (a
  section over a source with no sequence named after the source's system: the author's finding, an error, never
  transient). The pending-length code stays a warning; the viewer shows it only after the source's entry leaves
  `pending` (a settle, the same rule as the status line's "rows have not arrived"). Cost: one more code.

**Recommend C.** The unnamed system is a different fault with a different fix (name the sequence), so it is its own code
with a hint saying so, raised by core at resolve (a spec fault, never data-dependent). The transient code is held back by
the viewer until its address's entry is not `pending` (data-store.md owns the state, viewer-models.md says the rule: "an
issue of a pending address is not shown"). The viewer says it as an entry in the issues state and in the status line the
stories already show; a host reads `viewer.layoutIssues` / the model's issues.

No configuration row. The inline `strand: null` crash is a plain bug (fixes.md) and unrelated. The new code goes in core's
issue catalog and the resolve section of the spec. Partly a fix (the missing report: a bug), partly a decision (the
viewer's hold-back rule).

## 7. "A ruler naming its chain stays labelled"

**Seen (read, not run).** A ruler track's label is empty by default; an against view with two numberings (Complex, Antibody
against) gives each ruler a label naming its chain, so the reader can tell the two apart, and the label stays drawn
(the lane and `labelText` of `layout/tracks.ts`). The sibling finding in the same line, an empty track label still drawn as
empty text, is a plain bug (a draw of "" costs a text op and a hit).

**Verdict.** Intended: with two rulers over one view the name is the only signal which chain's numbers are which; hiding
it would hide the numbering's owner. It is the spec that names the chain (a ruler's own label), not a default the viewer
adds. So: not a finding. Say it in track-header-parts.md or layouts.md in one line ("a ruler is unlabelled unless its
spec names it; with more than one in a view, name them"), and let the viewer's accessibility name use it. Ask the author
only if the intent was the opposite (labels dropped when the ruler is hidden under the section's title); then the rule
is that a ruler's label hides with the section's title, and the lane that owns tracks decides.

No configuration row. The empty-text draw stays a fix: skip a track whose label is empty.

## Decided by the author

2026-11-20, each as recommended:
- 1, the rows rail beside the vertical overview: "Both, configurable" (`navigator.rowsRail: always | unlessOverview |
  none`, `always` by default; the replace rule on the positions axis only).
- 2, a narrow frame: "As proposed" (the overview narrows 96 → 48 px, hides below a 160 px plot with the rail kept, then
  the label column to its floor and the tree to its toggle, then the frame scrolls; `navigator.overview.minPlotPx`,
  `minPx`).
- 4, ring groups: "Fan up to 12" (`rings.cluster.fanMax`, a larger group listed in the hover card).
- 6, an opened source's issues: "As proposed" (`seq.resolve.system.unnamed`; `length.unknown` held back while the address
  is pending).
3, 5 and 7 are fixes (the drafts' causes confirmed by a failing test first). All seven go to the next fix lane, FX2,
with navigators.md, layouts.md and the issue catalog edited by the orchestrator at its merge.
