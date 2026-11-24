# Track groups: a `group` node over tracks, folded as session state to one summary row

<!-- Agreed by the author, 2026-11-24 ("Agree all"): close-out CL6; features 159 (groups and summaries), 70 and 141
     (groups of tracks, EC-103), 185 (collapse and groups) and 128 (collapse). A language change, `0.10.0` → `0.11.0`,
     shared with values-marks.md's `reference` param: one bump, in the language lane. -->

## Intent

A writer puts a section's tracks under named groups ("Features", "Variants", "Signals"); a reader folds a group to
one row that still shows what its tracks hold, and opens it again, as often as they scroll. A section of 30 tracks
then reads as six rows, and a figure can open with the groups it cares about open. A ranges track can also be drawn
in one lane with its labels off (the Ranges row's "collapse"), the form a folded group draws its ranges in.

What each row says is missing, closed here:

| Row | Missing (features.md) | Closed by |
|---|---|---|
| Tracks layout (159) | groups and summaries (EC-103) | the `group` node, the summary row |
| Trees (70), Collapsed clades (141) | groups of tracks (EC-103) | the same (the rows layer's clades are unchanged; a group holding the rows track folds it to its consensus row) |
| Track and group editing (185) | collapse and groups | `groups.*` commands, the group header, hide and move of a group |
| Ranges (128) | collapse | the ranges layer's `lanes: one` |

Library work: the language (a kind), the resolved view (a track's group), the tracks layout (headers, the summary
band), the layer types' summary forms, the viewer's session state and commands. Nothing here is application-level:
no group editor, no "group by attribute" (Rejected alternatives).

## Precedents

- **The previous implementation** (`packages/sequences`, read only): a `group` node under a section (`label`,
  `collapsed`, `hidden`), groups nested; collapse runtime UI state (`ui.expanded`, a set of group refs, never in the
  document); the `group.toggle` action; a track with `summary: true` shown under its group's header only while it is
  collapsed (`_spec/format.md`, Tracks and layers; `viewer/src/layout.ts` `items`; EC-103 "V/union-summary"); a group
  header a row of its own. **Taken**: the node, its three params, collapse as session state from an opening value,
  the authored summary shown only while folded, a header row. **Not taken**: nested groups (one level in version 1);
  a folded group with no summary track showing nothing (here an automatic summary); the summary as a flag on a track
  (here the group names it: one place).
- **This codebase's clade folds** (grids-and-rows.md, "Opening and closing clades"; `viewer/src/behaviors/rows.ts`,
  `models/view-bands.ts`): folds are session state in the `ViewModel` beside the viewport, the spec's `collapsed`
  their opening set, `clades.keep` one patch writing the set shown, a patch or undo of `collapsed` resetting the
  session's set, a toggle a hand's move and bulk folds and reveals jumps of the view's history, a summary row labelled
  and computed from what it folds, a reveal opening what holds its target. **Taken whole**, so the reader has one
  model of folding.
- **ProtVista / Nightingale** (UniProt's feature viewer): categories ("Domains & sites", "PTM", "Variants") whose
  collapsed row draws every feature of its type tracks merged in one track; expanded, a track a type. **Taken**: the
  merged ranges lane as the folded row. **Not taken**: categories fixed by the data provider (here the writer's).
- **UCSC Genome Browser**: composite tracks (subtracks under one parent, a parent's on/off) and the display modes
  dense (one line, no labels), squish (half-height lanes, no labels), pack (lanes, labels), full. **Taken**: dense as
  the ranges layer's `lanes: one` (one lane, labels off); the composite's parent as a hide of the whole group.
  **Not taken**: a composite folded to each subtrack's own dense line (N rows; here one row); squish as a mode (the
  tracks layout's `track.overflow: compact` already is it, layouts.md "Track sizing"); a mode switched silently by
  count.
- **IGV**: collapsed / expanded / squished per track (collapsed = one row, features overlaid); tracks grouped by a
  sample attribute (session state, a header per value). **Taken**: collapsed as `lanes: one` with overlap drawn
  shorter-over-longer. **Not taken**: grouping by attribute (deferred: a command writing groups as a patch, when a
  dataset of many samples asks).
- **JBrowse 2**: track categories in the track selector (nested), a track minimised to its header bar. **Taken**: the
  category name as the group's label. **Not taken**: categories living only in the selector (here groups are drawn
  in the view); nesting.
- **Jalview**: sequence groups with a consensus annotation; hidden sequences represented by one row. **Taken**: a
  folded rows track as its consensus row (the clade summary over every shown row).

What a user can do there and not here: group by an attribute on the fly (IGV: deferred, above); nest categories
(JBrowse: refused for version 1, additive later); drag a track into another group (JBrowse, IGV: the writer's,
Decision 7); each subtrack's own dense line (UCSC: refused, one row).

## Decision

### 1. A group is a language node (a language change: `0.11.0`)

```
group: { parents: ['section'], ids: 'required',
  params: { label?: string,                    // its id when absent
            collapsed?: boolean = false,       // the opening state (the session's set overrides it)
            hidden?: boolean = false,          // the group and its tracks not drawn; its tracks keep their own `hidden`
            summary?: 'auto' | 'none' | { ref, kinds: ['track'], in: 'group' } = 'auto' } }   // what a folded group draws
track.parents: ['section', 'group']            // was ['section']
```

- **One level**: a group's children are tracks only (templates included); a group in a group is a check error
  (`repr.schema.parent.invalid`, as any parent the description refuses). Why: no example needs nesting, and the header, cursor
  and summary rules stay one level; `group.parents` gaining `'group'` later is additive.
- **Its tracks are its children**, so a track in a group has the address `/views/v/<section>/<group>/<track>`. Why a
  tree over a membership ref: the order is the spec's order and the membership is where the track is, one place to
  read and to edit (Rejected 2). The cost: moving a track in or out of a group changes its address (a `move` patch,
  `repr/patch`'s op with a new parent); session state keyed by a track address (a band's scroll) starts fresh.
- **`summary`** names one of its own tracks by its id (a `ref` resolved `in: 'group'`; checked:
  `seq.lang.group.summary`, a track of another parent or none is an error), shown only while the group is folded, never while open (EC-103); `auto` draws the automatic summary
  (Decision 5); `none` draws the header line alone.
- **The version moves `0.10.0` → `0.11.0`**, pre-release, no migration (core `_spec` › Versions). A `0.10` spec reads
  unchanged. A reader without groups (EC-244, "collapsed groups and summaries" among what another renderer may lack)
  draws a group's tracks as the section's, every group open, its `summary` track left out; the version-aware validate
  (a carried item) lists `group` as a kind a lenient older reader flattens. **Provisional**: the kind's name. `group`
  is the previous implementation's word and ux-foundations' ("collapsing a group"), but the codebase already says
  scroll groups, block groups (`ResolvedSection.groups`, EC-210), copy groups, a text's groups of ten and the theme
  language's letter `group`; `track_group` reads longer and clashes with none. Settled by the author's pick.
- **Codes and the ref's form, as built** (CL6, 2026-11-24): the codebase's own codes (`repr.schema.parent.invalid`
  for any refused parent, `seq.lang.*` for the language's checks) and a ref by bare id in the group, not the draft's
  `seq.check.*` and `ref<track>`. Why: one catalog of codes, and a bare id reads as `across` does.
- **The builder**: `section.group(id, params)` returning a builder with `.track(...)`, as `section.track` does.

### 2. The resolved view: tracks stay flat, each names its group

- `ResolvedSection.tracks` stays the flat list in spec order (a group's tracks in place, depth first), so every
  reader of it today works unchanged; `ResolvedTrack.group?: string` (the group's address); `ResolvedSection.trackGroups?:
  readonly ResolvedTrackGroup[]` (`{ address, label, collapsed, summary?: string, count }`, in order; not `groups`,
  which is the blocks' name). A hidden group's tracks are left out as a hidden track is; the `summary` track is
  resolved with `summaryOf: <group>` so the layout can show it only while folded.
- **One walker** for a section's tracks in the spec (`tracksOf(sectionNode)`, descending into groups): today four
  readers walk a section's children for `kind === 'track'` (`view/view.ts:760`, `render/headers.ts:552`,
  `render/navigator-picture.ts:149`, `viewer/src/models/hidden-tracks.ts`); each moves to it.

### 3. Folding is session state, the clade folds' model

- **Where it lives**: the section's viewport in the `ViewModel` gains `groups?: readonly string[]`, the ids of its
  folded groups (absent: each group's `collapsed`), compared by identity as a band's `folds`. Why session (as the
  clades, grids-and-rows.md Decision 1): a reader folds as often as they scroll, and an undo stack of folds buries the
  edits. Why the section's viewport and not a band's: a group is a section's, and the history already stores and
  compares the viewport (`view-history.ts`).
- **The opening state** is each group's `collapsed`; `groups.keep` writes the set shown into the section's groups'
  `collapsed` (one patch, one undo step). A patch or an undo that changes a group's `collapsed` resets the session's
  set to the spec's (an edit of the opening state is meant to show); other edits keep it.
- **A story step** stores the folded set in its view (a step's view is data). **An export** draws the set it is given
  (`renderSvg`'s `groups` option, the viewer passing its own as it passes `folds`), else `collapsed`.
- **History**: a single toggle is a hand's move (not a step); `groups.collapseAll`, `groups.expandAll` and a reveal
  are jumps, so Alt+Left restores the set and the window together (the clades' rule).
- **A toggle keeps the group's header line where it is on the page** (the section's `top` adjusted, clamped), so the
  pointer stays on what it clicked; heights animate as a clade toggle does (`view.motion`, `view.jumpMs`; a cut under
  reduced motion).

### 4. The header

- **Open**: a header line of `header.groupPx` (20) above the group's first shown track, in the label column: a
  chevron (▾), the label, and on hover its slots and "⋯" (track-header-parts.md's geometry: `header.slots`,
  `header.slotPx`). The group's tracks' labels are indented `header.groupIndentPx` (12). Nothing is drawn across the
  data in the header line. Gaps: `trackGap` (4) between a group's tracks, `groupGap` (8) before and after a group
  (ux-foundations section 3 already says "4px within a group, 8px between groups").
- **Folded**: one row. Its label column is the header (chevron ▸, the label, a quiet count "4 tracks"; slots on
  hover), its band the summary (Decision 5). With `summary: none`, the header line alone.
- **Hit**: target `group-header` (the line, or a folded row's label column, outside its slots); its slots
  `header-slot` as a track's. The summary band's marks hit as their own items (Decision 5).
- **Gestures**: a click on the header (chevron or label) toggles (`groups.toggle`); the `seq:hide` slot hides the
  group (`track.hide` on the group's address: its `hidden`, one patch); a secondary click or "⋯" opens its menu.
- **Its menu** (`headerMenu`, a new entry scope `on: 'group'` beside `track`, `section`, `element`): Collapse or
  Expand (Space); Collapse all groups, Expand all groups (`groups.collapseAll`, `groups.expandAll`, the section's);
  Keep these folds as the opening state (`groups.keep`); Hide group; Show hidden ▸ (its hidden tracks, `track.show`);
  Move up, Move down (`track.move` on the group among the section's children); a host's `on: 'group'` entries. The
  section's menu gains Collapse all and Expand all when the section has groups.
- **Keys** (the `headers` behavior's header cursor, `behaviors/header-cursor.ts`): the group header is a stop in ↑ ↓
  between the section's line and its tracks' (a folded group one stop, its tracks skipped); Space toggles there (Space
  does nothing on a track header today: no conflict); Enter opens its menu (as on a track header). A cursor on a
  track's header whose group folds moves to the group's header.
- **Announced**: "group Variants, collapsed, 3 tracks" (the announcer; the canvas's accessible layer, `aria-expanded`
  on its item, planned with accessibility as the clade node is).
- **Export**: the header line and a folded row's label and count are exported as drawn; the chevron and the slots
  are chrome and are not (track-header-parts.md, Export).

### 5. The summary row: each layer kind in its summary form, composited

A folded group with `summary: auto` draws one band: its shown tracks' layers, each in the form its type declares
(`LayerType.collapse?: 'lane' | 'strip' | 'row'`, optional, additive; not `summary`, which is the ARIA description),
in the group's order:

| Kind | Form | Drawn |
|---|---|---|
| `ranges`, `secondary`, `variants` | `lane` | every item of every member in **one shared lane** (`lanes: one`, Decision 6), labels off; across members later over earlier, within a member shorter over longer; a device column holding more than one item takes the layer's density form (`ranges-density.ts`); variants as their ticks without heads |
| `pairs` | `lane` | each pair's two ends as 2 px ticks in the shared lane, in the pair's colour (arcs need height a row does not have) |
| `values` | `strip` | a strip a member, `group.summary.stripPx` (4) each, below the lane, the values layer's thin look (its colour by value or its series tint; layouts.md "The values look"); at most `group.summary.maxStrips` (4), the rest counted in the header ("+2 not shown") |
| `rows` | `row` | the alignment's consensus row: the rows layer's `summary` form (`consensus`, `logo`, `bar`) over every shown row, the root clade's summary row (grids-and-rows.md, "The summary row"), cached by (layer, entry version, theme) |
| `letters`, `logo`, `ruler`, `matrix`, `grid`, `frame`, `bands`; an extension's type without `collapse` | none | left out, counted in the header's hover card ("Sequence: not in the summary") |

- **Height**: the lane (one ranges lane, 22 px) when any member has a `lane` form, plus the strips, plus a row's
  height; at least one lane. A summary with nothing to draw is the empty-well rule's note row ("no items",
  `track.emptyPx`); pending members draw the one-lane placeholder (EC-169).
- **Hit and hover**: a summary mark hits as its own item (its layer's and track's address, so the hover card names
  the member: "R248Q, in ClinVar"), and is hovered, selected and highlighted as that item; where marks overlap the
  topmost answers, and its card says "+n more here". The band's empty space is empty space (the grammar's: a click
  clears, a double-click fits).
- **The layout lays a folded member's layers in their summary form only** (cache key `<layer>#summary`): a folded
  group packs no lanes and lays no labels, so folding makes a frame cheaper. Its members' preps (tier 1) are shared
  with the open form, so a toggle never prepares again.
- **An authored summary** (`summary: <track>`) is that track, drawn as any track, only while folded. Why keep it
  beside `auto`: a writer may want a density of all variants rather than their ticks (EC-103's case).

### 6. The ranges layer's collapse: `lanes: one`

- A layer param of `ranges`, `secondary` and `variants`: `lanes: 'packed' | 'one'` (`packed` when absent, as now).
  `one`: every item in one lane, labels off, shorter over longer, a column of several items its density form; the
  track one lane high (`lanesHeight(1)`, 22 px), so `track.overflow` never applies. A spec param (the author,
  grids-and-rows.md Decision 8: figure-meaning settings are spec params), set from the track's menu ("Lanes ▸ Packed,
  One lane", `seq:lanes`, menu only; `layer.lanes`, one patch, one undo step).
- **How it relates to a group's fold**: the same code path. A folded group's `lane` form is `lanes: one` over every
  member at once; `lanes: one` on a track is that form on one open track, its header and label kept. They are
  independent: a track's `lanes` survives its group folding and opening.
- **Provisional**: whether a track's lanes form is also session state (IGV's display mode is a session's). Recommended
  as a patch: it changes what a figure shows, and a reader who wants one row for several tracks folds a group.
  Settled by the author.

### 7. Hidden tracks, order and editing

- **A hidden track in a group** is left out as now (its row collapses, nothing takes its place); it is not in the
  summary. A group whose tracks are all hidden draws nothing, header included (as a section whose tracks are all
  hidden, track-header-parts.md). **A hidden group** (`hidden: true`) draws nothing; its tracks keep their own
  `hidden`.
- **The host's list** (`SeqViewer.hiddenTracks`, `mountHiddenTracks`): an entry gains `group?: { address, label }`;
  a hidden group is one entry (`kind: 'group'`, "Variants (3 tracks)"), its `track.show` showing the group. Hidden
  tracks stay the application's to list (the author, 2026-11-21).
- **Order**: `track.move` moves a track among its siblings (its group's, or the section's) and a group among the
  section's children; at a group's edge Move up or down is disabled ("first in its group"). Moving a track into
  another group, creating, renaming or removing a group: not built, the writer's (as a track's label and the track
  list are, track-header-parts.md), each a patch a host's entry may write.

### 8. What a fold does to the rest

- **The selection** is a location: a fold never changes it; selected items of folded tracks show selected in the
  summary.
- **Hover** projected from another section (a structure, a linked sequence) lights the item's mark in the summary.
  A header hover lights nothing (track-header-parts.md's refusal of EC-146 by default, the same for a group).
- **The keyboard cursor** on data is a position on the section's axis: a fold does not move it. On a rows track in a
  folded group it stands on the consensus row (the clades' rule, announced "in group Alignment"), back on its row when
  the group opens. The header cursor: Decision 4.
- **A reveal** (a focus, a link, a search step, a bookmark) landing on an item drawn only by a folded group's track
  opens the group, one jump; Alt+Left folds it again (the clades' rule, grids-and-rows.md Decision 3). A selection or
  a highlight does not open it (its marks show in the summary).
- **Track heights**: the sharing of heights among tracks of one kind ("only when close", layouts.md) reads every
  track, folded or not, so a fold never resizes another track. The summary band takes no part in the sharing.
- **Bridges** (`bridge` nodes naming a track) end on the summary band while the track's group is folded.
- **Navigators** whose picture names a track in a group read it as before (the walker, Decision 2).
- **Layouts**: groups fold in `tracks` (stacked blocks included: the set is the section's, applied in every block).
  In `text` and `rings` every group is drawn open with no header line and the folded set is kept, not applied
  (`seq.render.groups.layout`, info, once, when some are folded), as text unfolds every clade. Why: text repeats its
  bands every line and rings have no label column for a header; a summary there is not designed. `grid` draws its
  tracks' layers over the plot: groups are ignored there.

### 9. Commands

| Command | Kind | What |
|---|---|---|
| `groups.toggle { address }` | session (a hand's move) | fold or open one group |
| `groups.collapseAll { section? }`, `groups.expandAll { section? }` | session (a jump) | every group of the section, else of the view |
| `groups.keep { section }` | patch | the set shown into the groups' `collapsed` |
| `track.hide`, `track.show`, `track.move` on a group's address | patch | as on a track (the group's `hidden`, its place) |
| `layer.lanes { address, lanes }` | patch | a ranges layer's `lanes` |

The reveal opens groups through the view's existing reveal (`AddressModel`), not a command of its own.

## Defaults

The p53 protein section (`examples/src/p53`), written with three groups, the view at its fit:

```
 label column                  data
 p53                           MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWF…   (letters, no group)
 ▾ Features                                                                              (header line, 20 px)
     Domains                   [TAD ]   [PRD]        [ DNA-binding ............ ] [TET]
     Regions                   [ disordered ]                              [ disordered ]
     Motifs                         ▮              ▮▮                    ▮      ▮
 ▸ Variants  3 tracks          ▏ ▏▏  ▏   ▏▏▏ ▏▏  ▏▏▏▏▏▏ ▏▏▏ ▏▏▏ ▏▏▏▏  ▏ ▏  ▏▏  ▏ ▏  ▏      (one lane, every member's)
 ▸ Signals   3 tracks          ▁▂▃▅▇▇▆▅▃▂▁▁▂▃▅▆▇▇▇▇▇▇▇▆▅▄▃▂▂▃▄▅▆▆▅▄▃▂▁                      (conservation, 4 px)
                               ████████▓▓▓▓░░░░████████████████████▓▓▓▓░░░████          (pLDDT, 4 px)
                               ▁▁▁▁▂▂▃▃▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▂▃▅▆▇▇▆▅▃▂▁                (disorder, 4 px)
```

- **No setting**: every group open (`collapsed` false), `summary: auto`, labels indented 12 px, gaps 4 px within and
  8 px around a group. The header line's label in the label column's type (ux-foundations section 2), the count in
  `msui-text` at the quiet weight, the chevron 10 px in the label ink.
- **Folded rows**: the lane 22 px (a 12 px mark lane and padding, as a ranges lane), each strip 4 px with a 1 px gap,
  no labels; the summary's marks in their layers' own colours (a scheme's colour is the data, layouts.md "Colour by
  value": never tinted in a summary).
- **Empty**: a folded group whose tracks have no items in view draws the note row ("no items", 18 px); an open group
  whose tracks all collapse under `track.empty: collapse` keeps its header line with "no items".
- **Crowded** (100 tracks in 10 groups, all folded): ten rows of 22 to 38 px with their headers, about 400 px: inside
  `section.maxPx` (600), so the section does not scroll; open, its tracks scroll as layouts.md's section scroll says.
- **Widest view**: the shared lane at a genome's width is the ranges density form (stacked or majority, the layer's).
  **Closest**: items as marks, still no labels (a summary never labels; hover names).
- **Dark**: tokens only; the header line, chevron and count use the label column's tokens; nothing new.
- **Narrow label column**: the label cut with an ellipsis before the count, the count before the chevron; the
  chevron never dropped.

**Mock**: the p53 section above in the Storybook (`Seq* (next)/Design/Mocks/Track groups`), fake data: click a
chevron to fold and open (the lane and strips appearing in its place), Collapse all and Expand all from the section's
menu, a Variants group with an authored summary (a density values track) that shows only folded, "Lanes ▸ One lane"
on Domains, a hidden track in Signals restored from a hidden-tracks list, ↑ ↓ and Space at the header cursor, a
search hit inside a folded group opening it and Alt+Left folding it back. Whether it is built is the author's call.
**Precedent**: ProtVista's categories (the merged folded row).

## Configuration

| Item | Default | Why | Alternatives |
|---|---|---|---|
| a group's `collapsed` (spec) | false | a figure opens showing its tracks; the writer folds on purpose | true |
| a group's `summary` (spec) | `auto` | a folded group still shows where its data is | `none`, a track of the group (not built) |
| a ranges layer's `lanes` (spec) | `packed` | today's lanes and labels | `one` |
| `header.groupPx` (render config) | 20 | the section header's line (`header.sectionPx`) | any (not built) |
| `header.groupIndentPx` (render config) | 12 | membership read from the label column alone | 0 (not built) |
| `groupGap` (render config) | 8 | ux-foundations section 3 ("8px between groups") | any (not built) |
| `group.summary.stripPx` (render config) | 4 | the values layer's thin look reads from 3 px | 3 to 8 (not built) |
| `group.summary.maxStrips` (render config) | 4 | a folded row stays near one lane's height | any (not built) |
| the group header's click (bindings) | `groups.toggle` | the chevron and the label one target, as a tree view's row | `none` (the chevron only) (not built) |

## Rejected alternatives

1. **A group as a view over existing nodes** (by a track param or an attribute, IGV's grouping). No place for its
   label, its opening state or its summary track; a figure could not open folded; and every reader would recompute
   membership. Grouping by attribute returns later as a command that writes `group` nodes (one patch).
2. **A flat group**: a `group` node beside tracks and a `group: <id>` param on each track, keeping track addresses
   stable. Membership and order split in two places, a group's tracks could be apart (a rule needed for that), and
   the tree already says what holds what. The address change on a move is the price, and moves are rare.
3. **Folding as a patch** (one undo step a fold). Refused for clades by the author (grids-and-rows.md Decision 1),
   for the same reason here.
4. **A folded group as each member's own dense line** (UCSC's composite in dense). N rows for N tracks: it saves
   labels and lanes, not rows, and a 10-track group stays 10 rows.
5. **The section as the group** (today's answer in track-header-parts.md's precedents table). A section is an axis:
   each costs a ruler, a scroll and a label line, and groups on one axis would not share a window.
6. **Nested groups** (JBrowse's categories, UCSC's super-tracks). No example needs them; one level keeps the header
   cursor and the summary simple; additive later.
7. **Only an authored summary** (the previous implementation's: folded with no summary track shows nothing). A reader
   folding a group they did not write would lose sight of its data.

## Consequences

- **Amends**: building-blocks.md (View: a section's tracks may sit in groups); track-header-parts.md (the precedents
  row "collapses a group of tracks": built; the `on: 'group'` entry scope; the matrix rows below); ux-foundations.md
  section 8 (the matrix gains the rows below; "View state is not history" names the folded groups) and section 9;
  grids-and-rows.md (a rows track in a folded group shows its consensus row); layouts.md "Track sizing" (sharing reads
  folded tracks; `lanes: one` is one lane); changes.md (the session row: "the folded groups"); features.md rows 159,
  70, 141, 185 and 128 (moved by the lane, §features-move-with-merges) and the Storybook's Design pages (Features;
  a "Track groups" page).
- **Public API** (ask before changing, root AGENTS.md; each additive): the language (`group`, `track.parents`,
  `lanes`); `ResolvedTrack.group`, `ResolvedSection.trackGroups`; `LayerType.collapse`; the section `Viewport.groups`;
  `HiddenTrack.group` and its `kind`; the header entry scope `group`; `renderSvg`'s `groups` option; the commands in
  section 9; the hit target `group-header`.
- **Parts added**: the group header line and the folded row in the tracks layout's chrome (`tracks-headers.ts`,
  `headers-draw.ts`); the summary band in `tracks.ts`'s shown list (a band whose layers are its members' in summary
  form); the `lane`, `strip` and `row` forms in the ranges, variants, secondary, pairs, values and rows types; the
  `groups.*` commands in the `headers` behavior (or a `groups` behavior beside it, the lane's call); the walker
  `tracksOf`.
- **Edge cases**: EC-103 (groups collapse; a summary track shows only while collapsed) built, its test named;
  EC-244 (a reader lacking groups) recorded: the members drawn open, the summary track left out; EC-169 (pending is
  never empty) held in the summary; EC-146 (header hover) refused for groups as for tracks; EC-154 (the cursor) holds
  across folds (positions unmoved; the rows' consensus row a stop); EC-237 (a session as a spec: track order, hidden
  tracks) extended to groups' order and hidden groups; EC-140 (hide then show equals a fresh mount) extended to fold
  then open; EC-226 (lane order) unchanged by `lanes: one` (one lane has no order but its draw order, stated);
  EC-104, EC-105 (split rows) and EC-243 (collapse by an arbitrary field, with split) stay deferred with Ranges' split.
- **Performance representatives** (performance.md's form, each with its path and test):
  - 100 tracks in 10 groups (each 10 tracks: 6 ranges of p53-sized data, 4 values), at a 1,200 px width:
    `groups.collapseAll` relayout and draw under 16 ms, one toggle under 8 ms, open again under 16 ms
    (`core/src/render/_test/groups-time.test.ts`, heavy, `itWithin`); quick representative: 20 tracks in 4 groups,
    the same assertions on what is drawn without the bound.
  - A warm frame with groups folded allocates nothing in draw and hover (`core/src/render/_test/groups.quick.ts`, the
    `test:perf` tier, §performance-tested's hot loops).
  - A folded group lays no lanes: the lay count of a folded member's layers in open form is zero (quick lane,
    `groups.test.ts`).
- **Matrix cells** (ux-foundations section 8), each built when a test names its gesture:

  | | select | extend | add / remove | focus | back | zoom in / out | pan | hide | edit |
  |---|---|---|---|---|---|---|---|---|---|
  | **a group** (its header, its folded row) | no: a group is not an element; its items by their marks in the summary | — | — | no: focus an item; a reveal opens the group | Alt+Left after a bulk fold or a reveal | as a track | as a track | the `seq:hide` slot (`track.hide`, patch); shown again by the host's list (`track.show`) | fold, session state: click its header, Space at the header cursor (`groups.toggle`); `groups.collapseAll`, `groups.expandAll`; order (`track.move`); the opening set (`groups.keep`, one patch) |
  | **a summary mark** | as an element (its item) | as an element | as an element | as an element, opening the group | as an element | — | — | no: hide its track or group | no |
  | **a track or layer** (amended) | | | | | | | | | `layer.lanes` (patch) beside the row's edits |

- **The language lane** comes first and alone (the author's word, the version bump): the kind, `track.parents`,
  `lanes`, the check's two issues, the builder, the resolved view's `group` and `trackGroups`, the walker, the
  language description's docs. The render and viewer lane builds on it.

## Tensions

- **Session folds and saved figures** (settled, the clades' answer): a fold is not undone and not saved; `groups.keep`
  and a step's view are how a fold becomes the spec's. A reader who expects undo to reopen a group gets Alt+Left only
  after a bulk fold or a reveal, not after a single toggle. Standing, as for clades.
- **Addresses by nesting** (settled for nesting): a track moved into a group is a new address; bookmarks to items are
  unaffected (they name data, not tracks), session bands reset.
- **The summary is merged, so it hides overlaps across tracks** (settled): later over earlier and shorter over longer,
  "+n more here" on hover. A reader who needs every track's own line opens the group.
- **`group`, the word** (Provisional in Decision 1): one more meaning beside scroll, block, copy and text groups.
- **A ranges track's lanes as a patch, a group's fold as session** (Provisional in Decision 6): two ways to get one
  row, one undoable. The recommendation keeps the line where the clades drew it: the figure's look is the spec's,
  reading is the session's.
- **Values in a folded row are strips, capped at four** (Provisional): a group of many signal tracks shows four and
  counts the rest. Settled by a dataset with more (the Genomes' signal tracks): then the cap, or a merged envelope.
- **Pairs as end ticks** (Provisional): arcs say which ends join; ticks only where. Settled by the author at the mock.

## Self-check

- **Claim**: a group is the writer's structure (a node) and its fold is the reader's (session state), its folded row
  drawn from what it holds.
- **Negated, "a group is the reader's"**: groups made on the fly by attribute (IGV). It explains a session-only
  grouping and no node. Planned as a command that writes groups (Rejected 1), not built.
- **Negated, "a fold is the writer's"**: a fold as a patch. It explains `collapsed` as the opening state and
  `groups.keep`: built.
- **Inverses**: fold and open (`groups.toggle`); collapse all and expand all; a reveal opens and Alt+Left folds;
  hide and show (`track.hide`, `track.show` on a group); `lanes: one` and `packed`. Built.
- **Symmetric cases**: a track leaving its group (not built: the writer's, Decision 7); the cursor going into and out
  of a folded group's headers (Decision 4: one stop, a track header's cursor moved to its group's); a hover from the
  summary to another section and from another section into the summary (section 8). Built.
- **Plurals**: every group of a section (bulk), of the view (`section?` absent); several hidden groups in the host's
  list (one entry each). Built.
- **Carried**: the kind's name, `lanes` as patch or session, the strips' cap, pairs as ticks (Tensions).

## In the specs

| Decision | Where |
|---|---|
| The `group` kind, its params, `track.parents`; one level | core › Groups of tracks |
| `0.11.0`; a reader without groups draws them open | core › Versions |
| The `group.summary` and nesting checks | core › Groups of tracks |
| `ResolvedTrack.group`, `ResolvedSection.trackGroups`; the walker | core › Resolved view › Groups |
| The header line, the folded row, indent and gaps, the hit target | core › Render › Track groups |
| The summary band, its forms by kind, its height and hit | core › Render › Track groups |
| `LayerType.collapse` | core › Render › Track groups |
| `lanes: one` | core › Render › Track groups |
| Folded tracks share heights; the summary takes no part | design-only: not built yet (the sharing reads the bands shown: a fold can move another track's shared height) |
| Groups open in text and rings (`seq.render.groups.layout`) | core › Render › Track groups |
| `renderSvg`'s `groups` | core › Render › Track groups |
| The folded set in the section's viewport; `collapsed` the opening state; resets on its patch | viewer › ViewModel › Groups |
| A toggle a hand's move; bulk and reveal jumps | viewer › Behaviors › Groups |
| The `groups.*` commands and their kinds; `layer.lanes` | viewer › Behaviors › Groups |
| The header cursor's group stop, Space, Enter | viewer › Behaviors › Groups at the header cursor |
| The menu's `on: 'group'` entries | design-only: not built yet (the section's Collapse all and Expand all groups are: viewer › Hidden groups) |
| Hidden groups in the host's list | viewer › Hidden groups |
| The render config items | core › Render › Configuration |
| The matrix cells | design-only: ux-foundations.md's matrix is the instance's (§interaction-matrix) |
| Grouping by attribute, nesting, moving between groups | design-only: not built (Rejected 1, 6; Decision 7) |
