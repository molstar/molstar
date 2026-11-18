# UX foundations

## Intent

The author, 2026-11-17: "i would like to have better visual fidelity and nicer UI in general. so we should write "ux
foundations" at start too, similar to the ones the packages/ui has .. but target to showing the sequence data, when to
squircles, etc. … basically i want this to also look really nice. and the default behavior should be good too ..
consistent logic for selections, toggling, etc."

`@molstar/ui`'s foundations (its `styles/_spec`: color, type, spacing, sizes, alignment, shape, motion) govern the
chrome: headers, menus, legends, panels. These foundations extend them to the data: letters, cells, marks, tracks, and
what the pointer and keyboard do to them. The chrome uses `--msui-*` as it is; the data's tokens (`--seq-*`) are
derived from them, so a theme switches both.

## Provisional

This is the first statement. The numbers (sizes, steps, durations) are settled by the foundation stories, which draw
each rule on real data beside its violation, as `@molstar/ui`'s do; the grammar in section 8 is settled by its
conformance tests.

## 0. Configurable, with rules that hold

The author, 2026-11-17: "for "ux foundations" and "Extensibility": should "layouting/spacing/etc" be
extensible/configurable?" Yes, by the rule that what can be configurable is ([extensibility.md](extensibility.md)):

- **Every number here is a setting**: the spacing unit and scale, row heights, gaps, the label column, insets, letter
  sizes, density thresholds, radii, motion durations, and every color through the data tokens (`--seq-*`, derived from
  `--msui-*`). They come in presets (compact, comfortable, figure for export, print) that a viewer configuration or a
  section's style picks and overrides field by field.
- **Layouts share one internal interface**: tracks, text, rings and dot plots take a section of the resolved view, the
  settings and a size, and give geometry, hit-testing and draw lists; the interface opens as an extension point when an
  application needs one ([extensibility.md](extensibility.md)).
- **The rules hold whatever the values**, checked by one set of tests over every layout and preset: positions on the configured unit's grid; one label column a view; columns
  aligned to the device pixel across sections that share an axis; density steps that never move content; the emphasis
  order; square continuous data and squircle discrete marks; the grammar's fixed rules (section 8).

## 1. Color

- **Data and chrome never share a color role.** The accent (`--msui-accent`) is for actions and what is picked:
  selection and the focus cursor. Data colors (palettes, residue schemes, values) never use it, and selection never
  uses a data color. Why: a picked residue must read as picked whatever scheme colors it.
- **Palettes are tokens per theme.** Every categorical palette, ramp and residue scheme has light and dark variants,
  checked for distinctness (neighbors in a palette stay apart for the common color-vision deficiencies) and for
  contrast against the track background.
- **Ink follows its fill.** A letter on a colored cell takes dark or light ink by the cell's luminance (the canvas-ink
  rule), at 4.5:1 or better. Why: a scheme that colors cells must not cost the letters.
- **Missing is not zero.** Missing, unobserved and unmapped each have one neutral treatment (a gray family derived from
  the theme), distinct from the lowest value of any ramp.
- **Conflict is one color**, the theme's warning hue, for a mismatch, a conflicting map element, a failed load.

## 2. Type

- **Letters** are set in one monospaced data face with tabular metrics, centered in their cells on the cap height, in
  sizes on a short scale tied to the cell (section 3). Two ink weights only: regular, and strong for a reference row or
  an emphasized letter.
- **Labels and numbers** in the chrome use `@molstar/ui`'s type scale; rulers and numbering use tabular figures and
  never overlap (a label that does not fit is skipped, not shrunk).
- **One letter, one size per row**: a row shows letters for all its positions or none (EC-116).

## 3. Grid

- **The cell is the unit of data.** A position is a cell whose width is set by the zoom; letters show when the cell is at
  least the smallest legible size, then cells only.
- **Rows are on the 4px grid**: a letters row 20px by default (16 compact, 24 comfortable), feature rows 8 to 16px,
  value tracks in multiples of 4. Track gaps 4px within a group, 8px between groups.
- **One label column**: track and row labels sit in a column on `@molstar/ui`'s inset and label width, so the data
  starts at one x for every track of a view, in both layouts (tracks and text).
- **Columns stay aligned across sections** that share an axis (scroll groups), to the device pixel.

## 4. Density

What a cell shows steps down as it narrows, at fixed thresholds, without moving anything:

1. letters (with their colored cells);
2. cells only, colored;
3. dots or ticks per position (sequence present, letter hidden);
4. bars per run (a feature, a gap run);
5. density (a heat strip: coverage, conservation, variant counts).

Each step keeps the previous step's positions and colors, so zooming reads as resolution changing, never as the
picture changing. The thresholds are shared by every layer type. Each step's unit is a config item (`densityUnits`); by default letters are
measured in CSS px (10), the other steps in device px (cells 3, dots 1, bars 1/64) (default; to confirm with the author).
Why: a letter's legibility is in CSS pixels (a size that needs 14 device px at DPR 1 needs only 7 at DPR 2), while what
bounds the draw lists is a mark a device pixel. The letters step never starts below the width the
smallest letter size needs, so cells never show without letters.

## 5. Shape: when to squircle

- **Continuous data is square.** Letter cells, alignment grids, matrices, dot plots, value bars: adjacency carries
  meaning (a run reads as one run), so no corner is rounded.
- **Discrete things are squircles.** A feature, a domain, a ligand chip, a bookmark, a label pill, a selection's
  outline, a hover card: rounded with the superellipse `@molstar/ui` uses for surfaces (`--msui-corner-shape`), its
  radius proportional to its height (a third), capped by the surface radius.
- **Small things are square.** Below 6px in either dimension a discrete mark drops its rounding: a corner that can't be
  seen is noise, and adjacent items must not leave gaps.
- **Cut ends are square.** An item cut by the view's edge, an intron collapse or a block boundary has square ends on the
  cut side and rounded ends where it really ends (the frame rule, EC-100).
- **Outlines are concentric.** A selection or hover outline around a squircle mark takes the mark's radius plus the
  outline's offset, so the gap is even.
- **One path, every backend.** Squircles are drawn as a path in the draw lists, never left to a backend's rounded
  rectangle, so Canvas2D and SVG draw the same curve.

## 6. Marks and emphasis

- **Layers of emphasis, in a fixed order** from below to above: data; highlight channels (named, each a tint of its
  color under the data); selection (an accent tint under the data and an accent outline over it); hover (an outline at
  3:1 or better against its surroundings, EC-138); the keyboard cursor (an accent outline, thicker). A higher one never
  hides a lower one entirely. **Current pass order** (`selection` after `data`): the selection's tint is drawn
  translucent over the data, with its outline, not under it (default; to confirm with the author: the alternative is a
  separate under-data pass). Why over: a pass under the data is hidden by opaque cells.
- **The look as built** (the author's look at the selection, 2026-11-18, lane P): the selection is a soft tint (0.18)
  with one thin outline (1 px, 0.6), one band a section across its tracks and gaps (`emphasis.selectionShape`:
  `column`, or a box a track, `tracks`); the hover's outline is concentric, offset outward (1 px) and drawn above the
  selection; the keyboard cursor is a 1 px neutral line (`msui-text`), shown only after a key, not the thicker accent
  outline the list above names. A highlight channel is a tint (0.3) of its color; what is projected onto a linked
  section takes its own look at half strength (0.5).
- **Hover is the same everywhere**: one outline style for a cell, a range, a row, a pair; projected hovers on linked
  sections take the same style, lighter.
- **Selection is one look** in tracks, text, rings and dot plots.
- **Boundaries** (an insertion point, an indel, a zero-length item) draw as a hairline with a pinned caret, never as a
  cell (EC-127).

## 7. Motion

- **Navigation that jumps animates**: focus, reveal, a snapshot's view, zoom by a step (base duration, eased in and
  out); the eye follows the content to its new place.
- **Navigation that follows the hand does not**: drag to pan, wheel and pinch zoom track the input directly.
- **Data does not animate**: arriving regions draw in place (no fade on every tile); a density step is instant.
- **Reduced motion**: every animation is a cut.

## 8. Interaction grammar

One grammar for every layout, layer type and target. A gesture acts on **the element under the pointer at the
granularity of what was hit**: a cell (a residue, a column of an alignment), an item (a feature, a pair), a row (by its
header), a column (by the ruler), a track (by its header).

**The grammar is configurable** (the author, 2026-11-17: "the controls grammar should be configurable"). It is data: a
table of bindings from a gesture (with its modifiers and what it hits) to an action, and every action is a command
(`@molstar/ui/commands`). The table below is the default. An application composing the viewer replaces bindings, adds
its own, or picks a preset (for example one where a drag on data pans, as genome browsers do); a user's changes are a
settings layer over it. The rules that stay fixed whatever the bindings:

- **One configuration per viewer, applied to every view in it**: a gesture means the same thing in tracks, text, rings
  and dot plots of one viewer. Why: consistency is the point of a grammar; configurability must not bring back
  per-view behavior.
- **Granularity comes from the hit, not the binding**: a binding says "select", the hit says a codon or an item.
- **Modifiers are platform-aware** (Cmd on macOS, Ctrl elsewhere), never browser zoom keys or Tab (EC-155).
- **Every action stays reachable** by a command and from the keyboard, whatever the pointer bindings.
- **Conflicts are reported**: two bindings for one gesture is an issue when the configuration is built, not a race at
  run time.

| Gesture | Effect |
|---|---|
| hover | the element is hovered here and on every linked section; nothing else changes |
| click | the selection becomes the element; clicking the only selected element clears it |
| Shift+click, Shift+drag | extends the range the anchor is in to the element, keeping the other ranges (a range in one system; a run of rows) |
| Cmd/Ctrl+click | toggles the element in or out of the selection, keeping the rest; the anchor moves there |
| Cmd/Ctrl+drag | adds the range to the selection, keeping the rest |
| Alt/Option+drag | removes the range from the selection |
| click on empty space, Esc | clears the selection (Esc first cancels a gesture in progress) |
| drag on data | selects a range, snapped to the granularity hit (whole codons over nucleotides, whole items over items) |
| drag on an axis, navigator or with Space held | pans |
| wheel, two-finger scroll | pans; Cmd/Ctrl+wheel or pinch zooms about the pointer |
| double-click | focuses the element (zooms to it); on the element already focused, goes back; on empty space, fits the view |
| Shift+double-click, -, =, 0 | zooms out about the pointer; zooms out and in about the cursor or the centre; fits |
| Alt+Left, Alt+Right, the mouse's back button | back and forward through the view's jumps (focus, fit, reveal, a bookmark or link): session state, not undo |
| arrows | move the keyboard cursor; Shift+arrows extend the anchor's range; Cmd/Ctrl+Shift+arrows start a new range; Enter focuses; Space toggles |
| drag, wheel or pinch on a follower | moves its leader through the follow, so the follower moves as asked (`follow.followerMoves`: `lead`) |

- **Selection is a location**, held per system and projected: what is selected in a protein shows on its gene, its
  structure and its alignment rows.
- **Toggles are uniform**: hiding a track, collapsing a group, pinning a row, showing a layer: one gesture, one state
  shown in the chrome (an icon, never only a color), one undo step when it changes the spec.
- **View state is not history**: hover, selection, the cursor, the viewport are transient; edits to the spec (a
  bookmark, a hidden track, a theme) are patches, one undo step each, with a preview while a gesture is in progress that
  ends in its patch or a cancel.
- **Every gesture has a keyboard way and a command**, and the context menu lists what applies to the element under the
  pointer, in the same words as the commands.

### The interaction matrix

What a user acts on, against the verbs (the principles, "Interaction is a matrix", from the second slice's fold-back):
each cell is the gesture or key and the command, **planned** with the slice that builds it, or **no** with why. A slice
plan names the cells it fills; its gate checks them; a design that adds an action fills its row and column, and its
self-check's negation (an inverse, a symmetric case, a plural) names the cells to look at. As built after the second
slice:

| | select | extend | add / remove | focus | back | zoom in / out | pan | hide | edit |
|---|---|---|---|---|---|---|---|---|---|
| **an element** (a cell, an item) | click; `selection.set` | Shift+click (the anchor's range) | Cmd/Ctrl+click toggles | double-click; `address.focus` | double-click on it when focused | — | — | no: hide a track or layer | no: data is not edited in the viewer |
| **a range** (positions) | drag on data | Shift+drag | Cmd/Ctrl+drag adds, Alt+drag removes; `selection.add`, `selection.subtract` | Enter at the cursor | Alt+Left | — | — | — | — |
| **the selection** | — | Shift+arrows | Cmd/Ctrl+Shift+arrows start a range | Enter | Alt+Left | — | — | — | a bookmark from it: planned (a command, third slice) |
| **a leading section** | — | — | — | — | — | Cmd/Ctrl+wheel, pinch, -, =, Shift+double-click | wheel, drag on the axis, Space+drag, Shift+Alt+arrows | no: hide its tracks | its params by patch (height) |
| **a follower section** | as a leading one | as a leading one | as a leading one | as a leading one | as a leading one | moves its leader (`lead`) | moves its leader (`lead`) | as a leading one | as a leading one |
| **the view** | — | — | — | 0 fits; double-click on empty space | Alt+Left / Alt+Right, the mouse's buttons; `view.back`, `view.forward` | — | — | — | undo and redo through the host's stack |
| **a track or layer** | planned (a row selection, third slice) | — | — | planned | — | — | — | `track.hide`, `layer.hidden` | height, color, order by patch; previews while dragged |
| **a navigator** (strip, miniature, ideogram, vertical overview: planned, third slice) | — | — | — | a click centers the view there | Alt+Left | drag the window's handles | drag the window | the host hides the part | — |
| **a bookmark or link** | — | — | `refs.add`; removing one: planned | runs its action | Alt+Left | — | — | — | its label: planned |

Planned beyond this slice, from the second slice's use gate: a section's direction (reversed, following its leader),
another system's letters on a section (one axis for protein and gene), collapsed introns; the history's toolbar.

## 9. Chrome

- Built from `@molstar/ui` components and tokens only: track headers, legends, the hover card, menus, toolbars,
  navigators' frames. Section headers are `@molstar/ui` headers (32px); track header rows are rows (24px).
- A track header names its data and shows its state (hidden layers, filters, a pending load) as icons in fixed slots.
- Legends sit with their track or in one legend panel, never over data.
- Loading regions show a placeholder in place (`@molstar/ui`'s rule), never as empty data (EC-169).
