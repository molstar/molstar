# The novice passes' 14 design questions: draft answers (for the author)

From the fifth slice's [novice-findings.md](../../fifth-slice/novice-findings.md). There are no mocks (`plan.md`,
Re-plans, 2026-11-23). Groups A to C are interactions: their gestures are agreed here before a lane. Groups D to G are
a sentence each. Items 8 and 26 are one question.

Each answer names where it lands: a design, and the package spec section that changes with it. **Check** marks what
the draft assumed from the code and a lane confirms first.

## A. A layer's legend: hiding it, and seeing that it is hidden (5, 13)

What exists:
- `layer.hide` and `layer.show` are commands (patches: one undo step each).
- The track header's menu has "Hidden layers" with a "Show <layer>" entry for each.
- A click on a legend entry already **selects** that entry's items (`legend.select`, viewer › Legends). So hiding
  cannot be the entry click: one click with two meanings is the drift this answer must avoid.

Proposed:
- A legend row for a layer gets an **eye button** after its title: an `@molstar/ui` icon button labelled
  "Hide Contacts". A legend merged over several layers gets one button that opens a menu of its layers.
- **Hidden, the row stays.** The title and swatches are dimmed, the button turns to eye-off, "Show Contacts", and its
  entries have no action (no hover light, no select). This also answers 13: a hidden layer shows where its key was.
- Look 24 joins this lane: the Contacts entry gets its swatch, and its label stops doubling ("Contacts: by chromosome").

| Gesture | What it does |
|---|---|
| Click the eye | `layer.hide`; the row dims; the button reads "Show Contacts" |
| Click the eye-off | `layer.show`; the row back as it was |
| Tab to the button, Enter or Space | as a click |
| Click an entry of a shown layer | unchanged: selects its items (the click rule) |
| Click an entry of a hidden layer | nothing; the pointer says so (not-allowed) |
| Undo | the hide or show undone, one step |
| Hover the eye | a tooltip with the button's name |

Lands in: `design/themes.md` (Legends); viewer spec › Legends. That spec gains "a layer's row keeps its place when
hidden". Its "Hiding has two owners" rule is about hiding a **legend**, which this does not do: hiding a layer is the
spec's existing `layer.hide`.

**Question**: should a host be able to turn the eye off (a `mountLegend` option, `layers: false`)? Recommended: yes,
on by default. Why: a figure's legend may be printed, and the button is chrome there.

## B. The stem menu: partners and its route (8, 26, 28)

What exists:
- The menu lists one row per ribbon (`partnersOf`: a partner stem, its count and strand). So a chromosome with several
  stem pieces, or links on both strands, gets several bare rows ("chr9 · 375 links", "chr9 · 65 links").
- A chromosome's links to itself list it as its own partner.
- The menu opens only pre-opened in its story.

Proposed:
- **The route**: a right-click on a stem opens the stem menu at the pointer. Over the ring's empty space, the
  browser's own menu still opens. Shift+F10 and the Menu key on the focused stem open it too; this also checks
  unverified item 29.
- **Partners by chromosome**: one row per partner chromosome, its total links, by count. A chromosome reached by more
  than one ribbon gets a submenu with a row per piece, named by its band range and strand ("q34.1 · 375 links",
  "q21–q22 · 24 links, reversed").
- **Links within**: the stem's own chromosome reads "Within chr22 · n links", listed last, never as a partner.
- **Only the stem's actions**: the stem menu holds what acts on the stem. Move up, Reverse and Hide track are the
  track's and stay in its header menu.

| Gesture | What it does |
|---|---|
| Right-click a stem | the stem menu at the pointer |
| Right-click the ring's empty space | the browser's menu |
| Focus a stem, Shift+F10 or the Menu key | the stem menu at the stem |
| Hover or Right arrow on a chromosome with pieces | its submenu |
| Click a partner (or a piece) | shows the pair, one history step (as now) |
| Escape | closes the submenu, then the menu |

Lands in: `design/ring-extension.md` (the stem menu); ext-rondo spec › the viewer.

**Check**: what a "stem piece" is in `RondoLay` (the stems a chromosome splits into). The band range names it only if
each piece is one contiguous run.

## C. The square's way home (6)

What exists:
- A drag moves the square's window, and a drag on an edge or corner resizes it, snapped to whole blocks.
- A double-click on an empty spot moved the window to chr12–chr16.
- 0 (fit) did nothing.

Proposed:
- **A hint line** under the square, muted, one line: "Drag to move · drag an edge to resize · double-click for the
  whole genome". It is part of the square's part, and a host turns it off with an option (`hint: false`).
- **Double-click anywhere on the square** gives the whole genome, replacing the move. **0 with the square focused**
  does the same. Each is one history step (as FX24 made the square's moves).
- **No button**. Recommended because the hint names the gesture and a button is chrome the host may not want.
  Alternative: a small "Whole genome" button in the square's corner.

| Gesture | What it does |
|---|---|
| Drag inside the window | moves it (as now) |
| Drag an edge or corner | resizes it (as now) |
| Double-click anywhere | the whole genome, one history step |
| Focus, 0 | the whole genome, one history step |
| Back (⌥←) | the window before |

Lands in: `design/genome-coordinates.md` (the square); viewer spec › Input and the square's part.

## D. Cards (21, 25, 33)

- **21, the stem card**: the first row names what it counts ("This stem: 24 links to chr9"). The per-layer totals sit
  under a heading of their own, "chr22, all links" (Contacts 638, Translocations 24).
- **25, a pair's ribbon card**: both ends, each with its position and band ("chr4 q21.3 · 85.1–85.6 Mb", then
  chr7's). The bins row reads "n bins: zoom in to see them".
- **33, a Hi-C card**: in bases ("chr1 88.5–89.0 Mb × chr1 86.5–87.0 Mb") with its value; the bin indices go. The
  square's rulers read in kb or Mb.

Lands in: `design/tooltips.md`; ext-rondo spec (the cards), viewer spec › Hover card.

## E. A legend and a key (7, 9)

- **7, synteny's legend**: the bridges' colour by orientation (same, inverted) as legend entries, and a note that a
  bridge's width is its block's length. An inversion's colour differs from a same-strand bridge's, so chr3's inversion
  and a translocation no longer read alike. **Check** whether the bridge layer's colour is a channel today or fixed in
  the draw; if fixed, it becomes a channel, so the legend derives as every other one does.
- **9, the next feature**: `]` and `[` move the cursor to the next and previous feature on the focused track, skipping
  empty cells. They are listed in the keys panel. Alternative: Alt+Right and Alt+Left. **Check** that `]` and `[` are
  unbound in the keymap and the host's.

Lands in: 7 in `design/across-sections.md` (bridges) and the viewer spec › Bridges; 9 in `design/ux-foundations.md`
(the keys) and the viewer spec › Bindings.

## F. What tracks name (10, 11)

- **10, a canonical transcript**: marked when the data says so (GFF3 or GTF `tag=Ensembl_canonical` or `MANE_Select`).
  It is drawn first in its gene's stack, its label bold, and its card says "Canonical" (or "MANE Select"). With no tag,
  nothing is marked: the viewer never guesses. **Question**: should a collapsed gene draw only its canonical
  transcript? Recommended: not in version 1. Why: it is a new display mode; record it in the TODO.
- **11, the VCF's sample names**: the sample rows get a names gutter from the VCF header, the same part the alignments'
  rows use. The card reads "Sample S11 · 0/1 (heterozygous)", not "State S11 / Value 0/1". The methylation matrix's
  "State" (look 34) is fixed by the same change.

Lands in: `design/genome-formats.md`; core spec (the GFF3 and VCF readers' fields), viewer spec › Hover card.

## G. Empty states (12, 35)

- **12, "zoom in to read"**: the hatch names the gesture ("Zoom in to read: ⌘/Ctrl+scroll, or +"). A track that draws
  nothing at this zoom gets the same hatch, never an empty band. **Check** + and = first: the novice found they did
  nothing, perhaps the tool's keys.
- **35, links pending**: while the links load, the ring's hole says "Loading links…" (muted, centred), and the square
  shows the hatch with "Loading…". It is announced once, politely. On a failure, the issue's text replaces it.

Lands in: `design/region-reads.md` (the hatch) and `design/ring-extension.md` (the pending ring); viewer spec and ext-rondo spec (their pending states).

## Lanes this suggests

At most three run at once, never two in one file.

| Lane | Items | Model |
|---|---|---|
| Legend and layers | A (5, 13, look 24), E7 | strong: it changes the legend's rule |
| Ring | B (8, 26, 28, 29), D21, D25, G35's ring part | strong: the menu's structure |
| Square and keys | C (6), D33, E9, G12, G35's square part | Sonnet |
| Formats | F (10, 11, look 34) | Sonnet |
