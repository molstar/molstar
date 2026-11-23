# The novice pass (stage 5)

**Goal** (plan, stage 5; `loop.md`, Gates): a light agent plays a first-time user of Seq\* (next) in a visible browser,
one task script a view, and reports the friction it meets, in the shape of the third and fourth slices' passes
(`../../third-slice/briefs/novice-pass.md`). It changes no code.

**Where**: the Storybook of `packages-seq` at the port the orchestrator gives (already running; never start or stop a
server). Use the built-in browser's tools (`mcp__Claude_Browser__*`): navigate to a story's iframe
(`/iframe.html?id=<story id>&viewMode=story`), screenshot, hover, click, drag, type keys, read the page. The story ids are
in `/index.json` (those starting `seq-next-`). Light by default; for dark, `resize_window` with `colorScheme: 'dark'` and
`&globals=theme:dark` on the URL; set it back to light at the end. If `tabs_context` says the pane is hidden, stop and
say so: a hidden pane is never evidence.

**Who you are**: a biologist who knows what a chromosome, an assembly, a gene and its transcripts, a BED or VCF file, a
Hi-C contact map, synteny and a translocation are, and has never seen this viewer. Read only what a story shows on its
page (its description, legend, buttons, menus, the drawing). Do not read the code, the specs or the designs: if a task
needs them, that is a finding.

**Pointer paths** (B10): for every menu and submenu you meet, write the path you took to it (for example: right-click a
stem › "Open the pair" › …), whether you found it by hovering, right-clicking, a key or the page's text, and whether
Shift+F10 or the context-menu key opens it from the keyboard.

**The tasks**, one script a view: for each, try it with the page alone, note what you expected, what happened, and how
long it took you to find the gesture.

- **Find a gene by name** (Examples › Genomes › Views, and Features › Genome formats › GFF3 and GTF): find TP53 (or the
  name the page suggests) by typing it; read which chromosome and coordinates it sits at; open its transcripts and tell
  them apart; find its exons; move to the next gene with the keyboard.
- **Compare two assemblies across a bridge** (Views › Synteny, Liftover, Dot plot; Genome Ring › Two assemblies on one
  ring): find a block that is inverted between the assemblies; follow one gene from one assembly to its place on the
  other; tell what a bridge's colour and width mean from the legend.
- **Follow a contact on the ring** (Features › Genome ring: at the fit, One chromosome, A pair, Kinds of links and its
  variants, Hovering a stem, The stem menu): hover a link and read both its ends; tell a contact from a translocation
  and find the made-up chr9–chr22 fusion and its gene names; hide the contacts from the legend and bring them back;
  open a pair from a stem's menu; expand a folded stub; read the chromosome labels at the top of the ring (do any
  overlap?).
- **Switch chromosomes on the square** (The square beside, dragged, resized, clicked, scattered): move the window to
  another chromosome by dragging it; resize it from an edge and from a corner; pick a scattered set (chr4 and chr9 only)
  and read which cells say how many links; go back to the whole genome; undo.
- **Region reads** (Features › Region reads: Hi-C, matrix, pairs, zoom in): zoom from the whole genome to a region and
  read what the view says while data is loading; hover a contact cell and read both bins and the value.
- **Formats** (BED, GenBank, VCF): find a variant's genotype per sample; tell what a missing value looks like.
- **Everywhere you pass**: does hover say what you are on; does a selection look selected; can you undo; does the
  keyboard work once you click into the view; can a legend's "+n" be expanded and collapsed; is anything clipped,
  overlapping, unreadable in dark, or drawn twice.

**Report** (under 1,200 words): a findings table, one row a finding: the story id, what you tried, what you expected,
what happened, the kind (**bug**: it does something wrong; **missing design**: the page gives no way or no hint;
**look**: it works but reads badly), and a severity (blocks the task / slows it / cosmetic). Then the pointer paths
(one line a menu), the tasks you could not finish, and the three findings that would help a new user most. Attach no
more than five screenshots (the worst cases), saved in your scratchpad, their full paths in the report.

**Output**: the orchestrator writes `novice-findings.md` in this folder from the report, counted by kind, each finding
placed (`fixes.md`, a lane, `packages/seq/TODO.md`, a drift note).

**Model**: fast (sonnet). Read-only: no edits, no commits, no servers, nothing outside the browser and your scratchpad.

## Second short pass (2026-11-23)

The author, after the first pass's partial coverage: "Ok run the other short novice pass". It runs on the same brief,
on FX24's Storybook (the fixed cards and the square's history), over only what the first pass did not reach: Hovering a
stem; A pair; the scattered pick (chr4 and chr9 only) and the square's Back; the made-up chr9–chr22 fusion on Kinds of
links (its card, its genes); open a pair from a stem's menu (right-click a stem, and Shift+F10); expand a folded stub;
Views › Dot plot and Liftover; a Hi-C cell (Region reads › Hi-C, Matrix); GenBank and BED; the loading state (Links
pending); and dark on each. Report under 800 words, with the first pass's table shape and its numbering continued from 21.
