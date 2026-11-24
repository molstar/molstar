# The novice pass (close-out)

**Goal** (`loop.md`, Gates: after the last lane that draws; the author, 2026-11-24: "we can run the novice pass"): light
agents play a first-time user of Seq\* (next) in a visible browser, one task script a view, split by groups from the
start, and report the friction they meet, in the shape of the fifth slice's pass (`../fifth-slice/novice-findings.md`).
They change no code.

**Where**: the Storybook of `packages-seq` on port 6006 (already running; never start or stop a server). The built-in
browser's tools (`mcp__Claude_Browser__*`): **open a tab of your own** (`tabs_create`) and work only in it; navigate to
a story's iframe (`http://localhost:6006/iframe.html?id=<story id>&viewMode=story`); screenshot, hover, click, drag,
type keys, read the page. The story ids are in `/index.json` (those starting `seq-next-`). Light by default; for dark,
`&globals=theme:dark` on the URL. If `tabs_context` says the pane is hidden, stop and say so: a hidden pane is never
evidence. Close your tab at the end.

**Who you are**: a biologist who knows what a chromosome, an assembly, a gene and its transcripts, a protein's domains
and secondary structure, a GWAS and a ligand are, and has never seen this viewer. Read only what a story shows on its
page (its description, legend, buttons, menus, the drawing). Do not read the code, the specs or the designs: if a task
needs them, that is a finding.

**Pointer paths**: for every menu and submenu you meet, write the path you took to it (for example: right-click a
track's label › "Layers" › …), whether you found it by hovering, right-clicking, a key or the page's text, and whether
Shift+F10 or the context-menu key opens it from the keyboard.

**The groups** (one agent each, in parallel, each in its own tab):

- **A. Legends, layers and groups**: Features › Structure link › Focus And Back; Examples › Complex › Views (Ligands,
  As A Complex); Features › Legends; Features › Track groups (Opening, All Open, All Folded); Features › Layers. Tasks:
  hide one colour's layer from the legend and bring it back; hide every layer of a merged row at once and show them
  again; tell from the track's label that something is hidden; fold and open a group of tracks, read what a folded
  group shows; undo each step.
- **B. Finding things and reading values**: Examples › p53 › Search, Examples › Plasmid › Search; Features › Genome
  formats › GFF3 and GTF › Find By Name; Examples › Genomes › Views › Manhattan; Features › Tooltips › By Data. Tasks:
  find a residue, a motif or a gene by typing; step to the next and the previous match; go back to where you were;
  move to the next feature with the keyboard; read which variants pass the genome-wide line and what the dashed and
  dotted lines mean; read every transcript of a gene without scrolling inside the track.
- **C. The genome and the ring**: Features › Genome Ring (The Genome At The Fit, The Stem Menu, The Partner Keys, The
  Square Beside, The Square Clicked, Kinds Of Links); Features › Across sections (One Side Zoomed In, Two Genomes At The
  Fit). Tasks: open a chromosome on the ring and go back to the whole genome; open a pair from a stem's menu, by mouse
  and by keyboard; step through a stem's partners with the keyboard; pick a pair on the square and go back to all;
  zoom one side of a bridge and say whether the ribbons stay joined.
- **Everywhere you pass**: does hover say what you are on; does a selection look selected; can you undo; does the
  keyboard work once you click into the view; do menus look alike across the viewer; is anything clipped, overlapping,
  unreadable in dark, or drawn twice; do canvas and SVG twins match where a story shows both.

**Report** (each agent under 800 words): a findings table, one row a finding: the story id, what you tried, what you
expected, what happened, the kind (**bug**: it does something wrong; **missing design**: the page gives no way or no
hint; **look**: it works but reads badly), and a severity (blocks the task / slows it / cosmetic). Then the pointer
paths (one line a menu), the tasks you could not finish, and the three findings that would help a new user most. No
more than four screenshots (the worst cases), saved in your scratchpad, their full paths in the report.

**Output**: the orchestrator writes `novice-findings.md` in this folder from the reports, numbered across groups,
counted by kind, each finding reproduced in a visible pane before it is placed (`fixes.md`, a lane,
`packages/seq/TODO.md`, a drift note).

**Model**: Sonnet. Read-only: no edits, no commits, no servers, nothing outside the browser and your scratchpad.
