# The novice pass (stage 6)

**Goal** (plan, "Use gates"; `loop.md`, Gates): a light agent plays a first-time user of Seq* (next) in a visible
browser, one task script a view, and reports the friction it meets. It changes no code.

**Where**: the Storybook of `packages-seq` at http://localhost:6008 (already running; never start or stop a server).
Use the built-in browser's tools (`mcp__Claude_Browser__*`): navigate to a story's iframe
(`/iframe.html?id=<story id>&viewMode=story`), screenshot, hover, click, drag, type keys, read the page. The story ids
are in `/index.json` (those starting `seq-next-`). Light by default; for dark, `resize_window` with `colorScheme:
'dark'` and `&globals=theme:dark` on the URL; set it back to light at the end.

**Who you are**: a biologist who knows what a domain, an alignment, a plasmid map, a base pair, a PAE plot and Kabat
numbering are, and has never seen this viewer. Read only what a story shows on its page (its description, legend,
buttons, the drawing). Do not read the code, the specs or the designs: if a task needs them, that is a finding.

**The tasks**, one script a view: for each, try it with the page alone, note what you expected, what happened, and how
long it took you to find the gesture.

- **p53** (Examples › p53): find the DNA-binding domain and select it; zoom to residue 175 and read its letter; on the
  protein on its gene, find which exon holds residue 175; on the gene by region, move to the next region.
- **Globins**: in Family, find the most conserved column from the logo; fold a clade from the tree and unfold it; in
  Hemoglobin, find which residues contact the heme; split the copies and merge them again; in Lazy or Scale, move to
  row 5,000 with the keyboard.
- **Plasmid**: on the ring, find bla and read its strand; select the feature that crosses the origin and focus it; move
  the ring's window from the linear detail; find EcoRI's sites.
- **tRNA**: find the anticodon stem's pairs in the arcs; on the dot plot, hover a pair and select a rectangle of two
  ranges; tell a crossing pair from a nested one.
- **Complex**: find the interface residues of CDK2; on the PAE, hover a cell and read both residues and the error;
  drag a rectangle and read its mean; find which residues ATP contacts; hide the ligands past 2 Å and undo it.
- **Antibody**: read residue 100A's IMGT label; switch the rulers to IMGT and back; select CDR-H3 by Kabat and bookmark
  it, rename the bookmark, undo.
- **Features**: on Query, type an address that selects residues 100 to 120 of p53; on Themes, change a scheme and undo;
  on Navigators, drag the zoom strip's window, then release the mouse outside the canvas; on Highlights, tell two
  overlapping highlights apart; on Layers › Tracks in their wells, find the padding; on Files, find where a file's
  rows come from.
- **Everywhere you pass**: does hover say what you are on; does a selection look selected; can you undo; does the
  keyboard work once you click into the view; is anything clipped, overlapping, unreadable in dark, or drawn twice.

**Report** (under 1,200 words): a findings table, one row a finding: the story id, what you tried, what you expected,
what happened, the kind (**bug**: it does something wrong; **missing design**: the page gives no way or no hint;
**look**: it works but reads badly), and a severity (blocks the task / slows it / cosmetic). Then the tasks you could
not finish, and the three findings that would help a new user most. Attach no more than five screenshots (the worst
cases), saved in your scratchpad, their full paths in the report.

**Model**: fast (sonnet). Read-only: no edits, no commits, no servers, nothing outside the browser and your scratchpad.
