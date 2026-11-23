# Fixes: close-out

Open items for the next fix point or a running lane's fix round (the loop, Each step 7).

- **Fixed by CL2** (`a0397521c`). **Look 24, the Contacts legend entry** (the fifth slice's novice pass; CL1's brief, item 4, not landed: its code is in
  `ext-rondo/src/card.ts`'s `legendRondo`, CL2's files). The entry has no swatch, and its label doubles ("Contacts: by
  chromosome"). Fix: drop `<name>: ` when it repeats the legend's title, and give "by chromosome" a hue swatch. Kind:
  look. To CL2's fix round.
- **TP53-201 is first only in the SVG** (the orchestrator's look at CL4, Features › Genome formats › GFF3 › Tp53): the
  canonical transcript is bold in both copies, but on the canvas it is the third lane (under TP53-203 and TP53-227) and
  first in the SVG. Answer F10 says first in its gene's stack. Likely the TODO's "the canvas and the SVG of one spec
  lay out differently". Kind: bug. To a fix lane, with a twin test comparing the two copies' lane order.
- **A variant's card title is its raw id** ("20001:601:f95ba8e6>C"; the orchestrator's look, Examples › Genomes ›
  Views › Samples): before CL4, not its change. The title names the variant (its kind and position, or its ID). Kind:
  look. To a fix lane.
- **Stub labels drawn over each other** (the author, 2026-11-24, Features › Genome Ring › Kinds Of Links, the ring
  zoomed to chr2 from the square: "also a visual bug", a screenshot): two stubs at one angle draw their labels on the
  same spot ("chr11 · 335" and "chr1 · 3…" overprinted into one smear), and their wedges overlap. A stub's label needs
  the ring's label collision rule (stagger or drop the weaker with a "+n"), as the chromosome labels have. Kind: bug.
  Next fix point (a fix lane, with a test on the zoomed ring's label boxes not overlapping).
- **F9's pressed slot is not drawn on canvas headers** (the orchestrator's look at CL10, Features › Structure link ›
  Color By This Track; CL10's review): `pressed` reaches the DOM header parts (`aria-pressed`) but not the canvas
  header's slots (`headers-draw`), so the story the decision came from shows no pressed icon. Kind: bug. To a fix lane:
  the canvas slot drawn pressed at rest, one at a time, a test on the canvas draw list.
- **The pressed state refreshes only on a header run** (CL10's report, in the viewer spec): a host whose colouring
  changes outside a command is not seen until the next layout, hover change or header run. Kind: followup (TODO at
  the close).
