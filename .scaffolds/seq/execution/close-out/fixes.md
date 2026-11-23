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
