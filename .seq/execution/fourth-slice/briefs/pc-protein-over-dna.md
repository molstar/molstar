# Lane PC: a protein over its coding DNA

**Base**: `packages-seq` at the orchestrator's head (`git reset --hard packages-seq` first). Branch `seq/pc`. Model:
fast. The fix lane's rules ([fx-fix-lane.md](fx-fix-lane.md)) on headers, tests and hot files hold.

**Task** (the author, 2026-11-21: "are we able to reproduce this view in packages/seq? in appropriate storybook and with
our updated styling?"): the old Seq\*'s Case 14 (`packages/sequences/core/_test/format/cases.ts`, `case14`) in
`packages/seq`, in the new styling, as `Features/Coordinates/Protein Over Its Coding DNA` (with the other one-axis
stories; `packages/seq/storybook/src/_stories/features/coordinates/coordinates.stories.ts`), its spec an example view
in `@molstar/seqx-examples` beside `p53/axis.ts`'s `oneAxis` (a small file of its own, e.g. `coordinates/coding-dna.ts`,
or where the coordinates views live), registered in `views.ts`.

The view, one gene section on the plus strand, in genomic order with the whole fragment shown (UTRs and the intron not
collapsed: `pieces: 'whole'`), top to bottom: the ruler; the DNA's letters; the protein's letters a codon wide (through
a `codon` map with the two coding pieces as exons, as `case14`'s `cds`: 12–46 and 76–105, 0-based half-open; the
protein translated from them, the stop excluded); the exons as ranges (`Exon 1` 0–46, `Exon 2` 76–120), coloured by
exon; a domain (residues 5–14, 1-based) through the map, cut at the intron (if the design draws a connector across a
cut item, it shows; if not, say so in the report — do not invent one); a phosphosite at residue 17 (`pS17`) as a mark
(the `variants` layer, or whichever built-in draws a point feature). Labels: DNA, Protein, Exons, Domain, Phosphosite.

Check against the data at the fit view (every residue over its three nucleotides; the intron with no protein letters;
the domain's two pieces end where its residues' codons end); hover a residue: it and its codon light (emphasis on both
tracks); the card names the residue and its codon positions. A headless twin test (the drawn ranges against the
map), failing if the codon alignment is off by one. Look in light and dark on port 6011.

**Checks and report**: as [fx5-author-look.md](fx5-author-look.md) lists them; a report under 200 words with a
`.v6/changelogs/` entry's text.
