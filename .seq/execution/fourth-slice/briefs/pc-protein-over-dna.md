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
exon; a domain (residues 5–14, 1-based) through the map, cut at the intron and **joined by a connector** across it
(axis.md, "An item split by a path is joined by a connector", added by the author 2026-11-21: build it in the `ranges`
layer if it is not built, in the draw lists so the canvas and the SVG draw it alike; a hit on it is the item's; tests for
a split item, an item across a collapsed break (no connector) and an item in one run (none); **configurable** as the
library's other choices are (the author, 2026-11-21): the `ranges` layer's `connector` param (`line` default, `none`)
through the language (0.5.0's owner is the orchestrator: list the language change in the report if the bump needs
a hot file) and a `RenderConfig` preset for its width, ink and alpha, both rows in the core `_spec`'s Configuration
table (the config-keys test) and a test that `none` draws none); a phosphosite at residue 17 (`pS17`) as a mark
(the `variants` layer, or whichever built-in draws a point feature). Labels: DNA, Protein, Exons, Domain, Phosphosite.

Check against the data at the fit view (every residue over its three nucleotides; the intron with no protein letters;
the domain's two pieces end where its residues' codons end); hover a residue: it and its codon light (emphasis on both
tracks); the card names the residue and its codon positions. A headless twin test (the drawn ranges against the
map), failing if the codon alignment is off by one. Look in light and dark on port 6011.

**Checks and report**: as [fx5-author-look.md](fx5-author-look.md) lists them; a report under 200 words with a
`.v6/changelogs/` entry's text.
