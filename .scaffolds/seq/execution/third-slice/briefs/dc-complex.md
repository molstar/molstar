# Lane DC: the Complex dataset

Follow [datasets.md](datasets.md) for [examples](../../../design/examples.md) › Complex (CDK2 with cyclin A; ATP as a
one-position sequence drawn as a chip, molecules.md decision 14; AlphaFill's ligands; the interface and contacts as an
optional `pairs` layer the builder adds, never a default; the views of its table). Also **the grid over a complex**
([grids-and-rows](../../../design/grids-and-rows.md), "A grid over a complex"; the author, 2026-11-19: "can we have a
PAE plot version for complex whose tracks are rendered next to each other?"): a grid section whose x and y follow a
section of blocks (the chains side by side), an AF-shaped N×N PAE over the chains in order (synthetic, as the data
decisions say) read through FS-A2's matrix read (`shape: matrix`), chain boundaries drawn, a hover naming both residues
with their chains, a dragged rect split at chain boundaries, `gridSummary(view, rect)` (mean, min, max, count; a library
function in core, the agreed surface) and a story's readout panel ("Aligned A 360–821, Scored A 265–831, mean 10.1 Å").
Build on what is merged: blocks and chips (LB), grids and pairs, highlights (HL), the zoom strip. Replace the current
Complex stand-ins, keeping the existing stories (Blocks, PAE and contacts) working against the real anchors.

**Model**: strong (opus). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B seq/dc-complex
<that commit>`. **Storybook port**: 6015.
