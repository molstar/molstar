# Lane DT: the tRNA dataset

Follow [datasets.md](datasets.md) for [examples](../../../design/examples.md) › tRNA (PDB 1EHZ chain A with its
modified residues and pairs with kind and tertiary flag, in the `seq` format; Rfam RF00005's seed in Stockholm with
`SS_cons`, as published; the views: arcs nested and crossing, the base-pair matrix as a dot plot with crosshair hover
and rectangle selection, the covariance alignment with pairs). The modified residues are drawn as the section and
[molecules](../../../design/molecules.md) allow this slice (the parent letter; the `mark` is the fourth slice's: say so
if a view needs it). The Rfam seed's row count and size are "to verify when recording": record what it is.

**Model**: fast (sonnet). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B seq/dt-trna <that
commit>`. **Storybook port**: 6013.
