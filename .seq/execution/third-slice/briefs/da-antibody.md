# Lane DA: the Antibody dataset

Follow [datasets.md](datasets.md) for [examples](../../../design/examples.md) › Antibody (a Fab's heavy and light
chains numbered Kabat and IMGT, insertion codes, IMGT out of index order, the CDRs; a designed variant with indels
against the wild type; the views of its table: numbering schemes side by side, Kabat vs IMGT labels (EC-2 to 8), design
vs wild type with indel carets (EC-29), CDR selection and editing as bookmarks with undo (EC-162, 164)). The fourth
slice's Fc glycan (molecules.md) is not this lane's. Build on what is merged: the axis and numberings, rulers, the
text layout, bookmarks and links, row columns, chips. Replace the current Antibody stand-in if one exists.

**Model**: fast (sonnet). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B seq/da-antibody
<that commit>`. **Storybook port**: 6016.
