# Lane DA: the Antibody dataset

Follow [datasets.md](datasets.md) for [examples](../../../design/examples.md) › Antibody (a Fab's heavy and light
chains numbered Kabat and IMGT, insertion codes, IMGT out of index order, the CDRs; a designed variant with indels
against the wild type; the views of its table: numbering schemes side by side, Kabat vs IMGT labels (EC-2 to 8), design
vs wild type with indel carets (EC-29), CDR selection and editing as bookmarks with undo (EC-162, 164)). The fourth
slice's Fc glycan (molecules.md) is not this lane's. Build on what is merged: the axis and numberings, rulers, the
text layout, bookmarks and links, row columns, chips. Replace the current Antibody stand-in if one exists.

The numberings are declared with the `numbering` kind (NB, merged: `method: labels`, inline or a records source's
column, `scheme`); the first view that shows them is this one, so look at the ruler, hover and copied addresses in
both numberings in Storybook, light and dark. Records from SAbDab are CC BY 4.0: attribution and changes stated, as
Globins and Complex do (their `sources.json` and index). The shared mmCIF reader is `src/files/cif.ts`; the dataset
checks are `src/_test/sources.test.ts`.

**Model**: strong (opus; a strong slot was free, 2026-11-19). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B seq/da-antibody
<that commit>`. **Storybook port**: 6016.

**Checks**, **headers**, **hot files** and **stop when**: as [datasets.md](datasets.md) says. CY (source-backed axes)
may be merging beside you.
