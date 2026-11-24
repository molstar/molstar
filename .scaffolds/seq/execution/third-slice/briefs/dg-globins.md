# Lane DG: the Globins dataset

Follow [datasets.md](datasets.md) for [examples](../../../design/examples.md) › Globins (its section's entries, the
family alignment and tree, the 4HHB hemes as this slice allows: one-position sequences drawn as chips, molecules.md
decision 14; the views of its table). Build on what is merged: the rows layer, clades, the logo (saturated letters),
row columns (species and identity; RC), chips (LB), the zoom strip and navigators. The tree is read from the file as
published if FS-B is needed; else (the plan) inline Newick: say which and why. You own the `sources.json` test.
Replace the current Globins stand-ins with the recorded data, keeping the existing stories (Rows, Family, Folded,
Figure, Scale, ScaleWithNavigators, Lazy) working against the real anchors.

**Model**: strong (opus). **Base**: `packages-seq` at the commit this brief is in: `git checkout -B seq/dg-globins
<that commit>`. **Storybook port**: 6014.
