# Drift: clade names collide on repeated labels (CH's review)

**What**: a clade is named by its node's id, else its first and last leaves (grids-and-rows, Clades: "named by the
node's id, else their first and last leaves in tree order"). The rows model takes an inner node's label as its id, and
labels repeat: support values (`((A,B)95,(C,D)95)`). `byName` keeps the first node of a name, so hovering the second
"95" lights the first clade, and a fold or `clades.keep` of either names the same clade. The hit's extent (the rows) is
right; the name is not. Predates CH (folds have it), made visible by the clade hover.

**Where**: `packages/seq/core/src/render/rows-model.ts` (`buildClades`, names), `rows.ts` (`byName` lookups), the
`collapsed` param's addresses.

**Options**:
- (a) A node's id only where it is unique in the tree; a repeated or numeric-only label (a support value) falls back to
  `(clade :of [first last])`. Names stay stable across reorders; a writer's unique ids still name their clades.
- (b) Always `first..last`; ids only as labels. Simpler, but a writer's named clade (`Primates`) is no longer its name.
- (c) Refuse a tree with repeated inner ids (`seq.lang.tree.id`). Breaks common Newick files.

**Recommendation**: (a). Built as a small fix in GD (it touches the rows model's naming and its tests).
