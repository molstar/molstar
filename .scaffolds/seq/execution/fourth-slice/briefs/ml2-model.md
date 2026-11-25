# Lane ML2: the model: units, atoms, copies, the modified index, queries

**Task**: build the resolve and query half of [molecules](../../../design/molecules.md) (Decision: Complexes and units, Atoms to residues, Queries, Reading, Language, store, resolved model; Consequences `resolve/`, `query/`) and the units half of [maps](../../../design/maps.md): unit and atoms systems, the atoms group map applied at a path's ends only, held pairs with per-row units, copy groups derived equal for any unit order, the modified index, `ResolveConfig.maxUnits` (65,536), and the query forms (`component`, `modified` and the others the design lists; the forms `search`, `onto`, `expand`, `group` stay out, Gate 0, 11). Stay on resolve and query: no drawing (MR), no mmCIF (MC). `unitRows` is ST1's.

**Output**: the model with property tests (copies for any unit order; atoms collapsed then expanded cover what they started from; a path's result unchanged by a molecule's atoms form; a modified position's letter is its parent's) and a twin per consumer (follow, hover projection, selection, export read a unit as a system); `_spec` Resolve (Systems; Paths) and Queries (Forms); a quick representative of 200 merged copies.

**Model**: strong (opus).

**Where**: A worktree on a branch `seq/ml2-lane` from the base the orchestrator names (after the third slice's close). Hot files (plan.md) are the orchestrator's: list the changes you need in the report. Every modified file's header updated; a quick-lane test beside the code, failing on the parent for a fix; the look in a visible pane, light and dark, before the report.

**Report**: Under 300 words: commits; checks run and their result; tests shown failing on the parent; the public surface added (names, for the author's acceptance where the brief says so); drift notes; what you saw; the `.v6/changelogs/` entry's text.
