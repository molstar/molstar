# Lane CU: copy ends on units

**Task**: rewrite CE on units (Gate 0, 8): a pair's end on split copies (`copies: each`) resolves by unit, replacing grids-and-rows.md's "end on the first copy" and the caret (`seq.render.pairs.copy`). The parked branch `seq/cu-lane` (`189097018`) and its review are the cases, not a rebase; delete the branch after. The review's findings, each a test failing on the parent: (1) both ends on split copies drew nothing (the split hemoglobin lost every contact: it keeps them now); (2) "copies" counted as blocks holding a position, not the resolver's copies (one system, one range): the default merge mode with a zoom block dropped pairs and said a wrong issue, a pair within one split system vanished; (3) the issue counted filtered pairs and named only the first system a layer; (4) an end out of its block's range beside a copy end lost its mark and `pairs.systems`; (5) hits: a tie resolves to the first, arcs to the later, a Map iterator a layout; (6) tests by shape of the real example (Globins › Hemoglobin › Copies Split and merged). Read [grids-and-rows](../../../design/grids-and-rows.md) (the pairs layer's scope) and [molecules](../../../design/molecules.md) (copy groups, per-row units). Canvas and SVG alike.

**Output**: the pairs layer's ends on units, tests, `_spec` and grids-and-rows.md corrected (the design's change is the author's: list the text).

**Model**: strong (opus).

**Where**: A worktree on a branch `seq/cu-lane` from the base the orchestrator names (after the third slice's close). Hot files (plan.md) are the orchestrator's: list the changes you need in the report. Every modified file's header updated; a quick-lane test beside the code, failing on the parent for a fix; the look in a visible pane, light and dark, before the report.

**Report**: Under 300 words: commits; checks run and their result; tests shown failing on the parent; the public surface added (names, for the author's acceptance where the brief says so); drift notes; what you saw; the `.v6/changelogs/` entry's text.
