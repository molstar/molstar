# Lane NG: the resolved view's test consumer

**Task**: write the second consumer of the resolved view, test-only: a plain-SVG drawing of letters, ranges, values and variants from the resolved view and the store alone, no layout or draw list of core's (EC-183, [version-1](../../../design/version-1.md), Integration). Nightingale is dropped (the author, 2026-11-19): no dependency, no story. The field test: the same spec drawn by core and by the consumer agrees on what is where (positions, labels, colours by value) for a subset of the examples, and the consumer needs no private import. It lives in a `_test/` folder, not published.

**Output**: the consumer and its tests; a note in `_spec` Resolved view of what the consumer showed missing from the public view, if anything.

**Model**: fast (sonnet).

**Where**: A worktree on a branch `seq/ng-lane` from the base the orchestrator names (after the third slice's close). Hot files (plan.md) are the orchestrator's: list the changes you need in the report. Every modified file's header updated; a quick-lane test beside the code, failing on the parent for a fix; the look in a visible pane, light and dark, before the report.

**Report**: Under 300 words: commits; checks run and their result; tests shown failing on the parent; the public surface added (names, for the author's acceptance where the brief says so); drift notes; what you saw; the `.v6/changelogs/` entry's text.
