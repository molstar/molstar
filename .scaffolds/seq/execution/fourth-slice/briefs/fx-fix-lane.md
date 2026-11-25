# Lane FX: a fix lane over fixes.md

**Task**: take the open items of [fixes.md](../fixes.md) that need no design, at the fix point the orchestrator names (FX1 before stage 1, FX2 to FX4 after stages 1 to 3's gates, FX5 before M4). For each: reproduce it in the story or test named, find the cause, fix it with a test shown failing on the parent, look at it in light and dark. An item whose cause needs a design or a public API change is left open with the reason, not fixed. Read the `_spec` of the package it touches first.

**Output**: the fixes on one branch, one commit per item; fixes.md's Open lines moved to Fixed with the commit; no new feature.

**Model**: fast (sonnet), strong only when the orchestrator says a cause is in a model.

**Where**: A worktree on a branch `seq/fx-N (N the fix point)` from the base the orchestrator names (after the third slice's close). Hot files (plan.md) are the orchestrator's: list the changes you need in the report. Every modified file's header updated; a quick-lane test beside the code, failing on the parent for a fix; the look in a visible pane, light and dark, before the report.

**Report**: Under 300 words: commits; checks run and their result; tests shown failing on the parent; the public surface added (names, for the author's acceptance where the brief says so); drift notes; what you saw; the `.v6/changelogs/` entry's text.
