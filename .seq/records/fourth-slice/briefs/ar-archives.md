# Lane AR: archives over the assets port

**Task**: read and write zip with `@molstar/core/util/zip` ([version-1](../../../design/version-1.md), Integration; [security](../../../design/security.md): EC-91's bounds): entries resolved as assets by the assets port, a view saved with its data and opened again; every dataset saved and opened; the bounds (entry count, expanded size, ratio, names) by fuzz. The writer half is cut fourth if usage runs short (plan, Gate 0, 2).

**Output**: the archive reader and writer, fuzz tests, a story, `_spec` Data.

**Model**: fast (sonnet).

**Where**: A worktree on a branch `seq/ar-lane` from the base the orchestrator names (after the third slice's close). Hot files (plan.md) are the orchestrator's: list the changes you need in the report. Every modified file's header updated; a quick-lane test beside the code, failing on the parent for a fix; the look in a visible pane, light and dark, before the report.

**Report**: Under 300 words: commits; checks run and their result; tests shown failing on the parent; the public surface added (names, for the author's acceptance where the brief says so); drift notes; what you saw; the `.v6/changelogs/` entry's text.
