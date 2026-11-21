# Lane ST3: the fake structure viewer and its stories

**Task**: build the Storybook half of [structure-syncing](../../../design/structure-syncing.md) (The fake host's adapter; The fake structure viewer; Datasets): a component in `_stories/support/` importing nothing of Seq\*, on `@molstar/ui`, whose states are MolViewSpec validated with `@molstar/mvs-builder`; its own recorded chain table (from MX's 2OCJ and 1FIN records); rows of residue cells (unobserved as dashed gaps, ligands as chips, copies as rows), the state's `color` nodes, a hover outline, the selection, a camera frame that `focus` moves; its grammar (hover; click selects; Cmd/Ctrl+click toggles; Esc clears; double-click focuses; "select chain", "within 5 Å"), `latency`. Its adapter beside it translates on `structureFacts` and `unitRows` and gives `linkStructure` its `side`. Stories: Hover both ways, Selection both ways, Focus and back, Chain maps (SIFTS against an alignment by `priority` and `via`), Color by this track (a host's header entry through HP2), each with a headless twin on a fake clock (one with latency); the Design page. The twin tests: every lit cell in the fake in the image of Seq\*'s hover and back; the fake's colors equal `colorRuns` at each observed residue's preimage.

**Output**: the component, adapter, stories, twins and page.

**Model**: fast (sonnet).

**Where**: A worktree on a branch `seq/st3-lane` from the base the orchestrator names (after the third slice's close). Hot files (plan.md) are the orchestrator's: list the changes you need in the report. Every modified file's header updated; a quick-lane test beside the code, failing on the parent for a fix; the look in a visible pane, light and dark, before the report.

**Report**: Under 300 words: commits; checks run and their result; tests shown failing on the parent; the public surface added (names, for the author's acceptance where the brief says so); drift notes; what you saw; the `.v6/changelogs/` entry's text.
