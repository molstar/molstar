# Lane ST2: linkStructure and the viewer's four additions

**Task**: build the viewer half of [structure-syncing](../../../design/structure-syncing.md) (Wiring; "Color by this track"; What the viewer publishes for the wiring; Consequences `seqx-viewer`): `HoverModel.setLocation(loc, owner)`, `AddressModel.revealLocation(loc)`, `view.lastJump` and `viewer.retiredSignal` (accepted by the author as the design states; stop on anything past them), then `@molstar/seqx-viewer/link-structure`: `linkStructure` over a host's `side` (plain data), echoes never sent back, coalescing (one selection message a frame), a stale `null` ignored, focus with the jump's duration, clearing on either end and on retire, `via` per binding. No behavior, no port, nothing in the composition.

**Output**: the four additions and the entry with tests on a recording side (each path both ways, echoes, a stale null, one selection message a frame, focus on a fake clock, clearing at each of three ends); `_spec` the models' additions and a Structure link section.

**Model**: strong (opus).

**Where**: A worktree on a branch `seq/st2-lane` from the base the orchestrator names (after the third slice's close). Hot files (plan.md) are the orchestrator's: list the changes you need in the report. Every modified file's header updated; a quick-lane test beside the code, failing on the parent for a fix; the look in a visible pane, light and dark, before the report.

**Report**: Under 300 words: commits; checks run and their result; tests shown failing on the parent; the public surface added (names, for the author's acceptance where the brief says so); drift notes; what you saw; the `.v6/changelogs/` entry's text.
