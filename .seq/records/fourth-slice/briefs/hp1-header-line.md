# Lane HP1: core's header line

**Task**: build the core half of [track-header-parts](../../../design/track-header-parts.md) (Decision: Two headers, one entry; The track header's geometry; The section header; Chips; Hit, cursor, announced; Consequences `seqx-core`): the header line's slots as given (path data, the entry's name, enabled), chips (text, box), the section header under `section.header`, short-track marks, hits `header-slot`, `header-state`, `section-header`; the filter chip becomes the `seq:filtered` state; `TrackLayout.header` replaces `TrackLayout.filter`, `SectionLayout.header`. No registry kind, behavior or DOM parts (HP2).

**Output**: core changes with the property tests the design lists (slots never overlap each other, a chip or the "⋯"; a slot's x independent of the label; every hidden track counted once in its chip), drawing against the data for every example view, the records that change (a filter chip moves before the action columns), a quick representative (a pointer move over 200 tracks' label column allocates nothing); `_spec` Render (Tracks layout, Filters, Lanes and Hits).

**Model**: strong (opus).

**Where**: A worktree on a branch `seq/hp1-core's` from the base the orchestrator names (after the third slice's close). Hot files (plan.md) are the orchestrator's: list the changes you need in the report. Every modified file's header updated; a quick-lane test beside the code, failing on the parent for a fix; the look in a visible pane, light and dark, before the report.

**Report**: Under 300 words: commits; checks run and their result; tests shown failing on the parent; the public surface added (names, for the author's acceptance where the brief says so); drift notes; what you saw; the `.v6/changelogs/` entry's text.
