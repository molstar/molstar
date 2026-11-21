# Lane ST1: the structure library in core

**Task**: build `@molstar/seqx-core/structure` (a new entry; a library, not a port, no behavior: [structure-syncing](../../../design/structure-syncing.md), Decision: Who is who, Translation; Consequences `seqx-core`): `structureFacts`, `unitRows`, `toUnits` and `fromUnits` (the size gate and `seq.structure.*` issues), `colorRuns`; core's `SeqPorts.structure` stub removed. Seq\*'s locations only, no MolViewSpec types. Tests: both round trips, the gate monotone, order irrelevant, EC-25, EC-48, EC-65, EC-118, EC-1 to EC-8, a quick representative (one 2OCJ hover carried within a bound) and a heavy `maxUnits`-scale hover; a description round-tripped through DP's `structureNodes` and `structureFacts` gives the same unit facts. The size of the minimal viewer does not grow (its own entry).

**Output**: the entry, tests, `_spec` Structures (and Paths where `via` is read).

**Model**: strong (opus).

**Where**: A worktree on a branch `seq/st1-lane` from the base the orchestrator names (after the third slice's close). Hot files (plan.md) are the orchestrator's: list the changes you need in the report. Every modified file's header updated; a quick-lane test beside the code, failing on the parent for a fix; the look in a visible pane, light and dark, before the report.

**Report**: Under 300 words: commits; checks run and their result; tests shown failing on the parent; the public surface added (names, for the author's acceptance where the brief says so); drift notes; what you saw; the `.v6/changelogs/` entry's text.
