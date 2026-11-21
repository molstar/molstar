# Lane HD: host data: the fake host's adapter

**Task**: build the Storybook half of [seq-data-ports](../../../design/seq-data-ports.md) (The fake host's adapter; Consequences, Storybook): the fake viewer's adapter serving 2OCJ's nodes (and a pLDDT track by `values.field`) from its chain table through `withAdapter`, and the `Host data` story (the fake serving 2OCJ, "drop chain" followed by the sequences through `changes` and `reattach`) with its headless twin on a fake clock asserting no asset read and nodes equal to the default adapter's; the Extensibility page lists the data adapters and `withAdapter`.

**Output**: the adapter, story, twin and page text.

**Model**: fast (sonnet).

**Where**: A worktree on a branch `seq/hd-lane` from the base the orchestrator names (after the third slice's close). Hot files (plan.md) are the orchestrator's: list the changes you need in the report. Every modified file's header updated; a quick-lane test beside the code, failing on the parent for a fix; the look in a visible pane, light and dark, before the report.

**Report**: Under 300 words: commits; checks run and their result; tests shown failing on the parent; the public surface added (names, for the author's acceptance where the brief says so); drift notes; what you saw; the `.v6/changelogs/` entry's text.
