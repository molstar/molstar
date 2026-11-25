## Core: three-way merge of repr trees

- `@molstar/core/repr/patch`: `merge3(base, ours, theirs, options)`, a three-way merge of revisions by node. It returns
  the merged root with its net `changes` from ours, or the conflicts (`set`, `remove`, `insert`, `replace`, `order`,
  `invalid`; `move` reserved) with a validated `partial`. A node both sides inserted, placed two ways, is an `order`
  conflict. A stale patch is rebased as `merge3(base, current, applyPatch(base, p).root)`.
