## Fixes from the perf audit

- Fix: a visibility (or other state) flag set while a state builder is outstanding (made but not yet committed, or a
  commit waiting in the queue) is no longer lost from the state tree: saved sessions, undo and `setSnapshot` keep it.
  A builder can also be edited again after `getTree()` (it threw `Could not find node`)
