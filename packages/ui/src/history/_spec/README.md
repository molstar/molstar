# History: the specification

`packages/ui/src/history/`: generic undo and redo over an app's checkpoints, and a view of the steps. Code comments cite
its sections (`history §3`).

| File                  | What                                                                  |
| --------------------- | --------------------------------------------------------------------- |
| `history-model.ts`    | `HistoryModel`, `HistoryPort`, `HistoryStep`, `HistoryOrigin` (§1–§3) |
| `history-view.tsx`    | `HistoryView`, `historyRows`, `DefaultHistoryOrigins` (§4)            |
| `history-commands.ts` | `historyCommands`: `edit.undo` (⌘Z), `edit.redo` (⇧⌘Z) (§2)           |

## 1. Steps

- **DOM-free**, on core's states (no streams): the app gives checkpoints through a `HistoryPort` (`checkpoint()`
  synchronous and cheap, `restore(c)`, optionally `size(c)`); clocks injected (`now` for a step's time and the pause,
  `measure` for a checkpoint's cost).
- **A step is one undoable run, or one group**, with a checkpoint taken before it. `group` is atomic (all or nothing;
  with `onFailure: 'continue'` each run inside keeps its own checkpoint). `beginGroup`/`endGroup`: an origin's turn (an
  agent's), whose runs join one entry while it is the latest; another origin's call closes it.
- **Top-level work is serialized** (one lock): each run has its checkpoint, so a failed one is restored without touching
  another's work. Work that passes its context back (`within`) is part of that step: no lock, no step.
- A failed run restores its checkpoint and records nothing; one that changed nothing (`record` false) or nothing
  undoable (`undoable: false`) is no step.
- `invalidate(reason)`: a change outside the history makes every earlier step unsafe (and empties Redo).
- `capacity` (100) steps listed; only the latest `checkpoints` (10) keep theirs, older ones are listed as too old.

## 2. Undo and redo

- `undo()` restores the latest step; it never queues behind a run (blocked while one runs). It checkpoints the state it
  leaves, and `redo()` restores that, until a new step or `invalidate` clears the list.
- **Every action is a method** (for commands and scripts): `undo`, `redo`, `undoTo(id)` (that step and every later one,
  each as an Undo), `redoTo(id)` (up to that undone step), `moveTo(id | undefined)` (the position just after a step, or
  before every step: History's "Undo to here" and "Redo to here"), `undoTurn(turn, { discardLater })`, `clear`,
  `setCheckpointing`.
- States: `steps` (newest first), `redoable` (the next first), `canUndo`, `undoBlocked` (why not), `canRedo`, `running`,
  `checkpointing`, `lastMove` (what Undo or Redo last did: a new value each time), `checkpointsOff`.
- **Commands**: `historyCommands(history, category = 'Edit')` gives `edit.undo` and `edit.redo`, so the menu, a key and
  an agent do what History's buttons do. `when` is `undoBlocked`'s reason or "Nothing to undo"; a refused undo fails the
  command with its reason; a done one says "Undone: <step>".
- **A turn's runs join its entry in place**: each published step reads its own length of an append-only list, copied out
  when read, so a long turn is linear, not a copy per run (guarded: `_test/history.scaling.test.ts`, heavy lane).

## 3. Checkpoints that switch off

- Each checkpoint is timed (its size pass included) and sized. Three in a row over the time budget (50 ms), or one over
  the size budget (5 MB), turn per-step checkpoints off (`checkpointsOff` says so once): a step is then recorded unsafe
  ("Undo is off for this large scene"), while a group still takes its one checkpoint.
- Over the hard limit (500 ms or 50 MB), even a group's is dropped, and for 30 s none is taken at all (each try is a
  snapshot thrown away). `clear()` and `setCheckpointing(true)` try again at once.

## 4. The view

- `HistoryView` fills its parent (a flex column): a `ViewToolbar` with Undo and Redo (their tooltips: the step and the
  key, or why not; `aria-keyshortcuts`, `Control+Z` and `Control+Shift+Z` by default), the step count; the reason Undo
  is blocked as a muted line when there is something it would undo; then the rows.
- **The rows** (`historyRows`, pure): Start, the steps done oldest first (a group's runs under it, a level in), the
  position at the latest (`aria-current="step"`, the accent), then the undone steps Redo would put back, muted after it.
  A `VirtualList` of `ViewRow`s: the label, the origin's name and time as the description (the side tooltip).
- On a row (shown on hover or focus): "Undo to here" before the position, "Redo to here" after it, "Undo every step" on
  Start: all `moveTo`. An unsafe step shows a warning icon whose tooltip is the reason. Results go to `onResult` (the
  app says a failure, its status line).
