# History: the specification

`packages/ui/src/history/`: generic undo and redo over an app's checkpoints, and a view of the steps. Code comments cite
its sections (`history §3`).

| File                  | What                                                                                          |
| --------------------- | --------------------------------------------------------------------------------------------- |
| `history-model.ts`    | `HistoryModel`, `HistoryPort`, `HistoryStep`, `HistoryOrigin`, limits, labels (§1–§3, §5, §6) |
| `history-view.tsx`    | `HistoryView`, `historyRows`, `DefaultHistoryOrigins` (§4)                                    |
| `history-commands.ts` | `historyCommands`: `edit.undo` (⌘Z), `edit.redo` (⇧⌘Z) (§2)                                   |

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
- **Not recording** (`recording` option, `setRecording(false)`; or `recordIf(origin, label)` false for a top-level run
  or group, asked before it runs): it runs in order with no checkpoint and no step, and once it changed something (not
  `record` false) it is a change outside the history: `invalidate` with `notRecorded` (code `not-recorded`). Why: Undo
  past it would silently wipe it out. `undoable: false` differs: the run says it changed nothing the history restores,
  so earlier steps stay safe. Nested work follows its outer run.
- **An epoch** (bumped by `clear` and `invalidate`) stamps each checkpoint. A run whose checkpoint predates a `clear`
  records nothing; one that predates an `invalidate` records its step unsafe with that reason; neither restores its
  checkpoint on failure (it would undo the clear or the outside change too). Undo, Redo and `undoTurn` that a clear
  overtook while restoring leave the lists alone. Why: the run took its checkpoint before the change, and Undo would
  restore the old state.
- `clear()` forgets the steps and Redo, turns checkpoints back on and clears `checkpointsOff`; an open turn's next run
  starts a new entry under its label.

## 2. Undo and redo

- `undo()` restores the latest step; it never queues behind a run (blocked while one runs). It checkpoints the state it
  leaves, and `redo()` restores that, until a new step or `invalidate` clears the list.
- **Every action is a method** (for commands and scripts): `undo`, `redo`, `undoTo(id)` (that step and every later one,
  each as an Undo), `redoTo(id)` (up to that undone step), `moveTo(id | undefined)` (the position just after a step, or
  before every step: History's "Undo to <step>" and "Redo to <step>"), `undoTurn(origin | turnId, { discardLater })`,
  `clear`, `setCheckpointing`, `setRecording`, `setLimits`.
- `undoTurn` matches on the origin's group key (`historyGroupKey`: kind and turn), as turns group; a bare turn id is an
  `agent`'s. Why: two kinds with the same turn id are two turns.
- **A refusal has a code** (`HistoryRefusalCode`: `running`, `nothing-to-undo`, `nothing-to-redo`, `not-in-history`,
  `latest`, `turn-empty`, `later-steps`, `unsafe`, `no-checkpoint`) and its text from the labels (§6); an unsafe step
  has `unsafeCode` (`changed`, `not-recorded`, `too-old`, `checkpoints-off`, `too-large`). Why: an app says them its own
  way, and a view tells the transient `running` from the rest.
- States: `steps` (newest first), `redoable` (the next first), `canUndo`, `undoBlocked` (why not) and `undoBlockedCode`,
  `canRedo`, `running`, `checkpointing`, `recording`, `lastMove` (what Undo or Redo last did: a new value each time),
  `checkpointsOff`; `labels`, `limits`.
- **Commands**: `historyCommands(history, category = 'Edit', labels?)` gives `edit.undo` and `edit.redo`
  (`HistoryCommandIds`, `HistoryShortcuts`), so the menu, a key and an agent do what History's buttons do. `when` is
  `undoBlocked`'s reason, and Redo's is the model's `running` while a run is in flight, else "Nothing to redo"; a
  refused one fails the command with its reason; a done one says `undone(step)` / `redone(step)`.
- **A turn's runs join its entry in place**: each published step reads its own length of an append-only list, copied out
  when read, so a long turn is linear, not a copy per run (guarded: `_test/history.scaling.test.ts`, heavy lane).

## 3. Checkpoints that switch off

- Each checkpoint is timed (its size pass included) and sized. Three in a row over the time budget (50 ms), or one over
  the size budget (5 MB), turn per-step checkpoints off (`checkpointsOff` says so once): a step is then recorded unsafe
  (`checkpoints-off`), while a group still takes its one checkpoint.
- Over the hard limit (500 ms or 50 MB), even a group's is dropped (`too-large`), and for 30 s none is taken at all
  (each try is a snapshot thrown away). `clear()` and `setCheckpointing(true)` try again at once.

## 4. The view

- `HistoryView` fills its parent (a flex column): a `ViewToolbar` with Undo and Redo (their tooltips: the step and the
  key, or why not), the step count; the reason Undo is blocked as a muted line when there is something it would undo,
  but not for `running` (transient: the buttons are greyed; a note on every run flashed and named no control to stop);
  then the rows.
- **The keys** (`aria-keyshortcuts` and the tooltips): `keys` in ARIA's syntax when given; else the commands' keys now
  from the app's `commands` manager (after its keymap: a remap or a removal shows); else `HistoryShortcuts` for this
  platform (`KeyPlatform`: ⌘Z / `Meta+Z` on macOS, Ctrl+Z elsewhere; `ariaKeyShortcut`, `formatChord`). Why: the
  commands bind `mod`, and a fixed Control+Z was wrong on macOS and after a remap.
- **The rows** (`historyRows`, pure): Start, the steps done oldest first, the position at the latest
  (`aria-current="step"`, the accent), then the undone steps Redo would put back, muted after it. Under a group's step
  (a level in): its latest `runs` (10) runs and, above them, a row counting the earlier ones with "Show all" (pressed:
  every run the model keeps, at most `maxRuns`). Rows read `lastRuns`, never the whole list. Why: each run of a long
  turn publishes its step again; listing every run as rows was R² (12.5M rows for 5,000 runs), guarded
  (`_test/history-view.scaling.test.tsx`, heavy lane, mounted). A `VirtualList` of `ViewRow`s: the label, the origin's
  name and time as the description (the side tooltip).
- On a row (shown on hover or focus): "Undo to <step>" before the position, "Redo to <step>" after it (named after their
  step: a list of "Undo to here" told a screen reader nothing), "Undo every step" on Start: all `moveTo`. An unsafe step
  shows a warning icon whose tooltip is the reason. Results go to `onResult` (the app says a failure, its status line).

## 5. Limits

- `HistoryLimits` (options, `setLimits` at runtime, `limits`): `capacity` (100) undo steps kept, the oldest dropped;
  `checkpoints` (10, at most `capacity`) the latest keep theirs, older ones listed unsafe (`too-old`); `maxRedo` (100)
  redo steps kept, the farthest dropped (0: no Redo, and Undo takes no checkpoint of what it leaves); `maxRuns` (1000)
  runs listed in one step, the latest. Past any, what goes is thrown away, its checkpoints released. Why: an app bounds
  memory and a long session's list by what it can afford.
- A step past `maxRuns` keeps its checkpoint and its count (`runCount`, `droppedRuns`); its log keeps the latest
  `maxRuns`, cut to that once it doubles (a new list: the steps published before read the old one), so a long turn holds
  O(maxRuns) runs and copies each once. A lowered limit applies at once.

## 6. Words

- Every component takes `labels?: Partial<XLabels>` over an exported `DefaultXLabels` (English): strings, or functions
  for what is parametrized. The model takes `labels` in its options (`HistoryLabels`: the refusals' and unsafe steps'
  text), and `historyCommands` and `HistoryView` read the model's for the reasons they repeat. Defaults are
  domain-neutral ("Too large to undo", not a scene's). Why: base parts an app composes in its own language and domain.
