# Status: the specification

`packages/ui/src/status/`: what a status line shows (the latest message, running tasks, chips) and its views. Code
comments cite its sections (`status §2`). The app places the line (a stage's bottom-left slot) and the chips (its bottom
right); the views are standalone.

| File              | What                                                                    |
| ----------------- | ----------------------------------------------------------------------- |
| `status-model.ts` | `StatusModel`, `StatusMessage`, `StatusTask`, `StatusChip` (§1)         |
| `status-line.tsx` | `StatusLine`, `StatusPill`, `TaskPills`, `TaskPill`, `StatusChips` (§2) |

## 1. The model

- **DOM-free** (a `ReactiveModel`): its fade's timer is injected (`timers`, the globals by default), so Node tests and a
  headless host run it; `retire()` clears it.
- **Messages**: `post(text, kind, { detail })` shows a message (`current`), says it (`spoken`) and writes it into its
  log (`log`, a `LogModel`, given or its own: no source, so the log's Messages source shows them). It fades after
  `fadeMs` (4000); **an error never fades** (it stays until `dismiss`). A newer message restarts the fade; a stale fade
  never hides a newer message (`dismiss(id)`).
- **A message's scope** (`post(…, { scope })`, a plain string the app chooses: a pane's id, a document's): when the app
  shows several status lines (a stage a pane), the message shows in the line of its scope, not where the command came
  from (a shared toolbar acting on the active pane posts with that pane). One message on show at a time, wherever it
  shows. The status model knows no panes. Why: a result belongs next to what it changed.
- **`spoken`** is what the one polite live region says: the latest message posted or announced; the same words again get
  a trailing no-break space, so they are said again. `announce(text)` says without showing or logging (a selection's
  result).
- **Tasks from two writers**: `setTasks` mirrors one source whole (an app's task list); `putTask`/`removeTask` add one
  by one (an extension loading). `setTasks` never drops a put task. Overlay tasks (`overlay`) hold the whole app: in
  `tasks`, not the line's.
- **Kept by id**: `taskIds` (the line's, told only when the ids change, not on progress) and `task(id)` (one task's
  record, a `KeyedState` key). Why: a task's record is replaced on every progress event; a pill made anew restarted its
  spinner and replaced its Cancel under a press.
- `chips` (`setChips`): what the person works on, each with `activate` and `dismiss`.
- Guarded: a burst of n posts and n progress events grow about linearly (`_test/status.scaling.test.ts`, heavy lane).

## 2. The views

- **`StatusLine`**: a column growing upward: a visually hidden `role="status"` `aria-live="polite"` region (`spoken`),
  the message as a `StatusPill` (a button when `onOpenLog` is given: it sets the log's filter, `LogProblems` for an
  error, else `LogAll`, then calls the app), the task pills, then `after`. Ink on what is behind it, no frame, a halo of
  the background (`--msui-status-ink`, `--msui-status-halo`, which a stage may set). **`scope`** (several lines):
  `{ id }` shows that scope's messages; `{ rest: true }` the messages without a scope, the tasks and the live region
  (one line in the app takes the rest, so nothing is said twice); both on one line (the app's own stage). No scope:
  everything (one line in the app).
- **`StatusPill`**: an icon, its text (fading where cut), what follows, a ✕. With `autoDismiss` (ms) it goes on its own,
  a ring around its ✕ drawing the time left; the pointer over it or the focus in it holds the time (the timer is the
  scope's `restartable`: paused with hidden UI, gone with the line).
- **`TaskPill`**: made once per id (`For` over `taskIds`), one `watch` writes its label, its progress bar
  (`role="progressbar"`, hidden when indeterminate) and whether Cancel shows (hidden by its attribute; the click reads
  the latest `abort`).
- **`StatusChips`**: a raised chip each; a click activates it, its ✕ clears it.
- **Text is the app's**: `labels` (`StatusLabels`, `DefaultStatusLabels` English) on `StatusLine`, `StatusPill`,
  `TaskPills`, `TaskPill`, `StatusChips`, one prop on each (no per-string props); parametrized ones are functions. A
  task's Cancel is named after its task (`cancel(label)`, `Cancel Loading 1tqn`): two running tasks' buttons are told
  apart.
