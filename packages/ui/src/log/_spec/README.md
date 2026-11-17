# Log: the specification

`packages/ui/src/log/`: a generic log any app or view keeps (a console's output, an import's report, an agent's
actions), and its view. Code comments cite its sections (`log §2`).

| File           | What                                                                               |
| -------------- | ---------------------------------------------------------------------------------- |
| `log-model.ts` | `LogModel`, `LogEntry`, `LogFilter`, `boundLog`, `passesLog`, `isProblem` (§1, §2) |
| `log-view.tsx` | `LogView`, `LogKindIcons`, `logDescription` (§3)                                   |

## 1. The model

- **DOM-free** (headless): no document, no frame; a script, a test or a headless host writes and reads the same model.
  Its clock (`now`) and its scheduler (`schedule`, a microtask by default) are injected.
- **An entry**: `{ id, kind: 'info' | 'success' | 'warning' | 'error' | 'debug', text, time, source?, detail? }`, plain
  data (`save()` gives the entries; they survive JSON). `time` is ms since the epoch. No `source`: a message (what a
  status line showed); a source: what wrote it (a filter each).
- **`entries` are oldest first** (newest last). Why: a log reads down, and the view keeps its end in view.
- `sources` are the sources as first written (not sorted), told only when one is new: the filters do not move while
  entries come.
- `filter` is `{ kinds?, source? }`: the kinds shown (any of them; none: every kind) and the source (a name, `null` for
  the messages, absent for every source); the two combine (`LogAll`, `LogProblems` = warning and error).
  `toggleKind(kind)` (from every kind: that kind alone; else add or take out; the last out: every kind) and
  `setSource(source)` keep the other half. Why kinds and source apart: "errors from Transforms" is a question people
  ask, and a kind is what a person scans for. It and a text `query` (case-insensitive, over the text, the detail and the
  source) are the model's, one for every view of it, so a status line can open the log at its problems. `filtered` is
  derived from them.
- `clear()` forgets the entries and the sources; a filter that named a source goes back to every source.

## 2. Writing and the bound

- **`write` is batched**: an entry gets its id and time at once, and lands in `entries` at the end of the task with the
  rest of its burst: one change, one copy of the list. `flush()` (and `save()`) writes what waits now. Why: a load
  writes hundreds of entries in a task; one change a burst keeps the view and every reader at one pass.
- **One bound** (`size`, 300) for every source together: past it, the oldest entries go, warnings and errors last, so a
  burst never pushes a problem out while anything else is left to drop; only when the problems alone fill it do the
  oldest problems go (`boundLog`, one pass).
- Guarded: a burst of n and 4n entries, and n entries into a bounded log ten a flush, grow about linearly
  (`_test/log.scaling.test.ts`, the heavy lane; the quadratic control is caught).

## 3. The view

- `LogView` fills its parent (a flex column; styles §8.4: a view that virtualizes fills the body). A `ViewToolbar`: the
  filters as icon toggles, each named by its `aria-label`, said in its tooltip, `aria-pressed`: All, then a kind each
  (Info, Success, Warnings, Errors; debug shows under All), then, once something wrote with a source, Messages and a
  source each (pressed again: every source), then the text filter (`TextInput search`: Escape empties it). Then the
  entries.
- **The entries are a `RecycledList` of `ViewRow`s**: one line each, a kind icon (its color by kind), the text fading
  where cut, the detail, source and time as the description, the whole in the row's side tooltip (`data-tip-row`). Why:
  a long session's log costs a screenful of rows, and scrolling re-points them.
- **Kept at the bottom**: while the list ends in view, a change scrolls to its new end once a frame (`frames`, the
  scope's); scrolled up, it stays. The scroll is heard in the capture phase on the view's root (read only).
- `role="log"` with `aria-live="off"`: a log is live by default, but each message is said once, by a status line's
  polite region.
- Empty: "Nothing is logged yet." or "Nothing matches the filter." (`ViewEmpty`).
