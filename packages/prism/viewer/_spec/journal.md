# The session record and what the viewer stores

The session's record (an `@molstar/ai` thread), History, kept records, and light state in browser storage. The code is
the truth for how.

## 1. One record per viewer session

- **One record, one writer** (`journal/record.ts`, `createSessionRecord`): `run` per finished API call (any origin),
  `turn` per agent turn, `undo`/`redo` from the history, `session` for open, clear, save. One write chain: record order
  is the order things happened. Not recorded: reads, `record: false` operations (pointer-rate view state:
  `trajectory.frame`, `snapshot.seek`), `highlight.*`, `view.*`, `export.preview`, nested runs.
- **In memory by default.** A reload ends it and every checkpoint: a turn cannot be undone after a reload.
- **Kept on purpose**: "Keep this session's record" (History, `record.keep`) copies it into the spec's `persist.rows`
  ([viewer-spec.md](viewer-spec.md) §2.6), else IndexedDB `<persist.key>-records` (`prism-viewer-records` for the app's
  key `prism-viewer`; `prism-default-records` when the spec has no key, §6; `recordsDatabase`; loaded when asked), and
  appends there. Viewers on a page with different keys keep and list their records apart (README decision 29). **Only
  the record**: no state, files or checkpoints; state comes back only from a session file.
- **Keep runs in the write chain**: one at a time, same entry ids and positions, appends wait for it; a failed copy is
  removed and the record stays in memory. `unkeep` copies back into memory and forgets the kept rows.
- **Index of kept records** (threads have no listing by owner): an append-only thread in the same rows, folded when
  listed; touched on keep, at most once a minute while appended to, at a turn's end, on dispose.
- **The index stays small by folding** (`createRecordIndex`, `IndexFoldAt` = 64): past that many entries, and past twice
  its records, the thread is rewritten as one entry per kept record (newest title and time, the time first listed); the
  superseded touches and the forgotten records' tombstones go. Done in the background after a touch or a list. **Folding
  never drops a record**: a kept record idle for 90 days stays listed, marked in Home, and only the person forgets it
  (§7); the number of records is bounded by `keptRecordsMax` where the person keeps them, so the index is that many
  entries plus at most the fold bound.
- **Fenced per database, not per index** (`fenced`): one viewer makes several indexes on one database (History's Keep,
  Home's list, the Storage view) and unkeyed viewers share one (§6), so every list, touch, forget and fold on it goes
  one at a time: a Web Lock per database across tabs (`navigator.locks`), else a chain per database in this page
  (another page without Web Locks is not fenced: its lost touch is the next one's to write; a thread folded away under
  it is made again); the host's `persist.rows` by their object. Nothing is read or appended mid-rewrite. **A fold
  survives a closed page**: it writes the folded index to a copy sealed by its last entry before it rewrites the thread,
  and the next use finishes it (a sealed copy is written again; an unsealed one is let go, the thread untouched). The
  record lets go of its kept rows only after its index's fold is done.

## 2. Entries

- **Envelope** `{ codec: 'prism-viewer.record.v1', version: 1, value }`, strict JSON; kinds in `RecordEntry`.
- **A run** (`RunEntry`) holds its operation, summarized `values`, `origin`, result, and what it touched: `targets` (the
  text of its `structure`, `component`, `node` and `target` values, from the runner's checked values, cut at 80
  characters, URLs bare: `targetsOf`) and `changed` (the result's `changed`, small: at most 8 node refs, `more` counts
  the rest, and the `selection`, `view`, `sequence` flags: `changedOf`; a failure has targets, no `changed`). Fields a
  reader of an older record lacks are optional.
- **Values summarized** (`summarizeValues`): non-JSON said in words, long text cut, URLs without credentials, query and
  fragment (tokens stay out of kept records and `.molx`), oversized values replaced by a note.
- **Append-only**: an undo is a new entry naming the step and its actor (`by`, the origin of the Undo, Undo turn or Redo
  run, as a run's `origin`; absent when the history is driven outside the API). **An Undo turn is one entry** naming
  every step it took back (`step` is the newest, `steps` all of them when more than one), written when its run ends;
  History reads all of them as undone. A redo is an entry too.
- **A turn under way** is a begun streaming entry (`beginTurn`, from `agent/session.ts`), kept alive by a five-minute
  heartbeat against the thread's inactivity timeout; settled with reply, outcome, usage. An aborted turn reads as
  stopped. The entry has the turn's `model` (`provider:model`) from the prompt on.
- **The reply streams**: the agent's `onText` (the text of the model call under way) goes to `TurnWriter.reply`, which
  writes thread deltas coalesced to one per `ReplyMs` (300 ms): `{ text }` with what was appended since, `{ replace }`
  when the text is not an extension of what was sent (a retried call). A read of an open turn folds its deltas into
  `reply` (incrementally: a read costs what was added; reads of one stream go in turn, so two readers at once, History
  and the assistant, never fold a delta twice; the preview is let go once the turn is not under way); the final text of
  `settle` replaces them (a settle with none keeps what had streamed). Not previewed: a reply over 32,000 characters,
  more than 1,000 deltas, a stream that cannot be written (the record was kept meanwhile; the settle carries the text).
  An aborted turn (the page went away) shows no reply: the thread does not give an ended stream's deltas back.
- **Reads are windows**, newest first, at most 256 entries per page.

## 3. Shown back

- **History** (`ui/history.tsx`, Steps beside the Log): newest first; a run is a row, an agent turn one group, undone
  steps struck. Undo the latest step; Undo turn asks first when later steps would go too. A run's second line says what
  it was pointed at and changed ("1CBS · A:10-20 · 3 nodes, view", `metaOf`) and "Undone by the assistant" when it was
  not the person who undid it; a turn shows its reply as it streams, and its model's name (the provider in the tooltip).
- **`journal/history-model.ts`** reads only while the tab shows: a window, "Show earlier", then once per frame only what
  came since (from the oldest open turn). Loaded rows are bounded; a row that did not change keeps its object (a turn's
  reply compared by its text: a retried call can keep its length).
- **The assistant's history** (`AgentModel.history`) maps History's rows; the turn card holds Undo turn while the newest
  step is that turn's.
- **After a reload** only kept records exist: Home lists them while it shows (`journal/kept-records.ts`); one opens
  read-only in History in place of this session's steps (`openKeptRecord`), never undone or replayed.

## 4. The agent's session

- `createViewerAgent` ([agent.md](agent.md) §1) runs an in-memory `AgentSession` that lives as long as the page. A
  record that ends settles an open turn as stopped. Compaction of the agent's journal never shortens the record.

## 5. Session files and kept records

- **Opt-in** `session.save` `record` (prompts and replies are in it): `.molx` gets `record.jsonl` (an envelope a line,
  `recordLines`), `.molj` a `record` field. Opening one shows it read-only in History (`recordFromLines`), not kept.
- **Bound** `options.ui.keptRecordsMax` (20): History's Keep at the bound asks before forgetting the oldest; the
  `record.keep` command never forgets one.
- **Forget**: Home's list, History, the Storage view (§7).

## 6. Light state in the browser's storage

- **Names are per viewer** (README decision 29): the host builds each store with the viewer's `persist.key` as its
  `prefix` (`webStore`; the app: `prism-viewer`), the records' IndexedDB is `<key>-records` (§1; a viewer with no key:
  `prism-default-records`, never the app's `prism-viewer-records`, so an embedded viewer on the app's origin cannot
  list, open or forget the app's kept records; **viewers with no key on one origin share one kept-records list**, safe
  as the index is fenced per database, §1), and the assistant's preferences (model, own endpoints, allowed domains) are
  `prism-assistant@<key>.*` (the app's key, `prism-viewer`: the original `prism-assistant.*`; no key: `prism-default`).
  **The assistant's keys are the person's and are shared by every viewer in the browser** (`prism-assistant.key.*`,
  IndexedDB `prism-assistant-keys`).
- **Two stores**, never for data (`persist/stores.ts`, `webStore('local' | 'session', { prefix, version, budget })`):
  `persist.preferences` (localStorage, else sessionStorage, else nothing; ~64 KB: the theme, the layout, the person's
  keys) and `persist.restore` (sessionStorage, ~256 KB: this tab's light restore). A value of another version or over
  budget is dropped; a blocked store fails nothing. Settings are not kept.
- **The layout** (`ui/layout.ts`, key `layout`; `persist.layout`, viewer-spec.md §2.6): `viewState.save()` (docks,
  sizes, float, pin, hidden tabs) and `navigation.save()`'s docks (their open tab and collapse; not the card or the
  page, which are summoned), written 500 ms after a change that differs from what was written, again on `pagehide` and
  on dispose (before the tabs unregister). Restored at start before the tabs register, with the host's parts (mode,
  compact rail, buttons) left as the spec has them; the start rule still decides whether the main panel opens. Either
  part of another version, or not a layout, is ignored.
- **The light restore** (`sources/restore.ts`): a second after the last step (longer after a slow write), the plugin's
  snapshot without data, when every data cell comes from an address; with a local file or inline data nothing is kept,
  said once.
- **Runs only on the host's call**; the app calls it on a reload (navigation type `reload`) when the URL loads nothing.
  The history starts empty.
- **Opt-in**: the app passes `persist.restore` only with `?restore` (`apps/prism-viewer/src/spec.ts`).
- **Decided against**: autosave into IndexedDB and durable agent sessions (no second file system).

## 7. Versions, collection, and the Storage view

- **Every stored value has a version** (README decision 32): the record's codec, the sequence's snapshot state, the
  light stores, a stored layout (`StoredViewLayout`). Another version's value is migrated or dropped.
- **Collected**: light stores drop old versions when read or listed; kept records idle 90 days are marked "forget it?".
- **Deploys**: a lazy chunk that fails says "A newer version of the viewer is available: reload" and retries on next
  use.
- **The Storage view** (`storage/model.ts`, `viewer.storage`; `ui/storage.tsx`, end of Settings): sources `register`;
  listed only while shown; the origin, Forget per entry, Forget everything (asked). Viewer sources in
  `storage/sources.ts` (kept records, preferences, light restore); the assistant's in `assistant/src/viewer/storage.ts`
  (keys by provider and scope, never the key; allowed domains).
