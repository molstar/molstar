# The session record and what the viewer stores

The session's record (an `@molstar/ai` thread), History, kept records, and light state in browser storage. The code is
the truth for how.

## 1. One record per viewer session

- **One record, one writer** (`journal/record.ts`, `createSessionRecord`): `run` per finished API call (any origin),
  `turn` per agent turn, `undo`/`redo` from the history, `session` for open, clear, save. One write chain: record order
  is the order things happened. Not recorded: reads, `record: false` operations (pointer-rate view state),
  `highlight.*`, `view.*`, `export.preview`, nested runs.
- **In memory by default.** A reload ends it and every checkpoint: a turn cannot be undone after a reload.
- **Kept on purpose**: "Keep this session's record" (History, `record.keep`) copies it into the spec's `persist.rows`
  ([viewer-spec.md](viewer-spec.md) §2.6), else IndexedDB `prism-viewer-records` (loaded when asked), and appends there.
  **Only the record**: no state, files or checkpoints; state comes back only from a session file.
- **Keep runs in the write chain**: one at a time, same entry ids and positions, appends wait for it; a failed copy is
  removed and the record stays in memory. `unkeep` copies back into memory and forgets the kept rows.
- **Index of kept records** (threads have no listing by owner): an append-only thread in the same rows, folded when
  listed; touched on keep, at most once a minute while appended to, at a turn's end, on dispose.

## 2. Entries

- **Envelope** `{ codec: 'prism-viewer.record.v1', version: 1, value }`, strict JSON; kinds in `RecordEntry`.
- **Values summarized** (`summarizeValues`): non-JSON said in words, long text cut, URLs without credentials, query and
  fragment (tokens stay out of kept records and `.molx`), oversized values replaced by a note.
- **Append-only**: an undo is a new entry naming the step (one per step of an Undo turn).
- **A turn under way** is a begun streaming entry (`beginTurn`, from `agent/session.ts`), kept alive by a five-minute
  heartbeat against the thread's inactivity timeout; settled with reply, outcome, usage. An aborted turn reads as
  stopped.
- **Reads are windows**, newest first, at most 256 entries per page.

## 3. Shown back

- **History** (`ui/history.tsx`, Steps beside the Log): newest first; a run is a row, an agent turn one group, undone
  steps struck. Undo the latest step; Undo turn asks first when later steps would go too.
- **`journal/history-model.ts`** reads only while the tab shows: a window, "Show earlier", then once per frame only what
  came since (from the oldest open turn). Loaded rows are bounded.
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

- **Two stores**, never for data (`persist/stores.ts`, `webStore('local' | 'session', { prefix, version, budget })`):
  `persist.preferences` (localStorage, ~64 KB: the theme) and `persist.restore` (sessionStorage, ~256 KB: this tab's
  light restore). A value of another version or over budget is dropped; a blocked store fails nothing.
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
