# The session record: what happened, kept and shown back

Yes to both (David, 2026-10-15): the viewer runs the agent on an `AgentSession` (`@molstar/ai-session`), and it keeps a
record of the session in a fiber (`@molstar/ai-fibers`): what the person did, what the agent did, in order, shown back
in History. **In memory by default**; the person keeps a record on purpose, and a kept record is only a record (§1).
The record is the viewer's composition of the AI packages' building blocks, not something the library does (packages/ai/\_spec/session.md, "The session keeps no person-facing history").
Status: **accepted** (2026-10-15).

## 1. One record per viewer session

- **A viewer session** starts when the viewer opens empty or with a URL, and ends when the person clears it or opens
  another session file. Its id is the record's fiber id; the owner is the viewer (one owner id per browser profile).
- **Who writes**: the history manager (README §3.3) appends a `run` entry for every API call it records, whoever made
  it; the agent's composition appends a `turn` entry per turn from the session's hooks and `TurnResult`. Both go through
  one small writer in the viewer (`packages/prism/viewer/src/journal/`), so the order in the record is the order things
  happened.
- **In memory by default** (`createMemoryRows`): the record lives as long as the page. A reload ends it, and with it
  every checkpoint, so an agent's turn can no longer be undone after a reload; that is accepted.
- **Kept on purpose**: "Keep this session's record" (History's menu, and a command) copies the record into IndexedDB
  (`createIndexedDbRows('prism-viewer-records')`, `@molstar/ai-components`) and keeps appending there for the rest of
  the session. **Only the record**: no plugin state, no files, no checkpoints, no autosave, so IndexedDB never becomes a
  second file system; the state comes back only from a session file the person saved.
- The spec chooses the rows a kept record goes to (`persist.rows`, viewer-spec.md §2.6; the record's `keepIn`); without
  them "Keep" uses the browser's IndexedDB as above. The record itself starts in memory either way; headless hosts and
  evals use memory.
- **The copy is in the write order**: "Keep" runs in the record's one chain of writes, so what is appended while it copies
  waits and goes to the kept record; one keep runs at a time (the button is off while it does). A copy that fails is
  removed, and the record goes on in memory.
- **Listing kept records** (review F19): fibers have no listing by owner, so the viewer keeps a small index of kept
  records (id, title, first and last time; size not yet) in the same IndexedDB rows, written when a record is kept,
  forgotten, and while it is appended to at most once a minute and at the end of each turn (and when the viewer goes);
  each write is an entry, folded when listed (folding the index fiber itself as it grows: not yet). A fiber listing in
  `@molstar/ai-fibers` would replace it (README §4). Copying the memory record into IndexedDB reads its windows and
  writes them in order with the same entry ids and positions: settled entries settled, a turn under way begun, an aborted
  one aborted.

## 2. Entries

Envelopes `{ codec: 'prism-viewer.record.v1', version: 1, value }` of strict JSON:

| Kind      | What it holds                                                                                                         | Written                          |
| --------- | --------------------------------------------------------------------------------------------------------------------- | -------------------------------- |
| `run`     | the operation id and title, its values (targets as MolQL JSON, targets.md §2), origin (person, agent with its turn, script, url, link), the result (`did`, `changed`, or the failure), time and duration | settled, when the run ends       |
| `turn`    | the prompt, the focus it was sent with, the reply (streamed as deltas while it arrives, then settled), the ids of its runs, the outcome (done, stopped, failed), usage and model | `begin` at the prompt, `settle` at the end |
| `undo`    | which step was undone (a run, or a turn), by whom                                                                     | settled                          |
| `session` | opened (from what: empty, URL, a file), cleared, saved (to which file)                                                | settled                          |

As built (`journal/record.ts`), not yet everything in the table: a `run` holds its values as summarized JSON (what JSON
cannot carry said in words, a function as `(function)`, a URL without its query or fragment, long text cut; MolQL
targets and `changed` not yet). A `turn` holds the prompt, the focus, the reply only when it settles (no deltas yet: History
shows the turn's runs while it runs, its reply at the end), the outcome and usage; the ids of its runs are not in it (its
runs name the turn instead), and the model not yet (the agent does not pass it). An `undo` names the step and the turn it
was part of, not yet by whom; an Undo turn writes one `undo` per step of the turn (one entry naming the turn: not yet). A
`session` entry holds the event and the API's words for it.

- **Append-only**: an undo is a new entry that names the step; the History view shows the step as undone. Nothing is
  edited in place.
- **Bounded** by the fiber's limits (entry 1 MB, windows of at most 256); a run's values over a size are kept as a
  summary, and a target's "what it selected" is kept as counts and chains past a few hundred residues (critic, scale).
- **A turn waiting** for an approval or an `ask` keeps its streaming entry alive (a heartbeat every five minutes while it
  is open, against the fiber's hour), so the fiber's inactivity timeout never aborts it (review F52). A turn the fiber
  aborts anyway (the page went away mid-turn) reads as stopped, from what it was begun with. A record that ends (the
  viewer disposed) settles a turn still open as stopped.
- **Reads** are windows (newest first), so the History view is a `VirtualList` over windows, however long the session.
  As built (`journal/history-model.ts`): the newest window, "Show earlier" for the one before; after a change, once per
  frame and only while the tab is shown, only what came since the last read (`record.readAfter`, from the oldest turn
  still under way, so its settling shows), and the newest window again when more came than a window holds (G8 perf).
  A plain list until a session's rows call for virtualizing, contained, its rows `content-visibility: auto`. A read
  longer than the fiber's window (256) goes in windows.

## 3. Shown back

- **History** (a core tab, replacing Log; Log's messages stay in the status column): the record, newest first. A
  person's runs one by one, an agent's turn as one item with its prompt, reply and runs under it, undone steps marked.
  Undo on the latest step (the history manager's), Undo turn on the agent's latest turn; with later steps above it,
  the row says how many and offers to undo them too.
- **The assistant's history** (`AgentModel.history`) is History's rows (`rowsOf`) mapped to its items, read only while
  something watches it; its turn card holds Undo turn while the newest step is that turn's and can be restored.
- **Without the assistant** History still shows the person's runs; the assistant adds the turn's rendering (the turn card,
  the Inspector's Trace for a turn still in the harness's storage).
- **After a reload**: only kept records are there (not built yet: the list on Home). Home lists them (title, first and last entries, when); opening one
  shows its record, read-only, beside whatever is loaded now. Nothing in it can be undone or replayed into the current
  state; it is a record.
- **Evals and replay**: the eval's run log is this record (memory rows); replaying it in the app reads the same entries
  (README §3.3, Evals).

## 4. The agent's session

- **`AgentSession`** is made by the viewer's `createViewerAgent` (agent.md §1), over the harness's memory storage
  (`@molstar/ai-components`), with the viewer's toolkit (agent.md §3), its drive's policies (approvals, agent.md §4),
  and the account and focus each turn (agent.md §5); the assistant and headless hosts both call it. Like the record, it lives as long as the page: a
  reload ends a running turn, and a kept record shows it as stopped.
- The session's journal (what the model sees, compacted) and the record (what the person sees) are separate: compaction
  never shortens the record.

## 5. Session files and kept records

Not built yet: neither the record in session files nor the bound below. Forget is built in the Storage view (§7), not
yet in Home's list or History's menu.

- **In `.molx`, opt-in** (David, 2026-10-15): saving a session offers "Include the session's record"; the record goes in
  as `record.jsonl` beside the snapshot (one entry a line, the envelopes of §2). Opening such a file shows the record in
  History, read-only, as a kept record is (§3). `.molj` (one JSON file) carries it as a `record` field. Prompts and
  replies are in it, so the option is off by default and says so.
- **A bound on kept records** (David, 2026-10-15): at most a number of kept records (a setting, 20 by default); keeping
  one more offers to forget the oldest. Each kept record has "Forget" (Home's list, History's menu), which removes its
  fiber.

## 6. Light state in the browser's storage

Beside the record, two small stores for state that is cheap to keep (David, 2026-10-15), each under a size budget, never
a place for data:

| Store                       | Keeps                                                                                       | Budget  |
| --------------------------- | ------------------------------------------------------------------------------------------- | ------- |
| localStorage (`preferences`) | the layout, the theme, the person's settings, recent servers and IDs                       | ~64 KB  |
| sessionStorage (`restore`)  | this tab's light restore: after a reload, the scene comes back when it is light enough      | ~256 KB |

- **The light restore** is written after each step (debounced, and skipped while its last write took long): the
  plugin's snapshot without data (the state tree's transforms and params, camera, canvas; the sequence's state), stamped
  with the app's version, when every data cell comes from a URL or a known service and the snapshot fits the budget.
- **When it runs** (review F18): only on a reload of the same page (the navigation's type is `reload`), only after the
  viewer is `ready`, and only when the URL asks for nothing to load (URL parameters win: a `?pdb=` reload loads that,
  not the restore). A restore written by another version of the app is dropped. The data is fetched again from the
  network (the byte cache lives in memory and does not survive the reload); the history starts empty.
- **When it does not fit** (a large scene, a local file, a drop): only what can come back is kept, the loads by URL and
  the canvas, and the viewer says once what will not come back ("Local files are not kept across a reload; save a
  session to keep them").
- **Per tab and short-lived**: sessionStorage ends with the tab, so nothing accumulates; localStorage holds only
  preferences.
- Not a record and not undo: the checkpoints stay in memory (§1).
- **As built** (`sources/restore.ts`, `persist/stores.ts`): `webStore('local' | 'session', { prefix, version, budget })`
  stamps each value with the app's version and drops another's or one over the budget; the restore is written a second
  after the last step (longer after a slow write), is not kept for a scene with a file or data given by value (said
  once), and the app runs it on a reload when the URL asks for nothing. The preferences keep the theme; the layout and
  settings are not kept yet.

## 7. Versions, collection, and what is stored

- **Every stored value has a version** (README decision 32): the record's codec (`prism-viewer.record.v1`), the
  sequence's state in snapshots (`version: 1`), preferences (`{ version, … }`), the light restore (the app's version),
  scrambled keys, domain rules. A value of another version is migrated where a migration is written, else dropped and
  said once.
- **Collected**: kept records past the bound (§5) or older than a setting (90 days by default) are offered for
  removal; preferences for parts and commands that no longer exist are dropped when read; the light restore ends with
  its tab.
- **Deploys** (critic): a tab open across a deploy keeps working: hashed chunks of the last few versions stay on the
  site (the publish adds files, it does not replace the branch), and a chunk that fails to load anyway marks its
  feature failed with "A newer version of the viewer is available: reload" (as built: the assistant's chunk, in the
  status column).
- **The Storage view** shows what is kept and where: kept records (count, size, age), preferences, the light restore,
  stored keys (by provider and scope, never the key), domain rules, the origin that holds them, each with Forget, and
  Forget everything. It sits beside the assistant's Trace, Tools and Limits, and in Settings when there is no
  assistant.
- **As built** (`storage/`, `ui/storage.tsx`): a DOM-free `StorageModel` (`viewer.storage`: `register(source)`,
  `sources`, `refresh`, `forget`, `forgetEverything`, `watch`) with the viewer's sources (`storage/sources.ts`): kept
  records from the index (title, size of their rows, age; this session's marked; those not written for 90 days marked
  "forget it?"; Forget removes the rows and the index entry, and this session's kept record goes back into memory with
  `record.unkeep()`, so History keeps showing it), the preferences and the light restore (the stores list what they
  keep, `ListedStore.entries()`, and drop another version's values while listing). The assistant registers its stored
  keys (`assistant/src/viewer/storage.ts`: provider and scope, never the key; Forget goes through the setup's
  `removeKey`). The view is a folded "Storage" section at the end of the Settings tab, always (the assistant's Trace,
  Tools and Limits are not built); it says the origin, lists only while open and shown, and asks before Forget
  everything. IndexedDB is opened only when the records' database exists. Not built: domain rules (the assistant keeps
  none; "Allow this domain" lasts for the page), the assistant's own choices in localStorage (`prism-assistant.setup`,
  `.providers`) and its per-browser scrambling key (IndexedDB `prism-assistant-keys`), the bound on kept records (§5),
  removal offered on its own (the age is only marked), and a refresh while the view stays open (a Keep in History shows
  when the view is opened again).

Decided against (David, 2026-10-15): autosaving state into IndexedDB, and durable agent sessions; both would make
IndexedDB a second file system. Light state goes to localStorage and sessionStorage instead (§6).
