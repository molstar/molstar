# The session record: what happened, kept and shown back

Yes to both (David, 2026-10-15): the viewer runs the agent on an `AgentSession` (`@molstar/ai-session`), and it keeps a
record of the session in a fiber (`@molstar/ai-fibers`): what the person did, what the agent did, in order, shown back
in History. **In memory by default**; the person keeps a record on purpose, and a kept record is only a record (§1).
The record is the viewer's composition of the AI packages' building
blocks, not something the library does (packages/ai/\_spec/session.md, "The session keeps no person-facing history").
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
- The spec chooses the rows (`persist.rows`, viewer-spec.md §2.6); headless hosts and evals use memory.

## 2. Entries

Envelopes `{ codec: 'prism-viewer.record.v1', version: 1, value }` of strict JSON:

| Kind      | What it holds                                                                                                         | Written                          |
| --------- | --------------------------------------------------------------------------------------------------------------------- | -------------------------------- |
| `run`     | the operation id and title, its values (targets as MolQL JSON, targets.md §2), origin (person, agent with its turn, script, url, link), the result (`did`, `changed`, or the failure), time and duration | settled, when the run ends       |
| `turn`    | the prompt, the focus it was sent with, the reply (streamed as deltas while it arrives, then settled), the ids of its runs, the outcome (done, stopped, failed), usage and model | `begin` at the prompt, `settle` at the end |
| `undo`    | which step was undone (a run, or a turn), by whom                                                                     | settled                          |
| `session` | opened (from what: empty, URL, a file), cleared, saved (to which file)                                                | settled                          |

- **Append-only**: an undo is a new entry that names the step; the History view shows the step as undone. Nothing is
  edited in place.
- **Bounded** by the fiber's limits (entry 1 MB, windows of at most 256); a run's values over a size are kept as a
  summary.
- **Reads** are windows (newest first), so the History view is a `VirtualList` over windows, however long the session.

## 3. Shown back

- **History** (a core tab, replacing Log; Log's messages stay in the status column): the record, newest first. A
  person's runs one by one, an agent's turn as one item with its prompt, reply and runs under it, undone steps marked.
  Undo on the latest step (the history manager's), Undo turn on the agent's latest turn.
- **Without the assistant** History still shows the person's runs; the assistant adds the turn's rendering (the turn card,
  the Inspector's Trace for a turn still in the harness's storage).
- **After a reload**: only kept records are there. Home lists them (title, first and last entries, when); opening one
  shows its record, read-only, beside whatever is loaded now. Nothing in it can be undone or replayed into the current
  state; it is a record.
- **Evals and replay**: the eval's run log is this record (memory rows); replaying it in the app reads the same entries
  (README §3.3, Evals).

## 4. The agent's session

- **`AgentSession`** is created by the assistant's composition (`@molstar/prism-assistant`), over the harness's memory
  storage (`@molstar/ai-components`), with the viewer's toolkit (agent.md §3), its drive's policies (approvals,
  agent.md §4), and the account and focus each turn (agent.md §5). Like the record, it lives as long as the page: a
  reload ends a running turn, and a kept record shows it as stopped.
- The session's journal (what the model sees, compacted) and the record (what the person sees) are separate: compaction
  never shortens the record.

## 6. Light state in the browser's storage

Beside the record, two small stores for state that is cheap to keep (David, 2026-10-15), each under a size budget, never
a place for data:

| Store                       | Keeps                                                                                       | Budget  |
| --------------------------- | ------------------------------------------------------------------------------------------- | ------- |
| localStorage (`preferences`) | the layout, the theme, the person's settings, recent servers and IDs                       | ~64 KB  |
| sessionStorage (`restore`)  | this tab's light restore: after a reload, the scene comes back when it is light enough      | ~256 KB |

- **The light restore** is written after each step (debounced): the plugin's snapshot without data (the state tree's
  transforms and params, camera, canvas; the sequence's state), when every data cell comes from a URL or a known service
  and the snapshot fits the budget. After a reload the viewer reloads it: the data is fetched again (from the byte cache
  where it is still held), the history starts empty.
- **When it does not fit** (a large scene, a local file, a drop): only what can come back is kept, the loads by URL and
  the canvas, and the viewer says once what will not come back ("Local files are not kept across a reload; save a
  session to keep them").
- **Per tab and short-lived**: sessionStorage ends with the tab, so nothing accumulates; localStorage holds only
  preferences.
- Not a record and not undo: the checkpoints stay in memory (§1).

## 5. Session files and kept records

- **In `.molx`, opt-in** (David, 2026-10-15): saving a session offers "Include the session's record"; the record goes in
  as `record.jsonl` beside the snapshot (one entry a line, the envelopes of §2). Opening such a file shows the record in
  History, read-only, as a kept record is (§3). `.molj` (one JSON file) carries it as a `record` field. Prompts and
  replies are in it, so the option is off by default and says so.
- **A bound on kept records** (David, 2026-10-15): at most a number of kept records (a setting, 20 by default); keeping
  one more offers to forget the oldest. Each kept record has "Forget" (Home's list, History's menu), which removes its
  fiber.

Decided against (David, 2026-10-15): autosaving state into IndexedDB, and durable agent sessions; both would make
IndexedDB a second file system. Light state goes to localStorage and sessionStorage instead (§6).
