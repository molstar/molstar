# Agents: the Mol\* UI client API

An agent (an in-app assistant, or an outside one over MCP) controls Mol\* through the same commands people use, sees a
compact account of what is loaded and what the person points at, and every change it makes is visible, attributed, and
undoable in one step. The harness (`@molstar/ai-*`, [../../ai/](../../ai/)) runs the model loop and calls tools; it
knows nothing about UI. This document is the Mol\* side: what the harness is given and how the UI shows an agent at
work.

## 1. Principles

- **Tools are commands.** The command registry ([commands.md](commands.md)) is the tool registry: no second API to
  maintain, and an agent can do nothing a person could not type in the command menu. Panels with custom UI (§2.1 of
  commands.md) are driven by the commands and actions behind them, never by clicking the DOM.
- **The agent sees what the person sees.** Every prompt carries a short state account and the person's focus (§3);
  details are read on demand with read-only tools.
- **Visible and attributed.** Each run shows in the history with its source (you, the assistant, an outside agent, a
  script); the view moves while the agent works, and the person can stop it at any time.
- **One step back.** Before an agent changes anything, a checkpoint is taken; the whole turn undoes at once, back to the
  checkpoint (state, camera, selection), however many commands it ran. Only while nothing has come after it: the next
  action (the person's command, another turn) releases the checkpoint, so at most one snapshot is ever held (no growing
  memory, and no undo that would silently discard later work).
- **The person approves what leaves the page or cannot be undone** (downloads from a server, remote state upload,
  closing a session): such commands wait for Allow or Deny.

## 2. Commands as tools

| Command (`CommandDefinition`) | Tool (`ToolDefinition` of any-agent)                                                 |
| ----------------------------- | ------------------------------------------------------------------------------------ |
| `id`                          | `name` (`structure.component.add`)                                                   |
| `title`, `alias`, `category`  | `description` ("Add component… Also: create component. Structure.")                  |
| `params(ctx, values)` (PD)    | `inputSchema`: JSON Schema generated from PD (formats become `format` + description) |
| `when(ctx)`                   | listed only when available; the reason is returned if called anyway                  |
| `check(args, ctx)`            | validation before run (the harness does not validate)                                |
| `run(args, ctx, signal)`      | `invoke`: progress frames from the command's task, then the result                   |
| `CommandResult`               | result parts: what was done, what changed (nodes, selection, view), undoable         |
| `opens` (custom UI, planned)  | not a tool; its panel's actions are commands of their own                            |

Commands add two fields for agents:

- `read?: true`: no effect (reads run together and need no checkpoint): `scene.get`, `selection.get`,
  `structure.sequence.get`, `snapshot.list`.
- `approval?: true`: the run waits for the person (the harness's `permission` decision parks the turn).

**Params as input** (`paramsShape(params, { formats })`, `@molstar/prism-assistant/tools/params`): PD params become
`ai-tools` field shapes, every field optional (a param has a default), and the arguments become the params' values over
the defaults (`toValues`; `toArgs` reads values back). Select: a string enum, labels in the description; multi-select: a
list of them; boolean, text (multiline: verbatim), numeric (integer when its step is), interval, vec3, mat4 (number
lists); color: a name or hex; color list; url; group: an object; mapped: a variant on `name` with its `params`;
conditioned: a variant on `condition` with its `value`; converted: its displayed form; object list: a list of objects;
refs: a node ref; script: language and expression. Files, line graphs, opaque values and hidden params are left out and
keep their defaults.

Generated commands (presets, themes, queries, state actions) are tools too; there are hundreds, so they are grouped and
_deferred_ (declared by group, loaded when the agent asks for the group), as Stories loads its tools by skill.

## 3. What the agent sees

Each prompt starts with:

- **State**: the loaded structures and volumes as an outline (id, label, components with their representations and
  colors, hidden items marked), the snapshots, the animation, the canvas (background, illumination), capped like
  Stories' `stateBlock` (at most 40 entries; the rest by `scene.get`).
- **Focus**: what the person points at right now and never stored: the selection (as residues per chain and as a
  MolScript expression), the focused or picked node, the open panel and section, the hovered residue.

Details come from read tools; a repeated read of unchanged state is collapsed (Stories' `RepeatedReads`).

## 4. The UI client API

The surface the harness binds to (in-page for the built-in assistant; behind an MCP server for outside agents):

```ts
interface MolstarUIClient {
  /** The tools: commands as ToolDefinitions, grouped; deferred groups listed by name. */
  tools(filter?: { groups?: string[] }): readonly ToolDefinition[];
  /** Runs a command for an agent: attributed, checkpointed per turn, approval-aware; frames as any-harness expects. */
  invoke(call: { id: string; args: Json; turn: string; source: CommandSource }, signal: AbortSignal): AsyncIterable<ToolFrame>;
  /** The state account and the focus (§3). */
  state(): StateAccount;
  focus(): FocusAccount;
  /** Turns: a checkpoint before the first change; end closes the undo step; restore returns to the checkpoint. */
  beginTurn(turn: { id: string; prompt: string; source: CommandSource }): void;
  endTurn(id: string, outcome: 'done' | 'stopped' | 'failed'): void;
  restore(turnId: string): Promise<void>;
  /** What happens, for observers (the history, an outside agent's view): runs, outcomes, focus changes. */
  readonly events: Observable<ClientEvent>;
}

type CommandSource = { kind: 'person' } | { kind: 'assistant' } | { kind: 'agent'; name: string } | { kind: 'script'; name?: string };
```

`CommandSource` becomes core's `RunSource` (contexts.md §5: person, agent, mcp, behavior, link; `via: 'assistant'` for a
run the person confirms). The checkpoint is the mode's (`Checkpointer`, contexts.md §5): in the viewer a Mol\* state
snapshot (data, camera, selection, components), so restoring is exact even for changes the state tree's undo does not
cover; in the story mode the document's inverse operations.

The views bind to `AgentModel` (`@molstar/prism-assistant`, `model`), which the client will back: `turn` (the latest
`AgentTurn`: its prompt, status, streamed `text`, `runs`, `checkpoint` held or not, `undone`, `retries`, `usage`),
`view` (`'card' | 'line' | 'none'`), `history` (`HistoryItem`: a person's `CommandRun` or an `AgentTurn`), `ask`,
`stop`, `approve`, `undo`, `dismiss`, `hide`, `reopen`, and `settings` (`AISettingsModel`: providers and keys, the
model, refresh, own endpoints). A run's status (`RunStatus`) is `running` where the turn's (`TurnStatus`) is `thinking`
or `working`, and adds `denied`.

## 5. In the UI

- **The command menu asks** (`views/command-menu`, the `assistant` port, `CommandMenuAssistant`: `ready`, `propose`,
  `ask`, `setup`). Under the list, always there, an _Ask the assistant_ bar in the accent's soft wash: it shows what is
  typed (with nothing typed, _describe what you want, in your words_, and a click focuses the input) and takes it with
  ⌘↵ (Ctrl+↵; the footer lists it), or with ↵ when no command matches (the bar is then the default, outlined, and the
  empty list reads _No command matches: ask the assistant._). Words after a prefix (`@`, `#`, `:`, `?`) are not asked.
  The assistant (any-ui's `assist`: a one-off agent) proposes one command when the words fit it, read into its values:
  the menu opens the ordinary inline sentence filled in, over it _Proposed by the assistant from “…”: check it, then
  Enter to run_; nothing runs until Enter. Words that do not fit one command go to the agent as a turn (in the mock,
  _color it red_ proposes `selection.color`; _color everything red_ is a turn). Without a model, the bar is quiet (no
  tint, no shortcut) and reads _Set up an AI assistant_; it, and ⌘↵, open the Assistant tab.
- **An agent turn** (the harness's progress) shows in the status column at the bottom left of the canvas, under its
  tasks (the shell's `agent` slot), as a card that asking opens: the request and its state (Thinking…, Working…, Waiting
  for you, Done, Stopped, Failed), the model's words as they stream, a note when the model was retried ("the provider
  was busy"), and each command as it runs, in words (its title with its arguments, not the tool name): a spinner, its
  progress as a thin line when the command reports one, its outcome or error under it, its duration; reads quieter. Rows
  follow their run by id, so a spinner does not restart as the records change. Stop while it runs. A command waiting for
  approval shows what it does, why it asks ("Fetches data from www.ebi.ac.uk"), Deny and Allow; Deny marks it denied and
  the model is told: it goes on another way, or says it cannot and stops (the mock still ends the turn there). At the
  end, one line (commands that change something, failures, time, tokens in and out) and _Undo turn_ while the checkpoint
  is held (_Undone_ after). Esc (unless a field, list, menu or dialog has it) or ✕ closes the card to a line of the
  status column, like a task's (`StatusPill`: the request and its state, the spark pulsing while it thinks or works, the
  state in the accent while it waits; a running turn goes on): a click opens the card again, its ✕ takes the line away
  (the turn stays in History).
- **History** (with the agent module, the Log tab at the rail's foot becomes History): every command, newest first; the
  person's one per row with how long ago, an agent's grouped under its request (by the agent's name) with its answer,
  its commands listed while it runs and folded under a one-line summary once done (folds kept as records change), and
  the checkpoint: _Checkpoint before it_ with _Undo turn_, or _Undone: back to how it was before_, or _No longer
  undoable: something came after it_, or _Nothing to undo_ (only reads). Messages (the old log) one switch away.
- **The Assistant tab**, at the rail's foot by Settings (the rail says _Set up assistant_; as Stories untold sets up
  any-model-adapters): first the model in use (a select of the models of providers with a key; a refresh beside it asks
  those providers' `/models` for what they offer now; it applies from the next request), then _Keys_: one list of
  providers (Anthropic, Gemini, OpenAI, and one's own), one key each, for as many as one likes. Each shows its key's
  state (no key, until reload, for this session, kept in this browser, from the server) and Add key / Change (Cancel
  while open) / Remove; a key from the server has none of these. Entering a key also chooses how long it is kept (_Until
  reload_: memory only; _This session_; _This browser_, the default). With no key yet, a word on what the assistant is
  replaces the model and the first provider's key is open. Under the list, _Your own endpoint_ folds open: an
  OpenAI-compatible server (name, URL (https, or http on localhost), OpenAI chat completions or OpenAI Responses, an
  optional key and how long to keep it, its models typed or fetched from its `/models`; OpenRouter, Ollama, vLLM,
  LiteLLM; it must allow CORS, and a failed fetch says to check the endpoint, the key and CORS). Keys go only to their
  provider; kept ones are stored scrambled (encrypted with a key this browser keeps, as Stories does), so the storage
  read as it is gives no working key, which the tab says, along with what it does not guard against (code running in the
  page). With a server-side agent, keys come from its environment instead.
- **The Inspector** (`@molstar/prism-assistant/inspector`, a separate entry: the tab and the card do not pull it in), a
  page over the whole app (`ShellPage`) with three tabs over `@molstar/ai-trace`'s data:
  - _Trace_: an overview strip (a column per turn, a block per step; model and tool calls told apart by shape and lane,
    failed in danger, a retry hatched, a running step pulsing, the selected one outlined in the accent; Steps or Time),
    the rows (turns, requests, model and tool calls, nested runs indented; j/k move), and the detail: a turn's facts and
    its input per model call, a model call's _What it saw_ (the window's makeup, its share, what is new, why a cache
    prefix broke, the entries), _Answer_ and _Timing_ (a waterfall per attempt), a tool call's _Result_, _Arguments_,
    _Definition_ (a link to Tools) and _Timing_.
  - _Tools_: the groups (core, each skill, single; loaded, unavailable here and why) and their tools (reads, asks,
    disabled; loaded at which turn; size in tokens), a tool's facts and its input as fields or JSON Schema, a group's
    guidance as the model reads it.
  - _Limits_: the turn budget per kind (used against the value, up to its ceiling; editable, with the default and
    Reset), the model's context window, compaction trigger and thinking, credits as information (never a limit).

  The header has the tabs, working or idle, and the budget kind nearest its ceiling. Ways in: the commands _Agent trace_
  (⌥T), _Agent tools_, _Agent limits_; the Assistant tab's _Inspect_ block (Trace, Tools, Limits with their counts); the
  turn card's _Trace_ link (at that turn).

- While an agent runs, the person's own commands still run (they are attributed to the person); a command of theirs
  after the turn releases its checkpoint, and History says so.

**State (mocks, `@molstar/prism-assistant`: `views/`, `inspector/`, `testing/mock-agent.ts`,
`testing/mock-inspector.ts`; Shell › Assistant: Set up, Not set up, Trace, Trace running, Trace what it saw, Tools,
Tools at a tool, Limits, Inspect):** the Ask bar with proposed commands (`selection.color`, `data.download-structure`,
`focus.target`) and agent turns, the turn card and its line (streaming, steps with progress and durations, a retry note,
approval and denial, Stop, failure, summary, Undo turn, one checkpoint released by the next action), History with turns
and checkpoints, the Assistant tab (key needed; providers, model and refresh, own endpoint). Not yet: the `read` and
`approval` command fields (they exist only on `CommandRun`), the PD → JSON Schema generator, the state and focus
accounts, the client API, the harness wiring, an answer in a line from the menu's assistant (the port proposes a command
or asks the agent), keeping keys (the mock holds them in memory whatever the scope; encrypted storage is to come) and
keys from a server, the MCP server (a separate piece of work).

## 6. Opt-in: not in every build

The agent is a package an app includes, `@molstar/prism-assistant`, not part of the UI: prism-ui does not depend on it
and none of its files name it (`src/_test/opt-in-boundary.test.ts`), and the UI knows it only through ports: the command
menu's `assistant` (`CommandMenuAssistant`), the shell's `agent` slot (`ShellViewProps.agent`, under the status column's
tasks), and tabs and commands registered like any contribution (the agent's `mockAgentAddOn` adds the Assistant tab and
swaps Log for History in the mock shell). An app that does not include it ships none of it. Planned, not built: one that
does registers it as a contribution (`agentContribution({ … })`), with the harness and the model adapters
(`@molstar/ai-*`) as its own dependencies, and the spec can still start it hidden (`ai: { enabled: false }`) for apps
that include it but turn it on per user. Each mode's agent parts (its skills and guidance, its state and focus accounts,
its drive overrides) are an `agent` folder of the mode's package (`@molstar/prism-viewer/agent`,
`@molstar/prism-stories/agent`), the only code there that imports the assistant or `@molstar/ai-*`, imported only by
apps with the assistant; an extension brings its skills the same way
([../../ai/analysis/viewer-mode.md](../../ai/analysis/viewer-mode.md) §2.1, §8).

## 7. Skills

The viewer mode's tools and skills (selections first), the generic tool layer and session, the Inspector (trace, tools,
limits) and fibers are analysed in [../../ai/analysis/viewer-mode.md](../../ai/analysis/viewer-mode.md).
