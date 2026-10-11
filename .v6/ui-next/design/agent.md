# Agents: the Mol\* UI client API

An agent (an in-app assistant, or an outside one over MCP) controls Mol\* through the same commands people use, sees a
compact account of what is loaded and what the person points at, and every change it makes is visible, attributed, and
undoable in one step. The harness (any-harness, moving into the Mol\* repo) runs the model loop and calls tools; it
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

| Command (`UICommand`)        | Tool (`ToolDefinition` of any-agent)                                                 |
| ---------------------------- | ------------------------------------------------------------------------------------ |
| `id`                         | `name` (`structure.component.add`)                                                   |
| `title`, `alias`, `category` | `description` ("Add component… Also: create component. Structure.")                  |
| `params(ctx)` (PD)           | `inputSchema`: JSON Schema generated from PD (formats become `format` + description) |
| `when(ctx)`                  | listed only when available; the reason is returned if called anyway                  |
| `check(values)`              | validation before run (the harness does not validate)                                |
| `run(values, ctx, signal)`   | `invoke`: progress frames from the command's task, then the result                   |
| `CommandResult`              | result parts: what was done, what changed (nodes, selection, view), undoable         |
| `opens` (custom UI)          | not a tool; its panel's actions are commands of their own                            |

Commands add two fields for agents:

- `read?: true`: no effect (reads run together and need no checkpoint): `scene.get`, `selection.get`,
  `structure.sequence.get`, `snapshot.list`.
- `approval?: true`: the run waits for the person (the harness's `permission` decision parks the turn).

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

type CommandSource = { kind: 'person' } | { kind: 'assistant' } | { kind: 'agent'; name: string } | { kind: 'script' };
```

The checkpoint is a Mol\* state snapshot (data, camera, selection, components), so restoring is exact even for changes
the state tree's undo does not cover (camera, selection, focus).

## 5. In the UI

- **The command menu asks.** Under the list, always there, an _Ask the assistant_ bar in the accent's soft wash: it
  shows what is typed and takes it with ⌘↵ (Ctrl+↵), or with ↵ when no command matches (it is then the default, and the
  list says so). The assistant (any-ui's `assist`: a one-off agent) proposes a command with its values, and the menu
  opens the ordinary inline sentence, filled in; nothing runs until Enter. It can also answer in a line, or hand the
  request to the agent when it takes several steps. Without a model, the bar is quiet (no tint, no shortcut) and reads
  _Set up an AI assistant_; it opens the Assistant tab.
- **An agent turn** (the harness's progress) shows as a card in the status column at the bottom left of the canvas,
  under its tasks: the request and its state (Thinking…, Working…, Waiting for you, Done, Stopped, Failed), the model's
  words as they stream, a note when the model was retried, and each command as it runs, in words (its title with its
  arguments, not the tool name): a spinner, its progress as a thin line when the command reports one, its outcome or
  error under it, its duration; reads quieter. Stop while it runs. A command waiting for approval shows what it does,
  why it asks ("Fetches data from www.ebi.ac.uk"), Deny and Allow. At the end, one line (commands, failures, time,
  tokens in and out) and _Undo turn_ (back to the checkpoint). Esc or ✕ closes it to a line of the status column, like a
  task's (the request and its state, a pulsing spark while it runs; a running turn goes on): a click opens the card
  again, its ✕ takes the line away (the turn stays in History).
- **History** (the Log tab grows into it): every command, newest first; the person's one per row, an agent's grouped
  under its request with its answer, its commands folded under a one-line summary, and the checkpoint before it with
  _Undo turn_; Messages (the old log) one switch away.
- **The Assistant tab**, at the rail's foot by Settings (as Stories untold sets up any-model-adapters): one list of
  providers (Anthropic, Gemini, OpenAI, and one's own), each with its key (as many as one likes): its state (none, until
  reload, for this session, kept in this browser, from the server), Add key / Change / Remove; entering a key also
  chooses how long it is kept (_Until reload_: memory only; _This session_; _This browser_). With no key yet, a word on
  what the assistant is comes first and the first provider's key is open. Then the model in use (only providers with a
  key; a refresh beside it asks the providers' `/models` for what they offer now), and one's own OpenAI-compatible
  endpoint (name, URL, OpenAI chat completions or OpenAI Responses, an optional key and how long to keep it, its models
  typed or fetched from its `/models`; OpenRouter, Ollama, vLLM, LiteLLM; the server must allow CORS). Keys go only to
  their provider; kept ones are stored scrambled (encrypted with a key this browser keeps, as Stories does), so the
  storage read as it is gives no working key, which the tab says, along with what it does not guard against (code
  running in the page). With a server-side agent, keys come from its environment instead.
- While an agent runs, the person's own commands still run (they are attributed to the person); a command of theirs
  after the turn releases its checkpoint, and History says so ("No longer undoable: something came after it").

**State (mocks, `views/agent`, `_stories/mocks/agent.ts`, Shell › Assistant and Assistant Setup):** the Ask row with
proposed commands (color, download, focus) and agent turns, the turn card (streaming, steps with progress and durations,
a retry note, approval, Stop, failure, summary, Undo turn), History with turns and checkpoints, the Assistant tab (key
needed; providers, model, own endpoint). Not yet: the PD → JSON Schema generator, the state and focus accounts, the
client API, the harness wiring, the MCP server (a separate piece of work).

## 6. Opt-in: not in every build

The agent is a module an app includes, not part of the core: nothing in the package imports `views/agent` or
`model/agent` (a test enforces it), and the core knows it only through ports: the command menu's `assistant`
(`CommandMenuAssistant`), the canvas's `agent` slot, and tabs and commands registered like any contribution. An app that
does not include it ships none of it; one that does registers it as a contribution (`agentContribution({ … })`), with
the harness and the model adapters as its own dependencies (any-harness moves into the Mol\* repo as a package of its
own). The spec can still start it hidden (`ai: { enabled: false }`) for apps that include it but turn it on per user.

## 7. Later: skills

Stories untold's MolQL selection skills (how to write a selection from words: residues, chains, ligands, "around",
"within", secondary structure) come over as the agent's selection skills, with the selection commands as their tools. A
separate slice, after the UI is done.
