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
  checkpoint (state, camera, selection), however many commands it ran.
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

- **The command menu asks.** Text that matches no command (or starts with `?`) shows an _Ask_ row; the assistant
  (any-ui's `assist`: a one-off agent) proposes a command with its values, and the menu shows the ordinary inline
  sentence, filled in and marked as proposed; nothing runs until Enter. It can also answer in a line, or hand the
  request to the agent when it takes several steps.
- **An agent turn** shows as a card at the bottom of the canvas: the request, each command as it runs (with progress),
  Stop; at the end, its outcome and _Undo turn_ (back to the checkpoint). A command waiting for approval shows Allow /
  Deny in the card.
- **The history** (the Log tab grows into History): every run, grouped by turn, with its source; checkpoints before
  agent turns, each with _Restore_; a failed run with its error.
- While an agent runs, the person's own commands still run (they are attributed to the person); restoring a turn's
  checkpoint also discards what the person did after it, and says so.

**State (mocks):** the Ask row and a proposed command, an agent turn card (steps, Stop, approval, Undo turn), the
History tab with sources and checkpoints. Not yet: the PD → JSON Schema generator, the state and focus accounts, the
client API, the harness wiring, the MCP server.

## 6. Later: skills

Stories untold's MolQL selection skills (how to write a selection from words: residues, chains, ligands, "around",
"within", secondary structure) come over as the agent's selection skills, with the selection commands as their tools. A
separate slice, after the UI is done.
