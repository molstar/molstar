# The agent in the viewer

How an agent drives the viewer (`src/agent/`): tools over operations, skills, what it sees each turn, what it may do
without asking, and running without a page. Built on `@molstar/ai` (`packages/ai/_spec`). The person's side (setup,
keys, turn card, Inspector) is the assistant's: [assistant/\_spec](../../assistant/_spec/README.md). The code is the
truth for how.

## 1. Where it lives

- **Tools, skills, session: the viewer's `agent/` folder**, DOM-free, on `@molstar/ai`. The viewer itself never imports
  it (the AI SDK, zod); the assistant's chunk brings it ([viewer-spec.md](viewer-spec.md) §2.4).
- **The declarations are plain data** (`agent/skills.ts`: `CoreTools`, `Skills`, `withChange`; no imports), so a spec
  reads or extends them without the toolkit. The spec's `agent.core`, `agent.skills` change the defaults by id
  (`{ add, remove }`, an added id of a default's replaces it); `agent: false` means no tools. Extensions bring skills
  (`ViewerExtension.skills`, a lazy one's in its stub), never core tools.
- **PD params become tool fields** through `paramsShape` (`agent/params.ts`).
- **`createViewerAgent({ viewer, models, guidance?, rules? })`** (`agent/session.ts`) makes the `AgentSession` and keeps
  it: `turn`, `resume({ decisions })`, `stop`, `busy`, `parked`, `policy`, `lastTurn`, and `inspect` (client, budget,
  window: read-only, for the Inspector). The assistant, evals and an MCP server call it alike.
- **The assistant is a factory** (`AssistantFactory = (host: AssistantHost) => AssistantUI`, `src/spec.ts`): the host
  gives the viewer and `createAgent(models, { rules })`; models come from the person's keys, which the assistant keeps.
  The Inspector is the assistant's ([assistant/\_spec](../../assistant/_spec/README.md), "The Inspector").
- **A turn at a time**: `turn()` throws while one runs or is `parked` (waits for an approval); `stop()` ends either, a
  parked turn settled as stopped (record entry, history group).
- **One run path**: tools run `viewer.api.as({ kind: 'agent', turn })`; the turn is one history group (Undo turn undoes
  it at once; [sequence.md](sequence.md) §7.3) and one record entry ([journal.md](journal.md) §2), which has the turn's
  model (`provider:model`) from the prompt on and the reply as the model says it (`onText` to `TurnWriter.reply`).

## 2. Tools

- **Thin adapters over operations** (`operationTool`, `agent/tools.ts`): flat fields in the model's words, normalized
  and checked (`paramsShape`), mapped to the operation's values; the result in one order: did, changed, warnings, next.
- **Words are the viewer's, not the screen's**: structures, chains, residues, the selection, the view; never tabs or
  panes (`ViewerGuidance`). A view-only operation's tool says the structures stay as they are. Messages read by people
  and the agent alike are plain; the tools add their own words.
- **Registry values are the field's choices as they are now**; between turns `toolkit.refresh()` describes the tools
  again and replaces the parts whose definitions changed (a theme registered, a feature loaded or gone).
- **What a tool says it cannot do** is what this viewer cannot do at all (`supported`: "Not in this viewer: …"); what
  cannot be done now is the call's answer, so descriptions do not change with the scene.
- **Targets are strings** to the model (a `TargetParam` is a text field; [targets.md](targets.md)); the shared params
  (`target`, `structure`, `language`) are said once in the guidance, not per schema.
- **Tokens are kept small**: a select's names are its enum only; one over 16 options has no enum, a wrong name answers
  with all of them (`LongSelect`, `params.ts`).

## 3. The core and the skills

- **Core** (from the first turn): `disclose_skills`, `disclose_tools` (the toolkit's), `get_scene`, `undo_turn` (its
  previous turn, never the person's later steps), `ask`, and `CoreTools` (`get_selection`, `load_entry`, `load_url`,
  `select`, `clear_selection`, `focus`, `apply_preset`, `apply_style`, `color`, `set_background`, `reset_camera`).
- **Skills** (disclosed on request, `Skills`): selections, styles, canvas, camera, measure, superpose, snapshots,
  animation, export; extensions' skills over their own operations.
- **The core stays small**: a new tool goes into a skill unless nearly every turn needs it.
- **A tool whose operation the viewer lacks is left out**, and a skill with none is left out.
- **A lazy skill** (its `feature` a stub extension) is in the catalog by its summary with its tools' approval flags
  (`pending`, by `mayAsk`); disclosing it loads the extension, then its tools.
- **`ask`**: when the intent is unclear the agent asks; the tool tells it to end its turn, and the answer is the next
  request.
- **Read operations** (`def.read`) become read tools (`read: true`, idempotent) where a skill declares them.

## 4. Policies

As little in the way as is reasonable (README decision 26); decided at call time on the values the operation will run
with (`agent/policies.ts`, the tool's `approvedBy`).

- **Preauthorized**: changes to the viewer's state while the turn can be undone. With undo off (checkpointing off for a
  large scene, [sequence.md](sequence.md) §7.4) the first changing call asks, once for the turn.
- **Ask every time**: `destructive` (removes or replaces), `upload` (sends anything out), and `undoable: false` (what no
  undo restores: a snapshot added, a file written for the person, which the page then saves). A denied call is told to
  the model.
- **Network, per address**: a `network` operation says what a call reaches (`addresses(values)`; `[]`: only the known
  services or the viewer's configured servers). Allowed: https to `KnownServices`, hosts the person wrote in the request
  (for that turn), domains the person allowed. Anything else asks, and so does an operation without `addresses` and an
  address not https written whole. A waiting call with blocked https hosts offers Allow <host> (`blockedHosts`,
  `hostsOf`); the rules (`createDomainRules`, session or browser) are the assistant's to keep.
- **Operations that run others** (an extension's `script.run`): their inner calls join their step and are not asked
  again; such an operation declares its own approval and addresses for what it reaches, and refuses for the agent what
  would ask on its own (the command-language extension's _spec §9).
- **Untrusted text**: the model sees compact state, tool results and metadata, never structure files. Text from data
  (titles, labels, entity names, file names, the person's notes) is `quoted()` (`agent/account.ts`): in quotes, escaped,
  controls and bidi/zero-width marks removed, at most 80 characters.
- **Keys and models are the person's**: the agent gets a `ModelSource`; keys never reach the model, the record or a log
  ([assistant/\_spec](../../assistant/_spec/README.md), "Keys").

## 5. What the agent sees each turn

- **State**: the scene tree's outline (`account.outline`): each turn at `components` level (structures, components with
  representations and atoms in a line, volumes, measurements; at most `StateLines` = 24 lines, the rest via
  `get_scene`); `get_scene` takes `summary`, `components` or `full` (default: representations and hidden rows too).
- **Focus**: the selection's and the focus's summaries (`account.focus`); the record keeps it with the turn.
- **The person's notes**: their own changing runs since the agent's last turn, in the viewer's words, quoted.
- **Repeated reads** of unchanged state are collapsed (`@molstar/ai`'s `RepeatedReads`).

## 6. Without a page

- **Headless**: `createViewerAgent` over a page-less viewer (`@molstar/prism-viewer/headless`); tools that need the
  canvas or files say why they cannot run.
- **Evals** (`agent/eval.ts`; README decision 19): a task's prompts run on a headless viewer, checked on the scene, not
  the model's words; Undo turn of each turn must restore the outline and the data tree's transforms. Scripted models in
  unit tests, real ones in `*.live.ts`.
