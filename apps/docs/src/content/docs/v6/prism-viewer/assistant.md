---
title: The Prism assistant
description: "The optional AI assistant of the Prism viewer: your keys and how they are kept, what the agent may do without asking, undoing a turn, the session record, and running the agent headless."
---

The assistant lets you ask the viewer for things in words: load an entry, select a ligand, color a chain, change the
background. It is optional. The viewer works without it, and the default specs do not include it.

When the app includes the assistant, you get an **Assistant** tab, an **Ask** row in the command menu and a card under
the status line that shows the running turn. The assistant's code is a chunk of its own that loads after the viewer is
up.

<!-- prettier-ignore-start -->
:::caution[Work in progress]
The assistant is the first slice. The last section lists what is not built yet.
:::
<!-- prettier-ignore-end -->

## How it works

The assistant runs an agent over the viewer's [API](/docs/v6/prism-viewer/api/). Its tools are the viewer's operations
in the model's words, so what the agent does goes through the same path as a person's call: it is a history step, it is
in the record, and Undo restores it.

A turn is one request and everything done for it. One turn runs at a time: asking while one runs, or while one waits for
your approval, is refused until you answer or stop it.

Each turn starts with a short account of the scene: the outline of structures, components and representations (up to 40
lines; the agent reads more with `get_scene`), what you have selected and focused, and what you did yourself since the
agent's last turn. The agent never sees structure files. Text that comes from data, such as a title or a chain name, is
quoted as data in that account and is never treated as an instruction.

When your request is unclear, the agent asks with the `ask` tool and ends the turn. Your answer comes as the next
request.

### Tools and skills

The agent starts with a small set of tools: `get_scene`, `get_selection`, `load_entry`, `load_url`, `select`,
`clear_selection`, `focus`, `apply_preset`, `apply_style`, `color`, `set_background`, `reset_camera`, `undo_turn` and
`ask`, with `disclose_skills` and `disclose_tools` to load more. The rest are in skills that the agent loads when a
request needs them:

| Skill        | For                                                            |
| ------------ | -------------------------------------------------------------- |
| `selections` | Select around a target, expand, invert; the forms of a target. |
| `styles`     | Components, representations, transparency, hiding entities.    |
| `canvas`     | Lighting, illumination, camera mode, clipping.                 |
| `measure`    | Distances, angles, dihedrals, labels.                          |
| `superpose`  | Superposition by chains, TM-align, ligands.                    |
| `snapshots`  | Save views of the scene and step through them.                 |
| `animation`  | Spin, rock, wiggle, trajectory frames.                         |
| `export`     | Images and session files.                                      |

[Extensions](/docs/v6/prism-viewer/extensions/) add skills for their own operations. A skill of a lazy extension is
listed from the start, and loading it loads the extension.

## Your keys

The assistant calls a model provider directly from your browser, with your own key. A key goes only to its provider.

Providers are Anthropic, Google and OpenAI. You can add your own OpenAI-compatible endpoint, such as OpenRouter, Ollama
or vLLM. Its address must be https, or http on `localhost`, and must not contain a login. Add the key in its own field.

You choose how long a key is kept:

| Scope            | Kept in                                      |
| ---------------- | -------------------------------------------- |
| **Until reload** | Memory. It is gone when the page reloads.    |
| **This session** | `sessionStorage`: this tab, until it closes. |
| **This browser** | `localStorage`: until you remove it.         |

A kept key is always scrambled, in every scope, so reading the browser's storage as it is shows no key. How strong the
scrambling is depends on what the browser offers:

| The browser offers      | The key that scrambles it                                                                      |
| ----------------------- | ---------------------------------------------------------------------------------------------- |
| WebCrypto and IndexedDB | A per-browser AES-GCM key that cannot be exported. It is generated once and kept in IndexedDB. |
| WebCrypto only          | An AES-GCM key from the source code.                                                           |
| Neither                 | A reversible scramble (XOR and base64) with the key from the source code.                      |

Only the first protects a stored value against being copied to another browser. The other two only hide it from a casual
read. The Assistant tab says that kept keys are scrambled and that this does not guard against code running in the page.
No scheme does: every script on the page can use a key while it is in memory. A key is not sent to the model and is not
in the session record.

Every page served from the same origin can use stored keys and kept records. Serve the app from an origin that serves
nothing else.

The model you chose and your own endpoints are kept in `localStorage`, apart from the keys.

## What it may do without asking

The aim is that the assistant stays out of the way for what you can undo, and asks for the rest.

**Without asking:**

- Reads, and view changes: selecting, focusing, moving the camera, highlighting.
- Changes to the scene while the turn can be undone: loading, styling, coloring, adding components, changing the canvas.
- Network access to the known data services: RCSB PDB, PDBe and the EBI, PDBj, AlphaFold DB, ModelArchive, UniProt and
  Zenodo (https only, subdomains included), and to https addresses that you wrote in your request.

**Asking each time**, with Allow and Deny on the turn card:

- Removing or replacing: removing a structure, a component, a representation, a measurement, a snapshot or a state cell,
  clearing a snapshot list, opening a session in place of the scene.
- Sending something out, such as uploading a state.
- A change that no undo restores, for example adding a snapshot or making a file (an export, a saved session); the
  viewer saves a file the assistant made for you.
- Any other network address, and any network operation that does not say what it reaches. Only https goes without
  asking.
- The first change of a turn, once, when undo is off because the scene is large.

A denied call is told to the model, which carries on without it. Stop ends a turn that is running, or one waiting for an
answer.

## Undo turn

A turn is one history step, so one Undo takes back everything it changed. While the turn is the latest step, the turn
card has **Undo turn**. **History** has it as well, on the agent's latest turn. If you made steps after the turn,
History tells you how many and asks before it undoes those too. The agent can also call `undo_turn` to take back its own
previous turn, but never your later steps.

Like any undo, this restores the scene (structures, components, representations, canvas, the sequence's style) and not
the camera, the selection or the focus. See [History and undo](/docs/v6/prism-viewer/api/#history-and-undo).

## The session record

The viewer keeps a record of the session: what you did and what the agent did, in order. **History** shows it, newest
first. Your runs are listed one by one. An agent's turn is one item with its request, its reply and the runs under it,
and undone steps are marked. "Show earlier" reads further back. Without the assistant, History still shows your own
runs.

By default the record is in memory and ends when the page reloads. A reload also ends every checkpoint, so a turn cannot
be undone after it.

**Keep this record** in History copies the record into IndexedDB in your browser and goes on appending to it. This keeps
only the record, which holds your requests and the agent's replies. It does not keep the scene, files or undo. The scene
comes back only from a session file that you saved.

For a viewer you create yourself, `persist.rows` in the spec chooses where "Keep" writes, and `viewer.record` is the
record's object.

## For developers

### Adding the assistant to your viewer

```ts
import { assistant } from '@molstar/prism-assistant/spec';

const spec = {
  ...FullPrismViewerSpec(),
  agent: { assistant: assistant({ provider: 'google', model: 'gemini-3.8-flash' }) },
};
```

`assistant(options)` returns a loader: the assistant, the agent and the AI SDK are a chunk of their own, loaded when the
assistant is first used, never in the page's startup bundle. Its options are all optional:

| Option              | What it is                                                                                                                              |
| ------------------- | --------------------------------------------------------------------------------------------------------------------------------------- |
| `provider`, `model` | The model until the person chooses one, when its provider has a key. By default Anthropic's `claude-sonnet-5-5`.                        |
| `provided`          | Models the page brings with their providers (a server's), used without a key of the person's. A function loads them with the assistant. |
| `keys`, `prefs`     | Where the person's keys and their choice of model are kept (by default the browser's storage).                                          |

The keys are always the person's, typed into the Assistant tab: a page never passes a key.

`agent` can also change the viewer's skills and core tools. Each is a change to the defaults, by id, so you never
restate them: `add` appends (an entry with a default's id replaces it) and `remove` takes defaults out.

```ts
agent: {
  skills: { add: [mySkill], remove: ['export'] },
  core: { add: [{ id: 'get_camera', operation: 'camera.get' }], remove: ['load_url'] },
},
```

The declarations, `CoreTools` and `Skills` with their types, are plain data in `@molstar/prism-viewer/agent/skills`,
which imports nothing. With `agent: false` the viewer has no tools, and the agent keeps only `get_scene`, `ask` and
`undo_turn`.

### Running the agent headless

The agent, its tools and its policies do not need the assistant or a DOM. `createViewerAgent` makes one over any viewer,
headless included:

```ts
import { createAnthropicProvider } from '@molstar/ai/providers/anthropic';
import { DefaultPrismViewerSpec } from '@molstar/prism-viewer/default-spec';
import { createPrismViewerHeadless } from '@molstar/prism-viewer/headless';
import { createViewerAgent } from '@molstar/prism-viewer/agent/session';

const viewer = await createPrismViewerHeadless({ spec: DefaultPrismViewerSpec() });
await viewer.api.load.files({ files: [new File([bytes], '1cbs.cif')] });

const anthropic = createAnthropicProvider({ resolveApiKey: async () => process.env.ANTHROPIC_API_KEY! });
const agent = createViewerAgent({
  viewer,
  models: {
    providerFor: (id) => (id === 'anthropic' ? anthropic : undefined),
    default: { provider: 'anthropic', model: 'claude-sonnet-5-5' },
  },
});

let result = await agent.turn('Color the ligand red.');
while (result.status === 'awaiting-approval') {
  const decisions = Object.fromEntries((result.approval ?? []).map((p) => [p.callId, 'allow' as const]));
  const next = await agent.resume({ decisions });
  if (!next) break;                       // it was stopped meanwhile
  result = next;
}
console.log(result.status, result.text);

agent.dispose();
```

The options are `viewer`, `models` (a provider per provider id and a default model), and optionally `storage` for the
conversation (memory by default), `guidance` (more instructions after the viewer's) and `rules` (the domains a person
allowed). A key is the host's to provide: `resolveApiKey` here reads it from the environment.

| Member                          | What it is                                                                        |
| ------------------------------- | --------------------------------------------------------------------------------- |
| `turn(prompt, hooks?)`          | Runs a turn. Throws while one runs or waits.                                      |
| `resume({ decisions, hooks? })` | Continues a parked turn with `allow` or `deny` for each waiting call, by call id. |
| `stop()`                        | Stops the running turn, or ends a parked one.                                     |
| `busy`, `parked`, `lastTurn`    | Whether a turn runs; whether one waits for approval; the last finished turn's id. |
| `policy`                        | What goes without asking, and its domain rules (`policy.rules`).                  |
| `toolkit`                       | The tools.                                                                        |
| `inspect`                       | The conversation's harness, budget and context window, read-only.                 |

### Evaluations

`runEval` runs a task on a viewer: it loads the task's files, runs the agent through the prompts, and checks the scene,
not the model's words. A task can also require that Undo turn restores the scene.

```ts
import { readFile } from 'node:fs/promises';
import { evalReport, runEval, type EvalTask } from '@molstar/prism-viewer/agent/eval';

const task: EvalTask = {
  id: 'select-range',
  title: 'A selection by residue range',
  files: ['1cbs.cif'],
  prompts: ['Select residues 10 to 20 of chain A.'],
  undoes: true,                                     // Undo turn must bring the scene back
  check: async (viewer) => {
    const r = await viewer.api.select.get();
    return [{ name: 'something is selected', pass: r.ok && Array.isArray(r.value) && r.value.length > 0 }];
  },
};

const viewer = await createPrismViewerHeadless({ spec: DefaultPrismViewerSpec() });
const result = await runEval(viewer, task, {
  models,                                           // as above
  file: async (name) => new File([await readFile(name)], name),
});
console.log(evalReport(result));
```

An approval is answered with allow, or with `approve: 'deny'` on the task. The result has the `checks`, each turn's
tools and token use, and `runs`, the log of the API calls with their arguments and origins. The repository's unit tests
run the tasks with a scripted model, and `pnpm test:live` runs them on real models. That costs money and needs keys in
`.env.local`.

## Not built yet

- Allowing a domain from the turn card, and a list of the allowed domains. Policies support them, and the card offers
  Allow and Deny.
- The Storage view that shows what the viewer keeps in the browser, with Forget.
- Kept records listed on Home, and the record inside a saved `.molx` or `.molj`.
- Skills for the sequence view and for writing TypeScript against the API.
- Choices in a tool's description that follow the registries as extensions add to them.
- An MCP server over the headless viewer.
