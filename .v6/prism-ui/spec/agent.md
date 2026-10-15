# The agent in the viewer: tools, skills, account, policies

How an agent drives the Prism viewer: the tools over the API ([operations.md](operations.md)), how they are disclosed,
what the agent sees each turn, what it may do without asking, and how it runs without a page. Built on the AI packages
(`packages/ai/_spec`) and Stories Untold's agent (in that repository `designs/architecture/agent.md`, `harness.md`).
Status: **accepted** in shape (2026-10-15); the tool lists are a first pass.

## 1. Where it lives

- **Tools, skills and the session factory**: the viewer's `agent/` folder (`packages/prism/viewer/src/agent/`),
  DOM-free, on `@molstar/ai-tools` and `@molstar/ai-session`. Extensions bring theirs in `ViewerExtension.tools` and
  `ViewerExtension.skills` (viewer-spec.md §3). PD params become tool fields through `paramsShape`, which moves from
  `@molstar/prism-assistant` into `agent/params.ts` (the viewer does not depend on the assistant; review F12).
- **The session**: `createViewerAgent({ viewer, adapters, storage })` in `agent/session.ts` makes the `AgentSession`
  (the viewer's toolkit, drive policies, account and focus, the record's hooks; [journal.md](journal.md)). The
  assistant calls it, and so do headless evals and an MCP server, with no DOM (review F11).
- **The assistant** (the Assistant tab, the Ask bar, the turn card, the Inspector, the Storage view) is
  `@molstar/prism-assistant`, loaded lazily when the app passes it: `AssistantFactory = (viewer, agent) =>
  AssistantUI` (viewer-spec.md §2.4).
- **One run path**: a tool calls `api.run(id, values, { origin: { kind: 'agent', turn } })`, so the history, the
  checkpoint and the run log are the same as a person's.

## 2. Tools

- **Thin adapters over operations**: a tool takes flat fields (few, named in the model's terms), normalizes what small
  models send (Stories' rules: "3" a number, "A, B" a list, an enum in any case, "" absent), maps to an operation's
  values and renders its result in one order: what it did, what changed, warnings, what it read, the next step.
  Several tools may share one operation (`color` and `transparency` over `representation.*`).
- **Words are the viewer's, not the screen's**: structures, chains, residues, the selection, the view; not tabs or
  panes. View tools (open a panel, focus, select, frame) say they leave the structures as they are.
- **Registry values are open strings** (a preset id, a theme name), checked when the tool runs, with the current
  choices in the description, refreshed as the registry changes (review F2).
- **Targets are strings** ([targets.md](targets.md)); the selections skill's guidance teaches the forms.

## 3. The core and the skills

| Core (from the first turn)                                                                         | Skills (disclosed on request)                                                       |
| -------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------- |
| `list_skills`, `disclose_skills`, `disclose_tools` (always)                                         | `selections`: MolQL and SeqQL, PyMOL text, addresses; `select_around`, `expand`      |
| `get_scene` (the outline), `get_selection`                                                          | `styles`: presets, quick styles, colors and themes, transparency overlays            |
| `load_entry` (PDB, AFDB, EMDB, …), `load_url`                                                       | `measure`, `superpose`                                                              |
| `select`, `focus`, `clear_selection`                                                                | `snapshots`, `animation`, `trajectory`                                              |
| `apply_preset`, `apply_style`, `color`                                                             | `canvas`: lighting, illumination, clipping                                          |
| `set_background`, `reset_camera`                                                                   | `sequence`: style, colors, tracks                                                   |
| `undo_turn` (reads the history; the person's Undo is theirs)                                       | `export`; extensions' skills (`geo-export`, `assembly-symmetry`, `model-archive`, …) |
| `ask` (a question to the person: options, a fill-in, or both; the turn waits for the answer)       | `viewer-code`: TypeScript against the API in a sandbox (api.md §7)                  |

- **The first scope** (README §3.3): background and canvas, presets and styles, MolQL and SeqQL selections, fetching
  structures, focus and select, opening panels. The core stays small (Stories' is 24 tools, about 6.8k tokens); a new
  tool goes into a skill unless nearly every turn needs it.
- **Disclosure** is the toolkit's: `disclose_skills` loads guidance and tools, several in one call; extras and single
  tools load by name or when called; where the model allows, a load is appended so the cached prompt stays
  (packages/ai/\_spec/tools.md, "The toolkit"). `list_skills` is the viewer's tool over the toolkit's catalog.
- **A dynamic toolkit** (decision 25; today `createToolkit` is fixed at creation, review F2): tools and skills are added
  and removed between and during turns. An extension loaded adds its tools and skills; a lazy extension's skill is in
  the catalog by its stub's summary, and disclosing it loads the chunk, registers its tools and appends their
  definitions; a feature removed, or a registry changed, removes or replaces tools, and the removed definitions are an
  eviction hint at the next compaction. A mechanism in `@molstar/ai-tools` and `ai-session`; what is registered is the
  viewer's.
- **`ask`**: when the intent is unclear (load a second copy? which chain?), the agent asks instead of guessing; the turn
  parks as it does for an approval, the person picks an option or types, and the answer is the tool's result. As built
  so far (plan 4.1): `ask` tells the agent to end its turn, and the person's answer comes as the next request; the
  parked question comes with the assistant's turn card.

## 4. Policies (from Stories Untold)

As little in the way as is reasonable (decision 26):

- **Preauthorized**: operations that change only the viewer's state while the turn can be undone. When it cannot (undo
  off for a large scene, sequence.md §7.4, or steps made unsafe by a change outside the API), the turn card says so
  before the turn and the first changing call asks once for the turn (review F14).
- **Network, per address** (review F13): known data services (PDB, PDBe, RCSB, AlphaFold DB, EMDB, ModelArchive,
  UniProt) and URLs the person wrote in the request are allowed. Any other address asks when a call reaches it, whatever
  the operation (the approval is decided at call time from the URL, through the harness's `permission` decision):
  **Allow once**, **Allow this domain** (for the session, or in this browser), **Deny**. Domain rules are listed in the
  Storage view and in Settings, each removable. https only.
- **Ask every time**: sending anything out (`upload`: remote states) and closing or clearing a session (`destructive`).
  A denied call is told to the model ("The person did not allow this call.").
- **Untrusted text**: the model never sees structure files, only compact state, tool results and extracted metadata.
  Text that comes from data (a structure's title and labels, an entity's name, an MVS snapshot's description, a file
  name) is quoted as data wherever it reaches the model, in tool results and in the per-turn state block alike: in a
  fence, escaped, capped in length (review F15).
- **Keys** are the person's, kept by the assistant for as long as they chose (the mocks' `KeyScope`): until a reload
  (memory only), for this tab (sessionStorage) or in this browser (localStorage). **Stored keys are always scrambled**
  (David, 2026-10-15), building on Stories' `secrets.ts` (ported into `@molstar/prism-assistant`), so the storage
  read as it is shows no key, in every scope, whatever the browser allows. The page keeps what it needs **in memory**
  (the scrambling key, and a provider key while it is in use); which scrambling key is used depends on what the
  browser offers (IndexedDB is not always granted, WebCrypto needs a secure context):

  | The browser offers     | The scrambling key                                                                          | Against                                        |
  | ---------------------- | ------------------------------------------------------------------------------------------- | ---------------------------------------------- |
  | WebCrypto and IndexedDB | a per-browser AES-GCM key, non-extractable, generated once, kept in IndexedDB and in memory | reading a stored value, or copying it to another browser: no working key (copying the whole browser profile copies the key too) |
  | WebCrypto only         | an app key in the source (AES-GCM, imported from a constant)                                | a naive read: no key in plain sight            |
  | neither                | the same app key, as a simple reversible scramble (XOR, base64)                             | a naive read                                   |

  Every scope (until a reload, this tab, this browser) is offered in every row; in the last two the Assistant tab says
  that a stored key is scrambled, not protected. A value scrambled one way and read where another is in use (IndexedDB
  granted later, or withdrawn) is opened by trying the per-browser key, then the app key. Nothing here guards against
  code running in the page. A key is never in the model's context, the session record, a session file or the log.
- **The origin is shared** (review F17): every page served from the viewer's origin can use stored keys and read kept
  records. The app is deployed on an origin of its own (app.md §1), and the Storage view says which origin holds what.

## 5. What the agent sees each turn

- **State**: the scene's outline (structures with their components, representations and colors; volumes; snapshots;
  the canvas), capped (Stories' state block: at most 40 entries, the rest through `get_scene`); labels from data are
  quoted (§4).
- **Focus**: what the person points at now and never stored: Mol\*'s selection (as addresses per chain and as MolQL),
  the focus, the open panel, the hovered residue.
- A repeated read of unchanged state is collapsed (`RepeatedReads`).

## 6. Without a page

- **Headless** (`@molstar/prism-viewer/headless`): the same toolkit and `createViewerAgent` over the same API on a
  page-less viewer; tools that need the canvas or files declare those modules and say why they cannot run.
- **Evals** come later (README §2.19): a task, an agent on the headless viewer, a check on the state account; the run
  log replays in the app.
- **MCP**: an outside agent calls the same tools through an MCP server over the headless viewer, with the same policies;
  later.
