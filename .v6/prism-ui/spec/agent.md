# The agent in the viewer: tools, skills, account, policies

How an agent drives the Prism viewer: the tools over the API ([operations.md](operations.md)), how they are disclosed,
what the agent sees each turn, what it may do without asking, and how it runs without a page. Built on the AI packages
(`packages/ai/_spec`) and Stories Untold's agent (in that repository `designs/architecture/agent.md`, `harness.md`).
Status: **accepted** in shape (2026-10-15); the tool lists are a first pass.

## 1. Where it lives

- **Tools and skills**: the viewer's `agent/` folder (`packages/prism/viewer/src/agent/`), DOM-free, on
  `@molstar/ai-tools` (`defineTool`, `createToolkit`). Extensions bring their skills in `ViewerExtension.skills`.
- **The assistant** (the Assistant tab, the Ask bar, the turn card, the Inspector, the session over `@molstar/ai-session`)
  is `@molstar/prism-assistant`, loaded lazily when the app passes it (viewer-spec.md §2.4).
- **The session**: an `AgentSession` (`@molstar/ai-session`) made by the assistant's composition, its turns written to
  the session record ([journal.md](journal.md)).
- **One run path**: a tool calls `api.run(id, values, { origin: { kind: 'agent', turn } })`, so the history, the
  checkpoint and the run log are the same as a person's.

## 2. Tools

- **Thin adapters over operations**: a tool takes flat fields (few, named in the model's terms), normalizes what small
  models send (Stories' rules: "3" a number, "A, B" a list, an enum in any case, "" absent), maps to an operation's
  values and renders its result in one order: what it did, what changed, warnings, what it read, the next step.
  Several tools may share one operation (`color` and `transparency` over `representation.*`).
- **Words are the viewer's, not the screen's**: structures, chains, residues, the selection, the view; not tabs or
  panes. View tools (open a panel, focus, select) say they change only what the person sees.
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

- **The first scope** (README §3.3): background and canvas, presets and styles, MolQL and SeqQL selections, fetching
  structures, focus and select, opening panels. The core stays small (Stories' is 24 tools, about 6.8k tokens); a new
  tool goes into a skill unless nearly every turn needs it.
- **Disclosure** is the toolkit's: `disclose_skills` loads guidance and tools, several in one call; extras and single
  tools load by name or when called; every tool is declared from the start, and where the model allows, a load is
  appended so the cached prompt stays (packages/ai/\_spec/tools.md, "The toolkit").

## 4. Policies (from Stories Untold)

- **Preauthorized**: every operation that changes only the viewer's state (the history undoes the turn as one step).
- **Asks first** (`approval`): fetching from an address the person did not give or that is not a known data service
  (`network`), sending anything out (`upload`: remote states), and closing or clearing a session (`destructive`). A
  denied call is told to the model ("The person did not allow this call.").
- **Network**: loads from the known services (PDB, PDBe, RCSB, AlphaFold DB, EMDB, ModelArchive, UniProt) and from URLs
  the person wrote; https only.
- **Untrusted text**: the model never sees structure files, only compact state, tool results and extracted metadata;
  free text that arrives with data (an MVS snapshot's description, a file's title) is marked as data in the result, never
  as instructions.
- **Keys** are the person's, kept by the assistant for as long as they chose (the mocks' `KeyScope`): until a reload
  (memory only), for this tab (sessionStorage) or in this browser (localStorage). **Stored keys are always scrambled**
  (David, 2026-10-15), building on Stories' `secrets.ts` (ported into `@molstar/prism-assistant`), so the storage
  read as it is shows no key, in every scope, whatever the browser allows. The page keeps what it needs **in memory**
  (the scrambling key, and a provider key while it is in use); which scrambling key is used depends on what the
  browser offers (IndexedDB is not always granted, WebCrypto needs a secure context):

  | The browser offers     | The scrambling key                                                                          | Against                                        |
  | ---------------------- | ------------------------------------------------------------------------------------------- | ---------------------------------------------- |
  | WebCrypto and IndexedDB | a per-browser AES-GCM key, non-extractable, generated once, kept in IndexedDB and in memory | reading or copying storage: no working key     |
  | WebCrypto only         | an app key in the source (AES-GCM, imported from a constant)                                | a naive read: no key in plain sight            |
  | neither                | the same app key, as a simple reversible scramble (XOR, base64)                             | a naive read                                   |

  Every scope (until a reload, this tab, this browser) is offered in every row; in the last two the Assistant tab says
  that a stored key is scrambled, not protected. A value scrambled one way and read where another is in use (IndexedDB
  granted later, or withdrawn) is opened by trying the per-browser key, then the app key. Nothing here guards against
  code running in the page. A key is never in the model's context, the session record, a session file or the log.

## 5. What the agent sees each turn

- **State**: the scene's outline (structures with their components, representations and colors; volumes; snapshots;
  the canvas), capped (Stories' state block: at most 40 entries, the rest through `get_scene`).
- **Focus**: what the person points at now and never stored: Mol\*'s selection (as addresses per chain and as MolQL),
  the focus, the open panel, the hovered residue.
- A repeated read of unchanged state is collapsed (`RepeatedReads`).

## 6. Without a page

- **Headless** (`@molstar/prism-viewer/headless`): the same toolkit over the same API on a page-less viewer; tools that
  need the canvas or files declare those modules and say why they cannot run.
- **Evals** come later (README §2.19): a task, an agent on the headless viewer, a check on the state account; the run
  log replays in the app.
- **MCP**: an outside agent calls the same tools through an MCP server over the headless viewer, with the same policies;
  later.
