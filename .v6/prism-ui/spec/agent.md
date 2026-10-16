# The agent in the viewer: tools, skills, account, policies

How an agent drives the Prism viewer: the tools over the API ([operations.md](operations.md)), how they are disclosed,
what the agent sees each turn, what it may do without asking, and how it runs without a page. Built on the AI packages
(`packages/ai/_spec`) and Stories Untold's agent (in that repository `designs/architecture/agent.md`, `harness.md`).
Status: **accepted** in shape (2026-10-15); the tool lists are a first pass.

## 1. Where it lives

- **Tools, skills and the session factory**: the viewer's `agent/` folder (`packages/prism/viewer/src/agent/`),
  DOM-free, on `@molstar/ai` (`tools/`, `session/`). Extensions bring theirs as skills,
  `ViewerExtension.skills` (a lazy one's in its stub; never core tools); the spec's `agent.core` and `agent.skills`
  change the viewer's defaults by id (`{ add, remove }`; an added id of a default's replaces it), never restating them
  (viewer-spec.md §2.4, §3). The declarations (`CoreTools`, `Skills`, their types) are plain data in `agent/skills.ts`,
  which imports nothing, so a spec that reads or extends them does not bring the toolkit (`@molstar/ai` tools, zod) into the page
  (review B10). PD params become tool fields through `paramsShape`, which moves from
  `@molstar/prism-assistant` into `agent/params.ts` (the viewer does not depend on the assistant; review F12).
- **The session**: `createViewerAgent({ viewer, models, storage?, rules? })` in `agent/session.ts` makes the `AgentSession`
  (the viewer's toolkit, drive policies, account and focus, the record's hooks; [journal.md](journal.md)). The
  assistant calls it, and so do headless evals and an MCP server, with no DOM (review F11). Its `ViewerAgent` keeps the
  session to itself: `turn(prompt, hooks?)`, `resume({ decisions, hooks? })`, `stop()`, `busy`, `parked`, and
  `inspect` (the harness client, budget and window, read-only, for the Inspector; review B10).
- **The assistant** (the Assistant tab, the Ask bar, the turn card, the Inspector, the Storage view) is
  `@molstar/prism-assistant`, loaded lazily when the app passes it (`agent: { assistant: assistant({ provider, model })
  }`, from `@molstar/prism-assistant/spec`, types only and one dynamic import; the model is the default until the
  person chooses, in effect while its provider has a key): `AssistantFactory = (host: AssistantHost) =>
  AssistantUI`, the host giving the viewer and `createAgent(models)`, since the models come from the person's keys,
  which the assistant keeps (viewer-spec.md §2.4). As built: the Inspector is a page over the app
  (`assistant.page.inspector`, an extension view `pages` entry), opened by the Assistant tab's Inspect rows (Trace with
  its turns, Tools with its counts, Limits), the turn card's trace and the commands Assistant trace (⌥T), tools and
  limits (registered with the assistant, so only once it has loaded). Its data (`inspector-sources.ts`): the trace read
  from the session's harness (`agent.inspect.client()`, re-read on its `watch()`), the Tools view from the toolkit as it
  changes, the limits from `agent.inspect.budget()` and the model's window, read-only (the session takes them at creation and
  has no setter).
- **A turn at a time**: `turn()` refuses while one runs or waits for an approval (`parked`); `stop()` ends either, and a
  parked turn is then settled as stopped (its record entry, its history group).
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
  choices in the description, refreshed as the registry changes (review F2). As built: a select becomes the field's
  choices; between turns the toolkit describes its tools again and replaces the parts whose definitions changed (a
  theme an extension registered, a feature loaded or left out), so the next turn sees them.
- **What a tool says it cannot do** is what this viewer cannot do at all (a module it lacks, headless); what cannot be
  done now (nothing loaded) is the call's answer, so descriptions do not change with the scene.
- **Targets are strings** ([targets.md](targets.md)); the selections skill's guidance teaches the forms.

## 3. The core and the skills

| Core (from the first turn)                                                                         | Skills (disclosed on request)                                                       |
| -------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------- |
| `disclose_skills`, `disclose_tools` (always; the catalog is `disclose_skills`' description)         | `selections`: MolQL and SeqQL, PyMOL text, addresses; `select_around`, `expand`      |
| `get_scene` (the outline), `get_selection`                                                          | `styles`: presets, quick styles, colors and themes, transparency overlays            |
| `load_entry` (PDB, AFDB, EMDB, …), `load_url`                                                       | `measure`, `superpose`                                                              |
| `select`, `focus`, `clear_selection`                                                                | `snapshots`, `animation`, `trajectory`                                              |
| `apply_preset`, `apply_style`, `color`                                                             | `canvas`: lighting, illumination, clipping                                          |
| `set_background`, `reset_camera`                                                                   | `sequence`: style, colors, tracks                                                   |
| `undo_turn` (its previous turn, never the person's later steps)                                    | `export`; extensions' skills (`geo-export`, `assembly-symmetry`, `model-archive`, …) |
| `ask` (a question to the person: options, a fill-in, or both; the turn waits for the answer)       | `viewer-code`: TypeScript against the API in a sandbox (api.md §7)                  |

- **The first scope** (README §3.3): background and canvas, presets and styles, MolQL and SeqQL selections, fetching
  structures, focus and select, opening panels. The core stays small (Stories' is 24 tools, about 6.8k tokens); a new
  tool goes into a skill unless nearly every turn needs it.
- **Disclosure** is the toolkit's: `disclose_skills` loads guidance and tools, several in one call; extras and single
  tools load by name or when called; where the model allows, a load is appended so the cached prompt stays
  (packages/ai/\_spec/tools.md, "The toolkit").
- **As built** (plan 4.1, G6): the first scope's tools are there but for SeqQL selections (the sequence's own skill is
  not built) and opening panels (the agent speaks of no panes; a view tool may come with the assistant's UI).
- **A dynamic toolkit** (decision 25; today `createToolkit` is fixed at creation, review F2): tools and skills are added
  and removed between and during turns. An extension loaded adds its tools and skills; a lazy extension's skill is in
  the catalog by its stub's summary, and disclosing it loads the chunk, registers its tools and appends their
  definitions; a feature removed, or a registry changed, removes or replaces tools, and the removed definitions are an
  eviction hint at the next compaction. A mechanism in `@molstar/ai` (`tools/`, `session/`); what is registered is the
  viewer's.
- **`ask`**: when the intent is unclear (load a second copy? which chain?), the agent asks instead of guessing; the turn
  parks as it does for an approval, the person picks an option or types, and the answer is the tool's result. As built
  so far (plan 4.1): `ask` tells the agent to end its turn, and the person's answer comes as the next request; the
  parked question comes with the assistant's turn card.

## 4. Policies (from Stories Untold)

As little in the way as is reasonable (decision 26):

- **Preauthorized**: operations that change only the viewer's state while the turn can be undone. When it cannot (undo
  off for a large scene, sequence.md §7.4, or steps made unsafe by a change outside the API), the turn card says so
  before the turn and the first changing call asks once for the turn (review F14). A step no undo restores
  (`undoable: false`: a snapshot added, a file written) asks every time. As built: undo off is read from checkpointing;
  unsafe steps are not considered yet.
- **Network, per address** (review F13): known data services (PDB, PDBe, RCSB, AlphaFold DB, EMDB, ModelArchive,
  UniProt; Zenodo, plan 5.2) and URLs the person wrote in the request are allowed. Any other address asks when a call reaches it, whatever
  the operation (the approval is decided at call time from the URL, through the harness's `permission` decision):
  **Allow once**, **Allow this domain** (for the session, or in this browser), **Deny**. Domain rules are listed in the
  Storage view and in Settings, each removable. https only. As built (G6): an operation says what a call reaches,
  `addresses(values)` on the values it will run with (`[]`: the known services and the servers the viewer is configured
  with); one that does not say always asks, and so does an address that is not https written whole. The request's
  https hosts go for that turn. Allow this domain is built: a waiting call whose addresses leave one https host not
  allowed offers Deny, Allow once and Allow <host> (with Remember: in this browser, the assistant's prefs; else for
  the session); the policy names it (`blockedHosts`, `hostsOf`), the rules are the assistant's (`createDomainRules`,
  passed to `createAgent`) and the Storage view lists them with Forget. Not yet: a listing in Settings; a
  color theme or preset that fetches as it applies (validation reports, assembly symmetry) is reached through the core
  `color` and `apply_preset` tools, which carry no addresses: the ones registered today reach only the known services.
- **Ask every time**: sending anything out (`upload`: remote states) and closing or clearing a session (`destructive`).
  A denied call is told to the model ("The person did not allow this call.").
- **Untrusted text**: the model never sees structure files, only compact state, tool results and extracted metadata.
  Text that comes from data (a structure's title and labels, an entity's name, an MVS snapshot's description, a file
  name) is quoted as data wherever it reaches the model, in tool results and in the per-turn state block alike: in a
  fence, escaped, capped in length (review F15). As built: quoted (`"…"`, controls, separators, bidi and zero-width
  marks removed, 80 characters) in the state block (labels and descriptions), the focus, and the person's notes; tool
  results still name data unquoted ("Loaded 1CBS."), which is next.
- **A host's own models** (`assistant({ provided })`; a function makes them when the assistant loads): a provider the host brings with its adapter, listed
  first with its key "from the server" and in effect without a key of the person's; the stories bring a scripted model
  this way (keys and the choice kept in memory, never the person's).
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
  quoted (§4). As built: the scene tree's outline at a level (`account.outline`): each turn `components` (structures,
  their components with representations and atoms in a line, volumes, measurements; at most 24 lines), `get_scene` a
  `detail` of `summary`, `components` or `full` (the default: representations and hidden rows too); no colors,
  snapshots or canvas yet.
- **Tokens** (David, 2026-10-15): what every request carries is kept small. The shared params (`target`, `structure`,
  `language`) are said once in the guidance, not in each schema; a select's names are its enum only (no list in
  words), and a long one (over 16: color themes) has no enum, a wrong name answering with all of them; no `$schema`.
  Measured on 1CBS, a first request went from about 5,070 input tokens to 3,800 on Haiku 5.5, 3,140 to 2,510 on
  Gemini 3.8 Flash; the evals pass. A turn with a tool call is two requests, the second mostly cached (Anthropic).
- **Focus**: what the person points at now: Mol\*'s selection (as addresses per chain and as MolQL), the focus, the open
  panel, the hovered residue. As built: the selection's and the focus's summaries; the record keeps it with the turn
  (journal.md §2).
- **The person's notes**: what they did since the agent's last turn, in the viewer's words, quoted.
- A repeated read of unchanged state is collapsed (`RepeatedReads`). Not yet: no state version is passed.

## 6. Without a page

- **Headless** (`@molstar/prism-viewer/headless`): the same toolkit and `createViewerAgent` over the same API on a
  page-less viewer; tools that need the canvas or files declare those modules and say why they cannot run.
- **Evals** (README §2.19; plan 4.5, `agent/eval.ts`): a task, an agent on the headless viewer, checks on the scene,
  and Undo turn restoring it (the outline and the data tree's transforms); the run log replays in the app (not yet).
- **MCP**: an outside agent calls the same tools through an MCP server over the headless viewer, with the same policies;
  later.
