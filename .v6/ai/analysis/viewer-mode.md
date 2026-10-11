# The agent in the viewer mode: tools, skills, UI, fibers

Read-only sweeps (2026-10-11) of Stories Untold's agent (`stories-untold-redux/lib/stories/stories-agent`, `any-molql`,
`any-molstar`, `any-toolscript`, `mvs-code`), its agent UIs (`stories-ui/src/inspector`) and fibers
(`vendor/anystar/any-components/src/fibers` with its conformance, `lib/stories/stories-fibers`). What the viewer mode
needs, what becomes a generic `@molstar/ai-*` package, and what Prism writes itself.

## 1. In short

- **Stories' tool layer is about 60% generic** once five seams are cut: its document API (`StoryToolApi`), its result
  wording, the `Request:`/`Focus:` markers, the skill tables as globals, and the document's identity as the "state
  version". That part becomes `@molstar/ai-tools` and `@molstar/ai-session` (§4).
- **Every write tool is written against the story document** (operations applied all or nothing, desired state, an
  inverse for undo). The viewer has no such document: its tools are Prism commands over the live plugin, as
  [agent.md](../../prism-ui/design/agent.md) already says. 36 of Stories' 72 tools adapt, 24 are the story's only, 8 are
  generic (§2).
- **Selections are the best reuse**: the three-field input, the guidance, the generated PyMOL and MolScript references,
  `query_reference`, and the projection and bounds helpers. Most of `any-molql` is already `@molstar/query-language`;
  the rest is small (§3).
- **The know-how is in the guidance and the result text**, not the schemas: flat inputs, argument normalization for
  small models, results in a fixed order (did, changed, warnings, body, next), errors that name the valid values,
  repeated reads collapsed, reads scheduled together. All of it carries over.
- **The trace, the tool list and the limits are one Inspector** in Stories (tabs Trace, Tools, Limits; about 1.5k lines
  of React) over framework-free view models. Prism rebuilds the views on its runtime over the same data (§5).
- **Fibers know nothing of the harness**: a generic, append-only, user-facing transcript (one entry per turn) that
  survives compaction, which the harness's journal does not. The kernel ports as it is; its storage needs a small port
  of its own (§6).

## 2. Tools for the viewer mode

Tools are commands: one registry, a command's params as the tool's JSON Schema, `read`, `approval` and `ui` in its
metadata (contexts.md §4). Commands the person never types (`elements.get`, `bounds.get`, `query_reference`, the
disclosure tools) are commands with `ui: false`. A command that selects takes **one selection param type** with Stories'
three forms (`target`: keywords, `B:200-210`, `C:122@OG1`, a component or node, `selection` for the current one;
`items`: MVS component expressions; `query`: PyMOL text or a MolScript expression), so the command menu accepts them
too.

**Core (loaded from the first turn)**

| Group          | Commands                                                                 | Stories                                              |
| -------------- | ------------------------------------------------------------------------ | ---------------------------------------------------- |
| State          | `scene.get` (the state account, capped; one node by ref)                 | `get_story`, adapted                                 |
| Load           | `structure.load` (PDB, AlphaFold, URL, file; assembly; preset)           | `load_structure`, its result prose kept              |
| Selection      | `selection.set`, `selection.clear`, `selection.get`, `elements.get`      | `select_elements`, `clear_selection`, `get_elements` |
| Components     | `component.add`, `component.update`, `component.remove`, visibility      | `add_region`, `update_region`, `remove_region`       |
| Representation | `representation.set` (replace, add, restyle), `representation.remove`    | `set_look`, `remove_look`                            |
| Camera         | `camera.focus`, `camera.reset`                                           | `focus_view`; reset is new                           |
| Structures     | `structure.remove`                                                       | new                                                  |
| Generic        | `skills.list`, `skills.disclose`, `tools.disclose`, `remember`, `forget` | as they are                                          |

**Skills (loaded on request)**

| Skill         | Commands                                                                                                                                                                                                       | Stories                                            |
| ------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------- |
| selections    | `query_reference`, `bounds.get`                                                                                                                                                                                | ported nearly verbatim (§3)                        |
| style         | `preset.apply`, `background.set`, `assembly.set`, `model.set`, spin                                                                                                                                            | partly; presets redefined on Mol\*'s own           |
| volumes       | `volume.load` (box around a selection, detail to fit), `volume.set_look`, remove; volume streaming: turn it on for a structure, what it shows (around the selection or the camera, the whole cell), its levels | `load_volume`, `set_volume_look`; streaming is new |
| superposition | `structures.align`, `structure.transform`                                                                                                                                                                      | `align_structures`, `set_transform`                |
| clipping      | `clip.set`, `clip.remove`, `camera.set_view`                                                                                                                                                                   | `cut_open`, `remove_clip`                          |
| measurements  | `measurement.add` (distance, angle, dihedral, label), list, remove                                                                                                                                             | new (Stories draws distances as shapes)            |
| interactions  | `interactions.find`                                                                                                                                                                                            | `find_interactions`                                |
| snapshots     | `snapshot.add`, `list`, `apply`, `update`, `remove`                                                                                                                                                            | new (scenes are the story's)                       |
| export        | `image.export`, `state.export`, `structure.export`; an extension's own (`video.export` from mp4-export, §2.1)                                                                                                  | new                                                |
| code          | run code against the viewer (`@molstar/ai-toolscript`), MVS scripting (`mvs-code`)                                                                                                                             | `run_story_code`, the mvs-code skill, adapted      |

Not in the viewer mode: scenes, the story's text, animation timelines, narration, claims and articles, the helper
agents; these are the story mode's, which is Stories Untold with some changes, with its own tool set. Open for the
viewer mode: web and data (`fetch` with approval, reading a file, coloring by a table, MCP servers).

Code comes over for both modes: `any-toolscript` (a sandbox for code the agent writes, no Mol\* dependency) as
`@molstar/ai-toolscript`, and `mvs-code` (MVS scripting: the builder in the sandbox, convert, compare, its references)
next to the MVS packages, its imports moved to `@molstar/mvs-builder` and `@molstar/query-language`.

### 2.1 Extensions bring tools

A command is a tool by its metadata, so an extension's commands are tools wherever it is included: mp4-export's
`video.export`, with its params as the schema, needs nothing more. What an extension adds for the agent beyond that (a
skill: its guidance, which of its commands load together) is the agent's and ships only with it: an `agent` folder in
the extension (as in the modes, §8), which an app with the assistant registers with the extension. To design with the
contribution types (contexts.md §4).

Representation kinds, color themes and presets come from Mol\*'s registries when the tools are defined (Stories' enums
are hard-coded and would drift). The system prompt is rewritten from Stories' `guidance.ts`, keeping its general rules
(act through tools and batch calls, Focus is what the person points at, author numbering, read and fix failed calls,
never pass over a warning, decide, an honest last sentence).

## 3. Selections

The skill teaches the three forms, prefers `target` and `items` (exact) and keeps the query for distances and logic,
warns about the known traps ("X within 4 of Y" includes Y and nearby waters; Mol\*'s PyMOL reader fails on a group that
opens with two parentheses), and has every result say what it selected ("24 residues (180 atoms) in chain A: …") so the
model can check the count and chain. The references (PyMOL about 5k characters, MolScript about 5k, 160 symbol details
on demand) are generated from Mol\*'s symbol table and transpilers; the skill is about 3.5k tokens.

| Piece of `any-molql`                                                   | In Mol\* already                                               | To write                                   |
| ---------------------------------------------------------------------- | -------------------------------------------------------------- | ------------------------------------------ |
| text → MolQL (PyMOL, VMD, Jmol, MolScript), checking                   | `@molstar/query-language` (`compileScript`, validation issues) | the friendlier error text, as a wrapper    |
| MVS selector forms → MolQL                                             | `packages/mvs/runtime` (`rowsToExpression`)                    |                                            |
| keywords                                                               | `StructureSelectionQueries`                                    |                                            |
| `selectLoci`, `project` (counts and rows per level), `selectionBounds` |                                                                | small, over the live `Structure`           |
| `compactResidues` (`A:248,A:249` → `A:248-249`)                        |                                                                | a few lines                                |
| the generated references                                               |                                                                | a generator over `@molstar/query-language` |

The resolver (selection param → `Loci` on the live plugin) belongs to the viewer: commands and the command menu use it
too. The references and the guidance are the agent's and ship only with it.

## 4. The generic layer: `@molstar/ai-*`

| Package               | Holds (from Stories)                                                                                                                                                                                                                                                                                                                                | Lines |
| --------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----- |
| `@molstar/ai-tools`   | the handler pump (progress, cancel, failures), outcomes with `effect`, `createToolBindings`, nested tools, repeated reads (the state version injected), `toolSchedule`, argument normalization over a neutral field shape, a write hook; tool groups (skills): placement, keys, disclosure tools, catalog text; `./zod` for tools declared with zod | ~800  |
| `@molstar/ai-session` | `AgentSession` (host, client, runner, supervisor; turn, resume, stop, interrupted; structured hooks), the default drive policies, the turn budget and its extension, plain failures, the message composer, `runOnce` (the command menu's assistant), prices; `./compaction`                                                                         | ~1400 |
| `@molstar/ai-trace`   | `readTrace(reader, profile)` and its data types (turns, steps, attempts, the window a model saw, prefix breaks, usage); the profile gives the request marker, the window parts, what changed, prices                                                                                                                                                | ~650  |
| `@molstar/ai-fibers`  | §6                                                                                                                                                                                                                                                                                                                                                  | ~1700 |

None imports Mol\*, Prism, rxjs or React (the existing rules). zod validates, inside `ai-tools` only: a tool's input is
declared as a small field-shape interface (Prism maps PD to it) or, through `./zod`, as a zod schema; `ai-tools` turns
either into its validator and the JSON Schema, and normalizes arguments over the field shape.

Fixes on the way, from reading the code:

- **Effects.** Stories never sets a result's `effect`, so the harness records every failure as `unknown`, and its drive
  retries any transient failure. The layer carries `effect` through and the drive retries only what had no effect (or is
  idempotent).
- **Idempotency.** Stories is safe to re-run because its operations are desired state ("already so"). Prism commands are
  not: the handler keeps the run's key with its result and answers a `prior: possible` call from it.
- **Repair loop.** `not-enabled` and `no-handler` fall into "was not loaded; it is now" and would loop when a command's
  `when` turns it off; they get their own answers.
- **Changed, from data.** The trace finds what a call changed by parsing the result text; Prism's `CommandResult` gives
  it as data.

Prism writes: the commands with their metadata and the PD → JSON Schema converter, the selection param and resolver, the
state and focus accounts (and the request marker they share with compaction and the trace), the checkpointer, its drive
overrides (approval through the Assistant card; a denial is told to the model, which goes on or stops, §8), `AgentModel`
over `AgentSession`, keys and storage, and the policy for a turn interrupted by a reload (the viewer's state is gone, so
stop it rather than resume it).

## 5. The UIs: Inspector, with Trace, Tools and Limits

Stories opens it over the canvas (Alt+T). Prism opens it over the whole screen, a full-width panel (it needs the room),
with the three tabs, in `@molstar/prism-assistant` under a separate entry (`/inspector`) so the Assistant's tab and card
do not pull it in.

- **Trace**: an overview strip (a column per turn, a block per step; Steps or Time), a row list (turn, request, model
  call, tool call; nested runs indented) and a detail pane. A model call shows what it saw (the window's makeup:
  guidance, tools, state, conversation; its share of the window; what is new; why a cache prefix broke), its answer and
  a timing waterfall per attempt (waits, first output, streaming). A tool call shows its result, arguments, definition
  (a link to Tools) and timing. Re-read while a turn runs, not on every streamed token.
- **Tools**: groups (core, each skill, single) with their tools; per tool its description, input as fields and as JSON
  Schema, where it loads, its size in tokens, whether it is loaded now, why it is unavailable here. Prism adds `reads`,
  `asks` and enabled, which Stories has but does not show.
- **Limits**: rows with a label, what it bounds, a value, the default and a reset. Stories shows only its web-reading
  limits there and the turn budget as read-only text; Prism shows the turn budget per kind (used against the ceiling,
  editable), the context window and compaction trigger, and the model's thinking level. Credits are a display unit from
  an app-supplied price table, never a limit.
- Shared: one usage formatter (`k in (X% cached) · k out · credits`) for the card, History and the trace.

Style: no colored lanes; roles in quiet mono, state by icon and pattern (hatched retry, outlined block), danger for
failures, the accent only for the selected block and the primary action. `AgentModel` grows model-call steps with
attempts, the harness's full usage, and a parent run for nested runs.

## 6. Fibers

A fiber is an owner-scoped, append-only list of entries with idempotent commands: an entry streams byte-offset deltas,
then settles or aborts (interrupted, inactive, over a limit, closed); reads are windows by position and delta pages with
gap detection; everything is bounded by limits. The kernel (about 1,060 lines) is pure functions that decide a command
against the state the host gives it; Stories' host (650 lines) keeps a fiber per scope in a row store with a generation
check (memory, IndexedDB, files).

Stories keeps one fiber per story and writes one entry per turn itself (the prompt, the focus, the reply, what changed),
plus a small index for the history's tick rail; it does not stream deltas yet. That is what the harness lacks: its
journal is what the model sees and compaction drops turns from it, and its turn records hold no prompt or reply.

Port plan:

- `@molstar/ai-fibers`: the kernel as it is, the host, and its own storage port (the row-store shape plus a watch); not
  built on `HarnessStorage`, whose transactions are fixed to the harness's records. Imports nothing of the workspace.
- Its backends in `ai-components` (`src/fibers/`, next to the harness's in `src/harness/`): memory, IndexedDB and SQLite
  (`node:sqlite`, so a server or a headless host keeps conversations from the start). Stories' plain-files backend is
  not carried over; SQLite is the Node store.
- The conformance (8 cases, `node:assert` → Jest) as `@molstar/ai-conformance/fibers`, plus the cases Stories tests only
  in its own suite: an inactive stream aborted, durability across a restart.
- A turn observer in `ai-session` that writes an entry per turn, so any host gets a durable transcript; Prism's History
  reads it in windows.

## 7. Order

Reordered 2026-10-12: the generic packages and the self-contained ports first; the rest waits for a plan reconciled with
the Prism UI's and for the sequence viewer.

1. `@molstar/ai-tools` (done) and `@molstar/ai-session`, ported with their tests, the seams cut and the fixes above.
2. `@molstar/ai-trace` and the Inspector (Trace first, then Tools, then Limits), the views on mock trace data until a
   session feeds them.
3. `@molstar/ai-fibers`: the kernel, the host, the backends in `ai-components`, the conformance.
4. The PD → field shape converter (and so JSON Schema), for commands as tools.
5. **Pause**: reconcile this plan with the Prism UI's ([../../prism-ui/plan.md](../../prism-ui/plan.md)) and bring in
   the sequence viewer. Then, in an order decided there: Prism's side (command metadata, the selection param and
   resolver, the core commands of §2, the state and focus accounts, the checkpointer, `AgentModel` over `AgentSession`:
   a working assistant in the viewer mode), the selections skill and the viewer's system prompt, History over fibers,
   compaction in use, the other skills, `ai-toolscript` and `mvs-code`, the MCP server, evals.

## 8. Decided (2026-10-11)

1. **Tool schemas**: a neutral field shape (Prism maps PD to it) or a zod schema through `ai-tools/zod`; zod validates,
   isolated to `ai-tools`.
2. **Packages**: the generic layer as in §4. Each mode's agent parts (its tools' guidance and skills, its state and
   focus accounts, its drive overrides) live in an `agent` folder of the mode's package (`@molstar/prism-viewer/agent`,
   `@molstar/prism-stories/agent`), the only code there that imports `@molstar/ai-*` (an import rule checks it), and
   imported only by apps with the assistant. What both share (the Assistant's UI, the Inspector, `AgentModel` over
   `AgentSession`) is `@molstar/prism-assistant`. No package per mode's assistant.
3. **Selections**: one selection param type with the three forms, for commands and the command menu alike.
4. **Denial**: the model is told and goes on; when the denial means it cannot do what was asked, it says so and stops.
5. **Fibers' backends**: in `ai-components`, SQLite included.
6. **Code**: `@molstar/ai-toolscript` and `mvs-code` come over (the story mode needs both).
7. **Volume streaming** is a viewer tool; **extensions bring tools** through their commands, and skills through an
   `agent` folder (§2.1).
8. **Open**: web and data tools in the viewer mode.
