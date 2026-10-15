# The skeleton: contexts, modules, modes, commands, history

> **Reconciled (2026-10-15):** the Prism viewer's [spec](../spec/README.md) holds. No modes (§1, `PrismMode`), no
> module registry or contributions framework (§3–4); one mode-free viewer composed from `PrismViewerSpec`
> ([../spec/viewer-spec.md](../spec/viewer-spec.md)); commands and tools over one API, undo in the viewer's history
> manager ([../spec/sequence.md](../spec/sequence.md) §7). Still valid: origins (§5), the plugin boundary (§7), results
> and state (§8).

The contracts every other design builds on, frozen at the review gate ([../review.md](../review.md): B1, M1–M5, M16–M23)
for the first implementation. Where another document says otherwise, this one holds. Informed by how the agent harness
and the sequence viewer must be hosted ([../analysis/harness-and-seqstar.md](../analysis/harness-and-seqstar.md)).

> **Direction changed (David, 2026-10-14), to be reconciled with this document later.** Prism takes a library
> approach: prism-ui is the runtime for generic Mol\* UIs (the shell, the common controls, the canvas-ink abstraction,
> the agentic wrappers, commands), and the viewer and the story mode are opinionated compositions on it, each owning
> its features in its own code. There is no extension framework: where this document describes modules (`ModuleKey`,
> `PrismModule`, the registry, §3) and contributions (`UIContribution`, a mode's `contributions`, §4) as the way
> features are composed, the compositions own them directly instead, and prefer fewer packages. The rest (one mode per
> app, the context's lifetime, origins, the plugin boundary, results and state) still holds until the reconciliation.
> A composition may still take extensions, as the v5 Viewer does through its composition spec: the Prism viewer will,
> on its own terms; prism-ui, the library, takes no position on how.
> The first design under the new direction: [../../seq/prism-sequence.md](../../seq/prism-sequence.md) §2.

## 1. An app, its mode, its lifetime

- **One mode per app.** A Prism app runs one mode: the viewer (`@molstar/prism-viewer`, Mol\*'s own UI) or the story
  mode (`@molstar/prism-stories`). The mode is chosen before any context exists: a splash when the app opens (unless the
  host names one), Home's choice while nothing is loaded, or "Switch mode…", which says what goes (what is loaded;
  stories are kept in the browser) and then **disposes the app and creates a new one** in the other mode. There is no
  switching inside a running app, so no state, model or registration ever crosses from one mode to another.
- **Builds per mode.** As an app is composed of the mode it runs, a build can include one mode only (a story reader, a
  viewer without stories).
- **The remembered choice** (the last mode, the layout, the theme) goes to a host-given store (`PrismStore`, §8).

```ts
createPrism({ spec, mode }): Promise<PrismContext>              // headless: Node tests, the agent's tools, MCP
createPrismUI({ target, spec, mode? }): Promise<PrismApp>        // no mode: the splash first
interface PrismApp { readonly ctx: PrismUIContext; dispose(): Promise<void> }  // the caller owns it
```

Lifecycle of a context: `init()` (the plugin initialised; modules created in dependency order; contributions activated)
→ `attach(root, canvasHost)` (UI only: `plugin.mountAsync` into the canvas host; the canvas sources start) → `detach()`
→ `dispose()` (contributions deactivated, modules disposed in reverse order, the plugin disposed). Exactly one element
per context hosts the canvas.

**Embedding.** `<prism-viewer mode="viewer|story">` is the top-level element; it creates and owns its app, and has a
shadow root (host-page CSS stops at it; our `--msp-*` variables still theme it). Everything inside is light DOM; the
style registry adopts the sheets into that root; delegated events cross it (rendering.md §4.1). The embeddable parts
(`<msp-viewport>` and others, shell.md §7) stay light DOM, for hosts that style around them.

## 2. PrismContext: what every mode shares

```ts
class PrismContext {                                // @molstar/prism: no DOM, runs in Node
  readonly spec: PrismSpec;
  readonly mode: string;                              // the mode this app runs
  readonly plugin: PluginContext;                     // for sources, actions and composition code only (§7)
  readonly sources: PluginSources;                    // the frame-wide ones: tasks, log, canvas, hover, selection
  readonly commands: CommandManager<CommandContext>;  // the one registry (§4)
  readonly status: StatusModel;                       // messages, tasks, chips
  readonly history: HistoryModel;                     // runs, and what modules add (agent turns) (§5)
  readonly settings: SettingsModel;                   // the app's settings; modes and extensions add groups
  readonly modules: ModuleRegistry;                   // the mode's, the extensions', the agent's models (§3)
}

class PrismUIContext extends PrismContext {           // @molstar/prism-ui
  readonly spec: PrismUISpec;
  readonly shell: ShellModel;                         // Prism's own layout
  readonly keybindings: KeybindingService;
  readonly contributions: ContributionRegistry;       // the active UI contributions
  readonly theme: ThemeModel;
  readonly fullscreen: FullscreenModel;
  navigate(target: { tab?: string; overlay?: string; section?: string }): boolean;  // the port `opens` uses (§4)
}
```

Nothing mode-specific is a field: the scene, the selection model, snapshots, the state tree, a story document are
modules (§3). Views and commands get the models they work with, never the context; only composition code (a mode's or
contribution's factory, `createPrism`/`createPrismUI`) reads the context.

## 3. Modules: how a mode, an extension or the agent attaches models

```ts
interface ModuleKey<M> { readonly id: string }        // a typed key: `export const Scene = moduleKey<SceneModel>('scene')`

interface PrismModule<M> {
  key: ModuleKey<M>;
  requires?: readonly ModuleKey<any>[];               // created after these
  create(ctx: PrismContext): M;                       // sources, other modules (ctx.modules.get), the plugin via sources
  dispose?(model: M): void;
}

interface ModuleRegistry {
  get<M>(key: ModuleKey<M>): M;                       // throws, naming the key, when absent
  find<M>(key: ModuleKey<M>): M | undefined;          // an optional module (the agent)
}
```

- A module lives as long as the app: created at `init()`, disposed at `dispose()`. With one mode per app there is no
  other lifetime to manage (review Q1).
- The viewer's contracts and their plugin-backed implementations stay in `@molstar/prism-viewer`, which exports its
  modules (and a headless entry, `@molstar/prism-viewer/headless`, so `createPrism` tests can include them). The story
  mode's document and player are its modules; the agent's session is the agent's.
- A module that owns an outside context (the sequence viewer's `SeqViewerContext`) creates it in `create` and disposes
  it in `dispose`.
- Lazy parts (a panel's view model, an inspection per node) are created by the module that owns them and retired with
  it, or when what they follow goes (a keyed node); view-models.md §4.3.

## 4. Contributions, the mode, commands

```ts
interface PrismContribution {                         // headless: activated by PrismContext.init()
  id: string;
  modules?: readonly PrismModule<any>[];
  commands?(ctx: PrismContext): readonly CommandDefinition<CommandContext>[];  // close over the modules they use
  commandGroups?: readonly CommandGroup[];
  formats?: readonly [name: string, CommandFormat<CommandContext>][];
  commandContext?: CommandContextSlice<any>;          // what this contribution adds to the command context
  replaces?: readonly string[];                       // ids it replaces on purpose (§6)
}

interface UIContribution extends PrismContribution {  // activated by PrismUIContext: the headless part, then the views
  tabs?, overlays?, sections?, toolbar?: …;           // shell.md §5, with a reactive `when` (§6)
  paramEditors?, settings?, structureRows?: …;
}

interface PrismMode {
  id: string; title: string; icon: IconKey; description?: string;
  contributions: readonly (PrismContribution | UIContribution)[];   // headless use activates the headless parts only
  start?: string;                                     // the tab shown at start
  canvas?: CanvasParts;                               // content, overlays, playback, context, lead, aside (render functions)
  empty?(ctx: PrismContext): Observable<boolean>;     // nothing loaded (the wordmark, Home's mode choice)
  bindings?: readonly InputBindingGroup[];            // the Help tab's mouse and touch section
  account?(ctx: PrismContext): ModeAccount;           // the agent's state block and focus (§5)
  checkpointer?(ctx: PrismContext): Checkpointer;     // how the agent's turn is undone in this mode (§5)
}

interface PrismSpec { plugin: PluginSpec; contributions?: readonly PrismContribution[]; persist?: PersistSpec }
interface PrismUISpec extends PrismSpec { ui?: { contributions?: readonly UIContribution[]; layout?; theme?; … } }
```

The frame's own contributions (Log, Settings, Help, the theme, the command menu) are the spec's; the mode's come with
it. `DefaultPrismContributions()` (core) and `DefaultUIContributions()` (prism-ui) list the built-in ones.

**Commands (one registry).** `PrismContext.commands` is the only `CommandManager`, in `@molstar/prism`; `PluginCommands`
stay the plugin's internal bus and the classic UI's. The frame registers its commands once; the mode's and the
extensions' come with their contributions.

- **The command context is composed**: a frame base `{ busy, mode }` plus slices that contributions add, read by key
  (`ctx.slice(SceneContext)`). It has no `prism` field; commands close over their modules.
- **Running from code is the one run path**:
  `execute(id, args, { source, signal, turn? }): RunHandle { id; done: Promise<CommandResult>; cancel() }`. The menu,
  sentences, keybindings, the agent and MCP all go through it. Writes are queued in order; reads run alongside.
  Arguments are validated by rebuilding `params(ctx, args)` until the `reshapes` fields settle, then PD validation,
  formats and `check`.
- **Results** say what happened:
  `CommandResult { message?: string; effect: 'none' | 'applied' | 'unknown'; summary?: string }` (the summary is what a
  model reads). An assistant's answer is a separate `answered` outcome, never a run.
- **Metadata for tools**: `read` (no change; runs alongside, may be repeated), `approval` (asks first), `capabilities`
  (what it may touch; decides whether an agent outside the app may call it), `ui: true` (needs the UI: skipped headless
  and from the tool list). `opens` goes through `PrismUIContext.navigate`; without the UI, or when the target is not
  registered, `when` is false and says why.
- **As tools**: an id, a description, and the params as JSON Schema (a PD → JSON Schema converter; a `reshapes` field
  becomes `oneOf`, context-dependent selects become strings with a run-time format), plus a loose argument normalizer
  (as Stories does). The tool list is fixed for a turn; `when` is answered with its reason when a tool is called.

## 5. Origins, history, undo, the agent's account

- **Every change says who made it.**
  `RunSource = { kind: 'person' | 'agent' | 'mcp' | 'behavior' | 'link'; via?: 'assistant' }`: a run the person confirms
  after the assistant proposed it is the person's, `via: 'assistant'`. The origin travels from `execute` to the history,
  to the plugin transactions an agent run commits (tagged with its turn), and to the sources' events, so a link (the
  sequence viewer's) skips its own echoes.
- **History in core.** `RunRecord { id; commandId; title; args; source; turn?; status; startedAt; durationMs; result }`;
  `HistoryModel` holds runs, with an extension point for item kinds (the agent adds turns). Without the agent, the app
  shows Log over `StatusModel`; with it, the agent's contribution replaces Log with History (`replaces: ['log']`).
- **Undo.** The person's Undo undoes the last step. An agent turn is one step: the mode's `Checkpointer` takes one
  checkpoint before the turn (the viewer: a plugin state snapshot; the story mode: the document's inverse operations)
  and releases it at the next change not tagged with that turn (navigation, reads and answers do not release it). An
  `undoable` command commits exactly one transaction, which `execute` enforces.
- **While a turn runs**, the person's own commands still run (attributed to the person); in the story mode, document
  edits wait for the turn to end (as Stories does), so a turn stays one step.
- **The agent's account.** The mode gives the prompt a compact state block and the person's focus
  (`ModeAccount { state(): string; focus(): Focus }`): what is loaded, picked, selected, shown.

## 6. Registration and visibility rules

- **One rule for every registry** (tabs, overlays, sections, toolbar items, commands, modules): a duplicate id throws at
  activation, naming both owners, unless the newcomer declares `replaces`. A replacement keeps a shadow stack: undoing
  it restores the previous entry. `ui.commands.hide` and keybinding overrides apply in the manager, so the menu, keys,
  `execute` and the tool list all respect them.
- **One reactive `when`** for tabs, sections, overlays and toolbar items: `when?: (ctx) => Observable<boolean>`,
  subscribed once per item on activation (it follows the data, not only registrations).
- **At runtime, and lazily.** Contributions are not only the spec's: `PrismContext.contributions.add(c)` activates one
  in a running app (its modules, commands, tabs, sections, toolbar items, the agent's skills) and returns its undo;
  `remove(id)` retires it. A **lazy contribution** ships only a stub: an id, a title, what it shows before it is loaded
  (a command, an Export section, a toolbar item, a skill's summary for the agent) and `load: () => import(…)`. The first
  use (the command run, the section opened, the skill disclosed) loads it as a task in the status column, adds the real
  contribution in the stub's place (`replaces`), and replays what triggered it; a failed load says why and keeps the
  stub. Heavy features (video export, narration, the code skill) then cost nothing until someone uses them; the app's
  build has to split them into chunks (ESM with splitting, [../plan.md](../plan.md) Phase 8).
- **Render context**: a contribution's views get `{ shell, models, focus, subject, dock }`, where `models` are what its
  modules made; they are retired with the contribution.

## 7. The plugin boundary

- **By folder role, in any Prism package**: only `sources/` and `actions/` folders import `@molstar/plugin` (and
  `@molstar/model`/`@molstar/graphics` for the plugin's data); in the viewer, `api/` and a feature's
  `features/<name>/operations.ts` are its actions (operations.md §1). That holds in core, in a mode, and in an extension's
  prism entry; one reusable check runs on every package that has Prism modules.
- **Commands call actions** (never a `PluginCommand`); an extension's controllers that take a `PluginContext` (geometry
  export, MP4, volume streaming) are wrapped in actions.
- **Layout**: Prism never sets `plugin.layout`; an action calls `plugin.handleResize()` when the shell's canvas inset
  changes.
- **Not used**: `PluginLayout` and `leftPanelTabName`, `PluginCommands` for UI work, toasts, the UI-only maps,
  `customParamEditors`, behaviors that only drive the classic UI (shell.md §6).

## 8. Models and async work

- **Results and state are two channels.** A model method that does work returns `Promise<ActionResult>` (never rejects;
  also posted to status): the outcome for whoever asked (a command, the agent, a test). The state changes only through
  the reactive state, fed by the sources from the plugin's events, never set from a promise's result: one source of
  truth whoever made the change. Long operations also show as `busy` states.
- **Optimistic when the future is known.** A change is projected at once when its target state is fully determined by
  the request and the source will later report exactly that state (parameter edits that run a transform). Synchronous
  changes (a visibility flag) and changes whose result cannot be predicted are not projected.
- **Persistence**: `PrismStore { get(key); set(key, value) }` (localStorage in the viewer, IndexedDB for stories), given
  by the host; `persist: { store, key }` keeps `{ layout, mode, theme }`.

## 9. State

Not built: everything above is the contract for Phase 4–6. The mocks follow it in shape: one mode per app
(`PrismAppStory`, the splash), the frame's contributions, a mode's canvas parts and bindings, Help as a tab. They do not
yet have modules, `execute`, origins, or the core history.
