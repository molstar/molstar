# Commands, palette, and keybindings

A VS Code–style command layer: every user-facing operation is a named, searchable, context-aware command; the palette,
keybindings, menus, toolbar buttons, and scripts all go through it. This is what lets the new UI surface far less nested
UI: what used to be a card with a form becomes a command with arguments.

## 1. Relationship to `PluginCommands`

`PluginCommands` stay as they are: an internal, typed message bus with UUID ids that behaviors subscribe to. ui-next
commands are a separate, user-facing layer: string ids, titles, categories, arguments, enablement, keybindings. A UI
command's `run` may dispatch a `PluginCommand`, call an action function ([view-models.md](view-models.md) §6), or both.

The registry is framework-free and depends only on the plugin and the any-ui command core, so it can later move into
`@molstar/plugin` for headless use (scripting, automation, agents) without changes. That move is a separate decision.

## 2. Definition

Adapted from any-ui's `CommandDefinition<Ctx>`; the main change is that arguments are `ParamDefinition` params instead
of declarative definitions.

```ts
interface UICommand<P extends PD.Params = {}> {
  id: string;                         // 'structure.component.add', 'geo-export.export'
  title: string;                      // 'Add component…'
  category: string;                   // 'Structure', 'Selection', 'View', 'Export', …
  icon?: IconName;
  alias?: string | readonly string[]; // words that start the command with arguments: 'color by'
  shortcut?: KeyChord | KeyChord[];   // { key, mod?, alt?, shift? }; display derived per platform
  when?(ctx: CommandContext): true | string;     // reason string = listed, greyed, explained
  params?(ctx: CommandContext, values: Partial<PD.Values<P>>): P;  // arguments; recomputed when `reshapes` change
  reshapes?: readonly (keyof P)[];
  formats?: { [K in keyof P]?: ArgumentFormat }; // typed text input + completion for a field
  defaults?(ctx: CommandContext): Partial<PD.Values<P>>;
  presets?(ctx: CommandContext): readonly CommandPreset<P>[];     // search-only entry points
  preview?(values: PD.Values<P>, ctx: CommandContext): Preview | undefined;
  check?(values: PD.Values<P>, ctx: CommandContext): string | undefined;
  run(values: PD.Values<P>, ctx: CommandContext, signal: AbortSignal): CommandResult | Promise<CommandResult>;
  ui?: CommandUI<P>;                  // optional custom UI instead of the generated form
  group?: string;                     // parent group id ('select.query.structure-property'); top level if absent
  order?: number;                     // position within its group
}

/** A node in the command tree: a submenu in the palette, a menu, or a toolbar dropdown. */
interface CommandGroup {
  id: string;                         // 'select', 'select.query', 'select.expression'
  title: string;                      // 'Select', 'By query', 'Expression'
  parent?: string;                    // nesting; top level if absent
  icon?: IconName;
  order?: number;
  when?(ctx: CommandContext): true | string;
  /** Arguments shared by every command in the group (e.g. the selection mode); shown once in the group header. */
  shared?(ctx: CommandContext): PD.Params;
  /** Children computed from context or registries (entities, chains, ligands of the loaded structures). */
  children?(ctx: CommandContext): readonly (UICommand | CommandGroup)[];
}

interface CommandUI<P extends PD.Params> {
  placement?: 'inline' | 'dialog';    // under the palette input (default) or in a dialog
  render: ViewRenderer<CommandUIContext<P>>;   // dom, React, Preact, or an element (rendering.md §6)
}

interface CommandUIContext<P extends PD.Params> {
  invocation: CommandInvocation<P>;   // values, status, errors; set(), confirm(), cancel(), close()
  form: ParamsFormModel<P>;           // the generated form's model, to reuse fields inside a custom layout
  ctx: ReadonlyState<CommandContext>;
}
```

- **Arguments are PD params** (for now; the argument schema can be revisited once the palette is built). The same schema
  renders as an inline sentence in the palette when it is flat (select, number, text, boolean, color, and formatted text
  fields) and as a full form in a dialog when it is not (group, mapped, object-list…). Existing transformer and action
  params, extension params (`Mp4AnimationParams`, `GeometryParams`, …) become command arguments without rewriting.
- **Formats** are named argument types with `check(text) → value | message` and `complete(text) → options`. Built-in:
  `structure` (pick from the scene), `component`, `cell`, `selection` (residue list or MolScript, compiled with the
  existing `compileIdListSelection`/script support), `color` (names, hex, palette), `query` (registered selection
  queries), `snapshot`, `file`, `url`, `pdb-id`.
- **Custom UI** (`ui`) replaces the generated form for commands whose arguments are not a plain form: a multi-step flow
  (record ID → file list in Zenodo import; render → save in MP4 export), a picker with previews, a drawing tool. The
  custom UI drives the same invocation (setting values, confirming, showing progress and errors), so the command still
  has one id, appears in the palette, takes keybindings, and reports outcomes like any other. Without `ui`, the form is
  generated from `params`.
- **Results** (`message`, `detail`, `undoable`) become notifications; `undoable` offers Undo through
  `state.data.undo()`.
- **Groups** nest commands where there are too many to list flat (selection queries, presets, themes). Nesting is for
  browsing only: a search from the top level still matches every command at any depth and shows its path ("Select ›
  Structure property › Helix"), so nothing is harder to find. In the palette, Enter or → opens a group, Backspace on an
  empty input or ← goes back, and the input shows the path as a breadcrumb. A group can be opened directly
  (`palette.open-group` with the group id), which is how toolbar buttons and the selection bar reuse the palette scoped
  to a subtree. `shared` arguments are set once for the whole group (§10's selection mode) and passed to whichever child
  runs.

`CommandContext` is a `ReadonlyState` derived from view models:

```ts
interface CommandContext {
  plugin: PluginUIContext;                 // ui-next context (shell, commands, model, …)
  structures: readonly StructureNode[];    // hierarchy selection
  sceneTarget?: SceneNode;                 // node the user acted on (context menu, inspector)
  currentCell?: CellNode;
  selection: { elementCount: number; structureCount: number };
  focus?: FocusEntry;
  selectionMode: boolean;
  isBusy: boolean;
  activeView?: string;                     // shell focus, for scoped keybindings
}
```

## 3. Where commands come from

| Source                    | Examples                                                                                                                                                                                                                              | Lifetime                         |
| ------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------- |
| Built-in                  | Load file/URL/PDB/AlphaFold, undo, reset camera, screenshot, toggle selection mode, clear selection, add measurement, superpose, quick styles, snapshot add/next/previous/play, toggle view, change theme                             | Registered by `PluginUIContext`  |
| Contributed               | Extension contributions (geo-export, model-export, mp4 dialog, zenodo, assembly symmetry, …)                                                                                                                                          | While the contribution is active |
| Generated from registries | "Apply preset: …" (hierarchy and representation presets), "Color by: …" (themes), "Select: …" (queries), "Add representation: …", state actions ("Download structure", "Load G3D", "Download tunnels"), animations, background styles | Recomputed on `registry.version` |
| Contextual                | `actions.fromCell(currentCell)` as "Apply to <cell>: …"; node actions in scene/tree context menus                                                                                                                                     | While the cell/node is current   |
| View-scoped               | Tree navigation, palette navigation, dialog confirm                                                                                                                                                                                   | While the view is mounted        |

Generated commands are the largest win: today these choices are buried in nested menus inside cards. The palette lists
them flat, fuzzy-searchable, with `when` explaining why one is unavailable ("Select a structure first").

## 4. Registry

`CommandRegistry<CommandContext>` is any-ui's `CommandManager` with these changes:

- `register(commands, owner?)` returns an unregister function (unchanged); contributions and views register by mount
  scope.
- `query(text)` adds category and prefix filtering (§5) and ranks recents and contextual commands first.
- `open(id, tail?, preset?)`: a command with arguments opens an invocation (inline or dialog form); one without runs
  immediately. One run at a time, explicit cancel (unchanged).
- `execute(id, values)`: programmatic run without UI (scripts, tests, toolbar buttons with fixed arguments).
- `outcomes`: non-replaying event stream (unchanged), consumed by `NotificationsModel`.
- Recents persisted through an injected store (localStorage in apps, none in embeds).
- The `assist` port stays in the types with no implementation.

`CommandInvocation` keeps its state machine (`editing → running → done | failed | cancelled`), tail reading, and preview
sink; its value store is a `ParamsFormModel` instead of a declarative binding context.

## 5. Palette

- Opens with ⌘K / Ctrl+K and ⇧⌘P / Ctrl+Shift+P (both; the first is common on the web, the second familiar from VS
  Code), and from a toolbar button for touch.
- Modal overlay over the viewport; the canvas stays live (previews draw in 3D).
- Prefixes: none = commands; `@` = go to a scene node (structure, component, volume; focuses and selects it); `#` =
  snapshots; `:` = select residues by expression (`:A 10-20`, `:HEM`); `?` = help on prefixes and keybindings.
- Rows: title with matched letters, category, keybinding, disabled reason. Recents first on empty input.
- Typing an alias followed by text opens that command with the text read into its arguments ("color by chain-id", "show
  as cartoon").
- Flat arguments render inline as a sentence; Tab moves between fields, Enter runs, Esc closes. "Open as form" moves the
  invocation into a dialog for complex params. A command with `ui` shows its custom UI in the same place (inline or
  dialog, per `placement`).
- ARIA combobox: the input keeps focus, rows are `role="option"`.

## 6. Keybindings

A `KeybindingService` on the plugin root element:

- **Scopes**: `global` (anywhere in the plugin), `viewport` (canvas focused or hovered), `view:<id>` (a focused view),
  `palette`, `dialog`. The most specific scope wins; text-entry targets only receive ⌘/Ctrl chords (any-ui's rules for
  IME, repeat, and `e.code` fallback are kept).
- **Coexistence with canvas bindings.** Trackball keys (W/A/S/D/R/F/Q/E…) and behavior bindings (`I`, `O`, `T`, `G`,
  Ctrl+arrows) keep working through `InputObserver`. The service reads those bindings and reports conflicts in the Help
  view instead of overriding them. Later, behavior key bindings can be re-expressed as commands with `viewport` scope.
- **Platforms.** Chords are written once with `mod` (⌘ on macOS, Ctrl on Windows and Linux); `formatChord` in
  `commands/keys.ts` shows them per platform (`⌘K` / `Ctrl+K`, `⇧⌘P` / `Ctrl+Shift+P`).
- **Defaults** are part of command definitions; the UI spec can override them (`keybindings: [{ command, key }]`). User
  remapping is out of scope for the first version.
- Help is generated from the registry, behavior bindings, and trackball bindings, so it is always complete.

## 7. Menus and buttons are commands

A command can also just open UI: a panel, a tab, or a panel at a section (`Export animation (MP4)…` opens the Export
panel at its `animation` section). Features that live in a panel are thus still found by the palette, and extensions add
sections and commands rather than buttons.

Context menus (scene node, state-tree cell, viewport right-click), toolbar buttons, and view headers reference command
ids, optionally with fixed arguments:

```ts
{ command: 'structure.component.add', args: { structure: node.ref }, group: 'add' }
```

A button is enabled per the command's `when`, and its tooltip shows the title and keybinding. This keeps every action
discoverable in the palette and avoids duplicate logic.

## 8. Built-in command map (first pass)

| Today                              | Commands                                                                                                                                                                                                                                           |
| ---------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Home tab actions                   | `data.open-files`, `data.download-structure` (presets: PDB, AlphaFold DB, PDB-IHM, URL), `data.download-volume`, generated state actions                                                                                                           |
| Source card                        | `structure.apply-hierarchy-preset`, `structure.set-model-index`, `structure.set-type`, `structure.select-in-scene`                                                                                                                                 |
| Components card                    | `structure.apply-preset`, `structure.component.add`, `.add-representation`, `.set-color`, `.remove`, `.toggle-visibility`, `.rename`, `structure.component-options` (dialog)                                                                       |
| Quick styles                       | `view.style.default`, `view.style.illustrative`, `structure.preset.cartoon` / `spacefill` / `surface`                                                                                                                                              |
| Selection toolbar                  | The `select` and `selection` groups, §10                                                                                                                                                                                                           |
| Measurements                       | `measure.distance`, `.angle`, `.dihedral`, `.label`, `.orientation`, `.plane`, `measure.options`                                                                                                                                                   |
| Superposition                      | `superpose.chains`, `.tm-align`, `.ligands`, `.atoms`, `.db`                                                                                                                                                                                       |
| Procedural animation               | `animation.wiggle.dynamics`, `.uncertainty`, `.clear`                                                                                                                                                                                              |
| Volume / particles                 | `volume.add-representation`, `volume.set-isovalue`, `volume.streaming.enable`, `particles.*`                                                                                                                                                       |
| Viewport controls                  | `camera.reset`, `.reset-zoom`, `.orient-axes`, `.reset-axes`, `view.toggle-panels`, `view.expand`, `view.fullscreen`, `view.illumination.toggle`, `view.xr.toggle`, `view.settings`, `screenshot.copy`, `screenshot.download`, `screenshot.dialog` |
| Snapshots / animation / trajectory | `snapshot.add`, `.next`, `.previous`, `.play`, `.save` (molj/molx), `.open`, `animation.play`/`.stop`, `trajectory.next`/`.previous`/`.first`                                                                                                      |
| State tree                         | `state.undo`, `state.remove`, `state.toggle-visibility`, `state.set-current`, contextual actions                                                                                                                                                   |
| Focus                              | `focus.target` (ligands etc. as completions), `focus.clear`, `focus.recenter`                                                                                                                                                                      |
| Shell                              | `view.show:<view-id>`, `view.theme`, `palette.open`, `help.keybindings`                                                                                                                                                                            |

## 9. Extensibility

Extensions extend the command layer in four ways, all through their `UIContribution` ([shell.md](shell.md) §5) or the
plugin registry, never by patching built-in definitions:

| Mechanism                | What an extension adds                                                                                                         | Example                                                                                                                         |
| ------------------------ | ------------------------------------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------- |
| `commands`               | New commands, with `group` placing them anywhere in the tree, including built-in groups                                        | `geo-export.export` under `export`; `assembly-symmetry.enable` under `structure`                                                |
| `commandGroups`          | New groups, top level or nested under built-in ones                                                                            | `volseg` group with its segment commands; `select.validation` under `select.query`                                              |
| Registry entries         | Nothing UI-specific: registered selection queries, presets, themes, state actions, and animations appear as generated commands | rcsb's `hasClash` query under Select › By query › its category; dnatco presets under Apply preset; `DownloadTunnels` under Data |
| `formats`, `menus`, `ui` | Argument types, context-menu and toolbar entries that reference commands, custom invocation UI                                 | A `segment` format with completions; a toolbar button for the PAE plot; the MP4 export flow                                     |

- **Registration follows activation.** A contribution's commands and groups are registered when it is activated and
  unregistered when it is deactivated (a run of theirs is cancelled, their open form closed), the same rule as views.
- **Ids are namespaced** by the contribution (`<contribution>.<name>`); a duplicate id is an error at activation.
  Placing a command in another contribution's group by id is allowed and is how built-in groups are extended.
- **Apps decide what is visible.** The UI spec can hide commands or whole groups and rebind keys
  (`ui.commands: { hide: ['superpose.*'], keybindings: […] }`), so an app shows a subset without forking definitions.
- **Script languages** follow the plugin: importing a transpiler module (`@molstar/model/script/transpilers/pymol`)
  makes its language available (`Script.getAvailableLanguages()`), and the Expression group lists one command per
  available language. A custom language registered with `registerTranspiler` appears the same way.

## 10. Selection commands

Today's selection toolbar (selection mode only) has four mode buttons (union, subtract, intersect, set) that each open
the same nested query menu, a residue-list helper, a Theme panel whose form mixes target, action, and representation,
plus create/remove component, undo, help, and close. It is already command-like but each piece is its own UI. In ui-next
it becomes two command groups, used from the palette, from context menus, and from a slim selection bar.

### 10.1 `select`: changing what is selected

One group with a **shared `mode` argument** (`add` | `remove` | `intersect` | `replace`, default `add`) shown once as a
segmented control in the group header, instead of four buttons each opening a copy of the menu. Typed tails accept the
words too ("select helix remove"). Modifier keys on Enter choose the mode for one run: Enter adds, Shift+Enter removes,
Alt+Enter intersects, ⌘/Ctrl+Enter replaces.

```
Select ›                                  shared: mode = add | remove | intersect | replace
  All · Clear · Invert
  By query ›                              generated from the selection-query registry, grouped by `category`
    Structure ›            Polymer, Protein, Nucleic, Ligand, Water, Ion, Lipid, Carbohydrate, …
    Structure property ›   Trace, Backbone, Sidechain, Sidechain with trace, Helix, Beta strand/sheet, …
    Bond property ›        …
    <extension categories> e.g. Validation › Clashes; Confidence › pLDDT > 70; Membrane › Transmembrane
  In the loaded structures ›              computed per structure from the scene (dynamic children)
    Entities ›             polymer and carbohydrate entities (today's "Polymer/Carbohydrate Entities")
    Chains ›               (new)
    Ligands ›              ligands and non-standard residues
    Elements ›             element symbols present
  Residues…                               form: "A 10-20, B 5" with id type (auth/label), via compileIdListSelection
  Expression ›                            one command per available script language, text argument with syntax check
    PyMOL…  ·  VMD…  ·  Jmol…  ·  MolScript…
  Around the selection…                   (new) radius in Å, whole residues on/off
  Expand to ›                             (new) residues · chains · entities · structures
```

- **Preview while browsing.** Moving through rows highlights the row's selection in 3D through the preview sink
  (`lociHighlights`), so the user sees what "Sidechain" or "Chain B" means before choosing. Typing an expression
  previews it once it parses; parse errors show inline with the language's message.
- **Counts on demand.** The highlighted row shows how much it would select ("Helix: 412 residues"), computed for that
  row only, cancelled when the highlight moves.
- **Granularity** is not part of `select`; it is a picking setting (`selection.granularity`, below) shown in the
  selection bar.
- The `:` palette prefix is a shortcut into Select › Residues or Expression: `:A 10-20` reads as a residue list,
  `:pymol resn HEM` as a PyMOL expression.

### 10.2 `selection`: acting on the current selection

Commands whose target defaults to the current selection; `when` greys them out with "Nothing is selected".

```
Selection ›
  Focus  ·  Zoom to  ·  Clear  ·  Undo last change
  Color…                 color (named colors, recent colors, palette swatches, hex/RGB), representations (default all)
  Set transparency…      value, representations
  Set material…          material params, representations
  Clip…                  clipping params, representations
  Set emissive…          value, representations
  Reset color  ·  Reset material
  Create component…      label, representation (registered types)
  Remove from components
  Measure ›              Distance · Angle · Dihedral · Label · Orientation · Plane  (from the selection history)
  Granularity ›          Element · Residue · Chain · Entity · Model · Structure …  (radio; picking setting)
  Selection mode         toggle
```

- Today's Theme panel (`selection` + `action` mapped select + `representations`) splits into one command per action.
  Each has a small flat form that fits inline in the palette, and the search finds "transparency" directly instead of
  through an Action dropdown.
- **Color** is a `color` argument format whose completions are swatches (recent colors first, then the palette shown
  today) and which accepts names and hex. Lighten/darken and RGB editing live in the color field's popover. When the
  field is active, moving through swatches previews the color on the selection.
- The commands reuse `StructureComponentManager.applyTheme`, `add`, and `modifyByCurrentSelection` through
  `model/actions/selection.ts`, so the inspector and context menus run the same code.

### 10.3 The selection bar

In selection mode the viewport shows a slim bar instead of today's toolbar:

```
[ Residue ▾ ]  [ + − ∩ = ]  [ Select…  type a query, residues, or :expression        ]  [ 🖌 ] [ ⧉ ] [ ↶ ] [ ✕ ]
 granularity    mode         the `select` group as an inline palette                      color  comp.  undo  exit
```

- The input is the palette scoped to the `select` group, anchored to the bar: typing filters, Enter selects, the
  breadcrumb shows nesting. The mode control is the group's shared argument.
- The buttons open `selection.color`, `selection.create-component`, run undo, and leave selection mode; each shows its
  keybinding in the tooltip.
- Help moves to the Help dialog's selection section; the bar's `?` is dropped.
