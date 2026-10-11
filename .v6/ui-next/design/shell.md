# Shell, layout, and the UI spec

The shell is the frame around the viewport: the tab panel, optional docks, canvas overlays, the palette, and
notifications. It replaces `PluginContext.layout` (`PluginLayout`) and the fixed five-region `Layout` for ui-next, and
it is where contributions (built-in and extension UI) are placed.

## 1. Principles

- **Canvas first.** The viewport fills the shell; persistent chrome is the exception.
- **Summoned, not resident.** Short tasks (saving state, a screenshot, log history, help) open as temporary overlays on
  the canvas and close when done. Docked UI is for what the user keeps looking at while working.
- **One home for persistent UI.** A tabbed panel on the right holds Home, Structure, State tree, Settings, and extension
  tabs. Left, top, and bottom docks exist but are optional and empty by default; users and apps can move tabs there.
- **Appear when relevant.** Home shows when nothing is loaded; trajectory controls only for multi-model trajectories;
  snapshot controls only when there are snapshots; the selection bar only in selection mode; the log as one fading line.
- **Lean views.** Commands carry the operations ([commands.md](commands.md)); views show state and the common actions,
  not a form for everything.
- **Same spirit, not the same screens.** The workflow (load, style, select, measure, save) and the capabilities stay;
  the layout, density, and icons are redesigned.
- **Modular and embeddable.** Every tab, overlay, toolbar item, and section is a contribution; the same pieces work in
  the full app, a bare viewport, or single elements in a host page.

## 2. Structure

```
┌──┬──────────────────────────────────────────────────────┬───────────────────┬──┐
│  │ 1TQN ▾  Chain A ▾                                 ⌃  │                   │⌂ │
│  │ MALYGTHSHGLFKKLGIPGPTPLPFLGNILSYHKGFCMFDMECHKKYG…    │  right panel      │⧉ │
│  ├──────────────────────────────────────────────────────┤  (active tab)     │≡ │
│ l│ ┌ snapshots ─────────────┐                     [⌘]  │                   │⚙ │
│ e│ │ ‹ 3/12 ›  Title     ▶  │                     [⟲]  │  Home ·           │+ │
│ f│ │ markdown description…  │                     [📷] │  Structure ·      │  │
│ t│ └────────────────────────┘                     [💾] │  State tree ·     │  │
│  │                  canvas      ┌ hover label ┐    [◎]  │  Settings ·       │  │
│ (│                        ⌖ ────┤ GLU 213     │    [💡] │  extension tabs   │  │
│ o│                              └─────────────┘    [XR] │                   │  │
│ p│   ┌ overlay card (save state / screenshot / …) ┐ [⛶]  │                   │  │
│ t│   └────────────────────────────────────────────┘      │                   │  │
│ .│        ┌ ⧉ 4/20 |◀ ◀ ━━━━●━━━━━━ ▶ ▶| ⏯ ┐ timeline   │                   │  │
│ )│ ⓘ Loaded 1TQN · ⟳ Surface 40% ✕ · ● 12 residues       │                   │  │
│  ├──────────────────────────────────────────────────────┤                   │  │
│  │ bottom dock (optional, canvas width)                  │                   │  │
└──┴──────────────────────────────────────────────────────┴───────────────────┴──┘
  also: palette (modal), dialogs, notifications, drag-and-drop target
```

- **Floating docks.** The top dock (sequence) floats too: a glass card along the top of the canvas, between the left
  edge and what floats on the right, or a small glass pill when collapsed; the canvas runs under everything, so the
  sequence reads as part of the visualization. The snapshot panel sits below it (`--msp-canvas-top-inset`).
- **Right panel, floating.** By default the right dock floats over the canvas instead of taking space beside it: the
  canvas gets the full width, the active tab's content is a glass card (more opaque than the small controls, for
  reading), and a rail at the right edge carries the tab strip and the canvas toolbar as two glass groups. Collapsed,
  only the rail stays. Resizable; `setFloating(dock, false)` docks it beside the canvas instead; compact mode makes it a
  bottom sheet. The panel card ends above the bottom bar, so the bar keeps the full canvas width up to the rail
  (`--msp-canvas-rail-inset`, `--msp-canvas-bottom-inset`); the sequence card stops at the open panel
  (`--msp-canvas-right-inset`). If the bar still gets crowded, next steps are a compact timeline pill that expands on
  hover and messages/tasks condensing to icons.
- **Optional docks.** `left`, `top`, `bottom` take tabs moved there by the user or placed by the app spec (mesoscale
  explorer puts its entity tree on the left). The sequence is the only default occupant of `top`. Side docks run the
  full height; top and bottom docks span the canvas width only. A top or bottom dock with one tab has no tab strip: its
  content carries its own header (the sequence: structure and chain selectors, collapse), and collapsed it is a slim bar
  that expands on click.
- **Canvas toolbar.** A vertical strip of command buttons (palette, reset camera, screenshot, save state, selection
  mode, illumination, XR, expand/fullscreen, plus contributed items), each gated by its command's `when` and
  `PluginConfig.Viewport.*`. A button either runs its command or opens an overlay anchored to it.
- **Snapshot panel** (top-left of the canvas, §4.3): stepper and description. Nothing is stacked under it, so a
  description that grows or collapses moves nothing else.
- **Bottom of the canvas**: messages and tasks on the left, the current selection and focus as chips in the centre (the
  work in hand), the timeline (trajectory model index, animations; like a media player) on the right. The bar is its own
  size container: when what floats on the right narrows it, the chips take a row of their own above, and on a very
  narrow canvas everything stacks.
- **Glass and canvas theme.** Everything over the canvas uses flat glass (translucent fill, blur, hairline edge, no
  shadows) and takes its theme from the canvas background automatically, not from the UI theme, so it stays visible (as
  stories-untold's stage ink does). `model/canvas-theme.ts` computes it from the Canvas3D props by relative luminance
  (threshold 0.179, where light and dark text have equal contrast): the background color, a gradient's colors averaged
  over it by opacity; images and skyboxes count as the background color for now; a transparent canvas uses the host's
  fallback. The viewport model follows `events.canvas3d.settingsUpdated` and feeds `canvasTheme`.
- **Hover label** next to the cursor (§4.4).
- **Status line** (bottom of the canvas): the latest log message (fades after a few seconds; errors stay until
  dismissed), running tasks with progress and abort, and chips for the selection and focus (click to focus, ✕ to clear).
  Clicking the message opens the log overlay.
- **Overlays.** Cards over the canvas, anchored to the button that opened them, one at a time, closed with Esc or a
  click outside. Unlike dialogs they do not block the canvas.

## 3. Shell model

```ts
interface ShellLayout {
  docks: Record<DockId, { tabs: string[]; active?: string; size: number; state: 'open' | 'collapsed' | 'hidden' }>;
  canvas: { toolbar: boolean; statusLine: boolean; hoverLabel: 'cursor' | 'corner' | 'off' };
  mode: 'auto' | 'desktop' | 'compact' | 'embedded';
  expanded: boolean;          // fills the window (fixed position on the shell root, not <body> mutation)
  fullscreen: boolean;        // Fullscreen API on the shell root
}

class ShellModel extends ReactiveModel {
  readonly layout: ReadonlyState<ShellLayout>;
  readonly tabs: ReadonlyState<readonly TabEntry[]>;       // contributed tabs + where they are
  readonly overlay: ReadonlyState<OverlayEntry | undefined>;
  dock(id: DockId): DockModel;                              // per-dock states for granular rendering
  show(tabId: string, opts?: { reveal?: boolean }): void;
  hide(tabId: string): void;
  move(tabId: string, dock: DockId, index?: number): void;
  toggle(dock: DockId): void;
  openOverlay(id: string, anchor?: Element): void; closeOverlay(): void;
  openDialog(dialog: DialogSpec): DialogHandle;
  setMode(mode): void; setExpanded(v): void; setFullscreen(v): void;
  reset(): void;                                            // back to the spec's initial layout
}
```

- **Default layout**: right dock open with Home active while nothing is loaded and Structure once something is; top dock
  with the sequence (collapsed until a structure loads); left and bottom docks empty.
- **Modes** replace `controlsDisplay`. `desktop`: docks beside the viewport. `compact` (narrow or portrait, chosen by
  `auto` from the shell's own size via `ResizeObserver`, not window media queries): the right panel becomes a bottom
  sheet, overlays become sheets, the toolbar moves to the bottom edge. `embedded`: canvas, toolbar, palette,
  notifications only.
- **Persistence.** The layout is serialised to an injected store keyed by app id (localStorage in the viewer). It is not
  part of plugin state snapshots.
- **Ownership.** The shell model is a field of the new `PluginUIContext` (`ctx.shell`, §6), so the layout belongs to the
  UI context rather than to the core `PluginContext`.
- **Compatibility.** While old code still dispatches `PluginCommands.Layout.Update` or writes
  `behaviors.layout.leftPanelTabName`, a bridge in the new context maps them onto the shell (`showControls` → docks,
  `regionState` → dock states, `isExpanded` → `expanded`, tab names → tabs). New code uses `ctx.shell` and shell
  commands. The core `PluginLayout` and `leftPanelTabName` are deprecated once the old UI is retired; that is an
  API-ledger entry.

## 4. Default UI

### 4.1 Tabs (right panel)

| Tab                  | Content                                                                                                                                                                                                                                                                                                                                                                                                                                                         | Replaces                                                                  |
| -------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------- |
| **Home**             | Open files, download by ID (PDB, AlphaFold DB, PDB-IHM, URL as presets), examples, recently opened, a palette hint. Shown by default while nothing is loaded. Contributed import entries (Zenodo) appear here as commands.                                                                                                                                                                                                                                      | Home tab, custom import controls                                          |
| **Structure** (main) | Top: the scene tree (structures → components → representations; volumes; particles; measurements; generic entries), inline visibility, context menus. Header actions: apply preset, color by, add component (commands). Bottom: the inspector for the picked node (its parameters as a `ParamsFormModel`) followed by contributed sections for that node kind: assembly symmetry, volume streaming, kinemage groups on a structure; volume options on a volume. | Source, Components, Volume, Particles, Measurements list, extension cards |
| **State tree**       | The full state tree with the current node's editors and applicable actions.                                                                                                                                                                                                                                                                                                                                                                                     | Data tab                                                                  |
| **Settings**         | §4.5.                                                                                                                                                                                                                                                                                                                                                                                                                                                           | Settings tab, viewport settings panel                                     |
| Extension tabs       | Segmentation (volseg), PAE plot (model-archive; also movable to the bottom dock), app-specific tabs.                                                                                                                                                                                                                                                                                                                                                            | Large extension cards                                                     |

Quick styles, superposition, procedural animation, component options, measurement creation, and screenshot options move
to commands (with forms or custom UI where the arguments are rich).

### 4.2 Overlays

| Overlay        | Opened from                  | Content                                                                                                                                                                                                                                                                                      | Replaces                    |
| -------------- | ---------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------- |
| **Save state** | Toolbar 💾, `state.save`     | Name (auto-filled), three actions: **Snapshot**, **Download session** (`.molx`), **Download state** (`.molj`); the snapshot list below (apply, rename, reorder, delete, edit description). "More" reveals save options, key, markdown description, and remote sharing (upload, server list). | Plugin State tab            |
| **Screenshot** | Toolbar 📷                   | Preview with crop, size/format options, copy, download.                                                                                                                                                                                                                                      | Screenshot panel            |
| **Log**        | Status line message          | History with type filter and copy.                                                                                                                                                                                                                                                           | Bottom log region           |
| **Help**       | `help.*`, `?` palette prefix | Keybindings (generated), palette prefixes, mouse controls, selection help.                                                                                                                                                                                                                   | Help tab, viewport help     |
| **Animation**  | Playback cluster             | Choose an animation, its parameters, play/stop.                                                                                                                                                                                                                                              | Animation viewport controls |

Dialogs (modal) are kept for flows that must not be lost by a stray click: command forms opened "as form", MP4 export,
Zenodo import.

### 4.3 Snapshots and timeline

Each part appears only when it applies.

- **Snapshot panel** (top left, when snapshots exist): `‹ 3 / 12 ›`, the snapshot title, play/pause for the snapshot
  slideshow, and a progress bar during animated transitions. Keyboard: the existing Ctrl+←/→ bindings and commands
  `snapshot.next`/`.previous`/`.play`.
- **Snapshot description** (when the current snapshot has one): a markdown card in the snapshot panel with the existing
  markdown extensions (`#key` links apply snapshots, command links, color swatches, images, audio). Collapsible, text
  size adjustable, scrolls when long, remembers its collapsed state. Apps can replace it with their own renderer (a
  contribution), as mvs-stories and mesoscale-explorer do today.
- **Timeline** (bottom right, when the current structure comes from a multi-model trajectory or an animation runs):
  `4 / 20`, first, previous, scrubber, next, last, play/pause (model-index animation). Kept apart from the snapshot
  panel so neither moves the other.
- **Animation**: chosen and configured from the timeline (an overlay), played there.

### 4.4 Hover label

The label follows the cursor instead of sitting in a corner:

- Rows from `behaviors.labels.highlight` (HTML or markdown, as today), positioned from the hover event's page
  coordinates with an offset, flipped at the canvas edges, `pointer-events: none`.
- Positioned in `requestAnimationFrame` outside the reactive render path; hidden while dragging or rotating and when the
  pointer leaves the canvas.
- Setting `canvas.hoverLabel`: `cursor` (default), `corner` (today's behavior), or `off`. On touch devices it shows at
  the tap position for a moment.

### 4.5 Settings tab

Today's settings panel lists every canvas parameter at one level with truncated labels ("Camera Rese…", "Transparent
…"). The tab keeps the same parameters but organises them for use:

- **Search** at the top filters across all settings by label and description.
- **Groups by intent**, essentials first: Appearance (background, lighting/illumination, outline, occlusion, fog),
  Camera (mode, clipping, reset duration), Interaction (hover label placement, pick padding, trackball), Performance
  (resolution, multisample, transparency mode, Hi-Z), Layout and theme (UI theme, docks, status line), Advanced
  (renderer, handle, pointer, XR, behaviors tree). Groups come from a mapping over the existing PD params, so no param
  is lost and new ones land in Advanced.
- **Full labels**: labels sit above the control when they do not fit, with the PD description as a hint.
- **Reset** per group and for all.

### 4.6 Selection bar

Shown in selection mode at the bottom of the canvas: granularity, mode (add/remove/intersect/replace), the `select`
command group as an inline palette, color, create component, undo, exit ([commands.md](commands.md) §10.3).

### 4.7 Icons

Keep the vendored Material icons as the fallback set and move actions to one consistent modern line set (Lucide,
ISC-licensed; only the SVG path data that is used is copied into `components/icon/`, no runtime dependency), plus custom
glyphs for molecular concepts (cartoon, surface, ball-and-stick, ligand, volume, trajectory) that no generic set covers.
Every icon is chosen once per command, so buttons, menus, and palette rows agree.

Implementation: `scripts/icons.mjs` extracts the inner markup of the used icons from `lucide-static` 1.48.0 (a
devDependency; `pnpm icons`, `icons:check` in CI) into the generated `components/icon/_generated/lucide.ts`; the custom
molecular glyphs are hand-written in `components/icon/molecular.ts` on the same 24×24 stroke grid. `Icon` parses each
icon's markup once into a cached `<svg>` template and returns a clone with size, class, and label (`aria-hidden`
otherwise) set, so rendering an icon allocates one node tree and no parsing.

## 5. Contributions

A contribution is a plain, declarative object, listed explicitly in the spec (like `registry` entries):

```ts
interface UIContribution {
  id: string;
  commands?: readonly UICommand[];
  commandGroups?: readonly CommandGroup[];        // new groups, or nested under built-in ones (commands.md §9)
  tabs?: readonly TabContribution[];               // panel tabs (right by default, any dock)
  sections?: readonly SectionContribution[];       // Structure-tab inspector sections for a node kind
  overlays?: readonly OverlayContribution[];       // canvas overlay cards
  canvasItems?: readonly CanvasItem[];             // toolbar buttons, playback-cluster parts, other canvas overlays
  statusItems?: readonly StatusItem[];             // status-line chips
  menus?: readonly MenuItem[];                     // context menu entries (command refs)
  sceneEntries?: readonly SceneEntryContribution[];// scene-tree entries; replaces genericRepresentationControls
  paramEditors?: readonly ParamEditorContribution[]; // replaces customParamEditors (by transformer id or PD field)
  settings?: readonly SettingsGroupContribution[]; // groups or params added to the Settings tab
  formats?: readonly [name: string, ArgumentFormat][];
  models?(ui: PluginUIModel): ContributionModels;  // view models it owns, created on activation
}

interface TabContribution {
  id: string; title: string; icon: IconName;
  dock?: DockId;                                   // default 'right'
  order?: number;
  when?(ctx: CommandContext): boolean;            // e.g. only when a volseg entry exists
  render: ViewRenderer<ViewContext>;               // dom, React, Preact, or element (rendering.md §6)
}

interface SectionContribution {
  id: string; title: string; icon?: IconName;
  nodeKind: SceneNodeKind | readonly SceneNodeKind[]; // 'structure', 'component', 'volume', …
  when?(node: SceneNode, ctx: CommandContext): boolean; // e.g. assembly symmetry only when applicable
  render: ViewRenderer<SectionContext>;
}

interface OverlayContribution {
  id: string; title: string;
  render: ViewRenderer<OverlayContext>;            // opened by a command or a canvas item
}
```

- **Placement of extension UI.** Most extension panels become sections in the Structure tab (assembly symmetry,
  kinemage, volume streaming), commands (exports), or overlays/dialogs (MP4 export, Zenodo import). Only large,
  standalone views become their own tabs (segmentation, PAE plot).
- Contributions are activated after `plugin.init()`; each activation's registrations are undone on deactivation, as
  `plugin.register` does for registry entries.
- `when` replaces `isHidden` toggling inside cards; a tab, section, or item appears when its condition holds.
- An extension package exports its registry entry, its behavior, and (separately) its UI contribution, e.g.
  `@molstar/geo-export-extension/ui-next`. Apps list what they want; nothing registers UI as a side effect of a
  behavior.

## 6. The new `PluginUIContext` and UI spec

ui-next defines its own `PluginUIContext`, separate from `@molstar/plugin-ui`'s, carrying everything the new UI needs as
typed fields instead of untyped maps on core `PluginContext`:

```ts
class PluginUIContext extends PluginContext {
  readonly spec: PluginUISpec;
  readonly shell: ShellModel;                         // replaces plugin.layout + leftPanelTabName for the new UI
  readonly commands: CommandRegistry<CommandContext>; // user-facing commands (PluginCommands stay the internal bus)
  readonly keybindings: KeybindingService;
  readonly contributions: ContributionRegistry;       // active UIContributions; replaces customStructureControls,
                                                      // customImportControls, genericRepresentationControls,
                                                      // customParamEditors
  readonly sources: PluginSources;                    // created after init()
  readonly model: PluginUIModel;                      // root view model (scene, selection, viewport, tasks, …)
  readonly theme: ReadonlyState<Theme>;
}
```

The context owns the lifetimes: `init()` creates sources and the root model and activates the spec's contributions;
`dispose()` deactivates them and retires the model. Views and commands receive the context (or the parts they need),
never reach into `PluginContext.layout`.

The spec below is a starting point; it is expected to change as the views are built.

```ts
interface PluginUISpec extends PluginSpec {
  ui?: {
    contributions?: readonly UIContribution[];     // explicit, ordered; DefaultUIContributions() for the full UI
    layout?: Partial<ShellLayout> & { tabs?: Record<string, { dock?: DockId; hidden?: boolean }> };
    theme?: 'light' | 'dark' | 'auto' | ThemeTokens;
    commands?: { hide?: readonly string[] };        // command or group ids, globs allowed ('superpose.*')
    keybindings?: readonly { command: string; key: KeyChord | null }[];
    canvas?: { toolbar?: readonly string[] };       // command ids, order
    persistLayout?: false | { key: string };
  };
}
```

- `DefaultUIContributions()` in `@molstar/plugin-ui-next/default-spec` lists the built-in contributions; apps spread it
  and add extension contributions, mirroring `DefaultRegistry` and the v6 composition rule.
- Entry point: `createPluginUI({ target, spec })` builds the new `PluginUIContext`, initialises it, and mounts
  `<msp-plugin>`; no React render function is needed. `ui-next/react` offers `<PluginUI>` for React hosts.
- Today's spec fields map as follows: `components.controls.*: 'none'` → hide the tab or dock in `layout`;
  `components.viewport.view` → canvas items and overlays; `components.structureTools` → contributions;
  `remoteState: 'none'` → omit the remote-snapshots contribution; `customParamEditors` → `paramEditors`;
  `selectionTools.hide` → omit or reorder selection-bar items; `PluginConfig.Viewport.*` → toolbar commands' `when`
  (kept, so the viewer options keep working).

## 7. Embeddable elements

Top-level custom elements over the same models, for hosts that do not want the shell:

| Element                  | Content                                  | Used by                               |
| ------------------------ | ---------------------------------------- | ------------------------------------- |
| `<msp-plugin>`           | Full shell                               | Viewer, apps                          |
| `<msp-viewport>`         | Canvas + overlays + toolbar (+ palette)  | mvs-stories, docking viewer, examples |
| `<msp-sequence>`         | Sequence view                            | Hosts with their own layout           |
| `<msp-params-form>`      | PD form bound to a `ParamsFormModel`     | alpha-orbitals example, host settings |
| `<msp-snapshot-stepper>` | Snapshot i/N, previous/next, description | mvs-stories                           |
| `<msp-pae-plot>`         | PAE plot (model-archive contribution)    | alphafolddb-pae example               |

Each takes a `ui` or `model` property; several elements on a page can share one `PluginUIModel`.
