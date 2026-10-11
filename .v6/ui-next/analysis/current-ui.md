# Current UI: functional inventory

What `@molstar/plugin-ui` (`packages/plugin/ui/src`) does today, what it reads from the plugin, and how other code
extends it. This is the baseline ui-next must cover "in spirit". Paths: `P` = `packages/plugin/ui/src`, `C` =
`packages/plugin/core/src`.

## 1. Size and shape

57 TS/TSX files, 17,883 LOC; 3,662 LOC SCSS.

| Folder                         | LOC   | Content                                                            |
| ------------------------------ | ----- | ------------------------------------------------------------------ |
| `controls/` (incl. line-graph) | 5,752 | PD parameter editors, slider, color, icons, action menu, markdown  |
| `structure/`                   | 4,795 | Right-panel cards: source, components, measurements, superposition |
| root `*.tsx`                   | 2,476 | Layout, left panel, viewport overlays, sequence wrapper, toasts    |
| `state/`                       | 1,920 | State tree, transform editors, snapshots, animation                |
| `sequence/`                    | 1,018 | Sequence widget and wrappers                                       |
| `viewport/`                    | 821   | Canvas, screenshot, settings, help                                 |
| `custom/`                      | 501   | Volume-streaming custom param editor                               |

Hotspots: `controls/parameters.tsx` (2,213; ~30 classes), `structure/superposition.tsx` (785; algorithms + state surgery
in the view), `structure/selection.tsx` (664; 10 subscriptions, 21 setState/forceUpdate), `state/snapshots.tsx` (662),
`structure/components.tsx`/`volume.tsx`/`particles.tsx` (copy-paste hierarchy cards with data loading inside the view),
`state/tree.tsx` (subscriptions scale with tree size).

## 2. Infrastructure

- `PluginUIComponent` / `PurePluginUIComponent` (`P/base.tsx`, 71 LOC) are React class components that read the plugin
  from `PluginReactContext` and call `this.subscribe(obs, () => this.forceUpdate())`. There is no store or view model:
  ~70 classes read ~40 manager/behavior streams ad hoc and re-read `this.plugin.*` in `render()`.
- `useBehavior` (`P/hooks/use-behavior.ts`) exists but only 5 files use hooks.
- `CollapsableControls` is the card template for every right-panel tool and every extension panel (header, brand accent,
  `isHidden`, `renderControls()`); 19 external imports.
- `TransformControlBase` (`P/state/common.tsx`) is the apply/update shell for transforms and actions (auto-apply timer,
  busy tracking, reset to default). Subclasses: `UpdateTransformControl`, `ApplyActionControl`.

## 3. Spec, context, bootstrap

```ts
interface PluginUISpec extends PluginSpec {
  customParamEditors?: [StateAction | StateTransformer, StateTransformParameters.Class][];
  components?: {
    controls?: { top?; left?; right?; bottom? }; // component or 'none'
    remoteState?: 'none' | 'default';
    structureTools?; // right panel body
    viewport?: { view?; controls?; snapshotDescription? };
    sequenceViewer?: { view?; modeOptions?; defaultMode? };
    hideTaskOverlay?: boolean;
    disableDragOverlay?: boolean;
    selectionTools?: { controls?; granularityOptions?; hide?: {...} };
  };
}
```

`PluginUIContext extends PluginContext` adds `customParamEditors` and an untyped `customUIState`.
`createPluginUI({ target, render, spec })` creates the context, awaits `init()`, renders `<Plugin>` with an injected
React 18 `render`. Unused by the repo outside defaults: `sequenceViewer`, `selectionTools`, `hideTaskOverlay`,
`disableDragOverlay`; `PluginSpec.Action.customControl` is dead.

Extension UI enters through three untyped maps on core `PluginContext` (React components held by a non-UI package):
`customStructureControls` (right-panel cards), `customImportControls` (Home tab), `genericRepresentationControls`
(entries in the Components list). The plugin-composition design already lists moving these to a typed UI registry as a
follow-up.

## 4. Layout

`plugin.layout` (`C/layout.ts`, `PluginLayout`):
`{ isExpanded, showControls, regionState: { left: 'full'|'collapsed'|'hidden', top|right|bottom: 'full'|'hidden' }, controlsDisplay: 'outside'|'portrait'|'landscape'|'reactive', expandToFullscreen }`,
changed through `PluginCommands.Layout.Update`. `plugin.behaviors.layout.leftPanelTabName` is a separate subject that
extensions write directly (`sb-ncbr`, `volseg`, `g3d`). Expanded mode mutates `<body>` style and `<meta viewport>`.

Rendering (`P/plugin.tsx` `Layout`): five fixed regions (main, top = sequence, left = tabs, right = structure tools,
bottom = log). Geometry is fixed in Sass per `controlsDisplay` (left 330, right 300, top 100, bottom 70); no resizing,
no docking, no user-movable panels.

## 5. Functional inventory

**Left panel** (icon strip + body): Home (registry `StateAction`s as apply cards: download, open files, …; custom import
controls; remote snapshots), State tree (with dirty marker), Plugin state (snapshots), Help, Settings (Canvas3D params +
behavior tree).

**Right panel** (`DefaultStructureTools`, fixed order): Source (hierarchy picker, preset, model index, structure type,
conformation, focus row, selection stats), Measurements, Superposition (chains/TM-align/ligands/atoms/SIFTS), Quick
styles, Procedural animation, Components (preset/add/options/undo; groups with visibility, delete, actions menu,
per-representation editors; generic entries), Volume streaming, Volume, Particles, extension cards.

**Top**: sequence viewer (structure/mode/entity/chain/operator selectors; DOM residues; registers highlight/selection
providers). **Bottom**: log (last 10 entries).

**Viewport overlays**: camera reset (+ zoom, orient/reset axes), screenshot panel (preview, crop, copy, download, state
export/import, camera info), toggle controls, expand/fullscreen, settings panel (canvas params + region toggles + help),
illumination, XR, selection-mode toggle; top-left animation, trajectory, and snapshot stepper controls; snapshot
description (markdown); selection toolbar (granularity; union/subtract/intersect/set by query or residue list; apply
theme; create/remove component; undo; help); hover labels; toasts; task progress (abortable); drag-and-drop overlay.

**State tree**: recursive nodes with expand, current, delete, visibility, hover highlight; the current node shows
decorator editors and applicable actions.

**Snapshots**: options, add (name/key/markdown), local list (apply/edit/move/replace/remove), save/open `.molj`/`.molx`,
remote list/upload/fetch.

**Help**: static prose and binding tables generated from trackball and behavior bindings. **Keyboard**: none beyond
input-level Enter/Escape.

## 6. Parameter editors

`ParameterControls` maps PD types to controls via `controlFor(type)`: boolean, number (slider when bounded), text,
select (menu or cycle), multi-select, color, color-list (with offsets/presets), vec3, mat4, url, file, file-list,
value-ref, interval, group (presets, pivot, flat), mapped, converted, conditioned, line-graph, script, object-list.
Essentials vs advanced (`isEssential`), `category` groups, `hideIf(values)`. Per-transformer overrides through
`customParamEditors`; the only built-in one is the volume-streaming editor.

## 7. What reaches the plugin

The UI mixes `PluginCommands.*` (State, Interactivity, Layout, Toast, Camera, Canvas3D) with many direct manager calls
(`structure.selection.fromSelectionQuery`, `component.*`, `hierarchy.*`, `measurement.*`, `viewportScreenshot.*`) and
direct property/behavior writes (`plugin.selectionMode = …`, `leftPanelTabName.next(…)`, `snapshot.togglePlay()`). There
is no central, enumerable action registry, which is the main obstacle to a palette. The per-panel read/invoke map is in
[plugin-state.md](plugin-state.md).

## 8. Styling

Sass with build-time theme variants (`light`, `dark`, `blue`: separate CSS files from `@use … with (…)`), no CSS custom
properties, no runtime theme switching. Everything under `.msp-plugin`, ~229 `msp-*` classes, flat BEM-ish names, heavy
inline styles. Libraries ship compiled `lib/skin/*.css`; apps compile SCSS through `esbuild-sass-plugin`.

## 9. External usage

- Spec composition: viewer (layout flags, region state, `controls.*: 'none'`, ~20 `PluginConfig.Viewport.*`),
  docking-viewer (replaces `viewport.view`), mesoscale-explorer (replaces left and right panels, ~3k LOC of its own
  panels), mvs-stories (viewport only, inside custom elements), examples (hide all chrome; external React controls).
- Building blocks imported from plugin-ui by extensions/apps/examples: `context` 20, `controls/collapsable` 19,
  `default-spec` 14, `react18` 12, `spec` 11, `controls/parameters` 11, `controls/icons` 11, `base` 11,
  `hooks/use-behavior` 9, `controls/common` 9.
- `@molstar/plugin-extension` (`extensions/plugin`) has `PluginViewModel`/`PluginUIViewModel` (create/mount holders, not
  UI view models), `PluginCanvas`, and hooks; `examples/react` uses them.

See [extensions.md](extensions.md) for the per-extension inventory and port plan.
