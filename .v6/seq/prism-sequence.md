# Design: the sequence viewer in Prism (first slice)

Status: proposed 2026-10-14, for David's review; the first slice runs in the mocked viewer mode (David). Follows [plan.md §3](plan.md); facts from the code at the port
(`packages/sequences`, `packages/prism/{ui,viewer}`).

## 1. Scope

In the first slice:

- **One Seq\* viewer** where the mock's sequence strip is now (the floating top dock), **without the card around it**
  (no glass panel, no hairline): the viewer sits on the canvas.
- **Two presentations**, switched by a small control: **compact text** (letters in groups of 10, wrapping) and
  **tracks** (one linear track, zoomable).
- **Hover and selection, both ways** with the 3D view; highlight is the hover's 3D side.
- **Tooltip** (the hovered position) and a **legend** (what the colors mean), in Prism's own components.
- **The canvas's ink**: it follows the canvas background as the other canvas controls do.

Later, each when needed: secondary structure, annotations and other tracks (providers already exist in Seq\*:
`secondary_structure`, `plddt`, `pae`, `sifts_coverage`), track menus, search, export, the hover card's actions,
several viewers at once, the story mode.

## 2. Owned by the Prism viewer, no extension mechanism

prism-ui is the runtime for generic Mol\* UIs: the shell, the common controls, the ink abstraction, the agentic
wrappers, commands. The Prism viewer (`@molstar/prism-viewer`) is an opinionated composition on it, and so will the
story mode be (David, 2026-10-14). So the sequence viewer is not an extension composed through modules and
contributions: the viewer owns it directly, in its own code, and is free to be specific about it.

- **Where**: a feature folder in the viewer's package, `packages/prism/viewer/src/sequence/`. No new package (David:
  fewer packages, not more) and no new prism-ui abstraction. The story mode, when it wants a sequence view, composes
  Seq\* its own way.
- **Dependencies**: prism-viewer depends on `@molstar/seq-viewer`, `@molstar/seq-core` and `@molstar/seq-molstar-poc`
  (the structure provider, `StructureLink`, `linkViews`, `StructureView3D`); prism-ui stays free of Seq\*. This
  settles plan.md §3's "prism-ui depends on seq-viewer" for the viewer's package. molstar-poc is private and named
  for a proof of concept; its link and provider as a real package is a follow-up (§8).
- **Extensions stay possible**: the Prism viewer will take extensions through its own composition spec, as the v5
  Viewer does; prism-ui, the library, is not opinionated about how. The sequence view is not one of them: it is the
  viewer's own feature.
- **contexts.md**: its §3–4 (modules with `ModuleKey`, `UIContribution`s, a mode's `contributions`) describe the
  extension mechanism this design does without. contexts.md is to be revised for the "runtime plus opinionated
  compositions" split (§8); until then, this design follows David's direction, not those sections.

### 2.1 Who owns the `SeqViewerContext`: the viewer's model

The viewer's root model (the object the viewer mode creates at start and disposes at its end; the mock's
`viewer-mode.tsx` today) creates the sequence model and owns it; the sequence model owns every Seq\* context.

```ts
// packages/prism/viewer/src/sequence/model.ts
export interface SequenceViewers {
  /** The viewers by id; the first slice creates only 'default'. */
  get(id?: string): SequenceViewer;
  readonly ids: ReadonlyState<readonly string[]>;
  dispose(): void; // disposes every viewer's context
}

export interface SequenceViewer {
  readonly id: string;
  readonly seq: SeqViewerContext; // createSeqViewer(spec, ports)
  readonly style: ReadonlyState<'compact' | 'tracks'>;
  setStyle(style: 'compact' | 'tracks'): void; // ctx.ui.setStyle(section, ...), UI state only
  readonly structure: ReadonlyState<SequenceStructure | undefined>; // what it shows (§3)
  readonly theme: ReadonlyState<SequenceTheme>; // the color scheme the letters use, and later each track's (§6)
  setTheme(theme: SequenceTheme): void; // written into the document the model builds
}

export function createSequenceViewers(deps: { ink: ReadonlyState<CanvasInk>; /* the structure source (§3) */ }):
  SequenceViewers; // follows deps.ink (§5) and keeps each viewer's skin in step
```

- Keyed by id from the start, so a second viewer later is another entry and another place in the layout, not a
  redesign.
- The context is created with `rafScheduler` (Seq\*'s default scheduler never fires), the default registry plus
  molstar-poc's `MolstarSvsEntry`, and loaders that refuse the network: the structure comes from the host (§3).
- It takes what it follows as plain arguments (the ink, the structure source), not a context object, so it is
  testable headless.

### 2.2 The dock tab

The viewer registers the tab itself, as the mock registers the strip today: `{ id: 'sequence', title: 'Sequence',
icon: 'sequence', dock: 'top', render: domView(SequenceControl) }`, registered while a structure is shown, replacing
`SequenceStrip`. The "Sequence" canvas button and the `view.sequence` command keep working, since they address the
tab.

### 2.3 The link: subscribed while the control is rendered

`SequenceControl` subscribes the two-way link when it mounts and drops it on cleanup:

```ts
const unlink = linkViews(viewer.seq, view3d, link, options); // hover both ways, selection both ways, focus
onCleanup(() => { unlink(); link.detach(); });
```

- `view3d` is a `StructureView3D` (molstar-poc's Mol\*-agnostic interface: `highlight`, `select`, `focus`, `onHover`,
  `onSelect`). Over the real plugin it is molstar-poc's `pluginView3D(plugin, structure)`, which lives in the
  viewer's `sources/` folder (contexts.md §7: only `sources/` and `actions/` import `@molstar/plugin`). In the
  stories it is that one, over the stories' real plugin (§3.1); a headless test uses a fake.
- The replay rule makes remounting safe: on subscribe, the current selection is sent as a selection and never as a
  clear, and a reassembly's carry (`origin: 'graph'`) is no one's move (interaction.md, replays and carries).
- linkViews' size gates stay at their defaults (hover 2,000 items, selection 20,000; over that, whole chains).

## 3. What the viewer shows

The structure the 3D view shows: one view, one section over a `molstar:structure` complex. With no blocks, a section
takes every polymer unit (chain instance) in order, then ligands.

```ts
const b = createSeqViewSpec({ title });
b.data.complex({ ref: 'structure', provider: 'molstar:structure', source_ref: b.data.source({ ref: 'file', ... }) });
b.view({ ref: 'sequence' }).section({ ref: 'main', style: 'compact' /* or 'tracks' */ }).track().letters();
```

- The 3D view and Seq\* share one parse: the loaders' `parsed` hook returns the host's trajectory as a
  `StructureSource` (`StructureSource.of(trajectory, format)`), so the `Structure` that `StructureLink.attach` maps is
  the one the 3D view draws.
- In the stories, the real plugin of §3.1 loads the structure and the sequence model takes its trajectory; headless
  tests use molstar-poc's recorded fixtures (p53 structures), so no new fixture is committed.
- The structure picker in the strip's header (structure, chain/polymer mode, instance) is not part of the first slice.

### 3.1 A minimal real plugin in the stories

The viewer-mode stories get a real `PluginContext` in the shell's canvas area (David, 2026-10-14), not yet hooked to
the rest of the mock UI:

- Composed as molstar-poc's `pluginSpec()` is (v6's explicit composition: highlight, select and camera-focus
  behaviors, MolViewSpec), mounted on the canvas host the shell already has (`CanvasPlaceholder`'s place), disposed
  with the story.
- It loads the mock's structure, 1TQN, from RCSB at run time (the stories may use the network; tests never do) and
  shows it with the default representation.
- **The background**: a control in the story changes the canvas background (`canvas3d.setProps({ renderer: {
  backgroundColor } })`); the canvas props drive `canvasInk(props)`, which the shell's canvas layer and the sequence
  model both follow, so the area behind the sequence viewer and its skin change with the canvas.
- The plugin code lives in the viewer's `sources/` folder (contexts.md §7), the mock UI stays as it is (structure rows,
  selection bar and the rest are still mocks).

## 4. The control

```
┌ top dock, no card ──────────────────────────────────────────────────────────────┐
│ MKLAVL... (compact: letters in groups of 10, wrapping)          [Aa | ≡]  [◐]   │
└──────────────────────────────────────────────────────────────────────────────────┘
```

- **The surface**: `mountSection(seq, 'main', host)` from `@molstar/seq-viewer/canvas`, one section, so no loop over
  sections and none of Seq\*'s DOM or React UI. The engine draws the letters, ruler, hover, selection and focus itself
  and handles wheel, drag, click and keys.
- **The switch**: prism-ui's `Choices`, two icon choices (compact text, tracks), calling `viewer.setStyle`. Compact
  is Seq\*'s `compact` preset (wrapped, groups of 10, no zoom, vertical scroll); tracks is the `tracks` preset (one
  track, the ruler rail, zoom and pan).
- **Height**: tracks: the dock follows `seq.layout.contentHeight('main')` (navigator, ruler and one row). Compact:
  the dock's height (`--msp-dock-size`), the text scrolling inside it.
- **No frame**: the tab asks for no card (`ShellTab.frame: 'none'`; no `--msp-glass-panel-bg`, shadow or filter);
  the canvas shows through. The surface's canvases clear rather than fill, so the viewer draws on the canvas
  directly. Seq\* draws its text halos in the skin's surface color (§5), which keeps letters legible on any background.

## 5. The ink: prism-ui's abstraction, the viewer's state

The canvas ink (the canvas background, and whether chrome over it is light or dark) is prism-ui's abstraction; the
viewer holds its state, and the sequence model subscribes to it and updates its own. Nothing in the sequence code
reads the DOM or CSS to find it.

- **The abstraction**: prism-ui exports `CanvasInk = { theme: 'light' | 'dark'; color: string }` and
  `canvasInk(props: ReadonlyState<CanvasBackgroundProps>)`, a disposable `ReadonlyState<CanvasInk>` (its owner
  disposes it; `fixedInk(color)` for the mocks), which is the `canvasTheme(props$)` it already computes
  (`ui/src/model/canvas-theme.ts`, WCAG luminance at 0.179) plus the base color, today handed to the shell only as
  `canvasTheme` / `canvasColor` props; the shell takes the ink state instead. The viewer creates the one ink state and
  gives it to the shell and to `createSequenceViewers`; in the mocks
  it comes from the mock's background, over the plugin from the canvas props source (view-models.md §3, `canvas`).
- **The subscription is the sequence model's**, not the control's: each viewer follows the ink while it exists, so its state
  is already right when the control mounts. On a change it builds Seq\*'s `Skin` from the ink and sets it
  (`seq.setSkin(skinFor(ink))`): the surface and dark flag from the ink, ink and accent colors from Prism's palette
  for that theme, the selection as the mock strip draws it (the accent at 25%). Seq\*'s data themes have dark
  variants; a skin change reassembles the view with them.
- **Not used**: Seq\*'s skin reader (it reads `--seq-*` from computed style and watches `data-msp-theme` attributes).
  An explicit skin keeps one source of truth and works the same in a shadow root, a story and a headless test. The
  control still sets `--seq-font-*` from `--msp-font-*`; fonts are not part of the ink.

## 6. Theme, tooltip and legend

The sequence model manages its themes; Seq\* colors. Which color scheme the letters use (and, later, each track) is the
model's state (`SequenceViewer.theme`), next to the ink it already follows (§5). The model writes the choice into
the document it builds; Seq\*'s themes and renderers draw the colors from it, so nothing about coloring is
reimplemented in Prism.

- **Legend**: a small prism-ui `Legend` component (items: swatch or glyph and label; ramps: gradient, low and high),
  opened from a button on the control. Its data is the model's own: the scheme it chose, read from the scheme's
  public definition (`KernelColorSchemes`, `ColorSchemeEntry` from `@molstar/seq-core/format`). Seq\*'s legend model
  (`legendGroups` in `viewer/src/chrome.ts`) stays internal to its DOM and React UIs: no new Seq\* entry.
- **Tooltip**: the sequence's hover feeds the shell's one cursor-following `HoverLabel` (shell.md: the sequence's
  hover feeds the same hover label as the 3D canvas). The control maps `seq.labels.current` (`{ title, lines }`,
  already public) to the label model's lines while the pointer is over the surface. No per-move rendering: lines
  change only when the label does.
- Track headers stay off (`headerWidth: 0`), so header tooltips, track menus and the hover card are not in play.
- **Later, data tracks**: a legend whose bands or range come from the data or the user's edits (pLDDT, a continuous
  annotation) is computed by Seq\*'s renderers today. With fixed bands (pLDDT's) the model knows them from the
  theme it set; for the rest, Seq\*'s per-track legend (`trackLegends`) gets an export then (§8).

## 7. Work, in order

| #   | What                                                                                                                                                                                         | Who            |
| --- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------- |
| S2  | prism-ui: the `CanvasInk` abstraction (`canvasInk(props$)`, the shell reading it), `FloatingDock`'s frameless variant, the `Legend` component with a story                                                                 | Sonnet, medium |
| S3a | the stories' real plugin (§3.1): composition, mount on the canvas host, 1TQN from RCSB, the background control and `canvasInk` over its canvas props | Sonnet, high   |
| S3  | `prism-viewer/src/sequence`: the model (the ink subscription and `skinFor`, the theme state and the legend data from it), document, control, link over `StructureView3D`; `pluginView3D` over the stories' plugin, a fake `StructureView3D` for headless tests                                         | Sonnet, high   |
| S4  | the viewer-mode mock creates the sequence model and its tab instead of `SequenceStrip`; the shell stories (including `CanvasBackground` for the ink); the subscription guard extended to prism-viewer | Sonnet, medium |
| S5  | Visual check in Storybook                                                                                                                                                                    | David          |

S3a and S3 need S2; S4 needs both. Over the real plugin, the same code gets its
`pluginView3D` in `sources/` once Prism's plugin sources exist (view-models.md §3); nothing in the control changes.

## 8. Follow-ups

| Item                                                                                                     | Trigger                              |
| -------------------------------------------------------------------------------------------------------- | ------------------------------------ |
| molstar-poc's link and structure provider as a real package (`@molstar/seq-molstar`), out of the PoC      | the sequence view leaves the mocks   |
| Several viewers (a second id, its own tab or dock)                                                       | a feature that needs two             |
| The strip header's structure picker (structure, chain or polymer mode, instance)                         | more than one structure loaded       |
| Tracks beyond letters (secondary structure, pLDDT, annotations) and their legend; Seq\*'s `trackLegends` exported for data-driven ranges | asked for in Prism                   |
| contexts.md revised for prism-ui as a runtime and the viewer and story mode as opinionated compositions (§2): its modules and contributions (§3–4) | before the next design that leans on contexts.md |
| The `--msp-ink` tokens (defined, unused) for chrome outside the canvas layer                             | a canvas control outside that layer  |
| todo.md's "the Mol\* hover sometimes needs a second move after the pointer was over the 3D canvas"       | the link over the real plugin        |
