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
icon: 'sequence', dock: 'top', frame: 'none', render: domView(SequenceControl) }`, registered while a structure is
shown, replacing `SequenceStrip`. The canvas stack's "Sequence" button (always shown while a structure is, with the menu of
§4) and the `view.sequence` command address the tab.

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
- **It shows what the user wants by default** (David, 2026-10-14): switching should be rare. The sequence follows what
  the user is working with, in order: the focused structure (and chain), the structure of the latest selection, the
  one in view, the first loaded. The same idea carries the views to come: a contact map or another track that follows
  the focus (a dynamic contact map around the focused residues, for one).
- **No picker in the strip** (one way of doing things across the app). When the default is not what the user wants,
  the structure tree chooses: the structure (and later the chain or instance) chosen there. The strip has at most one
  icon, the structure's name, which runs a command that opens the structure tree with that structure revealed; the
  tree's rows and the command menu ("Show sequence of…") are where another one is chosen. Commands rather than a
  bespoke picker also give the agent the same way in. The mock's mode select (chain, polymers, everything) goes too:
  tracks show every polymer chain as its own block, compact text every chain in turn with its name.

### 3.1 A minimal real plugin in the stories

The viewer-mode stories get a real `PluginContext` in the shell's canvas area (David, 2026-10-14), not yet hooked to
the rest of the mock UI:

- Composed as molstar-poc's `pluginSpec()` is (v6's explicit composition: highlight, select and camera-focus
  behaviors, MolViewSpec), mounted on the canvas host the shell already has (`CanvasPlaceholder`'s place), disposed
  with the story.
- It loads the mock's structure, 1TQN, from RCSB at run time (the stories may use the network; tests never do) and
  shows it with the default representation.
- **The background**: the mock's Canvas settings › Background (color, transparent) drives the plugin's canvas props and
  follows them; the canvas props drive `canvasInk(props)`, which the shell's canvas layer and the sequence model both
  follow, so the band behind the sequence viewer and its skin change with the canvas.
- The plugin code lives in the viewer's `sources/` folder (contexts.md §7), the mock UI stays as it is (structure rows,
  selection bar and the rest are still mocks).

## 4. The control

```
┌ band: the canvas' color behind the lines, fading out below them ───────────────┐
│ MKLAVL... (text: groups of 10, centered, wrapping; scrolls with edge fades)     │
│ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ │
└ (the structure may pass under the faded part; the camera frames it below) ──────┘
```

- **The surface**: `mountSection(seq, 'main', host)` from `@molstar/seq-viewer/canvas`, one section, no loop over
  sections, none of Seq\*'s DOM or React UI. The engine draws the letters, ruler, hover, selection and focus itself and
  handles wheel, drag, click and keys. **The control has no controls of its own** (round 2): the version, the colors and
  the legend are the canvas stack's **Sequence** button (§2.2), always shown while a structure is, with a menu: one single choice of **No sequence, Text, Tracks** (each with its own icon, the current one in the accent as the picking
  levels' list marks its own; No sequence hides the band), then the letters' colors (a single choice) and Legend…. Each item runs a command (`sequence.text`, `sequence.tracks`,
  `sequence.colors.<scheme>`, `sequence.legend`, and `view.sequence`), so the command menu has the same, one way of doing
  things. The tracks' rail (`sequence.rail.<ruler|strip|band|none>`) is a command and `SequenceViewer.setRail`, not a
  menu item: the ruler's rail is what Tracks uses.
- **Tracks is the default** (the model's initial style); Text is one menu click away.
- **Text** is Seq\*'s `compact` preset with the model's render props: whole blocks per line over the full width,
  centered (`wrappedAlign: 'center'`, the margins split on both sides; `wrappedMaxWidth` limits a line, default none),
  letters at most the panel text size (`--msp-font-size-sm`, read once from the tokens), in the layer's colors
  (`wrappedColors`), unobserved positions in a quiet gray (`unobservedGray`, the skin's `unobserved`); a focus scrolls its
  lines into view (`focusViewport`'s `reveal`); the scroll strip shows only while the pointer is over the band, the text
  scrolls or a drag is on, and the cut edges fade (the control drives a mask from the viewport's scroll state, written
  when a side turns on or off, never per step). **Tracks** is the model's own preset `prism_tracks`: no header column
  (`headerWidth` 0), the ruler and the letters, the rail a registered host preset (`setRail`).
- **The band** (a frameless top tab, `ShellTab.frame: 'none'`): part of the canvas area, over the canvas' top. It spans
  the canvas area (the left edge to the rail, or to a pinned panel's border; a floating panel lies over it as over the
  canvas), in the canvas' own color (the ink) behind its lines, **fading out** over `BandDissolvePx` before and
  `BandFadePx` after the last line: a gradient, no `backdrop-filter`. The 3D canvas keeps its full height under it and
  the same width source (`--msp-canvas-end`); the band takes the pointer over the lines (the surface), the faded part
  lets it through to the canvas (`pointer-events`). The structure may pass under the faded part; the camera's default
  framing (the first one, a reset) keeps the band clear (`paddedFocus`: a snapshot provider for `requestCameraReset`),
  and follows the band's height only while the camera is as the framing left it.
- **Height**: the dock size is the lines' height plus the air above (`--msp-space-3`) and `BandFadePx`: Tracks follow
  `seq.layout.contentHeight('main')`, Text has a fixed height (two or three lines; the rest scrolls). `ShellModel.setSize`
  has its own minimum for a top or bottom dock.

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

- **Letters' colors**: the model's own scheme `prism:residues` (Zappo's residue classes, colder and softer, with a dark
  variant), registered as a Seq\* color scheme from the viewer's side (a registry entry's `colorSchemes`; names are
  namespaced) and the default; Zappo, Taylor and Hydrophobicity are in the menu too. Seq\* judges a scheme's fit by the
  section's first block, so the scheme lists every polymer alphabet (a nucleic acid first would otherwise show the
  letters' own colors); its letters are colored by class where they are the same letters (A, C, G, T).
- **Legend**: a small prism-ui `Legend` component (items: swatch or glyph and label; ramps: gradient, low and high),
  in an overlay card opened from the Sequence menu. Its data is the model's own: the scheme it chose, read from the scheme's
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
| Focus-driven sequence views (a dynamic contact map on focus, tracks around the focused residues)          | the focus source over the real plugin |
| Choosing what the sequence shows through the structure tree and commands (§3): the strip's structure icon, the tree rows' and the command menu's "Show sequence", instances | more than one structure loaded       |
| Tracks beyond letters (secondary structure, pLDDT, annotations) and their legend; Seq\*'s `trackLegends` exported for data-driven ranges | asked for in Prism                   |
| contexts.md revised for prism-ui as a runtime and the viewer and story mode as opinionated compositions (§2): its modules and contributions (§3–4) | before the next design that leans on contexts.md |
| The `--msp-ink` tokens (defined, unused) for chrome outside the canvas layer                             | a canvas control outside that layer  |
| todo.md's "the Mol\* hover sometimes needs a second move after the pointer was over the 3D canvas"       | the link over the real plugin        |

## 9. As built (S3a, S3, S4)

Where the code differs from, or fills in, the sections above:

- **Files**: `packages/prism/viewer/src/sequence/` (`model.ts`, `document.ts`, `skin.ts`, `legend.ts`, `link.ts`,
  `sequence-control.tsx`), `src/sources/` (`viewer-plugin.ts`: the plugin on a host element with its canvas props, 3D
  hover label and structure as plain states; `hover-lines.ts`), the stories' `testing/real-plugin-story.tsx` and the
  `Shell/Viewer mode › Real plugin` story. `mockViewerMode({ plugin })` takes the real plugin; without it the mock keeps
  `SequenceStrip` (the other stories have no structure to give the model).
- **The model's input** is `HostStructure { key, title, source: StructureSource, structure }`; the document names its
  file `structure:<key>`, which the loaders' `parsed` hook answers with the host's source (a new key rebuilds the view).
  `SequenceViewer.structure` is set only once the document is in; `legend` and `setFonts` were added (the fonts are
  read once by the control from `--msp-font-*` into the skin: an explicit skin has no `--seq-font-*`).
- **The theme** is `{ letters: scheme name | null }` (default `zappo`, which has labeled groups for the legend); a change
  rebuilds the document (the set stays the same object, so the link survives) and re-applies the section's style.
- **Height**: `ShellModel.setSize` clamped every dock to 180 px, above the top dock's own size (88); a top or bottom dock
  now has its own minimum (40), as it is sized by its view, not a splitter. Tracks follow
  `layout.contentHeight('main')` plus the bar's height (from a `ResizeObserver` entry); compact restores the dock's size.
- **Hover label**: the shell's `HoverLabel` listened to the canvas element only, which never sees the pointer over the
  floating sequence; it now listens on the canvas's parent (the viewport), and the viewer mock merges three sources (the
  sequence's label, the plugin's loci label, the mock's) with the first non-empty winning, since Seq*'s hover also
  highlights 3D, which makes the plugin report a label of its own.
- **Ink from Mol***: the canvas props have no change stream; `PluginContext.events.canvas3d.settingsUpdated` fires for
  settings commands only, so `setBackground` (molstar-poc) announces itself and `viewer-plugin.ts` re-reads the props
  after an MVS load (which sets the background from the state, `background` is therefore passed into it).
- **Subscription guard** (`prism-ui/src/_test/subscriptions.test.ts`) now also scans the viewer's components (`testing`
  excluded).

### Round 2 (as built)

- **Shift** with a click or a drag in the sequence adds or removes each stretch (`select`'s `shiftPositions: 'toggle'`,
  set in the model's viewer spec); Seq\*'s own default stays `extend`.
- **A ligand** in the link: a position of a non-polymer unit (heavy atoms) lights and selects its whole residue in 3D
  (`StructureLink.toLoci` extends to whole residues for non-polymer chains).

- **No bar on the control**: the Sequence button's menu (`sequence/commands.ts`: `sequenceCommands`, `sequenceMenu`) and
  the `LegendCard` overlay replace it; the viewer mode registers the commands, the overlay and the button's `menu` over a
  real plugin. The dock's tab is the band alone.
- **Canvas background** is Canvas settings › Background: `registerMockCanvasSettings(shell, binding)` reads the plugin's
  canvas props and writes `setBackground` / `setTransparent` (the plugin announces its own change, so the ink follows);
  the story's own background bar is gone, the `background` arg only starts it. The story starts with no panel open
  (`ShellStoryProps.noPanel`), and without the mock's tasks and selection bar.
- **The shell**: a top tab with `frame: 'none'` is a band (`FloatingDock`'s `band` class). `ShellView` sets
  `--msp-canvas-end` (the rail, or the pinned panel with it) for the canvas and the band, which end together; the band's
  own `.canvasLayer` already ends at a pinned panel (a first version subtracted it twice: the band was half as wide
  pinned). `--msp-band-fade` and `--msp-band-dissolve` come from `BandFadePx` / `BandDissolvePx`.
- **Seq\* changes** (all default off or unchanged): render props `wrappedAlign`, `wrappedMaxWidth`, `wrappedColors`,
  `unobservedGray`; `focusViewport`'s `reveal` and `revealPad` with `LayoutManager.isWrapped` / `revealY`; the surface's
  canvases carry `data-pass`. Tests in `wrapped-spacing`, `wrapped-reveal` and `draws`.
- **Tooltip**: prism-ui's `installTooltips` showed a tip for an anchor a closing panel had already removed (the timer
  fired late, no `pointerout` came); it now checks the anchor (connected, rendered) when the timer fires and on the next
  pointer event, and watches it with a `ResizeObserver` while shown; a row's card beside its panel is not shown where it
  would lie over the panel.
- **What is not done**: a nucleic acid's letters take the residue classes' colors where the letters coincide (A, C, G, T);
  a ligand or a short unit in the text is a small box at a line's start (Seq\*'s todo.md); a zoomed-out cell class by
  letter ignores `unobservedGray`.

### Round 2, later (as built)

- **A predicted structure** (`Shell/Viewer mode › RealPluginPredicted`): AlphaFill's model of CYP3A4
  (`https://alphafill.eu/v1/aff/P08684`, mmCIF text, AlphaFold's AF-P08684-F1 with ligands transplanted from PDB
  homologs, pLDDT per residue), loaded by `ViewerPlugin.load({ predicted: true, ligands })`: MVS with the polymer's
  color theme `plddt-confidence` (the model-archive extension's provider, registered for the state), and only the heme
  and erythromycin shown in 3D (the transplanted ligands overlap in space; the others stay in the sequence). The sequence
  model builds a **pLDDT track** from `molstar:plddt` as Seq\*'s case 8 does (`values`, `encoding: 'line'`, domain 0–100,
  two rows high) in the four AlphaFold bands softened to the palette (`PrismPlddtBands`), and the legend lists the bands.
  In Text the letters keep the residue palette; a pLDDT tint for them would be one more `wrappedColors` source (not done).
- **A secondary structure track** (both stories, Tracks): `molstar:secondary_structure` (new in molstar-poc: helices,
  strands and turns from the model's own assignment, mmCIF's `struct_conf` and sheets; coil is no range), drawn as Seq\*'s
  `secondary` renderer does, in softened rose, amber and teal, with its legend block.
- **Track context menus act on the 3D view.** A secondary press on a track is Seq\*'s new `contextMenu` behavior
  (`events.menu`, origin `pointer`, with the pointer's place; regions param; the keyboard's menu key was there already).
  What a track offers is the model's table (`TrackMenus` in `model.ts`, by track ref): letters: Apply default coloring;
  secondary structure: Color by secondary structure; pLDDT: Color by pLDDT. The control opens prism-ui's
  `openContextMenu` at the pointer (flips at the window's edges, closes on Escape or a press elsewhere) and the browser's
  own menu is suppressed only over a track that offers one. The actions run through `HostActions`
  (`ViewerPlugin.actions.colorStructure`: one undoable Mol* state update, `updateRepresentationsTheme`, of the polymer's
  representations to `plddt-confidence`, `secondary-structure` or the default).
- **The selection bar is real** over the plugin (`sources/selection.ts`: Mol*'s own label of the selection's statistics,
  one thing named, more counted in lower case; Focus and Clear are wired, the mock's color, theme and measure actions are
  left out) and sits to the left of the toolbar. A ligand's position names its whole residue in 3D
  (`StructureLink.toLoci`).
