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

## 2. A "Sequence viewer" extension, composed into the viewer mode

No new package (David, 2026-10-14: fewer packages, not more): an `extensions/` folder in the viewer mode's package,
`packages/prism/viewer/src/extensions/sequence/`, exported as `@molstar/prism-viewer/extensions/sequence` and composed
into the viewer mode as contexts.md §3–4 composes anything: a **module** that owns the Seq\* contexts and a
**contribution** (the dock tab).

- prism-ui stays free of `@molstar/seq-viewer`; prism-viewer gains the Seq\* dependencies, and only an app that
  composes the extension bundles Seq\* (the folder is its own export, imported by nothing else in prism-viewer). This
  settles plan.md §3's "prism-ui depends on seq-viewer" in favour of the viewer mode's package, its own alternative.
- The extension imports only prism-ui's public exports and prism-viewer's own model types, never the mock's
  internals, so it can move to a package later without rewriting.
- prism-viewer then depends on `@molstar/seq-viewer`, `@molstar/seq-core` and
  `@molstar/seq-molstar-poc` (the structure provider, `StructureLink`, `linkViews`, `StructureView3D`). molstar-poc is
  private and named for a proof of concept; turning its link and provider into a real package is a follow-up (§7).

### 2.1 The module: who owns the `SeqViewerContext`

The extension's module is created with the app and owns every Seq\* context (contexts.md §3: "a module that owns an
outside context creates it in `create` and disposes it in `dispose`").

```ts
export const SequenceViewers = ModuleKey<SequenceViewers>('sequence-viewers');

interface SequenceViewers {
  /** The viewers by id; the first slice creates only 'default'. */
  get(id?: string): SequenceViewer;
  readonly ids: ReadonlyState<readonly string[]>;
}

interface SequenceViewer {
  readonly id: string;
  readonly seq: SeqViewerContext; // createSeqViewer(spec, ports), disposed with the module
  readonly style: ReadonlyState<'compact' | 'tracks'>;
  setStyle(style: 'compact' | 'tracks'): void; // ctx.ui.setStyle(section, ...), UI state only
  readonly structure: ReadonlyState<SequenceStructure | undefined>; // what it shows (§3)
  readonly theme: ReadonlyState<SequenceTheme>; // the color scheme the letters use, and later each track's (§6)
  setTheme(theme: SequenceTheme): void; // written into the document the module builds
}
// create(ctx): subscribes to ctx.canvasInk (§5) and keeps each viewer's skin in step
```

- Keyed by id from the start, so a second viewer later is another entry plus another contribution, not a redesign.
- The context is created with `rafScheduler` (Seq\*'s default scheduler never fires), the default registry plus
  molstar-poc's `MolstarSvsEntry`, and loaders that refuse the network: the structure comes from the host (§3).
- The mocks have no module slot yet (contexts.md §8–9, "not built"). The first slice adds the smallest one the viewer
  mode needs (a `modules` registry created and disposed with the mode), as contexts.md specifies it, not a parallel
  mechanism.

### 2.2 The contribution: the dock tab

`{ id: 'sequence', title: 'Sequence', icon: 'sequence', dock: 'top', when: <a structure is shown>, render:
domView(SequenceControl) }`, replacing the mock's `SequenceStrip` registration in the viewer mode. The existing
"Sequence" canvas button and the `view.sequence` command keep working, since they address the tab.

### 2.3 The link: subscribed while the control is rendered

`SequenceControl` subscribes the two-way link when it mounts and drops it on cleanup:

```ts
const unlink = linkViews(viewer.seq, view3d, link, options); // hover both ways, selection both ways, focus
onCleanup(() => { unlink(); link.detach(); });
```

- `view3d` is a `StructureView3D` (molstar-poc's Mol\*-agnostic interface: `highlight`, `select`, `focus`, `onHover`,
  `onSelect`). Over the real plugin it is molstar-poc's `pluginView3D(plugin, structure)`, which lives in the
  extension's `sources/` folder (contexts.md §7: only `sources/` and `actions/` import `@molstar/plugin`). In the mocks
  it is a mock that feeds the mock's hover label and status line.
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
- In the mocks there is no plugin; the extension parses a fixture structure itself (the mock's 1TQN, added as a small
  `.bcif` fixture, so the letters match the mock's 3D labels).
- The structure picker in the strip's header (structure, chain/polymer mode, instance) is not part of the first slice.

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
- **No frame**: `FloatingDock` gets a frameless variant (no `--msp-glass-panel-bg`, shadow or filter on its content);
  the canvas shows through. The surface's canvases clear rather than fill, so the viewer draws on the canvas
  directly. Seq\* draws its text halos in the skin's surface color (§5), which keeps letters legible on any background.

## 5. The ink: a behavior the context provides

The canvas ink (the canvas background, and whether chrome over it is light or dark) is a behavior the Prism context
provides; the extension subscribes to it and updates its own state. Nothing in the extension reads the DOM or CSS to
find it.

- **The behavior**: the viewer mode's context exposes `canvasInk: ReadonlyState<CanvasInk>`, `CanvasInk = { theme:
  'light' | 'dark'; color: string }`. It is the `canvasTheme(props$)` that prism-ui already computes
  (`ui/src/model/canvas-theme.ts`, WCAG luminance at 0.179) plus the base color, today only handed to the shell as
  `canvasTheme` / `canvasColor` props. The shell's canvas layer and the extension read the same state; in the mocks
  it comes from the mock's background, over the plugin from the canvas props source (view-models.md §3, `canvas`).
- **The subscription is the module's**, not the control's: each viewer follows the ink while it exists, so its state
  is already right when the control mounts. On a change it builds Seq\*'s `Skin` from the ink and sets it
  (`seq.setSkin(skinFor(ink))`): the surface and dark flag from the ink, ink and accent colors from Prism's palette
  for that theme, the selection as the mock strip draws it (the accent at 25%). Seq\*'s data themes have dark
  variants; a skin change reassembles the view with them.
- **Not used**: Seq\*'s skin reader (it reads `--seq-*` from computed style and watches `data-msp-theme` attributes).
  An explicit skin keeps one source of truth and works the same in a shadow root, a story and a headless test. The
  control still sets `--seq-font-*` from `--msp-font-*`; fonts are not part of the ink.

## 6. Theme, tooltip and legend

The extension manages its themes; Seq\* colors. Which color scheme the letters use (and, later, each track) is the
module's state (`SequenceViewer.theme`), next to the ink it already follows (§5). The module writes the choice into
the document it builds; Seq\*'s themes and renderers draw the colors from it, so nothing about coloring is
reimplemented in Prism.

- **Legend**: a small prism-ui `Legend` component (items: swatch or glyph and label; ramps: gradient, low and high),
  opened from a button on the control. Its data is the extension's own: the scheme it chose, read from the scheme's
  public definition (`KernelColorSchemes`, `ColorSchemeEntry` from `@molstar/seq-core/format`). Seq\*'s legend model
  (`legendGroups` in `viewer/src/chrome.ts`) stays internal to its DOM and React UIs: no new Seq\* entry.
- **Tooltip**: the sequence's hover feeds the shell's one cursor-following `HoverLabel` (shell.md: the sequence's
  hover feeds the same hover label as the 3D canvas). The control maps `seq.labels.current` (`{ title, lines }`,
  already public) to the label model's lines while the pointer is over the surface. No per-move rendering: lines
  change only when the label does.
- Track headers stay off (`headerWidth: 0`), so header tooltips, track menus and the hover card are not in play.
- **Later, data tracks**: a legend whose bands or range come from the data or the user's edits (pLDDT, a continuous
  annotation) is computed by Seq\*'s renderers today. With fixed bands (pLDDT's) the extension knows them from the
  theme it set; for the rest, Seq\*'s per-track legend (`trackLegends`) gets an export then (§8).

## 7. Work, in order

| #   | What                                                                                                                                                                                         | Who            |
| --- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------- |
| S2  | prism-ui: the module slot (contexts.md §3, minimal), the `canvasInk` behavior on the viewer mode's context (the shell reads it too), `FloatingDock`'s frameless variant, the `Legend` component with a story                                                                 | Sonnet, medium |
| S3  | `prism-viewer/src/extensions/sequence`: module (the ink subscription and `skinFor`, the theme state and the legend data from it), document, control, link over `StructureView3D`; the mock `StructureView3D` and the 1TQN fixture                                         | Sonnet, high   |
| S4  | the viewer-mode mock composes the extension instead of `SequenceStrip`; the shell stories (including `CanvasBackground` for the ink); the subscription guard extended to the new package | Sonnet, medium |
| S5  | Visual check in Storybook                                                                                                                                                                    | David          |

S3 needs S2; S4 needs S3. Over the real plugin, the same extension gets its
`pluginView3D` in `sources/` once Prism's plugin sources exist (view-models.md §3); nothing in the control changes.

## 8. Follow-ups

| Item                                                                                                     | Trigger                              |
| -------------------------------------------------------------------------------------------------------- | ------------------------------------ |
| molstar-poc's link and structure provider as a real package (`@molstar/seq-molstar`), out of the PoC      | the extension leaves the mocks       |
| Several viewers (a second id, its own tab or dock)                                                       | a feature that needs two             |
| The strip header's structure picker (structure, chain or polymer mode, instance)                         | more than one structure loaded       |
| Tracks beyond letters (secondary structure, pLDDT, annotations) and their legend; Seq\*'s `trackLegends` exported for data-driven ranges | asked for in Prism                   |
| The `--msp-ink` tokens (defined, unused) for chrome outside the canvas layer                             | a canvas control outside that layer  |
| todo.md's "the Mol\* hover sometimes needs a second move after the pointer was over the 3D canvas"       | the link over the real plugin        |
