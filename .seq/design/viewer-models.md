# The viewer's models

## Intent

The second review found the viewer's models unlisted, session state deferred though the first slice builds it, the
viewport and follow undesigned though "p53's protein on its gene" needs them, hover and selection without a model, the
gesture bindings without a type (`@molstar/ui`'s keymap knows key chords only), hover as a command allocating on every
pointer move, and no field for a command's change kind. This design lists the models and settles those.

## Decision

### The models

All are `ReactiveModel`s on `@molstar/core/reactive`, DOM-free, tested in Node ([sources/molstar-ui](../sources/molstar-ui/README.md)):

| Model | Holds | Change kind ([changes.md](changes.md)) |
|---|---|---|
| `SpecModel` | the spec's root (immutable), its revision (the root's hash), the preview overlay; applies patches; the `HistoryPort` checkpoint is the root | spec patch, preview |
| `ViewModel` | the view shown; per section its viewport; scroll groups and follows | session |
| `HoverModel` | the hover: a Location, an element key, its owner, its system | session |
| `SelectionModel` | the selection: element keys and a Location, the anchor, the cursor | session |
| `HighlightModel` | named highlight channels, each a Location or element keys | session |
| `StoreModel` | the data store ([data-store.md](data-store.md)) | data |
| `ConfigModel` | the configuration, read-only after composition | configuration |

- **Session state** is the `ViewModel`, `HoverModel`, `SelectionModel` and `HighlightModel`'s state, expressed in
  addresses and Locations; a host may save it beside the spec (`viewer.session()` and `viewer.restore(session)`), never in
  it.

### Viewport, scroll groups and follow

- **A viewport per section**: the visible range on its axis (in positions, fractional) and the vertical scroll; zoom is
  the range's width against the section's width.
- **A scroll group** shares one viewport among sections on equal axes (EC-113).
- **A follow** is a section whose viewport is derived from another's through a map: the leader's range projected along
  the path between their systems, at scale 1, 3 or 1/3 (a codon step); followers never lead (EC-114). For "p53's protein
  on its gene", the gene section follows the protein section through the codon map.
- **Viewports survive updates**: an equal axis keeps its viewport; a changed one keeps center and zoom when the center
  is still shown, else fits (EC-111); a move before the section has a width is kept as an intent (EC-112).
- Collapsed introns (layout maps) come after the first slice; the slice's gene section shows the gene's full range.

### Hover and selection

- **Hover** is a method on `HoverModel`, `hover.set(hit, owner)` and `hover.clear(owner)`, never a command: it runs per
  pointer move and allocates nothing (reused records). One owner at a time; a stale owner's clear does nothing (EC-141).
- **Projection**: the hover's Location is projected onto every linked section's system when the resolved view links
  them, along their paths (static or overlay), lazily per frame for the sections drawn, under the copy rule (EC-25); a
  nucleotide lifts to its codon (closure, EC-30, EC-144).
- **Selection** holds element keys (items, rows, chains) and a Location (positions); the grammar's toggle acts on
  element keys, extend on positions from the anchor in one system ([ux-foundations.md](ux-foundations.md)).
  Projections of the selection are drawn on linked sections as the hover's are.

### Behaviors (the first list)

`hover` (pointer to `HoverModel`), `select` (clicks and drags to `SelectionModel`), `navigate` (drag, wheel, pinch,
keys to viewports), `follow` (keeps followers in step), `cursor` (the keyboard cursor), `focus` (focus and reveal), each
with its params as config items. Later: `hoverCard`, `contextMenu`, `snapshots`, `reveal-on-highlight`.

### Gestures, bindings and commands

- **A gesture binding** is `{ gesture, modifiers, target, action }`: gestures are `click`, `double-click`, `drag`,
  `wheel`, `pinch`, `key`; targets are the hit's kind (cell, item, row, column, track, empty, axis); actions name either a
  model method (hover, drag-select, pan: continuous, allocation-free) or a command (discrete). Key gestures go through
  `@molstar/ui`'s `Keymap`; pointer gestures through this table, which the grammar's presets fill
  ([ux-foundations.md](ux-foundations.md), section 8).
- **Commands** are `@molstar/ui/commands` definitions with Seq\*'s metadata: `changeKind: 'patch' | 'session' |
  'navigate'`, kept beside the definition (`@molstar/ui`'s `CommandDefinition` has no such field; Seq\* keeps a table by
  command id, and proposes the field upstream if it proves general).

### Configuration items of the first slice

`foundations.preset` (`comfortable`), `density.letterMinPx`, `density.dotMinPx`, `bindings.preset` (`default`),
`store.cacheBytes` (64 MB), `store.inFlight` (6), `store.tileBins` (1,024), `hover.projectLinked` (true). Each a
`ParamDefinition` with a default; the slice may read them from a default object of this shape.

## Consequences

- The first slice builds `SpecModel`, `ViewModel` with scroll groups and one follow, `HoverModel` with projection,
  `StoreModel`, and the `hover` and `navigate` behaviors; selection and the rest follow it.
- Tests in Node: viewport rules EC-111 to EC-114, hover ownership EC-141, projection through the codon map, the
  bindings' table resolving gestures by target.
