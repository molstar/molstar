# @molstar/seqx-viewer

The viewer: it is composed from an explicit spec (registry entries, behaviors, bindings, configuration, ports), keeps
the interaction state in DOM-free models on `@molstar/core/reactive`, and draws the core's draw lists on a Canvas2D. It
imports `@molstar/seqx-core`, `@molstar/core` and `@molstar/ui`.

This spec holds the contracts of what is built, each a rule and its reason; code comments cite its sections. A section
is added in the change that builds what it describes.

## The entry

- `ViewerLanguageVersion` is the core's `LanguageVersion`: the viewer draws the language its core reads.

## Models

DOM-free `ReactiveModel`s on `@molstar/core/reactive`, tested in Node. Session state (the view shown, viewports, hover)
lives in them, never in the spec or its history. Time and frames are `@molstar/ui/dom`'s `DomHost`, patched in tests.

### SpecModel

`models/spec.ts`: the spec's immutable `root`, its `revision` (`hashNode(root)`), a transient `preview`, `shown` (the
preview, else the root) and `model`, the static model of what is shown by a resolver kept between specs.

- **It is the history's port**: a checkpoint is the root, `restore` sets it back and drops a preview. Why: structural
  sharing makes a root a cheap checkpoint, and undo restores exactly the spec.
- `apply(change, { against })` sets the root; an equal tree is no change; a change made against another revision is
  refused (`stale-revision`, EC-81). Run it in `history.run(…, SpecRunOptions)`: one step a change, none for no change,
  a failure when refused.
- A preview is never in the history: `cancelPreview` drops it, `commitPreview` ends it in one change.

### ViewModel

`models/view.ts`, `models/axis.ts`: the view shown (`show(address)`, the first view by default) and its resolved view; a
viewport per section, told per section (`viewport(section)`, a keyed state).

- **A viewport** is `{ start, end, top, width, pending, leader? }`: `[start, end)` in offsets along the axis in display
  order (0 at the first position shown), fractional. Why offsets: a section's `range` and a cyclic range over the origin
  read like a whole system; `positionAt` and `offsetOf` convert.
- **Scroll groups** (EC-113): sections on equal axes (system, length, pieces) share one range; each keeps its own `top`
  and width (`viewport.groupEqualAxes`, true).
- **Follow** (EC-113, EC-114): a follower's range is its leader's carried along the resolved follow path: the exact hull
  of the whole positions' image (`followRange`; over an intron, the hull), its ends moved in by the leader's fractions
  times the scale (1, 3 or 1/3), so it moves smoothly within a run. Follows chain (leaders computed first). Followers
  never lead: a move asked of one is ignored, or detaches it (`follow.followerMoves`: `ignore`; `detach` gives it its
  own range until a fit or `attach`).
- **Placed** only with a width and a known axis length, else `pending`: a move then is kept as the intent and placed
  when it can be (EC-112). Placing clamps the range within the axis and at most `viewport.maxPxPerPosition` px a
  position. A viewport that shows the whole axis is fitted and stays so across resizes and updates.
- **Resizes** keep the zoom about the center (`viewport.onResize`: `keep-zoom`; `keep-range`); a width of 0 (hidden)
  keeps everything, so a section shown again comes back as it was (EC-112).
- **Updates** (EC-111): an equal axis keeps its range; a changed one keeps its center (as a position) and zoom when the
  new axis shows the center, else fits; a new section fits; a detached follower stays detached while its axis is equal.

### HoverModel

`models/hover.ts`: the hover, a method and never a command (it runs per pointer move).

- **A hit** is `{ layer?, section?, system, position?, element?, extent?, target? }`; the section is read from the
  layer's address when absent, the target (`hitTarget`) from what was hit: an element is an item, a position a cell,
  neither empty space.
- `set(hit, owner)` copies the hit into one reused record (`current`) and moves `version`: no allocation, and nothing
  told for the same hit. One owner at a time: another owner's `set` takes it; `clear(owner)` by a stale owner does
  nothing (EC-141).
- **Projection** is lazy: `location()` (the extent, else the position) and `projected(section)` are computed when asked
  and kept for the version (and the view and model). The location is closed first (`hover.closeCodons`: a nucleotide
  gains its residue and lifts to its codon, every piece over an intron, EC-30, EC-144), then each part is carried to the
  section's system along the view's link (paths found under the copy rule, EC-25); `hover.projectLinked` false keeps
  only the parts on the section's own system.

### StoreModel

`models/store.ts`: the viewer's one `DataStore`, made with `store.cacheBytes`, `store.inFlight`, `store.tileBins` and
timed by `DomHost` (`DomHostClock`). `sync(model, assets)` opens the static model's sources (each bound to a runtime,
and each inline access plan) and attaches them under their addresses; a source whose entry is the object it was opened
from is kept (the resolver keeps entries by hash), a changed one is opened again, one gone is detached. A source that
cannot be opened is an issue at its address (`issues`). Retired with its owner, the store with it.

## Behaviors

Small policies over the models (`behaviors/`), each reading the bindings.

- `hover`: a move over a hit whose hover gesture is bound to `hover` sets the hover for its owner token; a move over
  nothing or a leave clears what that owner holds only. A host's other hover sources use their own tokens.
- `navigate`: a drag (`pan`), a wheel (`pan`, `zoom` about the pointer) and a pinch (`zoom`) through the bindings; wheel
  deltas within one frame are applied once at the frame (`navigate.coalesceWheel`, through
  `DomHost.requestAnimationFrame`). The view's moves are also commands (`view.fit`, `view.zoom-in`, `view.zoom-out`,
  `view.pan-left`, `view.pan-right`) with default keys for the keymap: every move reachable by a command.

## Bindings

`bindings.ts`: the interaction grammar as data, `{ gesture, modifiers, target, action }`.

- **Pointer gestures** (`hover`, `click`, `double-click`, `drag`, `wheel`, `pinch`); keys go through `@molstar/ui`'s
  keymap, and a `key` binding here is reported. Modifiers are `shift`, `mod` (Cmd on macOS, Ctrl elsewhere,
  `modifiersOf`), `alt`, `space`.
- **Targets** are a hit's kind or `any`; a binding for the kind beats `any`'s. The action is a model method (`hover`,
  `pan`, `zoom`, `drag-select`) or a command id.
- **Presets**: `default` (the grammar's table) and `genome-browser` (a drag on data pans, Shift+drag selects).
- `bindingTable(preset, overrides, { commands })`: each override layer replaces bindings of the same gesture, modifiers
  and target (`action: null` removes one); two of one key within a layer, an action that is neither a method nor a
  command with a change kind, and key gestures are issues (`seq.viewer.binding.*`) when the table is built. Resolving is
  map lookups, no allocation.

## Commands

`commands.ts`: `@molstar/ui`'s `CommandDefinition` has no change kind, so the viewer keeps a table by command id:
`patch` (undoable, one step), `session` or `navigate` (`CommandChangeKinds`, `changeKindOf(id, extra)`). Every command
the presets bind has a kind (tested).

## Configuration

`config.ts`: `DefaultViewerConfig`, one object of the shape the models read (`viewerConfig(over)` for one group's items
over the defaults): `foundations.preset`, `density.letterMinPx`, `density.dotMinPx`, `bindings.preset`,
`store.cacheBytes`, `store.inFlight`, `store.tileBins`, `hover.projectLinked`, `hover.closeCodons`,
`viewport.maxPxPerPosition`, `viewport.groupEqualAxes`, `viewport.onResize`, `follow.followerMoves`,
`navigate.wheelZoomStep`, `navigate.keyPanFraction`, `navigate.coalesceWheel`. Each a choice with its recommended
default; `ParamDefinition`s come with the settings form.
