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
- `apply(change, { against, context })` sets the root; an equal tree is no change; a change made against another
  revision is refused (`stale-revision`, EC-81). Run it in `history.run(…, SpecRunOptions)` with the run's context: one
  step a change, none for no change, a failure when refused.
- **Changes outside the history**: with a history attached (`attachHistory`), a change without an active run's context
  invalidates its earlier steps. Why: undo restores a whole root, and would silently drop the outside change.
- A preview is never in the history: `cancelPreview` drops it, `commitPreview` ends it in one change.

### ViewModel

`models/view.ts`, `models/axis.ts`: the view shown (`show(address)`, the first view by default) and its resolved view; a
viewport per section, told per section (`viewport(section)`, a keyed state). Every change is written in one batch,
leaders before followers, so a listener of one section reads the others consistent; a viewport equal to the last makes
nothing and tells nothing.

- **A viewport** is `{ start, end, top, width, pending, leader?, pieces? }`: `[start, end)` in offsets along the axis in
  display order (0 at the first position shown), fractional. Why offsets: a section's `range` and a cyclic range over
  the origin read like a whole system; `positionAt` and `offsetOf` convert. A saved session (`session()`, to come)
  stores positions, never offsets: offsets change meaning when the axis does.
- **Scroll groups** (EC-113): sections on equal axes (system, length, pieces) share one range; each keeps its own `top`
  and width (`viewport.groupEqualAxes`, true). An axis whose length is pending has no pieces: its group key holds the
  section's own `range`, so only equal ranges group.
- **Follow** (EC-113, EC-114): a follower's range is its leader's carried along the resolved follow path at its scale
  (1, 3 or 1/3). Each end is the continuous image of the leader's (`point` at the start; at a whole end, the far edge of
  the last element's last place), so it moves smoothly within a run and is exact at element edges, over an intron too;
  where an end maps nowhere (a UTR) it is the end of the hull of the whole positions' image (`followRange`); ends that
  pass each other on a linear follower (a leader range over the coding sequence's origin) give the hull. On a cyclic
  follower a range over the origin wraps: `end` passes the total and `pieces` is `[start, total, 0, end − total]`
  (`follow.cyclicArc`: `wrap`; `hull` takes the image's hull). A follower is clamped by `maxPxPerPosition` as a leader
  is. Follows chain. Followers never lead: a move asked of one is ignored, or detaches it (`follow.followerMoves`:
  `ignore`; `detach` gives it its own range until a fit or `attach`); a move with invalid arguments detaches nothing.
- **Placed** only with a width and a known axis length, else `pending`: a move then is kept as the intent and placed
  when it can be (EC-112). Placing clamps the range within the axis and at most `viewport.maxPxPerPosition` px a
  position. A viewport that shows the whole axis is fitted and stays so across resizes and updates.
- **Resizes** keep the zoom about the center (`viewport.onResize`: `keep-zoom`; `keep-range`); a width of 0 (hidden)
  keeps everything, so a section shown again comes back as it was (EC-112).
- **Updates** (EC-111): axes are compared by value. An equal axis keeps its range; a changed one keeps its center (as a
  position) and zoom when the new axis shows the center, else fits; an axis that was pending carries its intent when the
  length arrives and the new axis starts at 0 (its offsets were positions from 0); a new section fits; a detached
  follower stays detached while its axis is equal.

### HoverModel

`models/hover.ts`: the hover, a method and never a command (it runs per pointer move).

- **A hit** is `{ layer?, section?, system, position?, element?, extent?, target? }`; the section is read from the
  layer's address when absent, the target (`hitTarget`) from what was hit: an element is an item, a position a cell,
  neither empty space.
- `set(hit, owner)` copies the hit into one reused record (`current`) and moves `version`: no allocation, and nothing
  told for the same hit. An element is the same hover by its key, whatever extent object comes with it. One owner at a
  time: another owner's `set` takes it; `clear(owner)` by a stale owner does nothing (EC-141); an undefined owner is
  refused.
- **Projection** is lazy: `location()` (the extent, else the position) and `projected(section)` are computed when asked
  and kept for the version (and the view and model). The location is closed first (`hover.closeCodons`: a nucleotide
  gains its residue and lifts to its codon, every piece over an intron, EC-30, EC-144), then each part is carried to the
  section's system along the view's link (paths found under the copy rule, EC-25); `hover.projectLinked` false keeps
  only the parts on the section's own system.

### StoreModel

`models/store.ts`: the viewer's one `DataStore`, made with `store.cacheBytes`, `store.inFlight`, `store.tileBins` and
timed by `DomHost` (`DomHostClock`). `sync(model, assets)` opens the static model's sources (each bound to a runtime,
and each inline access plan) and attaches them under their addresses; a source whose entry is the object it was opened
from is kept (the resolver keeps entries by hash), a changed one is opened again, one gone is detached and its issue
cleared. A source that cannot be opened is an issue at its address (`issues`); a failed opening, or a coverage the store
refuses (its own report), is forgotten so the next sync tries again. Retired with its owner, the store with it.

## Behaviors

Small policies over the models (`behaviors/`), each reading the bindings.

- `hover`: a move over a hit whose hover gesture is bound to `hover` sets the hover for its owner token; a move over
  nothing or a leave clears what that owner holds only. A host's other hover sources use their own tokens.
- `navigate`: a drag (`pan`), a wheel (`pan` by the larger delta, so Shift with a vertical wheel pans sideways; `zoom`
  about the pointer) and a pinch (`zoom`) through the bindings; wheel deltas within one frame are applied once at the
  frame (`navigate.coalesceWheel`, through `DomHost.requestAnimationFrame`); a section gone by the frame is dropped;
  `dispose` drops what is pending and later frames do nothing. The view's moves are commands (`view.fit`, `view.focus`,
  `view.zoom-in`, `view.zoom-out`, `view.pan-left`, `view.pan-right`) with default keys for the keymap, single keys
  scoped to the viewer (`scope: 'canvas'`); `when` says why one cannot run (no section, a follower that does not move,
  nothing to focus) instead of a throw.

## Bindings

`bindings.ts`: the interaction grammar as data, `{ gesture, modifiers, target, action }`.

- **Pointer gestures** (`hover`, `click`, `double-click`, `drag`, `wheel`, `pinch`); keys go through `@molstar/ui`'s
  keymap, and a `key` binding here is reported. Modifiers are `shift`, `mod` (Cmd on macOS, Ctrl elsewhere,
  `modifiersOf`, which fills a reused record), `alt`, `space`; exact, or `'*'` for whatever is held, tried after the
  exact ones. Why: a modifier held for something else must not kill the hover.
- **Targets** are a hit's kind or `any`; a binding for the kind beats `any`'s.
- **Actions must exist**: a model method (`hover`, `pan`, `zoom`) or a command the viewer defines (`DefaultActions`, or
  the `actions` given). A binding to another is reported (`seq.viewer.binding.action`) and left out, so the gesture does
  nothing. The selection's part of the grammar (`SelectionBindings`) waits for the selection model and its commands, and
  is in no preset until then.
- **Presets**: `default` (hover with any modifiers; drag on an axis, empty space or with Space pans; wheel and
  Shift+wheel pan, Cmd/Ctrl+wheel and pinch zoom; double-click focuses an element, fits on empty space; Alt is unbound
  on purpose, left to the platform and the host) and `genome-browser` (a drag on data pans too).
- `bindingTable(preset, overrides, { actions, kinds })`: each override layer replaces bindings of the same gesture,
  modifiers and target (`action: null` removes one); two of one key within a layer conflict. Resolving is map lookups,
  no allocation.

## Commands

`commands.ts`: `@molstar/ui`'s `CommandDefinition` has no change kind, so the viewer keeps a table by command id:
`patch` (undoable, one step), `session` or `navigate` (`CommandChangeKinds`, `changeKindOf(id, extra)`). Every command
the presets bind has a kind (tested).

## Configuration

`config.ts`: `DefaultViewerConfig`, one object of the shape the models read (`viewerConfig(over)` for one group's items
over the defaults): `foundations.preset`, `density.letterMinPx`, `density.dotMinPx`, `bindings.preset`,
`store.cacheBytes`, `store.inFlight`, `store.tileBins`, `hover.projectLinked`, `hover.closeCodons`,
`viewport.maxPxPerPosition`, `viewport.groupEqualAxes`, `viewport.onResize`, `follow.followerMoves`, `follow.cyclicArc`,
`navigate.wheelZoomStep`, `navigate.keyPanFraction`, `navigate.coalesceWheel`. Each a choice with its recommended
default; `ParamDefinition`s come with the settings form.
