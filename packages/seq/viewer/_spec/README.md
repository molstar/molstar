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

`models/spec.ts`: revisions are the truth, patches derived. It holds the head commit (`SpecCommit`: the spec's immutable
`root`, its `revision`, `hashNode(root)`, and its parents' revisions), a `preview` (a patch), the story `step` shown,
`shown` and `model`, the static model of what is shown, by a resolver that keeps its entries by node hash, so a patch
recompiles the entries whose nodes (or the nodes they read) it changed and a view edit keeps every data entry.

- **Patches** (`applyPatch(patch, options)`): `repr/patch`'s `applyPatch` with `SeqLanguage`, its slots
  (`SeqValidateOptions`) and the checks of a change (`SeqPatchValidator`, the core's `seq.lang` checks of what the
  patch's `Changes` reach): all or nothing, one commit when it changed something (its parent the old head), nothing for
  no change. `requireBase`, `allowDangling` and `allowRebind` pass through; agents and remote channels set
  `requireBase`. Why commits: a recorded root is exact by construction, and an undo restores it.
- **A stale patch** (its `base` not the head's revision, and its `result` not the head's either) is refused
  (`repr.patch.stale`, EC-81), or, when the host opts in (`rebase`, or `edit.rebaseStale`), rebased as
  `merge3(base, current, applyPatch(base, p))` with the same options on both, the base root from `rootOf`. Conflicts
  refuse it with the conflicts and the partial merge (`seq.viewer.rebase.conflict`); an unknown base is
  `seq.viewer.rebase.base`. Why opt-in: a merge the writer did not see is a decision, not a default.
- **Recent roots**: `rootOf(revision)` is the last `revisions.keepRecent` (64) roots made or restored, then the
  history's `root(revision)`.
- **The session follows**: `beforeCommit(f)` is told each commit's `Changes` before its root is set; `AddressModel`
  rewrites the session's addresses through `changes.moved` there.
- **Whole roots** (`apply(change, { against, context })`): a document loaded or a host's own setting, one commit whose
  changes are the whole root; an equal tree is no change; another revision than `against` is `stale-revision`.
- **Previews are patches**: `setPreview(p)` applies `p` over the root with the language's normalization and validation,
  not its own checks (interactive rates); a later call replaces it (made from the root). When the root changes it is
  applied again over the new one (a preview survives an agent's patch mid-drag), and dropped when it no longer applies
  (`seq.viewer.preview.dropped` in `issues`, or quietly: `preview.dropped`); the gesture then ends in nothing.
  `cancelPreview` (Esc) drops it; `commitPreview` applies it in full: one commit or nothing. Never in the history but by
  its commit.
- **Story steps** (`enterStep(id)`, `leaveStep()`): navigation, no commit and no history step. What is shown is
  `step(applyPatch(root, preview))`: the step's stored patch over the root and the preview (the language, refs
  followed). Edits patch the root under the step and the step is applied again over it; one naming a node only the step
  makes is refused (`repr.patch.target`, its hint naming the step); one that sets a key the step sets too lands in the
  root and is reported masked (`seq.viewer.edit.masked`, an info in the result and `issues`). A step that no longer
  applies is shown as not applying (`step.applies` false, `seq.viewer.step.stale`), the root shown without it; an
  unknown id is `seq.viewer.step.unknown`. An agent reads `root`, with `step` beside it.
- **The history** (Ports, History): `attachHistory(port)` binds the host's history with the model as its participant
  (`checkpoint()` the head, `restore(c)` a move of the head that drops a preview and keeps the step, `changed` the last
  commit an edit made). `edit(label, work, { origin })` runs `work` as one host step: every patch it applies (`apply`,
  `commitPreview`) is part of it, a refusal fails it and the host restores its checkpoint; without a history it runs at
  once and a refusal moves the head back. `patch(label, patches)` is an edit of patches in order.
  `beginTurn(label, origin)` opens an agent's turn: its edits join one step. A patch applied outside an edit's context,
  with a history bound, is reported to the host (`outside`). Why the port: the host keeps one stack for everything, and
  the spec knows nothing of it.

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
- **A source-backed system's length** is fed in when its source's coverage gives it (`setLengthSource`, `refreshLengths`
  on each store change): the resolved view's axes (and the follows' leader axes) get their pieces, so a section on a
  streamed sequence is placed and its kept intent applies.
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

### SelectionModel and HighlightModel

`models/selection.ts`: session state as addresses and Locations, never in the spec. This version holds what a command on
an address sets; the gestures, the anchor, the cursor and drawing the highlights are the selection's own work.

- **A value** (`Emphasized`): element keys, a Location, what is drawn (the Location and the elements' extents) and the
  address it was set from; an equal value tells nothing.
- `SelectionModel.value`: the selection; the frame's selection pass draws its `shown`, and `renderSvg` passes it (its
  address when it has one, which the export resolves itself; else its Location), so the SVG shows the selection the
  canvas shows.
- `HighlightModel`: named channels (`set`, `clear`, `channel(name)` a keyed state, `names` in the order first set). Not
  drawn yet.

### AddressModel

`models/address.ts` (`viewer.addresses`): addresses in the spec shown and the store's nodes now (the core's Addressing
over `overlayQueryModel` of the static model and `StoreModel.currentOverlay`, the last frame's builder's overlay without
syncing sources), with the `address.*` settings as bounds.

- **Actions**: `focus` (each leading, movable section whose system the address has positions on, or reaches by a static
  path, shows their hull; one move a scroll group; followers follow), `reveal` (nothing when shown; a pan to its center
  when it fits; else a focus), `select`, `highlightAddress` (a channel, `links.highlightChannel` by default),
  `hoverAddress` (the links' own channel, `links.hoverChannel`, EC-182: the pointer's hover is left alone), `unhover`.
  Each answers an outcome: done, the resolution's status, a message, its issues.
- **A focus that waits** (an address of a node not held, data not loaded) is kept and runs when the address resolves:
  retried on the refs' state, the store's version and its node version; another focus replaces it. A focus whose section
  is not placed yet (no width, no length) is kept as the view's intent and reported not done
  (`Kept until the view is placed`), never `Focused`.
- **The session follows the spec**: `followMoves(changes.moved)`, called before a patch's root is set (the spec model's
  `beforeCommit`, so every commit's moves reach it), rewrites the session's addresses (selection, highlights, the
  waiting and the last focus) through moves and renames as the spec's refs are; on each new root they are resolved
  again, and one that names nothing is dropped with `seq.address.missing` in `issues`.
- **Refs**: `refs`, the spec's refs with what each names now (status, label, color, targets), a new list only when an
  entry changed; `addRef` adds one as one patch (`addRefPatch`, made against the revision just read) through
  `SpecModel.patch`: one step of the history (a promise while the history runs it); a ref that would dangle is refused.
- **Links**: `link(action, address, channel?)`, `parseLink(text)` with the `links.*` settings; `copyLink()`, the link
  (`copyLink.action`) of the selection's address, else of what the pointer is over (its element, else its position).
- **The session as addresses**: `selectionAddress()` (the address it was set from, else the shortest of what it holds),
  `setSelection(address)`; `focusAddress()` (the address last focused while no viewport moved since, else
  `shownAddress()`), `setFocus(address)`; `shownAddress()`, the region the first leading section shows as ids and forms
  (what a bookmark of the view stores, so bookmarks never chain).

### StoreModel

`models/store.ts`: the viewer's one `DataStore`, made with `store.cacheBytes`, `store.inFlight`, `store.tileBins` and
timed by `DomHost` (`DomHostClock`). `sync(model, assets)` opens the static model's sources (each bound to a runtime,
and each inline access plan) and attaches them under their addresses; a source whose entry is the object it was opened
from is kept (the resolver keeps entries by hash), a changed one is opened again, one gone is detached and its issue
cleared. A source that cannot be opened is an issue at its address (`issues`); a failed opening, or a coverage the store
refuses (its own report), is forgotten so the next sync tries again. Retired with its owner, the store with it.

- **The overlay** (tier 2): `overlay(model, view, { runtimes, assets })` gives a frame the overlay of the store's
  dynamic nodes (`OverlayBuilder`: one a node version and model, the model's configuration read through it, another
  builder when the runtimes change; its version is part of the frame's key) and, when it or the view changed, syncs its
  sources (`OverlaySources`: a dynamic annotation's rows, the nested sources the view's templates reach), whose failures
  are issues here while their node is held (cleared when it opens or goes). A rebuilt frame's shown nodes are held
  (`holdShown`, the layout's `shown`) until a later frame no longer shows them. Why here: the store and what is attached
  to it have one owner.

## Composition

`viewer-spec.ts`, `viewer.ts`, `default-spec.ts`: a viewer is composed from an explicit spec, as a Mol* plugin is.

- `SeqViewerSpec`: `registry` (the core's entries), `behaviors` (a provider and its params), `bindings` (a preset,
  override layers, extra action ids), `config` (`[item, value]` pairs), `ports` (factories with ownership), `ui` (the
  parts; none for a headless viewer), `history` (the history port). `DefaultSeqViewerSpec()`
  (`@molstar/seqx-viewer/default-spec`) lists the `seq` runtime, the `hover`, `navigate` and `addressing` behaviors, the
  default ports, `uiHistory()` and the canvas; an application spreads it (a headless composition drops `ui` and
  `history`).
- `createSeqViewer(spec, target?)` makes, in this order: the core's context (registry, configuration, ports), the models
  (spec, view, hover, store, selection, highlight, addresses), the history (`spec.history`, bound to the spec), the
  binding table, the behaviors (one of each name; a behavior's `commands` join the viewer's), then the UI into `target`
  when there is one and the spec has `ui.canvas`. Without a target it is headless: models, commands, frames, `hit`, SVG.
  Why the order: behaviors read the models and the bindings, the UI the behaviors.
- The document is not part of the spec: `load(root, { view? })` sets the root and opens the sources; `settle()` waits
  for the picture (sources open, wanted data arrived, the view at rest).
- **The frame** (`draw`): the renderer lays the view out when anything it reads changed, else redraws the emphasis
  passes alone (a hover); the layout's wants go to the store through the viewer's want set; the backend draws. A frame
  is asked for by `schedule`, one frame callback (`@molstar/ui/dom`, on `DomHost`) however many changes asked; the
  wheel's own frame is asked for first, so its deltas apply in the frame that draws them. Changes that ask: the model,
  the resolved view, a viewport, the hover, the store's version and settledness, a resize.
- **Settled**: the view counts as settled `canvas.settleMs` after its last move and while the store has nothing pending
  (open domains computed, lanes packed).
- **Size**: the width is the target's; the data box is the width less the layout's inset, label column and gap
  (`dataX`), set as every section's width. The height is the layout's content (`canvas.height: 'content'`: laid out at
  the last content height, again when it changed, at most `canvas.maxHeight`) or fixed. The renderer takes one window
  for all sections: the first section that leads gives it.
- **The picture as SVG** (`renderSvg`): the context's export of the shown spec, size and window through the viewer's
  store, at the last frame's device pixel ratio (the density steps are in device pixels, so it is the picture the canvas
  shows); the context's `renderSvg` of the same spec, size and window gives the same string.
- **Links** (`runLink(text)`): a link's actions, in order, each through the model method its command runs (the
  `addressing` behavior's); undefined when the text is no link or the behavior is not listed.
- **Retire** (`retire()`): the viewer is marked retired first (its disposer is registered last, so it runs first), so
  the teardown that follows asks for no frame; the frame callback and the settle timer are cancelled, then the mounted
  UI, the behaviors, the want set, the models and the context's owned ports go. After it no frame is pending and input
  is ignored. A spec that fails to compose (a behavior listed twice) retires what it made.
- **Behaviors by provider**: the input finds a behavior by its provider (`viewer.behavior(HoverProvider)`), whatever the
  spec calls it; one the spec does not list is a `seq.viewer.behavior.missing` issue (`behaviorIssues`), once, and the
  gesture does nothing.
- `rested` moves once each time the picture comes to rest: when the frame laid out settled is drawn, `canvas.settleMs`
  after the view's last move; or `canvas.settleMs` after the last layout that changed what is drawn while the view was
  still (data arriving, the content's height). A move before the settled frame is drawn cancels that rest. What follows
  the picture (the SVG beside it) redraws on it. `load` opens the sources once.
- **Signals publish on change only**: a state derived each frame (`layoutIssues`, `behaviorIssues`, the frame's issues)
  is set only when it differs (lists compared item by item). Why: a listener that redraws on each set (an SVG beside the
  canvas) otherwise redraws every frame, and the churn hides bugs.
- **A signal named for the picture fires on every cause** that changes the picture: `rested` after a move, after data
  arriving and after a resize, each driven in its test. Why: a signal tied to one cause misses the picture that changed
  for another.

## Ports

`ports/index.ts`: the viewer's defaults for the core's port factories (core `_spec`, Context). Assets:
`browserAssetsPort(policy?)`, a `DefaultAssetManager` under the browser policy with the page as base (a host passes its
own manager with `disposeWithViewer: false`). Text metrics: a canvas where there is a DOM (`browserTextMetrics`, the
same CSS fonts the backend draws with), else the core's bundled table.

### History

`ports/history.ts`: Seq* owns no undo stack; the host does. `SeqViewerSpec.history` is a factory given the participant
(`SpecModel`), called once per viewer; none: edits apply at once and nothing undoes them.

- **The participant** (the host drives Seq*): `checkpoint()` a `SpecCommit`, `restore(c)`, `changed`.
- **The history** (Seq* asks the host): `run(step, work, { changed, failed })` an edit as one step (`EditStep`: its
  label, its kind, always `patch`, its origin), `beginTurn(label, origin)` an agent's turn, `isActive(context)`,
  `outside(reason)` a change outside a run (whether it invalidates the earlier steps is the host's rule),
  `root(revision)` a root the host's stack holds, and its ownership (`disposeWithViewer`, `dispose`).
- **Only spec patches enter it**: session, navigation and configuration never do (Commands).
- **Adapters** (`ports/ui-history.ts`, `@molstar/ui/history`), in the composition roots only; the models never import
  it. `uiHistory(options)`: a stack of the spec alone (`UiHistory.stack`), the default spec's; it invalidates on an
  outside change. `hostHistory(stack, participants, key)`: a host's stack over `HistoryParticipants`, one `HistoryPort`
  of several participants by key (Seq*'s spec beside a structure's state), the spec joining while the viewer lives. Both
  resolve `root(revision)` from the commits their stack checkpointed.

## Canvas2D backend

`canvas/`: `CanvasBackend` draws the core's draw lists on three stacked canvases at a device pixel ratio, the order of
the passes.

- The canvases stack as the passes do, bottom to top: the base (background, data), the overlay (highlight, selection,
  hover, cursor, which takes no events) and the labels, so the picture is the SVG's (a test compares the two orders). A
  hover redraws the overlay only: it is cleared and its passes drawn; an overlay with nothing to clear and nothing to
  draw is not touched. A rebuilt layout draws the base and the labels (a canvas cleared first) and then the overlay.
- The scale a canvas is drawn at is its backing size over its CSS size (the backing size is rounded), not the ratio
  asked for.
- A rect with a radius is the core's `squirclePath` replayed; a path is replayed from its commands; a clip is `save` +
  `clip`, a pass left open is closed; a text past its max width fades through a gradient clip (never ellipsized); an
  image goes through one pixel buffer a size, copied at each draw.
- **Reused**: colors are strings kept per value, fonts per shorthand (`FontCache`, a role written as the token table's
  stack, written again when the tokens' stacks change), glyphs per code; nothing is allocated per frame once warm (a
  faded text makes a gradient; a dashed path an array).
- **Tokens** (`canvas/css-tokens.ts`, `readCssTokens(root)`): a mounted viewer reads the token table from the computed
  style of its container: `--msui-*` (`@molstar/ui`'s theme, the chrome tokens and `--msui-font-sans|mono`) and
  `--seq-*` (data colors a host sets) as `createTokens` takes them. A custom property keeps its text (`color-mix(…)`),
  so a value that is no hex or `rgb()` is resolved by the browser through a probe element's `color`. A token the CSS
  does not define falls back to the built-in table of the theme the root shows (`data-msui-theme` on the nearest
  ancestor or the document element, `light` by default; a theme whose `color-scheme` is dark takes the dark table). Why
  not the config's theme: a host's theme switch is CSS, and the viewer must follow it without being told.
- **A theme switch**: a `MutationObserver` on `data-msui-theme` of the document element and of the nearest ancestor that
  has it reads the table again; `viewer.setTokens` keeps the current one when its version is equal (no redraw, no new
  layout key), else the renderer's tokens change, the base canvas is redrawn and `tokensVersion` moves (the SVG beside
  the canvas redraws on it). `retire()` disconnects the observer. `viewer.renderSvg` draws with `viewer.tokens`, so the
  SVG shows the theme the canvas shows.
- `Context2D` and `CanvasLike` are what the backend calls of a 2D context and a canvas: the tests record the calls with
  a fake context.
- `canvasTextMetrics(ctx)`: the text metrics port over a context's `measureText`, the font set only when it changes.

## Input

`input.ts` (`ViewerInput`, `viewer.input`): pointer, wheel and key events as calls in canvas coordinates, DOM-free; the
mounted UI (`mount.ts`) translates events into them and a test makes the same calls. Every gesture goes through the
binding table; per move one hit record and one modifiers record are reused.

- **Pointer**: a move drags when a drag is begun (`down` on a target whose `drag` binding exists) and always hovers
  through the hover behavior. `leave` clears what this pointer owns.
- **A resting pointer is hit again** once a frame in which the layout was rebuilt (EC-143, `pointer.rehit`): the hover
  follows the data under the pointer while the view moves, with no pointer event. Why once a frame: several moves of the
  view in one frame change what is under the pointer once. Not while dragging, nor after a leave.
- **Wheel**: the `wheel` gesture with the event's modifiers (pixels; lines and pages scaled). **Pinch**: a trackpad's
  pinch arrives as a wheel with Ctrl held and small pixel deltas (`pointer.pinchFromCtrlWheel`, at most
  `pointer.pinchMaxDelta`): it is the `pinch` gesture with no modifiers (the Ctrl is the gesture's), zoom factor
  `exp(−deltaY · pointer.pinchSpeed)`. On macOS every Ctrl+wheel is a pinch (Ctrl is no zoom modifier there); elsewhere
  a larger delta is a mouse wheel with Ctrl, which is `mod`: zoom. Line deltas are `pointer.lineHeight` px, a page the
  viewport's own width or height. Over the label column the zoom is about the section's first visible position. True
  when a binding took it, and the page does not scroll.
- **Double click**: the bound command with the section and what to focus (the element's extent, else the position's
  cell, in offsets).
- **Keys**: the commands' shortcuts, run on the section under the pointer (else the first); the `@molstar/ui` keymap
  replaces this when the commands join a command manager.
- **Mounted**: a focusable container with the three canvases, listeners on the container only (a drag ends on
  `pointerup`, `pointercancel` and `lostpointercapture`; Space, held for a drag that pans on data, does not scroll the
  page while the canvas has focus), a `ResizeObserver` for the width, a media query for the device pixel ratio now (made
  again after each change, a frame asked for), a web font arriving measures text again; the function `mountCanvas`
  returns takes away every listener, the observer and the nodes.

## Behaviors

Small policies over the models (`behaviors/`), each reading the bindings.

- `hover`: a move over a hit whose hover gesture is bound to `hover` sets the hover for its owner token; a move over
  nothing or a leave clears what that owner holds only. A host's other hover sources use their own tokens.
- `navigate`: a drag (`pan`), a wheel (`pan` by the larger delta, so Shift with a vertical wheel pans sideways; `zoom`
  about the pointer) and a pinch (`zoom`) through the bindings; wheel deltas within one frame are applied once at the
  frame (`navigate.coalesceWheel`, through a frame callback that `dispose` cancels); a section gone by the frame is
  dropped; `dispose` drops what is pending and later frames do nothing. The view's moves are commands (`view.fit`,
  `view.focus`, `view.zoom-in`, `view.zoom-out`, `view.pan-left`, `view.pan-right`) with default keys for the keymap,
  single keys scoped to the viewer (`scope: 'canvas'`); `when` says why one cannot run (no section, a follower that does
  not move, nothing to focus) instead of a throw.

- `addressing`: the actions on addresses as commands (`address.focus`, `address.reveal`, `address.select`,
  `address.highlight` with a `channel`, `address.hover`, `address.unhover`, `address.copy-link` whose message is the
  link, `refs.add` with a `ref`), so a link, a key, a menu, an agent and a test drive the view the same way; `when` says
  why one cannot run (no address, no ref, nothing to link). `runLink` runs a link.

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
`patch` (undoable, one step), `session`, `navigate` or `none` (it only reads) (`CommandChangeKinds`,
`changeKindOf(id, extra)`). Every command the presets bind has a kind, and every `addressing` command (tested): a focus
or a reveal navigates, a selection, a highlight or a link's hover is session state, `refs.add` a patch; copying a link
changes nothing (`none`).

## Configuration

`config.ts`: `DefaultViewerConfig`, one object of the shape the models read (`viewerConfig(over)` for one group's items
over the defaults): `foundations.preset`, `density.letterMinPx`, `density.dotMinPx`, `bindings.preset`,
`store.cacheBytes`, `store.inFlight`, `store.tileBins`, `hover.projectLinked`, `hover.closeCodons`,
`viewport.maxPxPerPosition`, `viewport.groupEqualAxes`, `viewport.onResize`, `follow.followerMoves`, `follow.cyclicArc`,
`navigate.wheelZoomStep`, `navigate.keyPanFraction`, `navigate.coalesceWheel`, `canvas.height` (`content`),
`theme.source` (`root`: a mounted viewer reads its root's CSS; `config`: the core's `theme` item; a headless viewer
always uses the config), `canvas.maxHeight` (2048), `canvas.dpr` (`device`), `canvas.settleMs` (150), `pointer.rehit`
(true), `pointer.pinchFromCtrlWheel` (true), `pointer.pinchSpeed` (0.01), `pointer.pinchMaxDelta` (50),
`pointer.lineHeight` (16), `address.limit` (10,000), `address.maxDepth` (64), `address.refDepth` (8), `links.prefix`
(`seq-`), `links.maxLength` (8,192), `links.maxActions` (16), `links.highlightChannel` (`link`), `links.hoverChannel`
(`link-hover`), `revisions.keepRecent` (64), `edit.rebaseStale` (false: a stale patch is refused; true: rebased through
`merge3`), `preview.dropped` (`issue`; `silent`), `copyLink.prefer` (`['ref', 'ids']`), `copyLink.numbering`
(`ordinal`), `copyLink.action` (`focus`). Each a choice with its recommended default; `ParamDefinition`s come with the
settings form. As items in a spec: `viewerItem('group.item')` (typed by key), read by `viewerConfigOf(pairs)`; the
core's items (`SeqCoreConfig`) are the same pairs.
