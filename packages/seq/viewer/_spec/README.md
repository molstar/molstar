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
  `cancelPreview` (Esc) drops it; `commitPreview` applies it in full, in one batch with the preview's end: one commit or
  nothing. A move applies the patch once: what is shown is reused from `setPreview`'s result until the root or the
  preview changes. Never in the history but by its commit.
- **Story steps** (`enterStep(id)`, `leaveStep()`): navigation, no commit and no history step. What is shown is
  `step(applyPatch(root, preview))`: the core's `stateOf` over the root and the preview (a full state, the step's patch
  over it, or a delta from the step it starts `from`; the language, refs followed), cached by the revision shown and the
  step; `step` carries its `form` and its `state` (view, timing). Edits patch the root under the step and the step is
  applied again over it; one naming a node only the step makes is refused (`repr.patch.target`, its hint naming the step
  only when the step makes the target); one that sets a key the step sets too lands in the root and is reported masked
  (`seq.viewer.edit.masked`, an info in the result and `issues`). A step that no longer applies whole is shown as not
  applying (`step.applies` false, `seq.viewer.step.stale`): the ops that still apply are shown (the design: a step's
  state is what still applies, flagged), and a start it cannot reach is the core's `seq.story.step.from`; an unknown id
  is `seq.viewer.step.unknown`. **While a full-state step is shown, an edit is refused** (`seq.viewer.step.frozen`, its
  hint to edit the step), a patch, a preview and a whole root alike (`apply`: `frozen-step`): the step does not follow
  the spec, so the edit would show nothing. `load` leaves the step first: a new document.
- **The frame slot** (`setFrame(root)`, `frame`): an animation's frame shown over what is composed, `frame(step(...))`:
  no patch, no revision, no history step; the spec's hash never moves during an animation. `undefined` leaves it. The
  viewer ends a frame that changes the spec when the head moves under it (an edit, an undo, a load): the frame was made
  from the old root and would hide the change. An agent reads `root`, with `step` beside it.
- **The history** (Ports, History): `attachHistory(port)` binds the host's history with the model as its participant (a
  history bound before is detached first) (`checkpoint()` the head, `restore(c)` a move of the head that drops a preview
  and keeps the step, `changed` the last commit an edit made, told once the edit ends well: never one it took back).
  `edit(label, work, { origin, within })` runs `work` as one host step: every patch it applies (`apply`,
  `commitPreview`) is part of it, a refusal fails it and the host restores its checkpoint; without a history it runs at
  once and a refusal moves the head back. `within`, a host's run in progress (its work runs a Seq* command), makes the
  edit part of that run: no step and no wait of its own (else the two would wait on each other). `patch(label, patches)`
  is an edit of patches in order; `patches` may be a function of the root and revision the edit runs on, so a command's
  patch is made when it runs and commands issued together each apply. `beginTurn(label, origin)` opens an agent's turn:
  its edits join one step. A patch applied outside an edit's context, with a history bound, is reported to the host
  (`outside`). Why the port: the host keeps one stack for everything, and the spec knows nothing of it.
- **The session follows every move of the head**: `beforeCommit(f)` is told a commit's changes, and on an undo, a redo,
  a restore or a rollback the moves between the two commits along the recent lineage (each recent commit's parent and
  moves: forward going to a child, inverted going back); a rebase tells theirs' moves, which `merge3` does not infer.
  Why: a session address renamed by an edit must come back with its undo.

### ViewModel

`models/view.ts`, `models/axis.ts`: the view shown (`show(address)`, the first view by default) and its resolved view; a
viewport per section, told per section (`viewport(section)`, a keyed state). Every change is written in one batch,
leaders before followers, so a listener of one section reads the others consistent; a viewport equal to the last makes
nothing and tells nothing.

- **A viewport** is `{ start, end, top, width, pending, leader?, pieces?, scale? }`: `[start, end)` in offsets along the
  axis in display order (0 at the first position shown), fractional. Why offsets: a section's `range`, a cyclic range
  over the origin, a reversed axis and collapsed pieces read like a whole system; the core's `positionAt`, `offsetOf`,
  `cellOffset` and `elementAt` convert through the pieces and the direction (`models/axis.ts` re-exports them: one
  implementation, the renderer's too). A saved session (`session()`, to come) stores positions, never offsets: offsets
  change meaning when the axis does. Navigation is in display order (a drag, the wheel, the pan keys, the cursor's
  arrows move what is on screen the way the hand or key goes, whatever the direction).
- **Layouts** (the core's "Section layouts", Navigation): each move asked of a section is its layout's arithmetic
  (`SectionNav`: `pan` by a scroll's px deltas, `zoom` about a unit, `fit`, `reveal` the least move, `focus`, `step`,
  `clamp`) over its placement in the last frame (`setLayoutSource`: the viewer's frames, and the one map of layouts the
  frames draw with, so a layout navigates as it is drawn; `placed()` after each frame reads its estimates again). The
  record is read the layout's way: tracks a window; text `start` its anchor (the unit at the box's top, fractional
  within its line), `scale` its letter size (none: the preset's), `end` the lines in view's estimate; `setRange` is a
  focus (a text brings the range's first line to the top, its size kept). Before a width or a length a move is the
  intent, a window of offsets. Why the layout's: the viewer's own arithmetic was the tracks layout's, wrong in text and
  on rings. A leader's window passes a cyclic origin only where its layout `wraps` (draws a window there: a ring's;
  `NavContext.cyclic`); tracks and text keep it within the axis (their transform ends at the total).
- **Scroll groups** (EC-113): sections on equal axes (system, length, pieces, direction) of one layout share one range;
  each keeps its own `top` and width (`viewport.groupEqualAxes`, true). Why one layout: a range is read its layout's way
  (a window, an anchor); across layouts, follow. An axis whose length is pending has no pieces: its group key holds the
  section's own `range`, so only equal ranges group.
- **Follow** (EC-113, EC-114): a follower's window is the image of its leader's window edges in display order
  (`followWindow`, the core's, the renderer's too): the leader's offsets read as positions (a break's sides as a
  window's start and end), carried along the follow path, read as offsets along the follower's own pieces in its
  direction. Each end is the continuous image of the leader's (`point` at the start; at a whole end, `projectEnd`, the
  far edge of the last place), place by place (each nucleotide of a split codon its own edge), so it moves smoothly and
  is exact at element edges; on a coding axis (`pieces: coding`, no `pad`, no `range` cutting the coding sequence) the
  follower's offsets are exactly three times the leader's positions; where an end maps nowhere (a UTR) it is the end of
  the hull of the whole positions' image (`followRange`); ends that pass each other on a linear follower (a leader range
  over the coding sequence's origin) give the hull. On a cyclic follower a range over the origin wraps: `end` passes the
  total and `pieces` is `[start, total, 0, end − total]` (`follow.cyclicArc`: `wrap`; `hull` takes the image's hull). A
  follower is clamped by `maxPxPerPosition` as a leader is. Follows chain.
- **A follower leads** (`follow.followerMoves`: `lead`): a pan, a zoom or a pinch asked of an attached follower moves
  its leader through the inverse of the forward map (each end to the leader place whose `point` image it is), so the
  follower shows what was asked and every section moves in step; the ends are read in positions (on a reversed axis an
  end moved right moved down), so with `follow` a pan of the gene moves the protein the same way on screen, and with an
  explicit `genomic` on a minus strand the other way. Where an end is in no place's image (an intron, a UTR) it takes
  the nearest image the way it moved, so a pan steps over a gap and never sticks. A pan moves the follower's image of
  the leader's window before its own clamp, so the leader keeps its width (the window a clamp widened, carried back,
  would zoom it out); one whose other end lies past every image moves the leader by the one end's step. A zoom or a
  range with an end moved outward past every image takes the leader's axis edge on that side, so a zoom out never
  stalls; inward, the leader's end stays. A zoom or a pinch acts on the window shown: past the follower's clamp, the
  leader re-centers on what the follower shows, as a leader's own zoom does. A chain leads up to its first leader,
  bounded as `canMove` is (16 follows: a cycle moves nothing). `canMove` is the leader's. Why the inverse of the forward
  map, not of the path: what the follower shows is the forward map's, so only its inverse brings back the range asked.
  `ignore` drops the move; `detach` gives the follower its own range until a fit or `attach`; a move with invalid
  arguments moves and detaches nothing. Focus and reveal by address move the leading sections as before, and the cursor
  on a follower reveals nothing unless it may detach (`reveal.cursor`).
- **Bands** (`bands`, `band(track)`, `setBand(track, band)`, `bandsOf(section)`): a band's own viewport by its track's
  address, a rows band's scroll (`bandTop`, the band's top in display rows below its pinned band, fractional, as a
  step's view declares it; the layout turns it into px), row zoom (`rowScale`, px a row) and folded clades (`folds`, by
  name: session state, the layer's `collapsed` when absent; `folds: undefined` drops them). One map, replaced on each
  change (compared by identity: a frame's key), read into each leading section's frame viewport (`viewports[s].bands`),
  so the export draws the viewer's folds and scroll. **During a story's motion** the shown frame's `bandTop` and
  `rowScale` are merged over the model's bands, the session's folds kept (`_bandsFor`); a step's bands are placed in the
  model when it is entered (top and zoom; never its folds), and the model's are where a step's motion starts. A step of
  the history keeps each band by its top row's key (`KeptBand`: the key and the fraction past it, with its zoom and
  folds), the snapshot's own field beside the groups, restored under its folds and zoom and clamped (a folded key at its
  summary row, a hidden one at the nearest row shown) through `BandRows`, what names the rows (the `rows` behavior gives
  it, `setBandRows`); with none, a band is kept as its `bandTop`. Why keys: a step back after rows were hidden, folded
  or reordered lands on the same rows. `rowsAxis(system)`: a rows system as an axis in its layer's display order
  (Shift+click extends along it).
- **History** (`jump`, `back`, `forward`, `history`): a jump (a focus, a fit, a reveal, a bookmark's or a link's move,
  the cursor's Enter) runs through `jump(move)`, which pushes the viewports it left (every group's range, a detached
  follower's among them) when it moved the view, and drops the steps forward; `same` when the view already showed it,
  `none` when nothing moved. A pan or a zoom by hand pushes nothing, so a jump after it pushes where the hand left the
  view; a step zoom (`jump(move, false)`) animates without a step and, like a move by hand, keeps the steps forward. A
  step keeps each group's offsets with its signature (its axis's system, pieces and direction, and its layout), its
  centre as a position, its letter size and its members; steps compare by offsets on equal signatures only; onto a
  changed axis or layout (a reverse, a collapse, a layout patch since) a step is restored by position, the centre and
  the zoom (a window's middle; a text's lines in view, so its anchor lands on the window's first position and a window
  on the text's lines), the letter size within one layout, the centre clamped to what is drawn (one now hidden goes to
  its break), never dropped. `back` returns to the last viewports left, every leading section's range together
  (followers follow; one detached since is attached again, one attached since is detached again at its range then; a
  section gone since is left out), `forward` goes again; a step to where the view is now is skipped. A new document
  (`load`) clears the history; an edit of the same spec keeps it. Bounded (`view.historyLimit`, 50; the oldest goes).
  Session state: never the spec's history (view state is not history), so undo leaves it. `jumps` moves on each jump,
  back and forward.
- **A source-backed system's length** is fed in when its source's coverage gives it (`setLengthSource`, `refreshLengths`
  on each store change): the resolved view's axes (and the follows' leader axes) get their pieces, so a section on a
  streamed sequence is placed and its kept intent applies.
- **Placed** only with a width and a known axis length, else `pending`: a move then is kept as the intent and placed
  when it can be (EC-112). Placing clamps the range within the axis and at most `viewport.maxPxPerPosition` px a
  position. A viewport that shows the whole axis is fitted and stays so across resizes and updates. **A cyclic axis's
  viewport wraps**, leader or follower, in tracks and rings (their layouts' `wraps`; a text's anchor does not): its
  start goes round the origin and its end may pass the total, never clamped there, so a ring's window and a linear view
  pan across the origin; a follower of it shows the window over the origin as two pieces. One rule,
  `view.wraps(section)`, for the navigation, the window's drags and the input; `view.arcs(section)` reads two ends as
  the shorter arc (a ring, or a linear window over the origin now).
- **A viewport as given** (`setViewport(section, { start, end, top?, scale? })`): a leading section takes the window,
  the letter size and the top exactly, clamped as any move; no history step. What a story step's view and a motion
  stopped where it is drawn use; a follower follows and takes none.
- **Resizes** keep the zoom about the center (`viewport.onResize`: `keep-zoom`; `keep-range`); a width of 0 (hidden)
  keeps everything, so a section shown again comes back as it was (EC-112).
- **Updates** (EC-111): axes are compared by value. An equal axis keeps its range; a changed one keeps its center (as a
  position) and zoom when the new axis shows the center (a center a collapse now hides goes to its break), else fits; an
  axis that was pending carries its intent when the length arrives and the new axis starts at 0 (its offsets were
  positions from 0); a new section fits; a detached follower stays detached while its axis is equal.

### HoverModel

`models/hover.ts`: the hover, a method and never a command (it runs per pointer move).

- **A hit** is `{ layer?, section?, system, position?, ySystem?, yPosition?, element?, extent?, target? }`; the section
  is read from the layer's address when absent, the target (`hitTarget`) from what was hit: an element is an item, a
  position a cell, neither empty space; a layer's own target wins (`row`: a rows header; `cell` with `ySystem` and
  `yPosition`: a rows or grid cell; `tree-node`; `summary-row`). The record keeps `ySystem` and `yPosition` (NaN when
  none) flat beside the rest: a hit of two axes allocates nothing.
- `set(hit, owner)` copies the hit into one reused record (`current`) and moves `version`: no allocation, and nothing
  told for the same hit. An element is the same hover by its key, whatever extent object comes with it. One owner at a
  time: another owner's `set` takes it; `clear(owner)` by a stale owner does nothing (EC-141); an undefined owner is
  refused.
- **Projection** is lazy: `location()` (the extent, else the position) and `projected(section)` are computed when asked
  and kept for the version (and the view and model). The location is closed first (`hover.closeCodons`: a nucleotide
  gains its residue and lifts to its codon, every piece over an intron, EC-30, EC-144), then each part is carried to the
  section's system along the view's link (paths found under the copy rule, EC-25); `hover.projectLinked` false keeps
  only the parts on the section's own system. `projectedAll()` is the same for every section at once (the core's
  `projectOntoView`, one Location keyed by system), what the frame draws (Emphasis).

### SelectionModel and HighlightModel

`models/selection.ts`: session state as addresses and Locations, never in the spec.

- **A value** (`Emphasized`): element keys, a Location, what is drawn (the Location and the elements' extents) and the
  address it was set from; an equal value tells nothing.
- `SelectionModel.value`: the selection; the frame's selection pass draws its `shown`, and `renderSvg` passes it (its
  address when it has one, which the export resolves itself; else its Location), so the SVG shows the selection the
  canvas shows.
- **The grammar's gestures are methods** (the interaction grammar), by the hit: an item is an element (its key, drawn by
  its extent), a cell a position. `select(hit)` sets it, and clears it when it was exactly that element; nothing hit
  clears it. `extend(hit)` makes the range the anchor is in reach from the anchor to the hit on the anchor's system, the
  other ranges and elements kept (extends from one anchor replace one another; with no anchor there, a select;
  `selection.extendScope`: `selection` makes the whole selection that range). `toggle(hit)` adds or removes the element
  or the cell, keeping the rest. `add(location, anchor?)` and `subtract(location)` add positions or take them out (an
  element they cut becomes the positions left of its extent; one they miss stays an element). An add with an anchor
  moves it there; one without (an address's positions) keeps the anchor, and what an extend keeps gains what was added.
  A subtraction keeps the anchor; what an extend keeps becomes what is left less the run the anchor is in, so a part
  split off stays. `beginDrag(hit, axis, mode)`, `dragTo(position, hit?)`, `endDrag()`, `cancelDrag()`: the range
  between the two ends, whole items at the ends when it began on one and the pointer is over one, by the drag's mode:
  `replace` (the selection becomes it), `add`, `subtract` (to or from the selection the drag began with), `extend` (the
  anchor's range reaches to it; with no anchor on its system, a replace); a cancel gives back the selection, the anchor
  and what an extend keeps. The range between two ends is taken in display offsets on the hit's axis, piece by piece in
  its direction: across a collapsed break, the positions shown between the ends (`selection.acrossBreaks`: `shown`;
  `span` adds the hidden stretches, the hull). Why methods: a binding, a key and a test act alike, and a drag step is a
  call, not a command.
- **The anchor** is the range of the last element or cell set, toggled, added or dragged from, on its system (an item's
  extent there, else its cell); `set(value, anchor?)` replaces it; `clear()` drops it. With it the model keeps what an
  extend keeps: the selection less the anchor's range when the anchor was set (none after a click or a plain drag; all
  before it after an add; all but what was hit after a toggle). Why kept, not found at the extend: an extend that grew
  over another range and then shrank must give that range back.
- **An element's extent** is the hit's, kept while it is selected; one selected without a hit (by an address, then
  toggled) is resolved by `extentOf`, which the viewer gives (`locationOfResult` over what is shown).
- **The cursor** (`cursor`): a position on a section's system, `{ section, system, position }`, drawn on every section
  of that system; `setCursor`. A click with a position moves it there.
- **Its nodes are held**: the viewer holds the dynamic nodes the selection names (its elements' owners, its systems)
  through `StoreModel.holdSession('selection', …)` while it is selected (`selection.holdNodes`). Why: a selected
  transcript must stay when its region scrolls away.
- `HighlightModel`: named channels (`set`, `clear`, `channel(name)` a keyed state, `names` in the order first set,
  `version` moving on every change). Drawn as tints (Emphasis).

### Emphasis

`models/projection.ts`, `viewer.ts`: what the emphasis passes draw, one reused record set before each frame.

- **The selection** (its `shown`) and **its projection** onto the linked sections (`Projector.of`: the core's
  `projectOntoView`, only the systems it has no part on), **the cursor**, **the highlight channels** in the order first
  set (each its Location with its projection, `own` kept; its color `highlight.colors[name]`, else `highlight.color`, a
  hex literal or `token:<name>`), and **the hover** with its projection (`HoverModel.projectedAll`).
- **Projections are cached** per Location by identity while the view and the model hold (`Projector`, at most 64); the
  highlights' list is made again only when a channel, the tokens, the view or the model changed: a redraw projects
  nothing and allocates nothing in the renderer (held by the allocation tier).
- **The cursor shows once the keyboard moved it** (focus-visible: `Cursor.keyboard`, set by the `cursor` behavior's
  moves; a click places it unseen). The hover's own section is passed by address (`projected.hoverSection`), the
  selection's element keys as `elements` (a zero-length item's hairline).
- `renderSvg` passes the highlights too; the export projects the selection itself. **The hover and the cursor are
  transient**: never in an export, which shows the selection and the highlights alone.

### AddressModel

`models/address.ts` (`viewer.addresses`): addresses in the spec shown and the store's nodes now (the core's Addressing
over `overlayQueryModel` of the static model and `StoreModel.currentOverlay`, the last frame's builder's overlay without
syncing sources), with the `address.*` settings as bounds.

- **Actions**: `focus` (each leading, movable section whose system the address has positions on, or reaches by a static
  path, shows their hull in display offsets, hidden stretches taking no width; one move a scroll group; followers
  follow; what lies wholly in a hidden stretch moves the section to its break, centered, the zoom kept, and the outcome
  says so and names `section.collapse`: `movedToBreak`), `reveal` (nothing when shown; a pan to its center when it fits;
  else a focus), `select`, `highlightAddress` (a channel, `links.highlightChannel` by default), `hoverAddress` (the
  links' own channel, `links.hoverChannel`, EC-182: the pointer's hover is left alone), `unhover`,
  `clearHighlight(channel?)` (one channel or all). Each answers an outcome: done, the resolution's status, a message,
  its issues.
- **A focus that waits** (an address of a node not held, data not loaded) is kept and runs when the address resolves:
  retried on the refs' state, the store's version and its node version; another focus replaces it. A focus whose section
  is not placed yet (no width, no length) is kept as the view's intent and reported not done
  (`Kept until the view is placed`), never `Focused`; the address is the focus (`focusAddress()`) once the kept move
  lands. Before any spec is set (`SpecModel.loaded` false: the empty spec a viewer starts with) every address waits, so
  a link followed before the document arrives is kept and runs when it is set (retried with the refs' state). The record
  waits on the sections the move moved alone, and ends at another jump (a focus, back, forward) or when a new spec
  leaves the address naming nothing (`seq.address.missing`).
- **The session follows the spec**: `followMoves(changes.moved)`, called before a patch's root is set (the spec model's
  `beforeCommit`, so every commit's moves reach it), rewrites the session's addresses (selection, highlights, the
  waiting and the last focus) through moves and renames as the spec's refs are; on each new root they are resolved
  again, and one that names nothing is dropped with `seq.address.missing` in `issues`.
- **Refs**: `refs`, the spec's refs with what each names now (status, label, color, targets), a new list only when an
  entry changed; `addRef` adds one as one patch (`addRefPatch`, made against the revision just read) through
  `SpecModel.patch`: one step of the history (a promise while the history runs it); a ref that would dangle is refused;
  `bookmarkSelection({ id?, label? })` adds one for the selection (its id and label made from the root the edit runs on:
  undefined with nothing selected).
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

- **Session holds** (`holdSession(name, model, addresses)`): the dynamic nodes a session state names, with what they
  name (`Overlay.holds`), held in the store under the name, those it no longer names let go; retired with the model.
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
  (`@molstar/seqx-viewer/default-spec`) lists the `seq` runtime, the `hover`, `navigate`, `select`, `cursor`,
  `addressing` and `editing` behaviors, the default ports and the UI (`canvas`, and `history`: a mounted viewer with no
  `history` given gets a `uiHistory()` stack of its own); an application spreads it. A headless spread (`ui: undefined`)
  has no history, so its edits apply at once (`refs.add` and the editing commands answer at once); a mounted one runs
  them in turn through its stack (`execute` waits for their message, `runCommand` gives `''`).
- `createSeqViewer(spec, target?)` makes, in this order: the core's context (registry, configuration, ports), the models
  (spec, view, hover, store, selection, highlight, addresses), the history (`spec.history`, else `uiHistory()` when
  mounted with `ui.history`; bound to the spec), the binding table, the behaviors (one of each name; a behavior's
  `commands` join the viewer's), then the UI into `target` when there is one: the toolbar the spec lists (`ui.toolbar`,
  a function the host gives: `mountToolbar` of `@molstar/seqx-viewer/mount-toolbar`), then the canvas when the spec has
  `ui.canvas`. Without a target it is headless: models, commands, frames, `hit`, SVG. Why the order: behaviors read the
  models and the bindings, the UI the behaviors.
- The document is not part of the spec: `load(root, { view? })` sets the root and opens the sources; `settle()` waits
  for the picture (sources open, wanted data arrived, the view at rest).
- **The frame** (`draw`): the renderer lays the view out when anything it reads changed, else redraws the emphasis
  passes alone (a hover); the layout's wants go to the store through the viewer's want set; the backend draws. Each
  leading section's window goes in the frame (and `renderSvg`) as its own viewport (`viewports` by address, the first as
  `viewport`): sections of different axes share no offsets. A frame is asked for by `schedule`, one frame callback
  (`@molstar/ui/dom`, on `DomHost`) however many changes asked; the wheel's own frame is asked for first, so its deltas
  apply in the frame that draws them. Changes that ask: the model, the resolved view, a viewport, the hover, the
  selection, the cursor, a highlight channel, the store's version and settledness, a resize. A frame whose base passes
  were drawn again (`Frame.base`: the selection's tint under the data) redraws the base canvas.
- **Settled**: the view counts as settled `canvas.settleMs` after its last move and while the store has nothing pending
  (open domains computed, lanes packed).
- **Size**: the width is the target's; the data box is the width less the layout's inset, label column and gap
  (`dataX`), set as every section's width. The height is the layout's content (`canvas.height: 'content'`: laid out at
  the last content height, again when it changed, at most `canvas.maxHeight`) or fixed. The renderer takes one window
  for all sections: the first section that leads gives it.
- **Motion** (`view.jumps`): a jump (a focus, a fit, a reveal, back, forward, a step zoom) is a plan (the core's `plan`)
  from every leading section's viewport drawn last (`drawnViewports`: its window, a text's anchor and letter size) to
  the model's, requested of the animation port over `view.jumpMs` (240): each section its own view channel, eased in and
  out, its width and letter size geometric and its centre along, a ring's window round its shortest arc (`cyclic` totals
  of the sections that wrap). The view is unsettled while it runs; the model is there at once, and hits read the frame
  drawn. Pointer input (a press, a click, a double-click, a wheel) first stops a motion, a jump's or a story's, or a
  host's frame, where it is drawn (`stopMotion`: each leading section takes the viewport drawn, `setViewport`, the
  state's spec is shown, and a frame is drawn), so what it hits, the pointer's offset and a drag read one picture; a key
  and a side button act on the model (a repeated `-` steps from where the last one went). A move by hand ends any motion
  and clears a host's frame. A jump is a cut under reduced motion (`reducedMotion`, the mounted UI's media query), with
  `view.motion: 'cut'` or `jumpMs` 0; `story.reducedMotion` (`always`, `never`) is the stories' only. `drawnViewport` is
  the first section's window the last frame drew; `viewport` and `renderSvg` read the model during a jump (the
  destination), a story's frame as shown. Why the drawn window, not the model: the model's state is the destination
  (tests, links, the history read it), and only the eye needs the way there.
- **Stories** (`playStep(id | undefined, { run })`): entering a step (or leaving one) is a transition: the step's state
  is shown at once in the model (the spec's step, its view's sections placed from positions by `setViewport`, its view
  shown), and a plan from what was drawn to it, timed by the step (`duration`, `easing`; else `story.transition`), runs
  on the animation port (`run`, the default) or is shown at its start for the application to drive
  (`animation.progress`, `animation.clear`). Each frame is shown in the frame slot: the frame's spec
  (`SpecModel.setFrame`, when the plan changes the spec) and each leading section's viewport (its `top` too, and a
  frame's `bands`; the model's viewport has no bands yet: the rows layer's). A token color blends through the viewer's
  token table (`resolveColor`). A section entering or a view switched during a transition is placed as the frame's spec
  gives it (the view switched at the start). A cut under `story.reducedMotion` (`platform`: the platform's reduced
  motion), its plan's duration 0. A step entered elsewhere (the `step.enter` and `step.leave` commands, a host's
  `spec.enterStep`) is the same transition, from the spec and viewports drawn last: a snapshot's view is a jump that
  animates.
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

### Animation

`ports/animation.ts`: Seq* owns no clock; the application does (the history port's twin). `SeqViewerSpec.animation` is a
factory given the participant, called once per viewer; none: the page's frames.

- **The participant** (the application drives Seq*, `viewer.animation`): `progress(plan, t)` (the plan sampled at `t`
  into its frame) or `show(plan, frame)`, `clear()` back to the state, `draw()` now (a headless clock), `drawn` (moves
  per frame drawn while a frame is shown), `durationOf(plan)`.
- **The animation** (Seq* asks the application): `request(plan, { duration, done })`, cancelled by its handle; the
  viewer's own motions (a jump, a story step) go through it, so a host that animates a structure's camera runs both on
  one clock, or cuts both. Ownership as the history's (`disposeWithViewer`, `dispose`).
- **Adapters** (`animation.ts`, the entry `@molstar/seqx-viewer/animation`), in the composition roots:
  `frameAnimation()` (the default spec's: `@molstar/ui/dom`'s frame callbacks, timed by `DomHost.now` from the request,
  its last frame at `t` 1); `headlessAnimation()` (`{ port, clock }`: a fake clock that moves only by `tick(ms)`, each
  tick showing every running plan at its time and drawing it; `runAll(stepMs)`), what a headless twin and an export run
  on.

### History

`ports/history.ts`: Seq* owns no undo stack; the host does. `SeqViewerSpec.history` is a factory given the participant
(`SpecModel`), called once per viewer; none: edits apply at once and nothing undoes them.

- **The participant** (the host drives Seq*): `checkpoint()` a `SpecCommit`, `restore(c)`, `changed`.
- **The history** (Seq* asks the host): `run(step, work, { changed, failed, within })` an edit as one step, or part of
  the host's run `within` (`EditStep`: its label, its kind, always `patch`, its origin), `beginTurn(label, origin)` an
  agent's turn, `isActive(context)`, `outside(reason)` a change outside a run (whether it invalidates the earlier steps
  is the host's rule), `root(revision)` a root the host's stack holds, and its ownership (`disposeWithViewer`,
  `dispose`).
- **Only spec patches enter it**: session, navigation and configuration never do (Commands).
- **Adapters** (`ui-history.ts`, the entry `@molstar/seqx-viewer/ui-history`, on `@molstar/ui/history`), in the
  composition roots only; the ports and the models never import it. `uiHistory(options)`: a stack of the spec alone
  (`UiHistory.stack`), a mounted viewer's; it invalidates on an outside change. `hostHistory(stack, participants, key)`:
  a host's stack over `HistoryParticipants`, one `HistoryPort` of several participants by key (Seq*'s spec beside a
  structure's state), the spec joining while the viewer lives, under a key of its own a viewer unless `key` is given.
  Both resolve `root(revision)` from the commits their stack checkpointed.

## Canvas2D backend

`canvas/`: `CanvasBackend` draws the core's draw lists on three stacked canvases at a device pixel ratio, the order of
the passes.

- The canvases stack as the passes do, bottom to top: the base (background, data), the overlay (window, highlight,
  selection, hover, cursor, which takes no events) and the labels, so the picture is the SVG's (a test compares the two
  orders). A hover redraws the overlay only: it is cleared and its passes drawn; an overlay with nothing to clear and
  nothing to draw is not touched. A rebuilt layout draws the base and the labels (a canvas cleared first) and then the
  overlay.
- The scale a canvas is drawn at is its backing size over its CSS size (the backing size is rounded), not the ratio
  asked for.
- A placed glyph run draws each letter at its place, turned and scaled about it, through the context's transform (the
  base transform composed with the letter's, then set back).
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

- **Round a ring** (a section whose window `wraps`), the placement read by angle: a zoom keeps the unit under the
  pointer's angle (the turn nearest the window, `RingsNav.zoom`); a pan drag (Space+drag) moves the window by the units
  the pointer swept round the ring, the shorter way; a drag-select takes the shorter arc between its ends (over the
  origin, two runs); a double-click focuses an item by the arc its pieces cover (the circle less their widest gap; its
  whole extent from `selection.extentOfHit`, a layer's hit being the piece hit). Only these read two ends as an arc: a
  range set (`setRange`, a grip, a step back) is kept as given, wider than half the ring too.
- **Pointer**: a move drags when a drag is begun (`down` on a target whose `drag` binding exists) and always hovers
  through the hover behavior. `leave` clears what this pointer owns. A drag bound to a drag-select (`drag-select`,
  `drag-add`, `drag-subtract`, `drag-extend`) begins on what was hit (nothing selectable: no drag) and selects once the
  pointer moved 3 px (either way); its range reaches the unit the placement has under the pointer (`offsetAt`: in a
  text's reading order across lines, the nearest between lines, a line's last past its end), past the data too. A
  drag-select held past a scrolling box's top or bottom (a text that overflows) scrolls it toward the pointer, 8 px a
  second for each px past (to 1,200), timed by the clock, never by frames, the selection extending to what comes under
  the held pointer; it stops when the pointer comes back, the content ends or the drag ends. A window's part (the hit's
  `target: 'window'` with its `part`) is dragged by the `window` behavior before the bindings (Space+drag still pans).
- **Click**: `click(x, y, e)` through the bindings (the `select` behavior); the click the platform sends after a drag
  that moved, and the second of a double-click (`detail` 2), are none.
- **A resting pointer is hit again** once a frame in which the layout was rebuilt (EC-143, `pointer.rehit`): the hover
  follows the data under the pointer while the view moves, with no pointer event. Why once a frame: several moves of the
  view in one frame change what is under the pointer once. Not while dragging, nor after a leave.
- **Wheel**: the `wheel` gesture with the event's modifiers (pixels; lines and pages scaled). **Pinch**: a trackpad's
  pinch arrives as a wheel with Ctrl held and small pixel deltas (`pointer.pinchFromCtrlWheel`, at most
  `pointer.pinchMaxDelta`): it is the `pinch` gesture with no modifiers (the Ctrl is the gesture's), zoom factor
  `exp(−deltaY · pointer.pinchSpeed)`. On macOS every Ctrl+wheel is a pinch (Ctrl is no zoom modifier there); elsewhere
  a larger delta is a mouse wheel with Ctrl, which is `mod`: zoom. Line deltas are `pointer.lineHeight` px, a page the
  viewport's own width or height. A zoom is about the unit the placement has under the pointer (over the label column
  the section's first visible). A section that scrolls inside its box (`Placement.clip`, a text) takes a wheel without
  Cmd/Ctrl only while its content overflows the box, scrolling it through its layout's `pan`, and gives it back to the
  page at its top and bottom (scroll chaining). True when a binding took it, and the page does not scroll.
- **Double click**: the bound command with the section, what to focus (the element's hull on the section, its extent
  carried along the layer's path, else the position's cell, in offsets; a Shift+click's extend before it is given back,
  as a clear one is), the unit the placement has under the pointer (`at`) and the double-click's zoom factor (`step`,
  `navigate.doubleClickZoomStep`).
- **Side buttons**: a mouse's back and forward buttons are the `back-button` and `forward-button` gestures
  (`sideButton`), taken only when the command can act (a step to go): then the mounted UI keeps the page from going back
  on their press and release; else they are the page's.
- **Keys**: the commands' shortcuts, run on the section under the pointer (else the first): the first whose shortcut it
  is (a command `view.keys` lists answers its chords there, a printable key whatever Shift it needs) and whose `when`
  allows it (the history's keys, Alt+Left and Alt+Right, are the page's when there is no step to go), a preview's
  command before the rest (Esc cancels a gesture first, a drag-select of any mode giving back the selection it began
  with, then clears the selection); none: the key is the page's. A Space tap (`spaceUp` with no drag begun while it was
  held) is the Space key's command (the cursor's toggle); a Space that panned is no tap. The `@molstar/ui` keymap
  replaces this when the commands join a command manager.
- **Mounted**: a focusable container with the three canvases, listeners on the container only (`click` among them; a
  drag ends on `pointerup`, `pointercancel` and `lostpointercapture`; Space, held for a drag that pans on data, does not
  scroll the page while the canvas has focus), a `ResizeObserver` for the width, a media query for the device pixel
  ratio now (made again after each change, a frame asked for), a web font arriving measures text again; the function
  `mountCanvas` returns takes away every listener, the observer and the nodes.
- **Toolbar** (`toolbar.ts`, `toolbarButtons(viewer, commands)`, DOM-free; `mount-toolbar.ts` mounts `@molstar/ui`'s
  toolbar over it, above the canvas): a part the host places (`ui.toolbar: mountToolbar`; the default spec lists none:
  `@molstar/ui`'s toolbar with its icons and menus is about 85 kB, the default viewer 508 kB without it), one button a
  command, `view.back` and `view.forward` known (`DefaultToolbar`). Its title is the label and the first chord of the
  command's `view.keys` as the platform writes it (`Back (Alt+←)`; no key when the application removed them), also its
  accessible name; it is disabled, with the reason as its tooltip, while its side of `view.history` is empty (the state
  is reactive: the bar follows jumps, back and forward, and edits that clear the history). A click runs the command on
  the first leading section. Why a model apart: the twin of the mounted bar reads the same states in Node.

- **The placement is the input's inverse** (`offsetAt(section, x, y)`, `viewer.placementOf(section)`): the unit under a
  point, a section's targets (a window with its part, a text gutter's `row`), never the input's own arithmetic.

## Behaviors

Small policies over the models (`behaviors/`), each reading the bindings. `rows` (an alignment's rows: Commands, the
rows' commands) reads the last frame's rows layers through the behavior context's `layout()`.

- **A section of blocks on several systems** (its axis on no system): a hit names its `block`, and its positions are
  read on that block's axis from its `offset`: a drag-select stays in the block it began on (clamped at its ends), a
  double-click and a focus move to the block's offsets, the keyboard cursor does not step along the section's units, a
  bookmark of the view names none. Why: a position on the section's units is no system's (core `_spec` "Resolved view",
  Blocks).

- `hover`: a move over a hit whose hover gesture is bound to `hover` sets the hover for its owner token; a move over
  nothing or a leave clears what that owner holds only. A host's other hover sources use their own tokens.
- `navigate`: a drag (`pan`), a wheel (`pan` by the larger delta, so Shift with a vertical wheel pans sideways; `zoom`
  about the pointer) and a pinch (`zoom`) through the bindings; wheel deltas within one frame are applied once at the
  frame (`navigate.coalesceWheel`, through a frame callback that `dispose` cancels); a section gone by the frame is
  dropped; `dispose` drops what is pending and later frames do nothing; a drag passes both deltas (a text scrolls its
  lines by Space+drag). The view's moves are commands (`view.fit`, `view.focus`, `view.back`, `view.forward`,
  `view.zoom-in`, `view.zoom-out`, `view.pan-left`, `view.pan-right`, `view.pan-up`, `view.pan-down`: a line in a layout
  of lines, Shift+Alt+↑↓) with their keys from `view.keys` (the first the shortcut), single keys scoped to the viewer
  (`scope: 'canvas'`); `when` says why one cannot run (no section, a follower that does not move, nothing to focus, no
  step back) instead of a throw. A fit and a focus are jumps; a focus on what the view already shows focused goes back
  (`view.backOnFocused`). A step zoom (the keys, Shift+double-click) is about the context's offset (the pointer's), else
  the keyboard cursor's cell on that section when in view, else the centre, by the context's `step` or
  `navigate.wheelZoomStep`.

- `select`: clicks and drags to the selection through the bindings (`selection.set`, `selection.extend`,
  `selection.toggle` with the context's `hit`, `selection.clear` with Esc; a drag bound to `drag-select`, `drag-add`,
  `drag-subtract` or `drag-extend` drags in that mode), so a key, a menu and a test select alike; `selection.add` and
  `selection.subtract` take a `range` (`{ system, start, end }`), an `address` (its positions, elements by their
  extents) or the `hit`'s positions, so the keyboard and a host reach what the drags do. `when` says why one cannot run
  (nothing hit, nothing selected, nothing given). A click with a position moves the cursor there. **Copy**
  (`selection.copy-letters`, `selection.copyLetters` in the design; Cmd/Ctrl+C in `view.keys`, scoped to the canvas):
  the selection's letters as a section shows them: the cursor's section, else the first the selection projects onto with
  a `letters` layer; in its reading order (display order, a reversed section's too), complemented as drawn, only the
  positions shown (a collapsed stretch's are not), a run of consecutive units a line, a section of several letters rows
  as FASTA lines (`>label`, then the row's letters); `stored: true` copies the stored letters of the selection's span on
  that system instead. Letters not held are loaded first (three rounds at most). The text goes through the behavior's
  `write` (`SelectProvider`'s param; the platform's clipboard when there is one); the message is a count ("Copied 393
  letters"): a message is shown and announced, never a sequence. Why what is shown: what a reader selected is what they
  saw (the author's decision).
- `cursor`: the keyboard cursor's commands, on the cursor's section (else the context's, else the first; the first key
  places the cursor at the section's first position shown): `cursor.prev`, `cursor.next` (a cell in display order,
  whatever the direction, clamped to the axis, across breaks: a step over one says what it crossed, "crossed a break,
  1,240 positions hidden", in its message), `cursor.up`, `cursor.down` (a line up or down in a layout of lines, the
  layout's `step`, the same column or a shorter line's last; past the first or the last line, and on a window, the
  section above or below, the position carried along the view's link, else its first shown), `cursor.extend-prev`,
  `cursor.extend-next`, `cursor.extend-up`, `cursor.extend-down` (Shift: extends the range the anchor is in, from the
  cursor's cell when there is no anchor on its system; up and down by a line), `cursor.add-prev`, `cursor.add-next`
  (Cmd/Ctrl+Shift: a new range at the cursor's cell, the others kept, which Shift then extends), `cursor.toggle`
  (Space), `cursor.focus` (Enter: the selection, else the cursor's cell, focused as an address is,
  `AddressModel.focusLocation`: the leading sections move, so a cursor on a follower focuses its leader). The keys are
  `cursor.keys`; a cursor moved out of view brings its section to it by the layout's least move (`reveal.cursor`).
- `window` (`WindowProvider`, in the default spec; a host may leave it out): a window's drags and its slider, one for
  every window (a ring's, a navigator's), over its section's viewport: a drag on its `body` pans it, on its `start` or
  `end` resizes it about the other end, never narrower than its floor (`viewport.maxPxPerPosition` px a position, an
  inward drag on a floored window keeping its far end); the pointer comes in as the unit the placement has under it, so
  the behavior knows no layout; on a cyclic axis a drag unwraps it against its last place, a grip's from its own turn (a
  step over the origin is a short one on, not a lap back: a ring's grip dragged over its origin); Esc puts the view back
  where the press found it. As a slider: `window.pan-back`, `window.pan-on` (a `navigate.keyPanFraction` step),
  `window.grow`, `window.shrink` (by `navigate.wheelZoomStep` about its middle), the focused window's arrows and
  Shift+arrows once a mounted slider binds them; its `value` (`aria-valuetext`: "21 to 30 of 100") is said when the
  window settles (a drag's end, a step), never a frame.
- **A linked selection is revealed** (`reveal.linked`, the viewer's): a selection made elsewhere (another section, a
  link, a structure) scrolls each text section (`all`: each section; `none`) that shows none of its projection to its
  first part, the layout's least move, as a move by hand (no history step); never the section it is being made in (a
  drag's or a click's, else the cursor's).
- `addressing`: the actions on addresses as commands (`address.focus`, `address.reveal`, `address.select`,
  `address.highlight` with a `channel`, `highlight.clear` (one channel, or all: `address.highlight`'s inverse, session
  state, no default key), `address.hover`, `address.unhover`, `address.copy-link` whose message is the link, `refs.add`
  with a `ref`, `selection.bookmark`: a ref to the selection's address, labeled by it unless a `label` is given, under
  the `id` given or the next free `bookmark-<n>`, one patch and one undo step; its key is `view.keys`'
  `selection.bookmark`, `Mod+d`, scoped to the canvas, and the viewer's even with nothing selected), so a link, a key, a
  menu, an agent and a test drive the view the same way; `when` says why one cannot run (no address, no ref, nothing to
  link). `runLink` runs a link.
- `editing`: edits of the spec as commands, each one patch made against the revision read (`SpecModel.patch`, one step):
  `layer.hide`, `layer.show` (the `hidden` key set or unset), `layer.color` (the layer's theme), `track.hide` and
  `track.show` (the track's own `hidden`: its layers keep theirs, so hide then show is a round trip), `track.move`
  (among its siblings, `after`), `track.add` (a `node` inserted under the section `address`, after `after`, last without
  it; the id taken is refused), `refs.remove`, `refs.label` (a bookmark's `label` param set, or removed with `''`);
  `step.enter` and `step.leave` (navigation); `preview.cancel` (Esc). A gesture's preview is a method, never a command
  (it runs per pointer move): `preview(patch)` replaces the last, `endPreview(label)` ends it in one step,
  `cancelPreview()` in nothing. `when` says why one cannot run (no such node, no step shown, nothing to cancel).

## Bindings

`bindings.ts`: the interaction grammar as data, `{ gesture, modifiers, target, action }`.

- **Pointer gestures** (`hover`, `click`, `double-click`, `drag`, `wheel`, `pinch`, `back-button`, `forward-button`);
  keys go through `@molstar/ui`'s keymap, and a `key` binding here is reported. Modifiers are `shift`, `mod` (Cmd on
  macOS, Ctrl elsewhere, `modifiersOf`, which fills a reused record), `alt`, `space`; exact, or `'*'` for whatever is
  held, tried after the exact ones. Why: a modifier held for something else must not kill the hover.
- **Targets** are a hit's kind or `any`; a binding for the kind beats `any`'s. A collapsed break is its own target
  (`break`, the render hit's `target`): no default binding but the hover, so a double-click on it does nothing (the
  pointer aimed at the break; showing its stretch is `section.collapse`, a spec change).
- **Actions must exist**: a model method (`hover`, `pan`, `zoom`, `drag-select`, `drag-add`, `drag-subtract`,
  `drag-extend`) or a command the viewer defines (`DefaultActions`, or the `actions` given). A binding to another is
  reported (`seq.viewer.binding.action`) and left out, so the gesture does nothing.
- **Presets**: `default` (hover with any modifiers; drag on an axis, empty space or with Space pans; wheel and
  Shift+wheel pan, Cmd/Ctrl+wheel and pinch zoom; double-click focuses an element (on the one focused: back), fits on
  empty space; Shift+double-click zooms out about the pointer; the side buttons go back and forward; the selection's
  part, `SelectionBindings`: click sets, Shift+click extends the anchor's range, Cmd/Ctrl+click toggles, a click on
  empty space clears; `SelectionDrags` on data: a drag selects a range, Cmd/Ctrl+drag adds it, Alt/Option+drag removes
  it, Shift+drag extends the anchor's range; Alt+wheel is unbound on purpose, left to the platform and the host) and
  `genome-browser` (a drag on data pans instead; the selection's clicks kept, its drags dropped).
- **A tree's nodes** (`CladeBindings`, in `default`): a click on a `tree-node` toggles its clade's fold
  (`clades.toggle`: decided, a node is chrome, as a group's chevron; a preset may bind it to select), Cmd/Ctrl+click
  toggles its rows in the selection; a `summary-row` is an element (a click selects the clade's rows, Shift+click on a
  row extends in display order), its double-click opens the clade (`clades.toggle`, again folds it).
- `bindingTable(preset, overrides, { actions, kinds, keys })`: each override layer replaces bindings of the same
  gesture, modifiers and target (`action: null` removes one); two of one key within a layer conflict. Resolving is map
  lookups, no allocation. **The keys check** (`keyConflicts`, the viewer passes `view.keys`): two commands whose chords
  one key event matches conflict (`seq.viewer.binding.conflict`): a printable chord without Shift matches whatever Shift
  types it, so a `Shift+=` chord is the `+` of `view.zoom-in`'s (Shift as the common layouts type it).

## Commands

`commands.ts`: `@molstar/ui`'s `CommandDefinition` has no change kind, so the viewer keeps a table by command id
(`CommandChangeKinds`, `changeKindOf(id, extra)`): `patch` (the spec, one step of the history), `preview` (a gesture's
patch: ends in one patch or nothing, never a step itself), `session`, `navigate` (what is shown without changing the
spec: the viewports, a story step) or `none` (it only reads). Menus show what Undo will undo; an agent tells an edit
from a look. A focus or a reveal navigates, as `view.back` and `view.forward` do; a selection, a highlight or a link's
hover (and clearing a highlight) is session state; `refs.add`, `selection.bookmark` and the editing commands are
patches; entering and leaving a step navigates; `preview.cancel` is a preview's; copying a link changes nothing. The
selection's gestures and the cursor's moves are session state; the cursor's focus (Enter) navigates.

- **Tested as the table says** (`_test/changes.test.ts`): every command of the default viewer, run with the default
  spec's history, changes exactly what its kind names (a patch: the revision and one history step; a preview: the
  preview; session; navigate; `none`: nothing), and undo restores the spec alone: the session, the viewports and the
  step shown stay. Why a test of every command: a kind declared and not held misleads the menu and the agent alike.
- **Running** (`viewer.runCommand(id, ctx)`, `viewer.execute(id, ctx)`): a command that edits through the host's history
  runs in turn; `execute` gives its message when it has run, `runCommand` only a message given at once.
- **The axis's commands** (the `editing` behavior's, patches: one undo step each): `section.reverse` flips what the
  section shows (`genomic` and `reversed` swap; `follow` or `along` becomes the explicit opposite of what it resolved
  to), `section.collapse` switches `pieces` between `whole` (unset) and the section's `collapsed` (`coding` when
  absent), with no memory of an earlier form. On the section named (`address`), else the keyboard cursor's, else the
  first leading one; their keys are `view.keys`' (Shift+R, Shift+C). Why patches: direction and pieces are what the
  figure is, so undo (not back) returns them and the export says the same.

- **The rows' commands** (the `rows` behavior's, `behaviors/rows.ts`, on a track's rows layer, the first when none is
  named): folds are session state (the band's `folds`: no undo step): `clades.toggle` (a hand's move, no history step),
  `clades.collapseAll` (every clade but the root), `clades.expandAll`, `clades.collapseTo { depth | height }` and
  `clades.reveal { rows }` (opens the clades holding them) are jumps, so back restores the folds and the rows' window
  together; `clades.keep` writes the folds shown into the layer's `collapsed` (one patch); `rows.pin` (toggles: a pinned
  row is unpinned), `rows.hide`, `rows.show` (the hidden-rows chip's) and `rows.order` are patches of the layer's
  params, read from the spec as it is now; `rows.zoom-in` and `rows.zoom-out` (`]`, `[` in `view.keys`) double and halve
  `rowScale` (1/64 to 40 px a row, fixed: 1/64 shows 10,000 rows in about 160 px, 40 is twice a letter row) about the
  band's middle row. A fold keeps the row at the band's top there (a fold above the window moves it up, never off). **An
  edit of `collapsed`** (a patch, an undo) resets the session's folds to it; other edits keep them. The wheel's vertical
  delta over a rows layer's own band (its cells, its header) scrolls its rows (`bandTop`, px over the row height,
  clamped), its horizontal pans the columns; over another layer's `row` hit it pans as any. Clades fold only in tree
  order (`when` says). `rows.filter` waits for a `filter` param of the rows layer.

## Legends

`legends.ts`: `LegendModel`, derived when the view settles (never per frame: `viewLegends` over the settled layout,
merged), the entry under the hover (`updateHover`: a mark's value), and the reader's substitution (`cvdSafe`, from
`theme.cvdSafe`, part of the frame's key). `legendRows(viewer, { of, maxItems })` is what a view of them shows: a row a
legend, its entries with accessible names (`<title>: <value>`), the entries past `maxItems` counted. Why derived on
settle: a legend that changes during a pan reads as noise, and its cost stays off the frame.

- **The legend behavior's commands** (`behaviors/theming.ts`, patches unless said): `layer.channel` (a channel's value
  or form set or unset, on one layer or several in one patch), `scale.set` (a scale's params; `null` reopens a domain's
  end), `theme.value-color` (one value's color in a field form's palette, the listed values keeping theirs; refused on a
  scheme, rules or any other form: never replaced), `layer.filter` (a value's items hidden or shown again: one patch of
  the layer's `data`, a `not` of the values hidden; no value: cleared, as a click on the track's "n filtered" chip),
  `legend.select` (session: the items whose field has the value; Cmd/Ctrl toggles), `legend.hover` (session: a highlight
  named `legend`), `theme.cvd-safe` (session). `fittingSchemes(layer)` is what a menu offers (the schemes whose subject
  and measures fit the layer); `preview(layer, channel, value)` shows one hovered (a preview, ended by `layer.channel`
  or Esc).

## Configuration

`config.ts`: `DefaultViewerConfig`, one object of the shape the models read (`viewerConfig(over)` for one group's items
over the defaults; `view.keys` merges by command id). Each item is a choice with its recommended default; all of it is
read when the viewer is made. As items in a spec: `viewerItem('group.item')` (typed by key), read by
`viewerConfigOf(pairs)`; the core's items (`SeqCoreConfig`) are the same pairs. `ParamDefinition`s come with the
settings form.

| Item                              | Default                                                                                                                                                                                                                                                                                                                                                                             | Why the default; the alternatives                                                                                                                                                                                                                                            |
| --------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `foundations.preset`              | `comfortable`                                                                                                                                                                                                                                                                                                                                                                       | the UX foundations' spacing; `compact`                                                                                                                                                                                                                                       |
| `density.letterMinPx`             | 7                                                                                                                                                                                                                                                                                                                                                                                   | the least CSS px a position needs for its letter (a number: no reason recorded)                                                                                                                                                                                              |
| `density.dotMinPx`                | 2                                                                                                                                                                                                                                                                                                                                                                                   | the least px a position needs for a dot (no reason recorded)                                                                                                                                                                                                                 |
| `bindings.preset`                 | `default`                                                                                                                                                                                                                                                                                                                                                                           | the interaction grammar (Bindings); `genome-browser`: a drag on data pans                                                                                                                                                                                                    |
| `store.cacheBytes`                | 64 MB                                                                                                                                                                                                                                                                                                                                                                               | the loaded bytes the store keeps, least recently used evicted, never what is in view (the size: no reason recorded)                                                                                                                                                          |
| `store.inFlight`                  | 6                                                                                                                                                                                                                                                                                                                                                                                   | requests in flight at most (no reason recorded)                                                                                                                                                                                                                              |
| `store.tileBins`                  | 1,024                                                                                                                                                                                                                                                                                                                                                                               | bins per tile (no reason recorded)                                                                                                                                                                                                                                           |
| `hover.projectLinked`             | true                                                                                                                                                                                                                                                                                                                                                                                | the hover shows on every linked section; false: only sections of the hovered system                                                                                                                                                                                          |
| `hover.closeCodons`               | true                                                                                                                                                                                                                                                                                                                                                                                | a nucleotide lights its whole codon and its residue (EC-30, EC-144); false: the position alone                                                                                                                                                                               |
| `selection.holdNodes`             | true                                                                                                                                                                                                                                                                                                                                                                                | a selection's dynamic nodes stay while selected, so a selected transcript survives a scroll; false: only while shown                                                                                                                                                         |
| `selection.extendScope`           | `range`                                                                                                                                                                                                                                                                                                                                                                             | an extend changes the range the anchor is in, the others kept; `selection`: the whole selection becomes that range                                                                                                                                                           |
| `selection.acrossBreaks`          | `shown`                                                                                                                                                                                                                                                                                                                                                                             | what a drag or an extend across a collapsed break selects: the positions shown (what is seen); `span`: the hidden stretches too (a gene's whole span to export)                                                                                                              |
| `highlight.color`                 | `token:cat-5`                                                                                                                                                                                                                                                                                                                                                                       | a channel's color when `highlight.colors` names none: a hex literal or `token:<name>`                                                                                                                                                                                        |
| `highlight.colors`                | `{ 'link-hover': 'token:hover' }`                                                                                                                                                                                                                                                                                                                                                   | colors by channel: a link's hover in the hover's ink                                                                                                                                                                                                                         |
| `cursor.keys`                     | `ArrowLeft`, `ArrowRight`, `ArrowUp`, `ArrowDown`                                                                                                                                                                                                                                                                                                                                   | back, on, the section above, below; Shift with back or on extends from it                                                                                                                                                                                                    |
| `reveal.cursor`                   | true                                                                                                                                                                                                                                                                                                                                                                                | a cursor moved out of view brings its section to it (the layout's least move); false: it moves unseen                                                                                                                                                                        |
| `reveal.linked`                   | `text`                                                                                                                                                                                                                                                                                                                                                                              | a selection made elsewhere scrolls a text section showing none of it to its first part (a text hides most of itself; a window shows the projection where it is); `all`, `none`                                                                                               |
| `viewport.maxPxPerPosition`       | 32                                                                                                                                                                                                                                                                                                                                                                                  | the zoom's upper bound (no reason recorded)                                                                                                                                                                                                                                  |
| `viewport.groupEqualAxes`         | true                                                                                                                                                                                                                                                                                                                                                                                | sections on equal axes share one viewport (EC-113); false: one each                                                                                                                                                                                                          |
| `viewport.onResize`               | `keep-zoom`                                                                                                                                                                                                                                                                                                                                                                         | a section hidden for lack of width comes back as it was (EC-112); `keep-range`                                                                                                                                                                                               |
| `follow.followerMoves`            | `lead`                                                                                                                                                                                                                                                                                                                                                                              | a move asked of a follower moves its leader through the follow, so every section stays in step (the author); `ignore`; `detach` (its own viewport until attached again)                                                                                                      |
| `follow.cyclicArc`                | `wrap`                                                                                                                                                                                                                                                                                                                                                                              | a follower's range over a cyclic origin kept as the wrap, two pieces; `hull`                                                                                                                                                                                                 |
| `view.historyLimit`               | 50                                                                                                                                                                                                                                                                                                                                                                                  | steps back the history of jumps keeps, the oldest dropped first (no reason recorded)                                                                                                                                                                                         |
| `view.backOnFocused`              | true                                                                                                                                                                                                                                                                                                                                                                                | a double-click on what is already focused goes back; false: nothing                                                                                                                                                                                                          |
| `view.jumpMs`                     | 240                                                                                                                                                                                                                                                                                                                                                                                 | a jump's animation (focus, fit, reveal, back, forward, a step zoom; the length: no reason recorded); 0: a cut                                                                                                                                                                |
| `view.motion`                     | `auto`                                                                                                                                                                                                                                                                                                                                                                              | animated unless the platform asks for reduced motion; `cut`: always a cut                                                                                                                                                                                                    |
| `story.transition`                | `{ duration: 240, easing: ease-in-out }`                                                                                                                                                                                                                                                                                                                                            | a step's transition unless the step says: the base duration of a jump (no reason recorded for its length)                                                                                                                                                                    |
| `story.linger`                    | 2,000 ms                                                                                                                                                                                                                                                                                                                                                                            | how long a played step is shown before the next (no reason recorded)                                                                                                                                                                                                         |
| `story.view`                      | `{ path: auto }`                                                                                                                                                                                                                                                                                                                                                                    | pan when the windows overlap or touch, else out, across and in; `pan`, `out-across-in`                                                                                                                                                                                       |
| `story.color`                     | `{ space: oklab }`                                                                                                                                                                                                                                                                                                                                                                  | perceptual: a blend reads evenly; `srgb`                                                                                                                                                                                                                                     |
| `story.cut`                       | `{ at: 0.5 }`                                                                                                                                                                                                                                                                                                                                                                       | a cut switches mid-window, where the change is least seen (no reason recorded)                                                                                                                                                                                               |
| `story.reducedMotion`             | `platform`                                                                                                                                                                                                                                                                                                                                                                          | follows the platform's setting (`reducedMotion`); `always`, `never`                                                                                                                                                                                                          |
| `story.stepStart`                 | `spec`                                                                                                                                                                                                                                                                                                                                                                              | a patch step applies over the spec: any step reached directly; `from` per step                                                                                                                                                                                               |
| `story.store`                     | `patch`                                                                                                                                                                                                                                                                                                                                                                             | the encoder stores a patch over the spec (it follows edits); `full`, `smaller`                                                                                                                                                                                               |
| `story.maxChain`                  | 8                                                                                                                                                                                                                                                                                                                                                                                   | deltas at most between a step and its keyframe: reaching a step costs at most that many applies                                                                                                                                                                              |
| `view.keys`                       | `view.back`: Alt+ArrowLeft, Backspace; `view.forward`: Alt+ArrowRight; `view.fit`: 0; `view.zoom-in`: = and +; `view.zoom-out`: -; `view.pan-left`, `view.pan-right`, `view.pan-up`, `view.pan-down`: Shift+Alt+Arrow; `selection.copy-letters`: Mod+c; `selection.bookmark`: Mod+d; `section.reverse`: Shift+R; `section.collapse`: Shift+C; `rows.zoom-in`: ]; `rows.zoom-out`: [ | the view's commands' keys by command id, the first the shortcut, and the one shown in the toolbar's tooltip (Cmd/Ctrl+D bookmarks the selection: the browser's own bookmark key, taken only while the canvas has focus); an override merges by id (`[]`: none)               |
| `ui.toolbar` (in `SeqViewerSpec`) | none                                                                                                                                                                                                                                                                                                                                                                                | a function mounting the bar above the canvas (`mountToolbar`: Back and Forward, their keys and disabled state); the host places it, as it would a legend: the bar is `@molstar/ui`'s and heavy                                                                               |
| `ui.legends` (in `SeqViewerSpec`) | `legendsMount`                                                                                                                                                                                                                                                                                                                                                                      | a function mounting the Legend part (`mountLegend(element, viewer, { of, maxItems })`: the legends of the view or of one track, at most 12 entries each, its rows options with accessible names, roving focus, a click selecting a value's items, a hover highlighting them) |
| `navigate.wheelZoomStep`          | 1.25                                                                                                                                                                                                                                                                                                                                                                                | the zoom of one wheel notch (100 px of delta) with the zoom modifier, and of a key zoom (no reason recorded)                                                                                                                                                                 |
| `navigate.doubleClickZoomStep`    | 2                                                                                                                                                                                                                                                                                                                                                                                   | Shift+double-click's zoom out about the pointer                                                                                                                                                                                                                              |
| `navigate.keyPanFraction`         | 0.1                                                                                                                                                                                                                                                                                                                                                                                 | the part of the viewport a key pans (no reason recorded)                                                                                                                                                                                                                     |
| `navigate.coalesceWheel`          | true                                                                                                                                                                                                                                                                                                                                                                                | wheel deltas within a frame applied once, at the frame; false: at each event                                                                                                                                                                                                 |
| `theme.source`                    | `root`                                                                                                                                                                                                                                                                                                                                                                              | a mounted viewer reads its root's CSS (`--msui-*`, `--seq-*`), redrawn when `data-msui-theme` changes; `config`: the core's `theme` item (a headless viewer always)                                                                                                          |
| `theme.cvdSafe`                   | false                                                                                                                                                                                                                                                                                                                                                                               | the reader's color-blind substitution: a scheme not safe under the three deficiencies is drawn as its safe substitute (`theme.cvd-safe` toggles it; a reader's setting, never the spec's)                                                                                    |
| `canvas.height`                   | `content`                                                                                                                                                                                                                                                                                                                                                                           | the layout's content, tracks below the page not laid out; a number of CSS px                                                                                                                                                                                                 |
| `canvas.maxHeight`                | 2,048                                                                                                                                                                                                                                                                                                                                                                               | the tallest the canvas grows with `content`, CSS px (no reason recorded)                                                                                                                                                                                                     |
| `canvas.dpr`                      | `device`                                                                                                                                                                                                                                                                                                                                                                            | the device's pixel ratio; a number for an export or a test                                                                                                                                                                                                                   |
| `canvas.settleMs`                 | 150                                                                                                                                                                                                                                                                                                                                                                                 | after the last move the view counts as settled: open domains computed, lanes packed (no reason recorded)                                                                                                                                                                     |
| `pointer.rehit`                   | true                                                                                                                                                                                                                                                                                                                                                                                | a resting pointer is hit again once a frame in which the layout was rebuilt, so the hover follows the data (EC-143); false: with the pointer only                                                                                                                            |
| `pointer.pinchFromCtrlWheel`      | true                                                                                                                                                                                                                                                                                                                                                                                | a trackpad's pinch arrives as a Ctrl+wheel with small deltas: the `pinch` gesture; false: a wheel with Ctrl                                                                                                                                                                  |
| `pointer.pinchSpeed`              | 0.01                                                                                                                                                                                                                                                                                                                                                                                | a pinch's zoom is `exp(-deltaY · pinchSpeed)`                                                                                                                                                                                                                                |
| `pointer.pinchMaxDelta`           | 50                                                                                                                                                                                                                                                                                                                                                                                  | a Ctrl+wheel with a larger delta (px) is a mouse wheel, not a pinch (not on macOS)                                                                                                                                                                                           |
| `pointer.lineHeight`              | 16                                                                                                                                                                                                                                                                                                                                                                                  | px a line of a wheel in line mode scrolls; a page is the viewport's own axis                                                                                                                                                                                                 |
| `address.limit`                   | 10,000                                                                                                                                                                                                                                                                                                                                                                              | the most nodes or elements an address names: more is an issue, never a truncation                                                                                                                                                                                            |
| `address.maxDepth`                | 64                                                                                                                                                                                                                                                                                                                                                                                  | how deep below the root resolving walks                                                                                                                                                                                                                                      |
| `address.refDepth`                | 8                                                                                                                                                                                                                                                                                                                                                                                   | how many refs a chain passes through                                                                                                                                                                                                                                         |
| `links.prefix`                    | `seq-`                                                                                                                                                                                                                                                                                                                                                                              | this viewer's keys in a link; other keys are another viewer's                                                                                                                                                                                                                |
| `links.maxLength`                 | 8,192                                                                                                                                                                                                                                                                                                                                                                               | the longest link read, in characters                                                                                                                                                                                                                                         |
| `links.maxActions`                | 16                                                                                                                                                                                                                                                                                                                                                                                  | the most actions one link runs                                                                                                                                                                                                                                               |
| `links.highlightChannel`          | `link`                                                                                                                                                                                                                                                                                                                                                                              | the channel a link without `channel` highlights on                                                                                                                                                                                                                           |
| `links.hoverChannel`              | `link-hover`                                                                                                                                                                                                                                                                                                                                                                        | a link's hover has its own channel, never the pointer's (EC-182)                                                                                                                                                                                                             |
| `copyLink.prefer`                 | `['ref', 'ids']`                                                                                                                                                                                                                                                                                                                                                                    | the forms an address is copied in: a ref naming exactly it, then ids and forms                                                                                                                                                                                               |
| `copyLink.numbering`              | `ordinal`                                                                                                                                                                                                                                                                                                                                                                           | labels are copied in it; one that does not read back falls back to it                                                                                                                                                                                                        |
| `copyLink.action`                 | `focus`                                                                                                                                                                                                                                                                                                                                                                             | the action a copied link carries; `reveal`, `select`, `highlight`                                                                                                                                                                                                            |
| `revisions.keepRecent`            | 64                                                                                                                                                                                                                                                                                                                                                                                  | the recent roots kept by revision, for a rebase of a stale patch                                                                                                                                                                                                             |
| `edit.rebaseStale`                | false                                                                                                                                                                                                                                                                                                                                                                               | a stale patch is refused (EC-81): a merge the writer did not see is a decision; true: rebased through `merge3`                                                                                                                                                               |
| `preview.dropped`                 | `issue`                                                                                                                                                                                                                                                                                                                                                                             | a preview that no longer applies is dropped with a warning in the spec's issues; `silent`                                                                                                                                                                                    |
| `bounds.maxBytes`                 | 64 MB                                                                                                                                                                                                                                                                                                                                                                               | the spec's size (Language, Bounds in core's `_spec`)                                                                                                                                                                                                                         |
| `bounds.maxDepth`                 | 256                                                                                                                                                                                                                                                                                                                                                                                 | levels of nodes and nested values (EC-76)                                                                                                                                                                                                                                    |
| `bounds.maxQueryDepth`            | 64                                                                                                                                                                                                                                                                                                                                                                                  | nested forms in an address string                                                                                                                                                                                                                                            |
| `bounds.maxColumnLength`          | 2^31 − 1                                                                                                                                                                                                                                                                                                                                                                            | a column's length                                                                                                                                                                                                                                                            |
| `bounds.maxNodeIssues`            | 100                                                                                                                                                                                                                                                                                                                                                                                 | issues kept a node, the rest counted (EC-86)                                                                                                                                                                                                                                 |
| `bounds.maxPlaces`                | 10                                                                                                                                                                                                                                                                                                                                                                                  | places kept a code, the rest counted (EC-86)                                                                                                                                                                                                                                 |
| `bounds.maxBlockValues`           | 2^24                                                                                                                                                                                                                                                                                                                                                                                | the values (rows times columns) a table of an answer may hold, held by the store                                                                                                                                                                                             |

The spec's bounds are held by the viewer's patches, previews, story steps and whole roots (`SpecModel.apply` refuses
with `'out-of-bounds'` and `boundsIssues`; `viewer.load` rejects with the issues); the table bound by its store. The
nested sources' depth is the resolve config's (`maxNestedDepth`, 4).
