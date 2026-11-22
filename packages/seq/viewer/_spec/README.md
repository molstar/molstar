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

- **A source-backed system's length** (the `lengths` state the viewer passes `SpecModel`, set from `sourceLengths` over
  `StoreModel.openedCoverage` after each sync and store change): the model is derived from what is shown and the
  lengths, resolved again with the lengths the opened sources' coverages give, so such a system is laid out, its rows
  wrapped over a cyclic origin and its grids read as over inline letters (the core's `ResolveOptions.lengths`). Why the
  model, not the view: what reads a length (an item across the origin, a codon map, a grid, a block) is the model's.
- **Issues of a pending address** (`viewer.modelIssues`, a state; `model-issues.ts`, `shownModelIssues`): the model's
  issues as the viewer says them, again when the model, the store's reports or the picture's rest move; an issue of a
  pending address is not shown. `seq.resolve.length.unknown` at an annotation or a map over a system read from a source
  is held back until the source's opening settles (opened from the model's own entry, `openedCoverage`, or reported by
  the store: failed); every other issue of the model is shown as it is, `seq.resolve.system.unnamed` too. The stories'
  status line reads it beside the export's issues. Why: the code is transient while a source opens (its length arrives,
  the model is resolved again, it goes), so said then it only flashes.

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
  is `seq.viewer.step.unknown`; a spec model with no story part (its `story` option, `StoryPart`, what the viewer gives
  when its spec lists `story`) enters none, `seq.viewer.step.none`. **While a full-state step is shown, an edit is
  refused** (`seq.viewer.step.frozen`, its hint to edit the step), a patch, a preview and a whole root alike (`apply`:
  `frozen-step`): the step does not follow the spec, so the edit would show nothing. `load` leaves the step first: a new
  document.
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
nothing and tells nothing. Beside `models/view.ts` (the moves and the publish): the scroll groups across updates
(`view-groups.ts`), the bands (`view-bands.ts`), the follows' arithmetic (`view-follow.ts`) and the history's steps
(`view-history.ts`).

- **The model it resolves over** (`model`): the spec's, or for a spec with a section of a source's complex
  (`opensComplexes`) the opened model over the store's overlay (`StoreModel.overlayOf`, the frame's one a node version;
  core's Overlay, Opened sources), keyed as the frame's overlay is: the node version and the runtimes
  (`StoreModel.runtimes`, the context's; a late extension's set them again). The view, the frame's renderer, the store's
  sync (a held unit's inline letters attached), the hover, the highlights, the describer (hover card, live region) and
  copy read it, so an update reflows the section once a version and a frame or a hover resolves nothing.
- **A viewport** is `{ start, end, top, width, pending, leader?, pieces?, scale?, y? }` (`y`: a grid's own y window):
  `[start, end)` in offsets along the axis in display order (0 at the first position shown), fractional. Why offsets: a
  section's `range`, a cyclic range over the origin, a reversed axis and collapsed pieces read like a whole system; the
  core's `positionAt`, `offsetOf`, `cellOffset` and `elementAt` convert through the pieces and the direction
  (`models/axis.ts` re-exports them: one implementation, the renderer's too). A saved session (`session()`, to come)
  stores positions, never offsets: offsets change meaning when the axis does. Navigation is in display order (a drag,
  the wheel, the pan keys, the cursor's arrows move what is on screen the way the hand or key goes, whatever the
  direction).
- **Grids**: a grid whose y moves apart from x (`movesY`: its own window, `follow_y: none`, or a section it follows)
  takes a vertical pan (`panPx`'s dy, px of its plot's height) and a zoom (about `anchorY`, an offset down y) on y too.
  A fit of a grid that follows (its moves leading) fits its leaders, x's and y's (`0` over the plot shows all of it). A
  pan stops at an axis's end, the window's width kept (a follower's pan too, when its window fits its axis): carried to
  a leader past its end, it narrowed the leader instead. Following a section, the move is carried to y's leader as a
  follower's is (`_lead`, through y's follow); its own window is kept on the section apart from x's group (`y`, absent
  when all of y shows), clamped within y, cleared by a fit, and passed to the frame (`viewports[section].y`, a
  follower's entry too); a history step (`snapshot`, `record`, `jump`) keeps it beside the groups (`ys`, by section) and
  restores it, so a y-only move is a step and back undoes it (a y that follows is kept by its leader's group). With
  `follow_y: x` y moves with x. Why: the PAE's cyclin A is a section of its own, so a scroll over the plot must move it,
  not CDK2.
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
  its break), never dropped. On a section of blocks a step keeps its ends as `{ system, position }` (core's
  `placeOfUnit`) and its signature gains the blocks' systems, directions, scales and offsets: a step over a changed
  block list (an `include`, an order, a scale) is restored by those ends (`unitOfPlace`), an end on a block no longer
  shown at the nearest shown edge; never in units. `back` returns to the last viewports left, every leading section's
  range together (followers follow; one detached since is attached again, one attached since is detached again at its
  range then; a section gone since is left out), `forward` goes again; a step to where the view is now is skipped. A new
  document (`load`) clears the history; an edit of the same spec keeps it. Bounded (`view.historyLimit`, 50; the oldest
  goes). Session state: never the spec's history (view state is not history), so undo leaves it. `jumps` moves on each
  jump, back and forward. `snapshot()` gives every group's range and band now, as a jump keeps them, and
  `record(snapshot)` pushes one taken earlier as a step (unless the view is there now): a gesture's one step, pushed
  when it ends (a navigator's far drag, of its window or its rows); `jump` pushes through it.
- **The last jump** (`lastJump`, set before `jumps` moves, so a listener of `jumps` reads it):
  `{ n, kind, location?, durationMs }`, `n` the count `jumps` moves to, `kind` `focus` (an address's or a Location's,
  the cursor's Enter, a double-click on an element or a cell), `reveal`, `fit`, `zoom` (a step zoom, a grid cell's
  double-click), `move` (a navigator's or a rows rail's move, a bulk fold; a jump that names no kind), `back` or
  `forward`; `location` the Location a focus or a reveal named (`AddressModel`'s, an address resolved; a double-click's
  element extent or position), and on a step back or forward the one the jump that left the view at that step named,
  kept with the step while no hand moved the view since (a step pushed after a hand's move names none); `durationMs` is
  `view.jumpMs`. Why: `jumps` counts but says neither where nor how long; a link to a structure viewer, a story's camera
  and an analytics host read the same record.
- **A source-backed system's length** is the model's once its source opened (SpecModel, its `lengths`): the resolved
  view's axes get their pieces then, so a section on a streamed sequence is placed and its kept intent applies.
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

### NavigatorModel

`navigator.ts` (`NavigatorModel`, `viewer.navigatorModel(id)`, `viewer.navigator(def)`): DOM-free, one a navigator (a
declared node by its id or address, or a host's definition of the same shape, never saved; its picture's own params in
`def.params`). **The contract of a navigator a host draws itself** (DOM, WebGL; public, so a rename is a recorded
break): `frame.windows` (each window's true and drawn ends, axes, positions and rows, after `draw`), `values`,
`readout`, `focused`, `cursor`, `version`, `issues`, the picture's geometry for the host's own marks
(`xRange(start, end)`, a range's x in display order whatever the direction, pieces or collapsed introns;
`positionAt(x)`; `ReadoutLinePx`, the readout's line above the picture), and the input methods `resize`, `hit`, `move`,
`down`, `up`, `click`, `doubleClick`, `wheel`, `key`, `cancel`, `leave`, `setFocus`, `focusNext`: the host maps its
events to them and draws what the frame says; every rule (geometry, floor, gestures, keys, idle, announcements) stays
the model's. Why: the model already holds every rule, tested headless; a host replacing the canvas keeps them by calling
it. Its picture is the composition's (`context.pictures`, resolved with `context.kinds`: core's "Render", Navigators,
Pictures); a host's `content` the composition lacks is the axis, said in `issues` (`seq.view.navigator.picture`, as a
declared one's is in the view's); the readout at a position is the picture's `readout` when it gives one. It holds its
picture (a core `NavigatorRenderer` over the viewer's store, theme, image cache and held domains: the view renderer's,
`viewer.renderer`), draws it over the viewer's last frame (`draw`, after each frame the viewer draws: `viewer.drawn`),
wants what it reads through a want set of its own, and holds the lit part, idle per window and axis (`windowIdle`: the
whole axis at the zoom limit, or the section cannot move; rows: every row shown), the gesture in progress, the focused
slider (`focused`) and the readout (`readout`: the window's region at rest and while it moves, else the position, the
row, a band's name under the pointer; never the viewer's hover). Its size is the host's (`resize`), each side the form's
box for the view where the host gives none. Its pointer (`move`, `down`, `up`, `click` with the platform's count,
`doubleClick`, `wheel`, `cancel`) and keys (`key`) are in its own px, carried to the bound section's display units
through the picture and run through the `navigator` behavior. A press on a window focuses its slider; a press that never
moves acts at its click (outside the window it centres there, a miniature's both axes); a miniature's body drags both
axes; a rows drag is clamped to the rows less those in view, and pushes the rows at its press when it went farther than
the rows in view (`navigator.dragStep`). Its sliders' values (`values`, `aria-valuetext` as positions on the navigator's
system and rows, `disabled` when idle) are set when the view rests (`viewer.rested`) after this part's own gesture or
key, or while the part has no focus (`setFocus`, the host's): a move made elsewhere announces nothing (the design:
noise), its value said when the focus leaves; never during a drag, never a frame. Keys while focused (one tab stop, the
windows a roving slider: PageDown and PageUp move the focus to the next and the previous window; why those: arrows, Home
and End are the slider's own, Tab the page's): arrows pan (Up and Down on rows, clamped), Shift+arrows resize about the
centre, Home and End go to the ends (jumps), Enter centres on the readout's position, `=`, `-`, `0` act on the focused
window's section; Alt+Left and Alt+Right walk the history even while the window is idle; Esc cancels a gesture. The
wheel acts as over the bound section, zooming about the window's centre, never the pointer. A navigator with a zoom
strip (core's "Render", Navigators, The zoom strip: `navigator.zoomStrip`) is one model: its height the box's and the
strip's (`zoomStrip.heightPx` more, set by the frame that brings or drops it), or a height the host fixes holding both
(the strip out of it, the context the rest: drawn again in the same call, nothing overflows); a side part (`left`,
`right`) has none; `viewer.renderSvg` passes `navigator.zoomStrip` and `navigator.direction` to the export (its
`navigator` option), so the SVG is what the parts show; its windows, grips and gestures are the zoom strip's; its
context is no target (a hit's `context`: nothing lit, the default cursor, a press no gesture, a double-click nothing),
its readout the band and position there, over the mark the region (a press there says it too, a tap's hover). **Over
blocks** (core's navigator over blocks: a genome's chromosomes): the readout names the block and the position there
(`chr17: 7,668,402`), the slider the block and the assembly (`chr17: 7,668,402–7,687,550 of GRCh38`, across blocks both
ends); the context above a zoom strip is the genome and a target: a click on another block switches to it, one on the
view's own block centres there (jumps; **Switching blocks** below); a double-click fits the block under the pointer,
from that fit the whole (`navigator.doubleClick: toggle-fit`, at the fit too, its window idle); PageDown and PageUp on
its slider are `block.next` and `block.previous` (ARIA's large step). **`navigator.switch: include`**: a click that
never moved runs `section.include` with the block under it on the bound section (on the window too; a drag keeps its
gestures; on what the section already shows alone, nothing), Home (and `showAll()`, a host's title) removes it. **The
square** (core's `SquareNavigator`: two axes over blocks × blocks): the same model draws it; a click on a cell off the
window moves each axis by the switching rule (one section on both axes: by the x block), one step; a double-click fits a
cell's pair (one section on both axes: the hull of its two blocks), again the whole; a click with a switch shows the
cell's two blocks. A drag from a cell off the window sweeps the x section (one section on both axes: both), the
rectangle its release zooms to drawn live, snapped to whole cells (the blocks under its ends, both included: what is
drawn is what the release fits; `SquareNavigator.snapSweep`, the behavior's `dragTo(at, from)`), tinted and outlined as
the window is (the same range on y for one section, else the y window as it is); Escape cancels it. Not built on the
square: its grips and drags per axis. **Switching blocks** (`blockSwitchRange`, `block-switch.ts`; the author,
2026-11-21, IGV's: a picked chromosome whole): from a view showing its block whole (within a pixel) the new block whole;
zoomed inside a block, the same width centred at the same relative place in the new block, clamped to its ends; a width
past the new block's length, the new block whole. One rule for a click on the navigator, `block.next` and
`block.previous`, and the square's click.

### HoverModel

`models/hover.ts`: the hover, a method and never a command (it runs per pointer move).

- **A hit** is `{ layer?, section?, system, position?, ySystem?, yPosition?, element?, extent?, target? }`; the section
  is read from the layer's address when absent, the target (`hitTarget`) from what was hit: an element is an item, a
  position a cell, neither empty space; a layer's own target wins (`row`: a rows header; `cell` with `ySystem` and
  `yPosition`: a rows or grid cell; `tree-node`; `summary-row`). The record keeps `ySystem` and `yPosition` (NaN when
  none) and `target` flat beside the rest: a hit of two axes allocates nothing. The frame's `emphasis.hover` carries
  them all (the rows layer lights a clade from a `tree-node` or `summary-row` target); what a linked section gets is the
  clade's rows as a Location (`projectedAll`), no target.
- `set(hit, owner)` copies the hit into one reused record (`current`) and moves `version`: no allocation, and nothing
  told for the same hit. An element is the same hover by its key, whatever extent object comes with it. One owner at a
  time: another owner's `set` takes it; `clear(owner)` by a stale owner does nothing (EC-141); an undefined owner is
  refused.
- **A Location** (`setLocation(loc, owner)`): a hover with no hit (a structure viewer's hover a link carries in): the
  record holds the Location as its extent with no layer, section or system, so it is drawn as a projected hover on every
  section it reaches (`projectedAll`), with no card, no words and no address; owned as a hit is (another owner's `set`
  takes it, a stale owner's clear does nothing, EC-141); an empty Location clears it; the same Location object for the
  same owner tells nothing. Why: without it a link falls back to a highlight channel's tint, which reads as a bookmark.
- **Projection** is lazy: `location()` (the extent, else the position, with the cell's y part when the hit has one: a
  grid cell is two positions, so linked sections light both) and `projected(section)` are computed when asked and kept
  for the version (and the view and model). The location is closed first (`hover.closeCodons`: a nucleotide gains its
  residue and lifts to its codon, every piece over an intron, EC-30, EC-144), then each part is carried to the section's
  system along the view's link (paths found under the copy rule, EC-25); `hover.projectLinked` false keeps only the
  parts on the section's own system. `projectedAll()` is the same for every section at once (the core's
  `projectOntoView`, one Location keyed by system), what the frame draws (Emphasis).

- **Molecules** (`describe-molecules.ts`): `molecule()` is the molecule the hover is on (a molecule's residue or atom, a
  unit of one, a chip): its address, name, class, `identifiers` and unit, computed when asked and kept for the version;
  undefined on a sequence. Why identifiers: a host's card draws a structure from a SMILES or a CCD code; no depiction
  port (one consumer, the host's card, already has what it needs). The hover's words (the card and the announcer) add a
  modified residue's component with its name (`Component: PSU (pseudouridine)`, the spec's `compound` node or the host's
  `componentNames`), its alternatives (`Components: H2U | U`), a molecule's residue by its code, an atom by its name and
  element with its residue (`Atom: O1G (O)`, `Residue: ATP 1`); a pair of a complex's links names each end's unit.

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
- **A cell of two axes** (a grid cell: `GridCell`, its x and y systems and positions) selects a rect
  (`(rect :a <x range> :b <y range>)`, core's Location): `selectCell(cell, mode)` sets a one-cell rect (`replace`; a
  click on exactly that rect clears it), the rect from the cell anchor (`extend`), or toggles the cell's rect (`toggle`:
  subtracted when inside the selection, else added); `beginRectDrag(cell, mode)`, `rectDragTo(cell)`,
  `endRectDrag(cancel?)` drag a rect in `replace`, `add`, `subtract` (rects subtract exactly: core's Location) or
  `extend` (from the cell anchor); a cancel gives back the selection. Why rects, not two ranges: a range of each axis is
  every row and column, not the block between them.
- **The anchor** is the range of the last element or cell set, toggled, added or dragged from, on its system (an item's
  extent there, else its cell); `set(value, anchor?)` replaces it; `clear()` drops it. With it the model keeps what an
  extend keeps: the selection less the anchor's range when the anchor was set (none after a click or a plain drag; all
  before it after an add; all but what was hit after a toggle). Why kept, not found at the extend: an extend that grew
  over another range and then shrank must give that range back.
- **An element's extent** is the hit's, kept while it is selected; one selected without a hit (by an address, then
  toggled) is resolved by `extentOf`, which the viewer gives (`locationOfResult` over what is shown).
- **The cursor** (`cursor`): a position on a section's system, `{ section, system, position }`, drawn on every section
  of that system; `setCursor`. A click with a position moves it there. On a section of two axes it has a row too:
  `ySystem`, `yPosition` (a data row of an alignment's rows, a grid's y position) and on rows its `part`, `cell` (the
  default) or `header` (the row's label; `position` keeps the column it left); drawn by that section's layer of two
  axes, its column marked on the other sections of the system (none while it stands on a header). A click on a cell of
  two axes sets the row (a summary row's cell its first row), on a row's label or a summary row's the header.
- **Its nodes are held**: the viewer holds the dynamic nodes the selection names (its elements' owners, its systems)
  through `StoreModel.holdSession('selection', …)` while it is selected (`selection.holdNodes`). Why: a selected
  transcript must stay when its region scrolls away.
- `HighlightModel`: named channels (`set`, `clear`, `channel(name)` a keyed state, `names` in the order first set, a
  channel cleared and set again keeping its place (a channel set on each hover never moves above the others: the draw
  order and the legend stay put), `version` moving on every change). Drawn as stripes where channels overlap (Emphasis).
  A channel's value has an optional `style` (`tint`, `band`, `outline`; `address.highlight`'s `style`), kept when the
  spec changes and carried with the channel onto every linked section; a changed style is a changed channel.

### Emphasis

`models/projection.ts`, `viewer-emphasis.ts`: what the emphasis passes draw, one reused record set before each frame.

- **The selection** (its `shown`) and **its projection** onto the linked sections (`Projector.of`: the core's
  `projectOntoView`, only the systems it has no part on), **the cursor**, **the highlight channels** in the order first
  set (each its Location with its projection, `own` kept, and its element keys, `EmphasisHighlight.elements`: a layer of
  keyed items draws those it holds, the pairs layer its arcs and dots; its color `highlight.colors[name]`, else
  `highlight.color`, a hex literal or `token:<name>`), and **the hover** with its projection
  (`HoverModel.projectedAll`), its `ySystem` and `yPosition` copied with the rest (a grid's crosshair reads the cell's
  row; the arcs ending at either position light).
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
syncing sources, with the labels of the numberings read from sources that the store holds: the core's
`withSourceNumberings`), with the `address.*` settings as bounds. The query model is kept while the model, the overlay
and the store's version are the same. A form naming a source's numbering whose labels are not in asks for its tiles
through the model's own wanter (after the resolve, never inside it), and so does a row form over an alignment's rows not
held yet (its `RowsSource`, one a model), so a link in a numbering no ruler shows, or to a row no layer read, waits and
then resolves.

- **Actions**: `focus` (each leading, movable section whose system the address has positions on, or reaches by a static
  path, shows their hull in display offsets, hidden stretches taking no width; on a cyclic axis over its whole system
  the shortest arc holding them, over the origin when that is shorter (EC-12); one move a scroll group; followers
  follow; what lies wholly in a hidden stretch moves the section to its break, centered, the zoom kept, and the outcome
  says so and names `section.collapse`: `movedToBreak`), `reveal` (nothing when shown; a pan to its center when it fits;
  else a focus), `select`, `highlightAddress` (a channel, `links.highlightChannel` by default), `hoverAddress` (the
  links' own channel, `links.hoverChannel`, EC-182: the pointer's hover is left alone), `unhover`,
  `clearHighlight(channel?)` (one channel or all). Each answers an outcome: done, the resolution's status, a message,
  its issues. `focusLocation(loc)` and `revealLocation(loc)` move to a Location as `focus` and `reveal` move to an
  address's (`moved`, `pending`, `none`); every focus and reveal names its Location to the view's last jump (ViewModel,
  History).
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
and each inline access plan, of a links' ends the ones built: core's `LinkEnds`) with the store's `limits` (its
`maxRows` and `maxBytes`) and attaches them under their addresses; a source whose entry is the object it was opened from
is kept (the resolver keeps entries by hash), a changed one is opened again, one gone is detached and its issue cleared.
`detachStale(model)` does that detaching alone, at once on a model change (the viewer, before its sync), so no frame
reads an edited source's old file; the viewer syncs again when another view of the same model resolves (the plans that
view built: a ranges layer's links' ends); `openedCoverage(model, address)` is a source's coverage only when it was
opened from the model's own entry (what `sourceLengths` reads). A source that cannot be opened is an issue at its
address (`issues`), and an opened source's own `issues` (the rows a file of records left out, each with its count,
places and hint) are reported there too; a failed opening, or a coverage the store refuses (its own report), is
forgotten so the next sync tries again. Retired with its owner, the store with it.

- **Host updates**: an opened source's `changes` (core's Data, Protocol) ask a reopen for the next frame (`DomHost`'s),
  so changes inside one frame reopen once; the source is opened again from its entry and reattached in place
  (`DataStore.reattach`: nodes kept by key, tiles read again); a reopen that fails is an issue at its address and the
  reading before it stays; an earlier reopen settling after a later one is dropped. A source gone or retired stops
  listening. Why: a host edits in bursts, and the sequences follow its edit without a second source.
  - **What referred to what went**: the addresses of the nodes a reopen removed, and of the children (units) a kept node
    lost, go to the store model's `removed` option; the viewer's takes them out of the selection
    (`SelectionModel.dropGone`: elements, parts and rects on or under them, the anchor and the cursor on one) and, when
    something was, one issue at the source (`seq.viewer.selection.dropped`, info). A bookmark or a ref naming one is as
    when a patch removes that node: it resolves to nothing, said by the addresses' own issue. Why: a host's edit is
    data, not an edit of the spec; the session state follows it, the spec is not patched.
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
`viewer.ts` composes and wires; each concern is a module beside it: the behaviors (`viewer-behaviors.ts`), the frame
clock (`viewer-clock.ts`), motion and the frame slot (`viewer-motion.ts`), the viewports a frame reads as the slot shows
them (`viewer-viewports.ts`), the rails and navigators (`viewer-rails.ts`), the emphasis (`viewer-emphasis.ts`), what a
selection moves (`viewer-selection.ts`) and the frame (`viewer-frame.ts`).

- **Extensions** (`SeqViewerSpec.extensions`, `SeqViewerExtension`: the core's descriptor, core's "Context", Extensions,
  and the viewer's kinds as the built-ins': `behaviors`, `commands`, `bindings`, `headerActions`), made by
  `defineViewerExtension` (`@molstar/seqx-viewer/extension`, so a viewer with none bundles none of it), which writes
  `composeViewer`: composed, its entries are the context's, its behaviors after the spec's (its commands a behavior's,
  `<id>:commands`), its bindings a layer over the spec's naming the viewer's actions or its own (`<id>:…`), its header
  entries beside the built-ins'; one with viewer parts and no `composeViewer` is refused (they would be dropped).
  `viewer.extend(ext)` registers one at run time and resolves to its undo: the context's checks (nothing changed on an
  error; a format runtime or a query entry is refused, since it changes what the open spec resolves to; so are header
  entries, settings and a behavior with a part over the canvas, `attach`, since nothing reads them late), its behaviors
  (`composed` called), commands and bindings; a failure takes away only what it made; the view is resolved again with
  the new layer kinds (`ViewModel.setKinds`), so the frame is a new renderer's and a layer drawn as its fallback is
  drawn as its type; the undo takes them away in reverse (behaviors disposed) and the fallback comes back. The same
  descriptor again is counted; one composed is the viewer's. A lazy chunk (`viewer-extend.ts`): registering late is
  rare, and a composition that never does bundles none of it.

- **Navigators** (`mount-navigator.ts`): parts the host places, like the legends. `mountNavigator(element, viewer, id)`
  draws one navigator (a declared node's id, or a host's `viewer.navigator(def)`) into an element the host owns: its own
  canvases sized by the element (its form's thickness where the element has no height), its listeners on it alone, one
  tab stop (`role="slider"`, its focused window's `aria-valuetext` and `aria-disabled` as the model says them, never a
  frame; its focus told to the model), its readout (navigators.md "The window's look", 4: the part's first line,
  `ReadoutLinePx` (16) of the container's padding above the canvases, right-aligned at its top right corner, never
  moving with the window or the pointer, so never over the control in use or the picture; a side part's picture is the
  form's box less the line, the part the view's height; one line (the system's name gives way, to an ellipsis and then
  nothing; a value wider than the part runs left; the value and `· row N` never; `readoutParts` is the model's split);
  the window's region at rest (the one last acted on; a rows-only overview the short form, `401–440 / 10,000`), the
  position and band under the pointer while it hovers, the window's region as drawn, every frame, while a window moves
  from any source until the pointer moves again; a position as `name: 1,234`), its clicks with their count; its disposer
  takes it away. `ui.navigators` (`navigatorsMount`, on in the default spec) places the view's shown navigators by
  `place` around the canvas (top before it, bottom after it, left and right beside it in a row) through the same call,
  again when they change; each side part's room is the row less the other side parts (`NavigatorModel.setRoom`: a
  vertical overview of its form's size narrows and hides by the plot it leaves, core's "Render", A narrow frame; one
  hidden so counts as hidden for the rows rail, `narrowed.ts`, so under `unlessOverview` the rail comes back), and the
  canvas keeps the layout's least width (`narrow.minWidth`, its `min-width`), the middle of the row scrolling only past
  it. What a navigator draws is its picture, a registry entry (`navigator-picture`: the default spec lists the core's
  four, a host adds its own `<id>:<name>`; core's "Render", Navigators, Pictures); a host that draws its navigator
  itself mounts none and calls its `NavigatorModel` (Models, NavigatorModel: the contract). `ui.railToggle`
  (`railToggleMount`; for a viewer embedded in a host's cell) puts a button in the canvas's corner that hides and shows
  the ruler rails (`setRailsShown`, `railsShown`: the reader's session, never the spec). The rails
  (`navigator.form: rail`, the default): every linear section, beside a declared navigator too
  (`navigator.positionsRail`: `always`; `unlessNavigator` the old rule: no shown navigator binds the section along its
  positions, the view's declared ones and a host's (`viewer.navigator(def)`) alike, a hidden one bringing the rail back;
  `none`) in the row `navigator.railPlace` says (`railSections`, `railRows`: core's "Render", Navigators); each active
  while the `navigator` behavior is listed and its window is not idle, found again after the layout placed the view (an
  idle flip drawn in the same frame, the window pass alone). Each rows layer's band has a rows rail (`rowsRailTracks`,
  by its track: beside a shown overview or miniature over its section, `navigator.rowsRail`), live while the behavior is
  listed, idle only while every row shows and would at the zoom-in limit too (the layout's: navigators' idle rule, so a
  rail zoomed out to the whole band keeps its grips). The viewer's export leaves the rails out with its navigators (its
  picture without its chrome: no slim row); with `navigators: include` they are drawn plain.
- `SeqViewerSpec`: `registry` (the core's entries: format runtimes, layer types, layouts; and the viewer's own kind,
  `header-action`, the `headers` behavior's, left out of the core's context; `viewer.registry` as listed), `behaviors`
  (a provider and its params), `bindings` (a preset, override layers, extra action ids), `config` (`[item, value]`
  pairs), `ports` (factories with ownership), `ui` (the parts; none for a headless viewer), `history` (the history
  port), `story` (the story part). `DefaultSeqViewerSpec()` (`@molstar/seqx-viewer/default-spec`) lists the `seq` and
  records runtimes, the core's layer types and layouts (`BuiltInRenderRegistry`), the header entries
  (`BuiltInHeaderEntries`), the `hover`, `navigate`, `select`, `cursor`, `addressing`, `editing`, `navigator`, `rows`,
  `hoverCard` (its port `hoverCardOverlay()`), `keysPanel`, `headers` (its menu `headerMenuOverlay()`,
  `contextClick: everywhere`) and `linkedPairs` behaviors, the story part (`StoryPart`), the default ports and the UI
  (`canvas`, and `history: uiHistory`: a mounted viewer with no `history` given gets a stack of its own from the
  factory); an application spreads it. A headless spread (`ui: undefined`) has no history, so its edits apply at once
  (`refs.add` and the editing commands answer at once); a mounted one runs them in turn through its stack (`execute`
  waits for their message, `runCommand` gives `''`).
- **Injected, opt in**: what the spec does not list is not in the viewer's bundle. `createSeqViewer` imports no layer
  type, layout, navigator, story, history or rows piece: the context is `composeSeqContext` of the spec's registry (the
  canvas, the navigation's layouts and `renderSvg` draw its layer types and layouts alone; a layer of another type or a
  section of another layout is said, `seq.render.layer.type`, `seq.render.layout`, and laid out as tracks); the
  navigators (their models, `navigatorModel` and `navigator(def)`, and their export) are the `navigators` part
  (`NavigatorsPart`: without it a declared navigator is not drawn, the ruler rail it would replace stays, and one issue
  says so, `seq.viewer.navigators.none`; `navigator(def)` throws; an export asked to include them says so,
  `seq.render.navigator.export`, and keeps the rails), their gestures the `navigator` behavior's; the rows rail's drags
  are the `rows` behavior's (`RowsBehavior.railGesture`); story steps and transitions are the `story` part (none: jumps
  and steps are cuts, `spec.enterStep` refuses, `seq.viewer.step.none`); the mounted history is `ui.history`, a factory.
  Why: one rule for every part, listed or absent, and the size a host pays follows what it lists; the size script's
  minimal entry checks it (none of the unlisted modules in its bundle). **Lazy chunks**: a listed part used rarely or
  late (the header entries' icons, their menu and DOM parts, the headers' late parts, a runtime such as mmCIF, the
  molecule letters) may load by a dynamic `import()` on first use, which every bundler splits; the size script lists
  each entry's chunks. The first frame does not wait for one: what it would draw is drawn plainly (a slot's letter for
  its icon, core's built-in states before the headers' late parts) and drawn again when it lands; `settle` waits for the
  chunks asked for.
- `createSeqViewer(spec, target?)` makes, in this order: the core's context (registry, configuration, ports), the models
  (spec, view, hover, store, selection, highlight, addresses), the history (`spec.history`, else `ui.history()` when
  mounted; bound to the spec), the binding table (the selection's and the rows' actions bound only when their behavior
  is listed), the behaviors (one of each name; a behavior's `commands` join the viewer's), then the UI into `target`
  when there is one: the toolbar the spec lists (`ui.toolbar`, a function the host gives: `mountToolbar` of
  `@molstar/seqx-viewer/mount-toolbar`), then the canvas when the spec has `ui.canvas`. Without a target it is headless:
  models, commands, frames, `hit`, SVG. Why the order: behaviors read the models and the bindings, the UI the behaviors.
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
  and clears a host's frame. Whether motion plays is one item for jumps and stories alike, `view.motion`: `platform` (a
  cut under reduced motion: `reducedMotion`, the mounted UI's media query), `cut`, `animate`; a jump is also a cut with
  `jumpMs` 0. Why one item: two that both decide whether motion plays overlap, and a reader asks it once.
  `drawnViewport` is the first section's window the last frame drew; `viewport` and `renderSvg` read the model during a
  jump (the destination), a story's frame as shown. Why the drawn window, not the model: the model's state is the
  destination (tests, links, the history read it), and only the eye needs the way there.
- **Stories** (`playStep(id | undefined, { run })`): entering a step (or leaving one) is a transition: the step's state
  is shown at once in the model (the spec's step, its view's sections placed from positions by `setViewport`, its view
  shown), and a plan from what was drawn to it, timed by the step (`duration`, `easing`; else `story.transition`), runs
  on the animation port (`run`, the default) or is shown at its start for the application to drive
  (`animation.progress`, `animation.clear`). Each frame is shown in the frame slot: the frame's spec
  (`SpecModel.setFrame`, when the plan changes the spec) and each leading section's viewport (its `top` too, and a
  frame's `bands`; the model's viewport has no bands yet: the rows layer's). A token color blends through the viewer's
  token table (`resolveColor`). A section entering or a view switched during a transition is placed as the frame's spec
  gives it (the view switched at the start). A cut under `view.motion` (as a jump), its plan's duration 0. A step
  entered elsewhere (the `step.enter` and `step.leave` commands, a host's `spec.enterStep`) is the same transition, from
  the spec and viewports drawn last: a snapshot's view is a jump that animates.
- **The picture as SVG** (`renderSvg`): the context's export of the shown spec, size and window through the viewer's
  store, at the last frame's device pixel ratio (the density steps are in device pixels, so it is the picture the canvas
  shows); the context's `renderSvg` of the same spec, size and window gives the same string.
- **Links** (`runLink(text)`): a link's actions, in order, each through the model method its command runs (the
  `addressing` behavior's); undefined when the text is no link or the behavior is not listed.
- **Retire** (`retire()`): the viewer is marked retired first (its disposer is registered last, so it runs first), so
  the teardown that follows asks for no frame; the frame callback and the settle timer are cancelled, then the mounted
  UI, the behaviors, the want set, the models and the context's owned ports go. After it no frame is pending and input
  is ignored. A spec that fails to compose (a behavior listed twice) retires what it made. `retiredSignal` (an
  `AbortSignal`) is aborted as the viewer is marked retired, before its models go: a helper over the viewer (a structure
  link) ends on it and clears what it set while the models stand. Why: `retired` is a getter and the dispose hooks are
  protected; any helper needs the end.
- **Composed** (`BehaviorHandle.composed(viewer)`): called once every behavior and the frame are made, before the UI;
  what reads the viewer's commands or registry (the header entries' checks). A throw retires what was made.
  `viewer.headers` is the `headers` behavior's geometry (undefined without it).
- **Behaviors by name**: the input finds a behavior by its provider's name (`viewer.behavior<HoverBehavior>('hover')`,
  or a provider), so it imports no provider; one the spec does not list is a `seq.viewer.behavior.missing` issue
  (`behaviorIssues`), once, and the gesture does nothing. A handle under a built-in's name is checked for the methods
  the viewer calls on it (`BehaviorShapes`) before it is handed out: a host's provider of that name with another shape
  is a `seq.viewer.behavior.shape` issue, once, and the gesture does nothing. Why: the name is the only link, so the
  shape is what the viewer can check.
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

### Hover card

`hover-card.ts` (the data and `placeCard`), `mount-hover-card.ts` (the default adapter): what the pointer rests on is
shown by the host (the author: "so apps can easily provide their own styles"). The port is the `hoverCard` behavior's
param (`port`, a `HoverCardPortFactory` given the canvas's element, none headless; the author decided, 2026-11-19), not
`spec.ports`: the core's ports are the core's.

- **Supplied, never imported**: the behavior draws nothing itself; the default spec gives `hoverCardOverlay()` as the
  param's value, as it gives `frameAnimation()` as the animation port. A host passes its own adapter there and bundles
  no default card; no `port`: the card's state moves (`shown`) and nothing is drawn. A viewer that leaves the behavior
  out bundles neither.
- **Made once a viewer**: the factory is called at the mount (with the canvas's element), or headless at the first card
  (with none); a mount after a headless card gives that port up first (disposed when the viewer owns it). The mount's
  removal hides the card and gives the port up.
- **Ownership** (`HoverCardPort.disposeWithViewer`, as the history's and the animation's): false, the host owns it (one
  adapter shared by several viewers) and the viewer never disposes it; else the viewer does when it lets it go.

- **The data** (`HoverCard`): a `title`, `rows` of `{ label, value }` (`item` on a row naming an item at the place: a
  feature under it, a group's member) and `text`, the announcer's sentence for the place (`describe.ts`, Announcements):
  one description for the live region and the card. The title is said once: no row repeats it (a row's card is titled by
  its label). A card of more than `hover.cardItems` (8) item rows goes to the port cut (`cutItems`): the first ones,
  then a row `+k more`, the rest in `more`; the behavior's `shown` keeps the whole card and the live region says every
  item. Why: a ring's position under 20 features filled the page; a host that wants the list whole draws its own card
  from `more` (the port).
- **Providers** (`card-provider.ts`, `card-providers.ts`): the card is composed by the application, as Mol\*'s plugin
  composes its labels: registry entries `{ kind: 'card-provider', name, order?, describe(at, out) }` in the viewer spec,
  beside `header-action`, run by `order`, then name (the header entries' rule; no order: after the ordered ones). `at`
  is the place (the hit, its section and layer in the last frame, the floored positions, the view resolved, the store,
  the `hover.cardItems` cap) and the viewer's lookups (`CardReads`); `out` is the parts so far (`CardParts`: `title`,
  set when empty; `rows`; `words`, the sentence's phrases, never one a row; `whole`, set by a provider that says all
  there is; `moreCount`, matches past the cap not made into rows). The viewer's own sections are `BuiltInCardProviders`
  (`seq:region` 10, `seq:pair` 20, `seq:group` 30, `seq:item` 40, `seq:columns-hidden` 50, `seq:row` 60, `seq:grid` 70,
  `seq:residue` 80, `seq:place` 90, `seq:tooltips` 95, `seq:block` 110), listed by the default spec; their rows and
  words are the card's and the live region's as they were (the hover sweep's digests hold them). A host adds an entry,
  leaves one out of the registry it spreads, or lists its own under a built-in's name (a later entry under a name
  replaces the earlier). A registry listing none keeps the built-ins. Providers change rows, never the hit or the
  emphasis: the hover contract's floor is `seq:item` and `seq:place`, and a composition without either is said
  (`seq.viewer.card.floor`, a warning). A provider that throws has what it added dropped (rows, words, title), the rest
  shown, said once a provider (`seq.viewer.card.provider`). A run of small blocks is said by `seq:block` alone. Why
  entries, not internal sections: the application composing the card is the outside user; one registry rule orders
  header entries and providers.
- **Tooltips** (`card-tooltips.ts`, `seq:tooltips`, order 95): the spec's tooltips over the core's `tooltipsAt` (core
  "Render", Tooltips). It answers only on drawn marks (an element under the pointer, or a layer that marks every
  position it shows: letters, values, grid, rows, matrix, logo), never on an empty band (EC-254); it ignores `whole` (a
  note shows over an item, a region cell or a pair said whole, EC-260). The place it reads: the pointer's position and a
  cell's y, a pair's two ends; an item by key (by data). Its rows come after the place's, share the `hover.cardItems`
  cap with the place's items (a single row on the mark always shows), the matches past it counted in
  `HoverCard.moreCount` (`cutItems`'s `+k more` says `more.length + moreCount`). A row carries its Markdown tree
  (`HoverCardRow.rich`) for an adapter that draws it; `value` (the label, a row with no label) is its plain text. The
  live region says the rows shown and the count (`…; and 47 more`), never a million texts. Addresses resolve through the
  address model, once resolved.
- **Export titles** (`cardTitles`): the viewer's `renderSvg` passes its composed card (title and rows, a line each, at
  the item's middle) for an export with `card` titles.
- **By region** (`describe-region.ts`): a pairs bin, a region grid's cell and a matrix cell are said in the core's words
  (`pairsBinHover`, `gridCellHover`, `matrixCellHover`: core "Render", Hover texts by region), their rows the card's;
  first of the sections, the whole description. Why the core's: an export's accessible layer and a host's own card say a
  bin as the viewer does.
- **A residue over its codon**: a hit on a letters layer whose data is on another system than its section's through a
  codon step (a protein's letters over its gene; an alignment row, an offset or a pairwise map keeps the place's card)
  is titled by the residue (`p 12`), with its `Index`, its letter and the codon's positions on the section's system, a
  codon split by an intron naming each part (`Codon: g 6–7, 13`). Why: the place's own rows named the nucleotide under
  the pointer and the first letters layer's letter, not what the pointer is on.
- **The place** (`HoverCardPlace`): the `anchor` (the hovered cell or item, the cursor's cell: never covered), the
  `pointer`, the `view` (the canvas) and the side and position `placeCard` gives a card of `CardEstimate`'s size; an
  adapter places its own card by its measured size with `placeCard(anchor, pointer, view, size)`: beside the pointer,
  above the anchor, else below it, else to its right or left, inside the view. Why boxes: the host's card has its own
  size; the viewer knows what it must not cover.
- **The port** (`HoverCardPort`): `show(card, place, options)` (one card: a show replaces the last), `hide()`,
  `dispose()`, `disposeWithViewer`. `options` (`HoverCardShowOptions`): `lines` (`hover.tooltipLines`, 3) and, for a
  card with a link (`cardHasLink`), `link(href)` (run through the viewer's `runLink`, the card hidden), `enter()` and
  `leave()` (the enterable card, Behaviors, hoverCard).
- **The default adapter** (`hoverCardOverlay()`): `@molstar/ui`'s primitives in the canvas's element (`msui-surface`,
  its tokens and squircle, a header, a form of rows; the rows' elements kept across shows), absolutely placed,
  `aria-hidden` (the live region says the same), never in the draw lists or an export. Compact: the title a text line,
  each row a text line (not a control row's height), the value next to its label (the label column as wide as its
  longest label, no control inset). Why: a card is read at a glance, beside the pointer. **Its width capped**
  (`CardMaxPx`, 320, or `maxWidth`), bounding the value column (`max-content minmax(0, 1fr)`): a long value wraps there,
  at `,` `;` `|` `:` first (a zero-width break after each, `breakable`) and anywhere in a token with none, to `lines`
  lines, then an ellipsis (a line clamp); the title to two. Why: a VCF's `CLNDISDB` ran past the card's edge. A row's
  Markdown (`rich`) is drawn as elements (strong, emphasis, code, a link an anchor whose click runs `link`), never HTML,
  in the value's place (the label's, a row with no value). A card with a link takes the pointer (`pointer-events: auto`,
  its `pointerenter` and `pointerleave` the options' `enter` and `leave`); one without never does, so it never steals a
  hover from the data.

### Header menu and Header parts

The `headers` behavior's params, as the hover card's port is its behavior's.

- **The menu** (`menu`, a `HeaderMenuPortFactory` given the canvas's element, none headless; made at the first menu
  headless, given up at the mount): `show(menu, place, behavior)` (`HeaderMenu`: its target and groups of
  `HeaderMenuItem`s: name, kind, label, icon, chip text, key, enabled or `why`, checked, slot, command and args, a
  submenu), `hide()`. A choice runs `behavior.runItem(item)` (`viewer.execute`), a hovered scheme `behavior.preview`, a
  menu closed without a choice `behavior.closeMenu()` (a preview ends in nothing); `loaded()`, its DOM's chunk while
  pending (`settle` waits for it). `headerMenu(viewer, target)` is the data, DOM-free. The default adapter,
  `headerMenuOverlay()` (`mount-header-menu.ts`; its DOM, `header-menu-dom.ts`, a lazy chunk loaded at the first open),
  is `@molstar/ui`'s context menu at the place. An item with a submenu opens a **flyout beside its row** (right of it,
  left of the menu where there is no room, inside the window), the parent staying open and its row `aria-expanded`:
  never in the parent's place. A click, →, Enter, Space or the pointer resting on the row (`submenuRestMs`,
  Configuration) opens it; the same rest on another row closes it, arriving in the flyout cancels that (a diagonal move
  keeps it: the tolerance is that time window, not a geometric triangle; leaving the menu stops any rest); Tab in a
  flyout closes the whole menu; the row `aria-controls` its flyout, the flyout is `aria-labelledby` its row and scrolls
  inside the window; ← or Esc in it closes it alone, the focus on its row (a preview ends: `behavior.cancelPreview`),
  Esc again the menu; a choice closes all, a press outside closes all. The flyouts are `@molstar/ui`'s
  `ContextMenuItem.items` (`openContextMenu`), each a dismiss layer of its own; its keyboard is the menu's, and closed
  the focus goes back to the canvas. The view resolving again (a source's answer, an edit, a scheme previewed in the
  menu, an undo) keeps the menu while its header is shown, and closes it when its section or track is gone or hidden.
  Why: a load landing just after the menu opened, or its own preview, closed it at once (the menu opened nothing, a
  submenu closed at its first hovered row).
- **The parts** (`parts`: `canvas`, the default, or a `HeaderPartsPortFactory` given the canvas's element and the
  behavior): the layout lays out every header as the canvas's (the same boxes); once the port's `loaded()` resolves (the
  DOM in place) the lines are `hosted` and the canvas draws nothing in their slot, "⋯" and chip boxes (until then it
  draws them; an export draws the states). Why: the parts are the DOM's alone, never a button over a drawn chip. The
  port is given `headerParts(viewer)` (`HeaderPartsState`: each shown header's geometry, whether its actions show, and
  the items `headerMenu` gives it) on change, and `focus(target)` when the keyboard cursor comes onto a header. A press
  runs `behavior.run(name, target)`, the slot click's path. `headerPartsDom()` (`mount-header-menu.ts`;
  `header-parts-dom.ts` a lazy chunk loaded at the first update: `@molstar/ui` buttons and icons): a part at each slot
  box of a header whose actions show and a chip at each state's box (its words, or a mark's squircle); one tab stop a
  header (← → between its parts, Enter or Space runs, Esc back to the canvas); Tab never captured (EC-155); its events
  stopped at its root (the canvas sees no press there).
- **Ownership** as the hover card's (`disposeWithViewer`).

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

`canvas/`: `CanvasBackend` draws the core's draw lists on four stacked canvases at a device pixel ratio, the order of
the passes.

- The canvases stack as the passes do, bottom to top: the base (background, data), the overlay (window, highlight,
  selection, hover, cursor, which takes no events), the labels and the tip (the hovered cell read whole over the labels,
  redrawn with the hover; a backend given no tip canvas draws the pass last on the overlay), so the picture is the SVG's
  (a test compares the two orders). A hover redraws the overlay only: it is cleared and its passes drawn; an overlay
  with nothing to clear and nothing to draw is not touched. A rebuilt layout draws the base and the labels (a canvas
  cleared first) and then the overlay.
- The scale a canvas is drawn at is its backing size over its CSS size (the backing size is rounded), not the ratio
  asked for.
- A placed glyph run draws each letter at its place, turned and scaled about it, through the context's transform (the
  base transform composed with the letter's, then set back).
- A rect with a radius is the core's `squirclePath` replayed; a path is replayed from its commands; a clip is `save` +
  `clip` (a squircle clip, a track's well, the core's `squirclePath` traced first), a pass left open is closed; a text
  past its max width fades through a gradient clip (never ellipsized); an image goes through one pixel buffer a size,
  copied at each draw.
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
- `canvasTextMetrics(ctx)`: the text metrics port over a context's `measureText`, the font set only when it changes; its
  `capHeight` the `actualBoundingBoxAscent` of `H` (a letter's size follows from it).

## Input

`input.ts` (`ViewerInput`, `viewer.input`): pointer, wheel and key events as calls in canvas coordinates, DOM-free; the
mounted UI (`mount.ts`) translates events into them and a test makes the same calls. Every gesture goes through the
binding table; per move one hit record and one modifiers record are reused.

- **A drag ends wherever the button is released**: besides pointer capture (which fails silently where it is refused: an
  iframe boundary, a pointer already up), every drag a mount starts (`mount.ts`, `mount-navigator.ts`, a story's own)
  runs a drag session (`drag-session.ts`) for as long as it lasts: `pointerup`, `pointercancel` and a `pointermove` with
  no button on the element's own window (`ownerDocument.defaultView`, an iframe's own), the window's `blur` and the
  document turning hidden end it, as do the element's `pointerup`, `pointercancel` and `lostpointercapture`. Each path
  ends it once (`up` is idempotent), the listeners are removed with it and on unmount, and a press ends a drag whose
  release never came before it begins. Why: a window listener alone misses a release outside the browser window; capture
  alone misses where it is refused.
- **Overlays own their events** (`overlay-events.ts`, `isViewEvent`): an event is the view's only when it happened on
  the container or one of its canvases; in any other element inside the container (the header menu, the header parts,
  the keys panel, a part a behavior attaches) it is the overlay's. The pointer moving onto an overlay is a `leave` for
  the view (its hover cleared, the card gone, once), except onto the header parts (they stop the move: the pointer on
  them is still on their header, which keeps its actions shown), and the overlay's presses, clicks, double-clicks,
  wheel, keys and secondary click (`header-gestures.ts`) never reach the view; a release and the losses of the pointer
  still end the view's drag. Why: the overlays live in the container (its tokens, its place) so their events bubble to
  its listeners; a press on a menu item over the data began the view's drag and captured the pointer, so its click went
  to the container (the item never ran) and selected what was under it.
- **A window's part under the pointer** (a ring's, a rail's) is lit (`windowHit`, drawn as the emphasis's hover), and
  the window being dragged grows its rail wherever the pointer goes (`dragHit`, `Emphasis.drag`); the cursor follows the
  part (`partCursor`, `cursor` set by the mounted canvas, `NavigatorModel.cursor` alike): a grip `ew-resize`
  (`ns-resize` on the rows rail, an overview, a miniature's rows grips), held through its drag; the body `grab`,
  `grabbing` dragged; outside `pointer`; a ring's window `grab`. A finger's pointer (`pointerType: touch`) hits a window
  wider (`navigator.window.touchHitPx`).
- **A ruler rail** (Composition, Navigators): a press on its grab zone goes to the `navigator` behavior at the rail's
  scale (`railUnitAt`): its body pans, its grips zoom, outside it a drag sweeps; a press that never moved acts at the
  click the platform sends after it (the rail's: no selection): outside the window it centres there, on a grip or the
  body nothing; a double-click on the rail is the navigator's (`toggle-fit`, `navigator.zoomPx` at the rail's scale). An
  idle rail is no target: the press falls through to the ruler (on an alignment's columns, a column pick above the grab
  zone).
- **A double-click on a block** of a row with nothing under the point (its label, its band's ground; where the binding
  is a fit) or only a position (a ruler, letters, a value: no item to focus) fits that block (`view.focus` with its
  units), and from that fit the whole (`toggle-fit`, one level deeper: never a zoom to one base of a chromosome). An
  item under the point is focused as anywhere. The hover card names a block's row by its system: `Sequence` for a
  nucleic sequence (an assembly's chromosome or contig, a plain DNA or RNA sequence), `Chain` for a complex's or a
  protein's.
- **The rows rail** (`RowsRailGesture`, `behaviors/rows-rail.ts`): a press on its grab zone, in display rows at the
  rail's scale (`railUnitAt` of the page's y): its body scrolls the rows (`bandTop`, clamped to the last rows), a grip
  zooms them about the other end, which stays (`rowScale` within the rows' bounds, core's `MinRowScale` to
  `MaxRowScale`, never past the rows there are), each by the rows the pointer travelled (a floored window reaches every
  row); a drag farther than the rows in view, or a zoom, is one step (`navigator.dragStep`); Esc puts the band back. A
  press that never moved acts at its click: outside the window the rows centre there (a scroll, no step), on the window
  nothing; a double-click on the window resets the rows' zoom (`resetZoom`: `rowScale` unset, the band's top row kept; a
  jump), outside it nothing more. A vertical overview's rows grips (its window's top and bottom bars) are this gesture
  over its window, never a scroll; a double-click on its window resets zoomed rows the same, else it is the positions'
  (one gesture, one step). Why: the rows rail and the overview are one window over the rows. `]` and `[` (the rows
  behavior's `rows.zoom-in`, `rows.zoom-out`) zoom the rows; PageDown and PageUp (`cursor.page-down`, `cursor.page-up`,
  with no cursor shown on rows: as `rows.page-down`, `rows.page-up`) move the window by the rows it shows less one, the
  page's own key at an end.
- **The one rail** (core's "Render", Section scroll): every box that scrolls answers one gesture (a drag on its window,
  a press beside it centring there, a scroll and never a history step, Esc back to where the press found it). A section
  past its cap is the viewer's own (`section-scroll.ts`: its top in px, a tracks section's viewport `top` through
  `setTop` (a follower's too), a text's anchor moved by px; `SectionRailGesture`; `SectionWheel`; the stored top kept
  within the content after each layout, `clampSectionTops`): every composition scrolls it, a minimal one too. A rows or
  a lanes band's is the rows behavior's (`RowsRailGesture`; the rows rail adds its grips). The keys' scrollers
  (`behaviors/scrollers.ts`: a rows band in display rows, a lanes band in lanes, a section in px). **The wheel's
  chain**: its vertical delta (no modifier, mostly vertical) scrolls a rows band or a lanes band under the pointer, at
  its end the section past its cap, at the section's end the page (false); over a section past its cap and none of its
  scrolling bands it is the section's at once, coalesced to the frame (allocating nothing, as the navigation's wheel).
  **Keys**: PageDown and PageUp (`cursor.page-*`) and Home and End (`cursor.scroll-start`, `cursor.scroll-end`,
  `view.keys`) move the innermost scroller under the pointer that can (a cursor shown on rows takes the page keys first;
  the page keys fall back to the first rows band); at the ends the key is the page's. The keyboard cursor on a band's
  row below the cap scrolls the section the least to show it. Why one rail: a scroll affordance that looks the same but
  answers differently is worse than none.

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
  the section's first visible). A text (`Placement.clip`) takes a wheel without Cmd/Ctrl only while its content
  overflows the box, scrolling it through its layout's `pan`, and gives it back to the page at its top and bottom
  (scroll chaining); a tracks section past its cap takes the wheel's vertical delta by the one rail's chain. True when a
  binding took it, and the page does not scroll. **The wheel over a plot is the plot's**: over a grid's data box (x's
  extent, the blocks' gaps included, and the well's height; not the label column nor the y blocks' labels in it) it is
  always true, so the page never scrolls from inside the plot, even where a pan stops at an axis's end; outside it, as
  anywhere. Why: a plot whose wheel fell through at an end moved the page under the pointer.
- **Double click**: the bound command with the section, what to focus (the element's hull on the section, its extent
  carried along the layer's path, else the position's cell, in offsets; a Shift+click's extend before it is given back,
  as a clear one is), the unit the placement has under the pointer (`at`) and the double-click's zoom factor (`step`,
  `navigate.doubleClickZoomStep`). On a grid's cell the context also says `cell` and the unit down y (`atY`): the focus
  zooms a step about the cell on both axes (below, `navigate`). Every point of a plot's data box is a cell's: a gap
  between blocks resolves the cell's binding and zooms about that point, never fits (`empty`'s). On a grid whose y
  follows x (one window) the step keeps the point under the pointer on x, y following.
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
- **Announcements** (`announcer.ts`, `Announcer`, DOM-free; the mount owns its one `role="status"` `aria-live="polite"`
  region, visually hidden, in the container, and shows `announcer.message` in it, cleared and set so the same words
  twice are heard twice). The words are `describe.ts`'s (`Describer`), the hover card's too: a place goes through the
  description's sections in order (a pair: its label, its two ends, its kind; an item: its label, extent and category; a
  rows layer's row: its label and the column cell named; a grid cell: its two positions, by their blocks' labels over a
  complex, and its value; a cell: the section, the system and the position, its labels in every declared numbering of
  the system (the section's shown one first: `kabat, h 13 (kabat 82B, imgt 107)`), the letter, each value track's value,
  the items there; last a block's chain and copy), each adding a row and its words. The sections are internal (chrome
  parts behind internal interfaces until a host needs one); the viewer has one describer (`describerFor`, by its spec
  model), the announcer's and the card's. A declared numbering with no label at the place (absent, or its labels not
  arrived) says `–`. A value says the unit its layer's theme declares (a `field` form's or a scale's `unit`: `9 Å`).
  What is said: the keyboard cursor's place after each move of its own (section, system, 1-based position, the letter a
  letters layer shows there, the items under it by their labels and extents: the same words as the exported SVG's
  titles, `accessibility.ts`; on rows first its row's label, a summary row's as drawn, `HBA +1`, a lazy row not yet read
  `pending`; on the header the row alone; on a grid the cell's two positions and value:
  `Describer.describe({ cursor })`); the hover's when it rests (`restMs`, 500 ms: a pointer crossing the picture says
  nothing; the item's label and extent, else the cell's place, with its label in the numbering the section's layers
  show, the first naming one: `kabat, h 13 (kabat 82B), letter Q`); a selection's change once
  (`Selected 6 positions, p 3–8`, `Selection cleared`; the same selection again says nothing). One rate limit for all
  (`minGapMs`, 700 ms): a message never follows another sooner, the latest one waiting its turn (a held arrow key says
  where it is each gap, not each key). A message is a token (a counter and its text), so equal words are said again: a
  cursor key that leaves the cursor in place (the end of its axis) is announced (`viewer.commandRan` moves after every
  command). A message waiting its turn is replaced only by one of the same or a higher priority (selection over cursor
  over hover), and the pointer leaving drops a waiting hover. While a gesture is dragged the cursor and the hover say
  nothing and the selection says its end once, after the drag. Why one region and a limit: a live region queues what it
  is given and a reader hears the queue after the cursor moved on. A navigator slider's value is its own
  (`aria-valuetext`), never this region's. The highlight legend's Clear buttons and the legend's entries are real
  buttons and options in the page's tab order, so they need nothing here.
- **Forced colors** (`forced-colors: active`, whatever the `theme.source`: the root's CSS tokens through
  `readCssTokens`, a configured theme's table through `forcedColorTokens` over `readSystemColors`; `css-tokens.ts`):
  `readCssTokens` reads the system colors (`Canvas`, `CanvasText`, `Highlight`, `GrayText`, resolved by a probe
  element's computed `color`; the platform's defaults where none resolves) into the chrome and emphasis tokens: the page
  and tracks `Canvas`; ink, the sequence, lines, breaks, the ruler, the window's grips and a letter cell's hairline
  `CanvasText`; muted text, missing and pending `GrayText`; the accent, the hover and the window `Highlight` (the wash
  and the muted accent a translucent `Highlight`). The data's colors stay (the meaning is in them: a residue class, a
  category, a ramp), and the legend's swatches opt out of forcing (`forced-color-adjust: none`) to match. The mount
  listens to the media query's change and reads the table again; the table keeps its theme (`light` or `dark`: schemes
  pick by it), only its colors and version differ. An export never uses forced colors.
- **Toolbar** (`toolbar.ts`, `toolbarButtons(viewer, commands)`, DOM-free; `mount-toolbar.ts` mounts `@molstar/ui`'s
  toolbar over it, above the canvas): a part the host places (`ui.toolbar: mountToolbar`; the default spec lists none:
  `@molstar/ui`'s toolbar with its icons and menus is about 85 kB, the default viewer 508 kB without it), one button a
  command, `view.back` and `view.forward` known (`DefaultToolbar`). Its title is the label and the first chord of the
  command's `view.keys` as the platform writes it (`Back (Alt+←)`; no key when the application removed them), also its
  accessible name; it is disabled, with the reason as its tooltip, while its side of `view.history` is empty (the state
  is reactive: the bar follows jumps, back and forward, and edits that clear the history). A click runs the command on
  the first leading section. Why a model apart: the twin of the mounted bar reads the same states in Node. `help.keys`
  is known too and in `DefaultToolbar` (the keyboard icon, "Keys and gestures (?)": where to find the keys, said once,
  quietly), shown only when the viewer has the command (the `keysPanel` behavior).

- **The placement is the input's inverse** (`offsetAt(section, x, y)`, `viewer.placementOf(section)`): the unit under a
  point, a section's targets (a window with its part, a text gutter's `row`), never the input's own arithmetic.

## Behaviors

Small policies over the models (`behaviors/`), each reading the bindings. `rows` (an alignment's rows: Commands, the
rows' commands) reads the last frame's rows layers through the behavior context's `layout()`.

- **A section of blocks on several systems** (its axis on no system): a hit names its `block`, and its positions are
  read on that block's axis from its `offset`: a drag-select stays in the block it began on (clamped at its ends), a
  double-click and a focus move to the block's offsets, a bookmark of the view names none. **The keyboard cursor crosses
  blocks** (EC-154, amended): a step goes along its block in display order, past the block's last position onto the next
  block's first (before its first, the previous block's last; a block with no positions skipped), at the section's ends
  it stays; the step's message names the block and the position (`Moved, chr18, 1`); revealed at its block's place. The
  cursor stays on the copy it is on (a position two blocks show: the block a step left it on, kept by the cursor), its
  next step, reveal and message too. A run of small blocks hits as `run` (its band with no strip: the block under the
  point, the run's count; under the ruler rail the rail's place carries its `count` and `state`), its card
  "chrUn_KI270302v1 · Run: one of 2,841 · not read at this zoom".
- **`navigate.acrossBlocks`** (`ViewModel`, a row of blocks): `cross` (the default) one viewport over the blocks end to
  end; `stop` keeps the window in one block: after each move it is clamped to the block it was kept in (its width at
  most the block's), a range set (a jump, a story step, a click on the navigator) moves to the block holding its centre,
  a fit is the block's. The kept block is remembered by its copy key (system, range, form) per section, so a block list
  changed under it (`section.include`, a reverse) keeps the same chromosome, never whichever block takes its index. Why:
  a position on the section's units is no system's (core `_spec` "Resolved view", Blocks).

- `hover`: a move over a hit whose hover gesture is bound to `hover` sets the hover for its owner token; a move over
  nothing or a leave clears what that owner holds only. A host's other hover sources use their own tokens.
- `navigate` on a grid whose y moves (`movesY`): a wheel's dy pans y and its dx x (Shift: as elsewhere), a zoom and a
  pinch zoom y about the pointer's y too (`atY`, from the placement's `strip.y`: `Input.offsetYAt`).
- `navigate`: a drag (`pan`), a wheel (`pan` by the larger delta, so Shift with a vertical wheel pans sideways; `zoom`
  about the pointer) and a pinch (`zoom`) through the bindings; wheel deltas within one frame are applied once at the
  frame (`navigate.coalesceWheel`, through a frame callback that `dispose` cancels); a section gone by the frame is
  dropped; `dispose` drops what is pending and later frames do nothing; a drag passes both deltas (a text scrolls its
  lines by Space+drag). The view's moves are commands (`view.fit`, `view.focus`, `view.back`, `view.forward`,
  `view.zoom-in`, `view.zoom-out`, `view.pan-left`, `view.pan-right`, `view.pan-up`, `view.pan-down`: a line in a layout
  of lines, Shift+Alt+↑↓) with their keys from `view.keys` (the first the shortcut), single keys scoped to the viewer
  (`scope: 'canvas'`); `when` says why one cannot run (no section, a follower that does not move, nothing to focus, no
  step back) instead of a throw. A fit and a focus are jumps; a focus on what the view already shows focused goes back
  (`view.backOnFocused`). A focus on a grid's cell (`cell`) is a step zoom about it on both axes instead (`at`, `atY`,
  by `step`), their ratio kept, a jump of its own (a step of the history: again zooms further, back undoes each); why: a
  jump to the cell's two positions re-centred each axis on its residue, so a second double-click shifted the view
  instead of zooming. A step zoom (the keys, Shift+double-click) is about the context's offset (the pointer's), else the
  keyboard cursor's cell on that section when in view, else the centre, by the context's `step` or
  `navigate.wheelZoomStep`.

- `select` on a grid: a hit with a `ySystem` that is no rows axis and no element is a cell of two axes, and the same
  commands and drag modes act on its rect (`selectCell`, the rect drag): click a one-cell rect, Shift+click the rect
  from the anchor, Cmd/Ctrl+click toggles, a drag the rect, Cmd/Ctrl+drag adds, Alt+drag subtracts.
- **A rect over blocks** (`select` on a grid over a section of blocks, a grid over a complex): the drag is the
  behavior's own (`behaviors/blocks-rect.ts`, internal: not re-exported): each cell's display units on the hit strip's
  own blocks (`Hit.block`, `Hit.yBlock`: a copy split is told from its twin), the span between the two cells on each
  axis split at the blocks' boundaries, a rect a pair of parts (one per pair of systems), set on the selection by the
  drag's mode from the selection as it began (`set`, `add`, `subtract`), Esc giving that back. The model's rect drag is
  unchanged. Why units: a drag across a chain boundary is one gesture over the plot, its selection a part a chain.
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
- **The cursor on track headers** (the `headers` behavior's, `behaviors/header-cursor.ts`: it takes the `cursor`
  behavior's moves over once composed, so a viewer without it carries none of it): `cursor.prev` at a section's first
  position goes onto its first track's header (`part: 'header'`, `track` its address; the section's line its section's
  address), drawn with its actions (`Emphasis.cursor.track`); ↑ ↓ step between the section's headers (its line first),
  then on to the next section's; → back to the data where it left; Space and the extends do nothing there; Enter opens
  its menu. After a layout that hid the cursor's track the cursor goes to the next track's header, else the previous
  one's, else the section's line. A rows track has no header line: its rows' header is the rows' own.
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
  `cursor.keys`; a cursor moved out of view brings its section to it by the layout's least move (`reveal.cursor`). **On
  rows** (a section with a rows layer; the cursor's row a data row): the first key places it on the first row shown too,
  and ↓ from the section's other tracks (a cursor with no row: a click on its ruler) enters its rows there; `cursor.up`,
  `cursor.down` step display rows (the pinned band, then the scrolled rows, where a pinned row is in its place too; ↓
  from the band's last row lands on the first scrolled row shown, the band keeping its scroll; a folded clade's summary
  row one stop, the cursor on its first shown row), past the first or the last to the section above or below (entering
  rows from above on the first row shown, from below on the last shown); `cursor.prev` on the first column goes onto the
  row's header, `cursor.next` back to the column it left; `cursor.page-down`, `cursor.page-up` (PageDown, PageUp, in
  `view.keys`) move a cursor shown on rows a page of rows and the band a page with it (`rows.page-*`'s page), and with
  none shown scroll the band a page as `rows.page-*` do (those take no key: one key, one command); at the band's end the
  key is the page's; `cursor.first-row`, `cursor.last-row` (Mod+ArrowUp, Mod+ArrowDown in `view.keys`; the page's keys
  when there is no cursor on rows) the first and the last display row. The cursor stands on the copy of a pinned row a
  step or a click left it on (the rows' last hit, `RowsLay.hitRow`), kept by the cursor; with none, its row's own
  (`RowsLay.cursorRow`: a folded row's summary row, a hidden or filtered row's nearest row shown, drawn, said and carded
  there). Shift with an arrow on a cell extends the block: the rect from the cell the extend began on (a selection
  changed since begins another) to the cursor's, the rows shown between them in display order as data rows (sliced from
  the rows' axis, `RowsLay.rowsAxis`: a million rows cost what ten do); Cmd/Ctrl+Shift with an arrow (`cursor.add-*`,
  `cursor.add-up` and `cursor.add-down` among them) starts a new block at the cursor's cell, the others kept, which
  Shift then grows. The block is the behavior's, not the model's cell anchor (rows in display order, the others kept).
  On the header Shift with ↑ ↓ extends the run of rows from the anchor, and with ← → does nothing (the column it left is
  not shown). Space on a cell toggles it (`SelectionModel.selectCell`; a summary row's: its rows added, or taken out
  when all are in), on the header folds or unfolds the row's clade (a summary row's, else the leaf's parent; no tree:
  the row toggled). Enter on a summary row opens it; on a row's header it focuses the row's extent (its first to its
  last letter, held or read from its tile; a lazy row not read yet, every column; with a selection, the selection as
  before). Moves are arithmetic on display rows: a row never loaded is stepped over by its index. Each move scrolls the
  band to show the row by the least move (`reveal.cursor`). **On a grid** `cursor.up`, `cursor.down` step its y in
  display order (past its ends to the section above or below), `cursor.first-row`, `cursor.last-row` go to its ends,
  Shift and Cmd/Ctrl+Shift with an arrow grow and start a block (`cellsRect`); its y is shown by the least move (a y of
  its own, or one following a section, pans by the px the cell is past the plot, `ViewModel.panPx`; a y showing x's
  window reveals as x does). The moves are in `behaviors/cursor-moves.ts` (internal: the provider gives them the last
  frame's layout); the exported `cursorBehavior` has no frame, so it moves on one axis.
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
- `navigator` (`NavigatorProvider`, in the default spec; a host may leave it out, its navigators and rails then drawn
  with no targets): a navigator's gestures on a `window` behavior of its own, in the bound section's display units. A
  drag on the body pans by hand and an end zooms about the other (a floored window's end dragged inward keeps the view,
  outward zooms out from the drawn end: `begin`'s `drawn`); a press outside the window that does not move is a click
  that centres there (`navigator.center`, a jump: animated), dragged it sweeps a region and zooms to it on release
  (`navigator.zoom-to`; `navigator.dragOutside: pan` centres and then pans); a double-click fits, or fitted zooms so the
  window is `navigator.zoomPx` navigator px wide (`navigator.doubleClick`); Esc cancels, the view back where the press
  found it. One gesture is one step: a click centres with its own step (a jump), and the double-click it begins (the
  platform's, by the click's count: its second click does nothing, its `dblclick` follows) moves within that step, as
  input's first click defers a clear (no timer: the platform knows a double-click); a drag pushes the viewports at its
  press on release when the view went farther than one window width (`navigator.dragStep: far`; `never`, `always`). A
  window bound to a follower moves its leader (`follow.followerMoves`). The ruler rails take the same gestures at the
  rail's scale (Input).
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
- `hoverCard` (`HoverCardProvider`, in the default spec; `hover.card` false turns it off; the port is its `port` param,
  Ports, Hover card): the card shows `restMs` (the announcer's 500) after the hover last changed, beside the pointer,
  for a pointer over the canvas (a host's hover through `address.hover` shows none); one timer waits for the rest
  however many moves come, the card made then, never per move (the input tells it the pointer in a reused array: a move
  allocates nothing). Hidden at a press, a wheel, Esc, while a drag runs, on leave and at the next hover's change (a
  re-hit of the same element, cell or target under a resting pointer, a zoom's or a pan's, is no change: a cursor's card
  stays); a press, a wheel, Esc and a leave also clear the rest waiting (no card follows them). A key that moves the
  keyboard cursor shows the card at its cell (on rows its row's cell or header, on a grid its cell) after the next frame
  (`framed`, the viewer's, each frame: the one that drew the cursor, a reveal's scroll included). No card outlives what
  it describes: a frame laid out again (a zoom, a pan, an edit, a load, data or a source's labels arriving) places the
  cursor's card again at its cell and sends a hover's away until the hover rests again. A rectangle dragged on a grid
  shows its summary while it is dragged (`gridSummary`: the cells counted, their mean, least and greatest value), at
  most once a 50 ms. **The enterable card** (EC-166's bridge): a card with a link (a tooltip's `!seq-…`) waits
  `BridgeMs` (300) after the hover leaves its anchor or the pointer leaves the canvas, so the pointer can cross the gap
  to it; entered (the adapter's `enter`), it stays while the pointer is over it, the hover frozen; its `leave` hides it
  and arms the next rest; the grace ended without an entry, it goes. Esc, a press and a drag hide it at once. Why: a
  link a reader cannot reach is no link; a card with none never takes the pointer. `hover.card` is read at composition
  (the configuration is read-only after it): off, no card is shown and no port made. **A variant allele's card** (the
  genome formats' matrix row "a variant allele"): an item whose source row holds `ref` and `alt` (a VCF's, core's Data,
  VCF) is titled by its id `pos:ref>alt` (an allele past 32 letters, hashed in the id, as its first letters and its
  length: `7674964:GAGGCT… 469 bp>G`, `alleleTitle`), never its record's ID, with its Category, its record's `ID`,
  `QUAL` and `FILTER`, and its INFO fields with a value (every field not in core's VCF `FixedFields`, the source's
  order, the first 8, then `INFO: +k more`; a value past 40 characters cut with an ellipsis; an integer written as it
  is: an id is never grouped). Its words are read once a hover (`alleleOf` keeps a list's last, a weak map: 0 B a held
  hover, the storybook's `allele-card.quick.ts`). A behavior's part over the canvas is
  `BehaviorHandle.attach(container, viewer)` (`viewer-spec.ts`): `mountCanvas` calls it for every behavior that has one,
  and the function it returns takes the part away.
- `keysPanel` (`KeysPanelProvider`, in the default spec): `help.keys` (`?` in `view.keys`; the toolbar's button) opens a
  `@molstar/ui` surface over the canvas's corner listing what the view answers to, `describeKeys(viewer)` (`keys.ts`,
  DOM-free, so a host draws its own): generated from the viewer's commands with their keys (`view.keys`, else the
  shortcut, as the input reads them) and from the binding table's gestures, grouped (move, zoom, select, edit, view), an
  entry a command or a method once with every key and gesture that runs it; `help.keys` again or Esc closes it (on the
  canvas a running drag's Esc cancels the drag first; then an open panel's Esc closes it and does nothing else). Open or
  closed is the reader's (`none`: nothing in the models changes). Non-modal: opened, it takes the focus (its Close
  button); closed, it gives the focus back to what opened it (the canvas, a toolbar's button), else to the canvas; the
  focus may leave it (no trap). Its events are its own: no pointer event or key in it reaches the canvas (an arrow moves
  no cursor, Enter on Close runs no command), and Esc on it closes it wherever the focus came from. Its DOM is a lazy
  chunk (`keys-panel-dom.ts`) asked for at the first open (`loaded`, which `settle` waits for). Why: a reader opens it
  rarely, and the default viewer's size has no room for what it does not draw.
- `linkedPairs` (`LinkedPairsProvider`, in the default spec, removed to opt out; `linked-pairs.ts`): on a hovered or
  selected range element (a feature, a CDR), the pairs with an end in its extent, on every pairs layer of the view
  (`layers`: those named by address, all when absent), written to its own highlight channel (`linked-pairs`, element
  keys alone): drawn in the highlight pass in the channel's color, never the selection's (related to the selection, not
  in it: a copy of the selection leaves them out); cleared when the hover or selection goes; `on` says which lights them
  (`hover`, `selection`, both by default). Built on the public API alone: the query language
  (`<pairs>/(items :overlapping <extent>)`, pairs answering `items`; the extent a hit's, else `(extent :of <element>)`;
  a selection of range elements alone by its `shown`), the highlight channel, the hover and selection models. The
  hover's query is cached per hovered element key while the view and the model hold (an element read from a source whose
  address still waits is asked again), and the channel's value per set of keys: a hover allocates nothing a frame past
  its first. A pair or a position hovered is the pairs layer's own rule (its partners light), never the channel. Why a
  behavior: a host chooses an analysis, and the same query serves the selection.
- `headers` (`HeadersProvider`, `behaviors/headers.ts`, in the default spec; left out: states only, no slots, no menu):
  the track headers' entries and gestures (core's Track headers draws and hits them). Its eager part holds the entries
  (their checks and `register`), the menu's state, the ports and a click's way in; its late parts (`header-extras.ts`, a
  lazy chunk asked for once composed: the facts and the layout's parts, the menu's data, its opening and the runs, the
  geometry and the parts' state, the secondary click and the menu key, the words, the cursor on headers) work once they
  land. Until then the headers draw core's built-in states, a click on a header part is kept and run once they land
  (one, the last: never lost), and nothing is published. Its lazy chunks (the late parts, the icons, the gestures, the
  ports' DOM) are its handle's `loaded()`, which `settle` waits for; one landing after the viewer is retired is ignored.
  Why: a settled picture is the whole one, and no import outlives a test's viewer. A header part's click runs its
  entry's command (one path, `Input`): through the behavior, or without it for a state listed as an entry with the
  layers it is about. Internal (`@internal`, no entry exports them): the late parts, the DOM ports' modules, the icons'
  reader, `headersBehavior` (the provider is the way in), and the handle's `partsOf`, `framed`, `click`, `context`,
  `loaded` and `preview`.
  - **Entries** (`header-entries.ts`): `HeaderEntry` (`name` namespaced `<id>:<name>`, `on` `track`, `section` or
    `element`, `kind` `action` or `state`, `icon` `@molstar/ui`'s, `label`, `when`, `command`, `args`, `text`, `order`,
    `slot`, `items` a submenu of other commands or args: a layer's checkboxes, a chip's list), listed in the registry as
    `{ kind: 'header-action', entry }` (`SeqViewerRegistryEntry`) or registered at run time (`register(...entries)`,
    returning its undo). Registration is atomic: every problem in one error before anything changes (a name twice, no
    namespace, a command the viewer lacks, `seq.viewer.header.command`); the registry's are checked once every behavior
    is made (`BehaviorHandle.composed`). The label and icon default to the command's. **A command whose one param is a
    select is a submenu** of its options (an entry with no `items`): each item the entry's args with the choice under
    the param's name, the param's default checked, previewed while hovered through the command's `preview`; a command
    reads the choice from its args (a palette's) or its context (`viewer.execute`'s). Why: the command owns its choices,
    so the menu, a palette and an agent offer the same. Why entries: one list for a slot, a menu item, a DOM part, a key
    and an agent; the registry's rules (atomic, undoable, namespaced) for free.
  - **Built-ins** (`BuiltInHeaderEntries`, in the default spec's registry; a host leaves one out to remove it):
    `seq:hide` (`eye-off`, `track.hide`), `seq:color` (menu only, "Color ▸": one place, never also a slot; the first
    shown layer with more than one fitting scheme: `layer.channel`'s choices a submenu, previewed while focused or
    hovered; its arrow at the row's end, the menu's), `seq:layers` (more than one addressable layer: a checked item a
    layer, `layer.hide` or `layer.show`), `seq:move-up`, `seq:move-down` (`track.move`, a sibling there), the states
    `seq:filtered` (its count the track's, every thinned layer's rows; `layer.filter` clears every thinned layer's
    filter, one patch), `seq:hidden-layers` (`layer.show` of every hidden layer, one patch; its list each) and the
    section's `seq:hidden-tracks` (menu only, "Show hidden ▸": a flyout of each hidden track, then "Show all",
    `track.show` of every hidden track in one patch), the section's `seq:reverse` and `seq:collapse`. `seq:pending` is
    not built (the TODO).
  - **`HeaderFacts`**: plain data from the resolved view, the spec and the last layout (the address, label, section,
    system, linear, the section's track ids, the track's layers with type, label, hidden, filtered rows, color and
    fitting schemes, the hidden tracks, the element). Entries' functions are pure over them, run once a layout a shown
    header (`partsOf`, the frame's `headers`, a new function when the entries change: a new layout) and when a menu
    opens.
  - **Shown or enabled**: an entry's `when` decides whether it shows, its command's `when` whether it runs (a disabled
    item muted, its reason its `why`). An action with `items` opens them (its slot opens the menu at it); a state's
    click runs it, its secondary click lists its items, then itself (the `hiddenChip: menu` param: the click lists
    them).
  - **Slots** are the applying actions with an icon (or `slot: true`), by `order` then name, the first
    `actionColumns(…).n`; every entry is in the menu. A slot's icon (the entry's `iconMarkup`, else a built-in's by
    name: `@molstar/ui`'s per-icon modules, never its map of every icon) is read once into path data (`header-icons.ts`:
    the Lucide markup's `path`, `line`, `circle`, `ellipse`, `polyline`, `polygon`, `rect`; arcs and quadratics as
    cubics), a lazy chunk, loaded once an entry has an icon: until it lands, and for an icon no path reads, the slot
    draws its label's first letter; then the headers are laid out again.
  - **Gestures**: a click on a slot runs it (the "⋯" opens the menu), on a chip runs its state, on a label nothing (on
    touch it shows the slots, as a hover there, so a second tap on a slot runs it); never a selection's. A secondary
    click (`contextmenu` on the mounted canvas) opens the menu under `contextClick` (`everywhere`, the default spec's:
    at what it hits, an element's entries first, then its track's and its section's; `header`, the label column only;
    `none`, the browser's); Shift+secondary click is always the browser's. A long press (touch, the `longPressMs`
    param, 500) opens it. The gestures are a lazy chunk loaded at the mount (`header-gestures.ts`). The menu key and
    Shift+F10 (`headers.menu`, `view.keys`) and Enter on a header open it at the cursor's header (its section's first
    track's on data), once per ask: the browser's secondary click that follows within a second is taken and opens
    nothing (EC-156).
  - **The cursor off what was hidden**: a frame whose layout no longer has the cursor's track moves it to the next
    track's header, else the previous one's, else the section's line.
  - **Published on change** (by name from the frame, `framed()`: no listener; a frame that kept its layout reads one
    field and makes nothing): `geometry` (`viewer.headers`, `HeaderGeometry[]`: target, line, label, slot and chip boxes
    with their names, the well, `shown`), compared structurally; `parts` (`HeaderPartsState`) only with a parts port,
    again when the header whose actions show changes.
  - **Words** (`describe.ts`, the describer's `headers`): a slot its label and key ("Hide track", why when disabled), a
    chip its text and what a click does ("12 filtered: click clear the filter"), a label its track's label and states
    (its layers for the card), the cursor on a header "…, Enter for actions". The card of a header part is a title
    alone. A header's hover lights none of the track's items (EC-146: refused by default).
- `links.hover` (`LinksHoverProvider`, `links-hover.ts`; not in the default spec: how events from different parts
  compose is the application's, the examples list it): hovering a molecule (a chip, a residue, an atom) lights the
  positions its links reach (a complex's `pairs` children, each end on its row's unit); hovering a linked site lights
  the molecule's whole unit (a link between two polymers, a disulfide: the other end). Written to its own highlight
  channel (`linked`), cleared when the hover goes; a merged block's hover is each of its copies', and an end on a merged
  copy is lit on the block it is merged into, a residue on an atoms-form block as its atoms (the atoms entry's offsets).
  A link is never a map: nothing is projected across it. Contacts get no such behavior (a layer the host adds).
- `editing`: edits of the spec as commands, each one patch made against the revision read (`SpecModel.patch`, one step):
  `layer.hide`, `layer.show` (the `hidden` key set or unset), `layer.color` (the layer's theme), `track.hide` and
  `track.show` (the track's own `hidden`: its layers keep theirs, so hide then show is a round trip; `layer.show` and
  `track.show` take `addresses` too, every one in one patch: a chip's "show all"), `track.move` (among its siblings,
  `after`), `track.add` (a `node` inserted under the section `address`, after `after`, last without it; the id taken is
  refused), `refs.remove`, `refs.label` (a bookmark's `label` param set, or removed with `''`); `step.enter` and
  `step.leave` (navigation); `preview.cancel` (Esc). A gesture's preview is a method, never a command (it runs per
  pointer move): `preview(patch)` replaces the last, `endPreview(label)` ends it in one step, `cancelPreview()` in
  nothing. `when` says why one cannot run (no such node, no step shown, nothing to cancel).

### Molecules: the interaction matrix's cells

The cells the molecules design names, each built or not and why (the core's "Resolved view", Blocks; "Render", Blocks):

- **A section of a source's complex** after an update: reflows (the view resolved again over the new overlay); no
  history entry, because it is data, not an edit of the spec.
- **A section of a complex**: zoom and pan as its layout; hidden by its tracks; `classes`, `copies`, `arrange` by a
  patch (`SpecModel.patch`): no `section.set` command yet, because nothing but a host's own menu would run it.
- **A chip** (`form: whole`): a click selects its block's every residue (its hit's `extent`), Shift+click and
  Cmd/Ctrl+click as an element; double-click focuses that extent. Not built: a sibling section zooming to the chip's
  residues or atoms through the atoms map, because the viewer's focus carries a location along static paths only and the
  atoms map is applied at a path's ends by the model's projection, not by the focus. Zoom and pan: no, a chip is fixed
  in px and the window skips it (one unit). Hidden: no, by `classes` or the blocks by patch. Its form: the `block.form`
  command (Commands), no context menu yet.
- **A modified residue**, **a molecule's residue or atom**: as an element; a modified residue's mark hidden by
  `modified: none` on the layer (by patch). Atoms collapse onto residues for other sections where the model's projection
  does it.
- **A link** (an arc, or a site's mark on a `ranges` layer over links): a click selects the link, both ends lit, each on
  its row's unit; Shift and Cmd/Ctrl as elements; double-click its extent, two parts on two systems; hidden by its
  layer's items (`layer.filter`).
- **A range over blocks with chips**: no drag across a chip in reading order, because a drag on a section of blocks
  stays in the block it began on (Behaviors, a section of blocks on several systems): a position on the section's units
  is no system's. Enter focuses the selection.
- **The selection**: `(filter :test (= (prop modified) true))` and the other forms through the query command; a bookmark
  by `selection.bookmark`. **The view**: 0 fits, chips keep their px. **A track or layer**: `track.hide`; `modified` by
  patch. **A bookmark**: `refs.add`, its action run. **A follower**: never a section of several systems
  (`seq.render.follow.blocks`). **A navigator**: chips skipped as the window skips them.

## Bindings

`bindings.ts`: the interaction grammar as data, `{ gesture, modifiers, target, action }`.

- **Pointer gestures** (`hover`, `click`, `double-click`, `drag`, `wheel`, `pinch`, `back-button`, `forward-button`);
  keys go through `@molstar/ui`'s keymap, and a `key` binding here is reported. Modifiers are `shift`, `mod` (Cmd on
  macOS, Ctrl elsewhere, `modifiersOf`, which fills a reused record), `alt`, `space`; exact, or `'*'` for whatever is
  held, tried after the exact ones. Why: a modifier held for something else must not kill the hover.
- **Targets** are a hit's kind or `any`; a binding for the kind beats `any`'s. A collapsed break is its own target
  (`break`, the render hit's `target`): no default binding but the hover, so a double-click on it does nothing (the
  pointer aimed at the break; showing its stretch is `section.collapse`, a spec change). A header's parts are targets
  too (`header-slot`, `header-state`, `header-label`, `section-header`: core's Track headers): never a selection's (a
  click there selects and clears nothing; with the `headers` behavior a slot or a chip runs its entry; without it a
  state listed as a `header-action` entry runs that entry's command over the layers it is about (core's state
  `addresses`: the "n filtered" chip, `FilteredHeaderEntry` in `filter-chip.ts`, clears every filtered layer in one
  patch, the minimal composition lists it), the other states are facts a click does nothing on. Why: the chip has
  cleared its filter since the third slice, with or without the menu's behavior; one command, the entry's, either way.
  Not the select behavior's: a viewer without it still takes a header's click), and a drag on a track's label line pans
  as its track (`header-label`, in `default`). A slot is hit only while its header's actions show.
- **A track's row keeps its header's actions shown** (`Emphasis.row`, `ViewerInput.row`: the track whose row the pointer
  is in, the nearest track's in the section's gaps, the rail under a ruler included): anywhere in the row, a gap between
  blocks too, they stay; they go when the pointer leaves the row. Why: a hover that hits nothing in a gap made them
  blink (the author, 2026-11-21).
- **An open menu keeps its header's actions shown** (`Emphasis.menu`: the address of the header whose menu or submenu is
  open), canvas and DOM parts alike, wherever the pointer is. Why: the pointer leaves the line on its way to the menu,
  and the menu's anchor must not go with the hover. It closes as the menu does (Esc, a press outside, a choice).
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
selection's gestures and the cursor's moves are session state; the cursor's focus (Enter) navigates. A navigator's
`navigator.center` and `navigator.zoom-to` navigate (jumps).

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
- **A block's commands** (genome coordinates): `block.next` and `block.previous` (the `navigate` behavior's, jumps: a
  step of the view's history) go to the next or the previous block of a row in display order from the block holding the
  window's centre, by the switching rule (whole from whole, else the zoom kept) (one with no positions skipped; at the
  ends they say so), keys `Alt+PageDown`, `Alt+PageUp` in `view.keys` (no Shift: Alt+Shift types a character on macOS,
  EC-155); `block.reverse` (the `editing` behavior's, a patch, one undo step; `address` the block, a node's or a derived
  block's keyed form) turns a block over: its node's `direction` set to `reversed` or unset, a derived block's override
  (core's `blockOverridesOf`: a `block` child whose `of` names its sequence, an alias too, or whose id is its name; the
  last, the one that applies) likewise or one inserted (its id from the sequence's name), the rest still derived;
  `section.include` (a patch) sets a section's `include` (`include`: names, patterns, `{ min_length }`), or removes it
  (`include` absent or null: every block). Reached from a block's menu (a secondary click on a block's label or band,
  Shift+F10 at the cursor: "Reverse chr2", "Show only chr2" on a section of assemblies, "Show all" where an `include` is
  set: built-in header entries `on: 'section'` reading `HeaderFacts.block`, `assemblies`, `included`).
- **A locus** (`goToLocus(viewer, text, { section })`, `locus.ts`: what a host's locus box calls): the text as a locus
  first (`(locus :text …)` under the section, the first section of assemblies when none is named:
  `chr17:7,668,402- 7,687,550`, `17:7668402`, a name or an alias), its Location focused (a jump); a text no sequence is
  named (`seq.query.locus.name`) falls back to the sources' names (core's `findNames` over `store.named()`: a GFF3's
  `Name`, `ID`, `Alias`, EC-228), the first hit a section shows focused, its system read through the assemblies' names
  rule; a hit only on a sequence no section shows is `hidden` with why, never focused unseen; an aborted search
  (`signal`) settles as `cancelled`, never an error; else `none` with why.
- **`block.form`** (the `editing` behavior's, a patch, one undo step; a command only, in no built-in menu: a host's
  `on: 'element'` header entry may name it; no default key): sets a block's `form` (`whole`, `residues`, `atoms`;
  `address` the block, `form`). A block a section of a complex derives from its units has no node: the patch inserts one
  `block` child of that unit's id (a merged copies block's: its first unit's) with the form, an override of that block
  alone (the core's Resolved view, A section of a complex); the others stay derived, units added later listed. Undo
  removes the child.

- **The rows' commands** (the `rows` behavior's, `behaviors/rows.ts`, on a track's rows layer, the first when none is
  named): folds are session state (the band's `folds`: no undo step): `clades.toggle` (a hand's move, no history step),
  `clades.collapseAll` (every clade but the root), `clades.expandAll`, `clades.collapseTo { depth | height }` and
  `clades.reveal { rows }` (opens the clades holding them) are jumps, so back restores the folds and the rows' window
  together; `clades.keep` writes the folds shown into the layer's `collapsed` (one patch); `rows.pin` (toggles: a pinned
  row is unpinned), `rows.hide`, `rows.show` (the hidden-rows chip's) and `rows.order` are patches of the layer's
  params, read from the spec as it is now; `rows.order-next` (`Shift+O` in `view.keys`, scoped to the canvas) is the
  keyboard's way to a column title's click: each press takes the next step down the columns in the order they are drawn
  (ascending, then descending, then the next column, round to the first) and says it ("Ordered by label, descending").
  `rows.tree` (`Shift+T`, and a click on the tree's toggle, `tree-toggle`) opens a tree a narrow frame collapsed to its
  toggle (the band's `tree`, session state: the tree whole, the frame scrolling) and folds it back again. Why a key and
  no focusable title: the titles are drawn on the canvas, and a key reaches them without a DOM overlay; `rows.zoom-in`
  and `rows.zoom-out` (`]`, `[` in `view.keys`) double and halve `rowScale` (1/64 to 40 px a row, fixed: 1/64 shows
  10,000 rows in about 160 px, 40 is twice a letter row) about the band's middle row. A fold keeps the row at the band's
  top there (a fold above the window moves it up, never off). **An edit of `collapsed`** (a patch, an undo) resets the
  session's folds to it; other edits keep them. The wheel's vertical delta over a rows layer's own band (its cells, its
  header) scrolls its rows (`bandTop`, px over the row height, clamped), its horizontal pans the columns; over another
  layer's `row` hit it pans as any. Clades fold only in tree order (`when` says). `rows.filter` waits for a `filter`
  param of the rows layer. **A row column's title** (a `column-title` hit, its `state` the field): a click runs
  `rows.order { field, ascending }` (one patch), a click on the title of the field the rows are ordered by ascending
  runs it `descending`, a third click on it **`rows.unsort`** (the order before the first sort by a column, remembered
  by the behavior: `tree`, `data`, or none; the author, 2026-11-21: "it doesn't seem to be possible to unsort"), each a
  patch (one undo step, so undo restores the sort); the order read from the spec as it is now, before the double-click's
  drop, so a fast second click steps again. **`viewer.rowSort`** (`ReadonlyState<RowSortEntry[]>`: layer, track, field,
  direction, read from the spec; its `listen` the change event) is for the host, which draws the control (the stories'
  toolbar has a "Clear sort", enabled while sorted): the view draws no clear-sort chip of its own, as with hidden
  tracks. Its cells hit as their row's header (select and hover as the label does), naming the value in the hit's
  `state` (the hover keeps it for the targets that read it, `row`, `column-title` and `columns-hidden`, never another's,
  a logo's letter, so a move within a position stays one hover: a cut value reads whole in the tip pass under the
  pointer alone, and the card names it). **The hidden-columns chip** (a `columns-hidden` hit): its card lists each
  hidden column and why ("no room at this width"), announced the same; a click does nothing.

## Legends

`legends.ts`: `LegendModel`, derived when the view settles (never per frame: `viewLegends` over the settled layout,
merged), the entry under the hover (`updateHover`: a mark's value), and the reader's substitution (`cvdSafe`, from
`theme.cvdSafe`, part of the frame's key). `legendRows(viewer, { of, omit, maxItems })` is what a view of them shows (a
size entry's `size` is the channel's value in the layer's own units; the part draws circles by ratio to the row's
largest; titles by `legendTitle`) (`of` an address or a list, then `omit` removes those; `selectLegends` is the one
rule): a row a legend, its entries with accessible names (`<title>: <value>`), the entries past `maxItems` counted. Why
derived on settle: a legend that changes during a pan reads as noise, and its cost stays off the frame.

- **The legend behavior's commands** (`behaviors/theming.ts`, patches unless said): `layer.channel` (a channel's value
  or form set or unset, on one layer or several in one patch; on one layer's color its choices are a select param,
  `scheme`, the layer's fitting schemes, its default the layer's scheme (or, on a field form, the scheme whose palette
  it has): read when no `value` is given; its `preview` the gesture's preview; a scheme chosen on a color that reads a
  field is that scheme over the layer's own field, its scale, domain, mapping and listed values kept (a field form with
  more than a name and a palette, or listed values, has the scheme's palette put in it), so the first entry brings the
  original colors back exactly, and a letters or logo layer of no known alphabet has no fitting scheme: the menu has no
  Color entry where no palette is read), `scale.set` (a scale's params; `null` reopens a domain's end),
  `theme.value-color` (one value's color in a field form's palette, the listed values keeping theirs; refused on a
  scheme, rules or any other form: never replaced), `layer.filter` (a value's items hidden or shown again: one patch of
  the layer's `data`, a trailing `(not …)` of the values hidden; no value: cleared, every layer of `addresses` in one
  patch, as a click on the track's "n filtered" chip). The user's filter is only that trailing `(not …)` (`dataFilter`):
  what leads it, the annotation and a leading run of item forms (`<annotation>/(item :category Domain)`, a layer a
  feature type), is the layer's base, kept when the filter is set, toggled or cleared, its rejected rows not counted as
  filtered (core's Render, Filters: `splitDataForms` is the one parser). A legend entry's items are the base's forms and
  `(item :<field> v)` in one query (`itemsAddress`: one item form holding every param, or an `and` when the base has
  `and` or `or`), `legend.select` (session: the items whose field has the value; Cmd/Ctrl toggles), `legend.hover`
  (session: a highlight named `legend`), `theme.cvd-safe` (session). **Both ways** (core's Emphasis, the hover
  contract): an item hovered outlines the entry it reads as (`updateHover`: its painted value), and the entry hovered
  lights its items in view; where its address names none the store holds (a source's rows: a VCF's alleles; a derived
  table), the items the layout drew in its value (`LegendModel.entryItems`: their hit marks, never items read and out of
  view) are lit and selected instead. A legend over several layers (one entry for a light and a heavy chain's CDRs)
  lights and selects the entry's items on every layer it lists. An entry **has an action** (`LegendEntryRow.action`,
  `LegendModel.names`) when an address names its items (`addressable`, the behavior's) or items in its value are drawn
  at the last settled frame (`LegendModel.drawn`: the part redraws when it moves, a pan or a tile): only then is it a
  pointer target with a click, a hover and Enter; a swatch alone keeps the default cursor. `names` is a lookup and
  `entryItems` is read once a settled frame (a hover again makes nothing). Swept by the hover sweep: every entry a
  hovered item reads as has its action and lights every item it names in view, over its positions; an entry without one
  lights nothing. `fittingSchemes(layer)` is what a menu offers (the schemes whose subject and measures fit the layer);
  `preview(layer, channel, value)` shows one hovered (a preview, ended by `layer.channel` or Esc).
- **The Legend part's groups run one a line** (a column, 12 px apart; each title muted, 600, 11 px: grouping is seen
  without reading) and **its channels row keeps its place**: the row is outside the listbox (a listbox holds options),
  one line (`nowrap`, the first two channels named, the rest a `+n` chip), and `channels: reserve` (the default) holds
  it whether or not a channel is active: empty, it holds a hidden copy of a channel's row, so the height is the real one
  at any text size; `collapse` has the row only while one is (an unknown value throws).
- **The Legend part's options** (`mountLegend(element, viewer, options)`): `of`, `omit`, `maxItems` as above, `form`
  (`full`, the default; `compact`), `orient` (`vertical`, the default: legends stacked; `horizontal`: in a row,
  wrapping) and `channels`; an unknown `form`, `orient` or `channels` throws. **The compact form** is one line a legend
  (`nowrap`): its title, then its swatches unlabelled (each entry's `title` and accessible name carry the label), a ramp
  a 40 px bar with its domain's ends as text, steps joined cells, `+n` past `maxItems` (its title "n more"; the full
  form writes "+n more"); the same entries, roving focus, gestures and hover outline, only the drawing differs.
  **Focus**: items have ids and the listbox's `aria-activedescendant` follows the roving focus, so a screen reader says
  the entry's name; in the compact form it is also shown beside the groups while the part has focus. An `of` naming no
  layer gives an empty part, nothing said; rows come in the view's order, never `of`'s. A compact group wraps in a
  horizontal band and is clipped in a vertical one. **Hiding has two owners**: the figure's, the layer's `legend: none`
  (core, Themes), and the host's, `of` and `omit`; no session command hides a legend, and editing stays the
  application's. Why: a legend derives from the spec, so a third owner would make "why is there no key" unanswerable.
- **Hidden tracks** (the author, 2026-11-21: "hiding a track moves the UI around … give this to control of the app
  layer"): the view draws nothing for a track the spec hides (`hidden: true`): no section line, no chip, no stub, so a
  hide moves only the tracks below it, by its height (core's Track headers); a section whose tracks are all hidden
  collapses fully (`collapsedSection`: no height, gap, line or rail; the navigator rails, the cursor and the toolbar's
  commands skip it), and a jump to what is only on hidden tracks (a track or its nodes, or an item key every layer
  drawing it lies on a hidden track of, resolved over the whole spec) says "On a hidden track: show it to see it" and
  moves nothing; a key also drawn on a shown track jumps there. The viewer gives the application what it needs:
  **`viewer.hiddenTracks`** (`ReadonlyState<HiddenTrackEntry[]>`: address, label, section, each section's in spec order,
  read from the spec, so a section with none left is listed; its `listen` is the change event, told when the list
  changed), the command **`track.show`** (`address`, `addresses` or `all: true`, with `all` limited to the section
  `address` names when it names one: every one a single patch, one undo step) and `track.hide`, the menu's **"Show
  hidden ▸"** (`seq:hidden-tracks`, a menu-only section entry: a flyout of each hidden track, then "Show all" when there
  is more than one) so the keyboard and the palette reach them, and **`mountHiddenTracks(viewer, element, options?)`**
  (`mount-hidden-tracks.ts`, beside `mountLegend`): one `@molstar/ui` toolbar button (eye-off) with a menu of each
  hidden track (a section's name over its group when more than one section has some) then "Show all"; always mounted,
  disabled with "No hidden tracks" while there are none (a host's layout never shifts on a hide); `menuSide` its option;
  taken away by what it returns. The host places it outside the view (a toolbar, a panel): the default viewer's
  Storybook composition mounts it in its toolbar (`toolbarWithHiddenTracks`), the minimal viewer none. Why: a library's
  view never decides where a host's controls go, and chrome that appears on a hide moves what the user was looking at; a
  host that mounts no list leaves its hidden tracks reachable only by the menu, the palette and the spec. Every entry is
  listed by a name (`hiddenTrackName`): its label, its id without one, and for an empty label (a ruler's, the rows': the
  author's) its layer's `label` or `title`, else its type as a word ("Ruler", "Rows"); the list and the menu share it.
  Why: the list has no header to read a label from.
- **The Legend part lists the active highlight channels** after the legends: a row a channel (`data-highlight`) with a
  swatch in the channel's color (`highlightColorOf`, beside the legend rows) shaped by its style, and a Clear, which
  runs `highlight.clear` with the channel's name. The rows are redrawn only when the channel names, a channel's style or
  its color change (a channel moved to other positions redraws nothing), and on a token change (the colors); a focused
  Clear keeps its focus across the redraw. A style that is not `tint`, `band` or `outline` is a tint, reported once
  (`seq.viewer.highlight.style`); an unset style is a tint, so setting `tint` is no change.

## Bridges

`bridge-view.ts`, the viewer's bridge chunk (across-sections, Moving two sections joined by a bridge; Hover, selection
and the keyboard): loaded when a spec holds a bridge (`models/bridge-hooks.ts` registers it with core's `withViewKind`,
so the model waits for it as for render's chunk) and set as `BridgeView.parts`; `ViewModel.get`, `panPx` and `zoom`, the
input's press and wheel, the cursor's crossing, the describer and the viewer's resolved view call it, one optional call
each. Why lazy: a viewer that never shows a bridge pays for none of it (eager: the hooks alone, about 0.5 kB).

- **One move a section a gesture** (`carriedSections`): a pan or a zoom asked of a section is carried through the view's
  `lock: px` bridges, each section moved once by its mover (up its attached follows to the first leader, then its scroll
  group): a lock between sections a follow or a group already moves moves nothing more, and a cycle ends where it began.
  A carried section pans by the same CSS px (`ownPanPx`, its own follow and group only) and zooms by the same factor
  about the same page x (`ownZoom`). Why px: two assemblies' scales and gaps differ; the reader's slant is kept. The
  carrying is kept per resolved view and section, found again only when a follower was detached or attached since: a
  drag's moves allocate what the sides' own pans do (the storybook's `across.quick.ts`).
- **The band's pointer**: a press or a wheel over a band (`bandAt`, between its sides) is the bridge's: a drag pans both
  sides (target `empty`: it never selects), a wheel pans or zooms both about the pointer's x, whatever the lock, through
  the same closure; `ViewModel.get(bridge)` is its upper side's viewport (the wheel's px and anchor).
- **`bridge.align`** (`navigate`, a jump with the location aligned; `Shift+A`, the command's shortcut, `view.keys`
  rebinds it; checked against every built-in key and binding): the side named (`section`: a key's is the pointer's),
  else the keyboard cursor's, else the leading section (the first in the stack a bridge touches): never none, so in a
  chain A–B–C it is one consistent step (never C to B's old window); each neighbour across a bridge touching it brought
  to the hull on the far system holding the most aligned length among the ribbons with their near end in the window
  (`alignRange`), all in one step; a bridge named (`bridge`, `to`) one way. Nothing in view: "Nothing aligned in view",
  nothing moves. A stub (its far end has no place) is never read by an align, select across or the cursor.
- **`bridge.selectAcross`** (`session`, no key): the far ends of the drawn ribbons with an end in the selection, across
  the bridge named or every bridge touching the side, added in one change.
- **`bridge.lock`** (`patch`, one undo step, no key): toggles the bridge's `lock` between `px` and unset.
- **The cursor across a band** (↓ and ↑ from a side): onto the image of its position through the topmost drawn ribbon
  covering it (a reversed ribbon read backwards); none: the view's link, else the far side's window, as between any
  sections.
- **A ribbon's card**: "peach chrX:101–200 ↔ grape chrY:201–300", its strand and its aligned length (`bp`, `kb`, `Mb`);
  the sides by their ids; a stub's far end ", not shown"; a bin pair's count and its bins' width ("37 pairs · 16 kb bins
  · zoom in to see them"), never an aligned length. Built from the hit ribbon alone (`BridgeLayout.ribbon(key)`), never
  the list. The commands' change kinds are the commands table's (`commands.ts`).
- **ARIA words** (`bridgeAria`): each bridge its `label` and the ribbons in view ("3 ribbons between peach and grape in
  view"; none: "no ribbons in view"); the mount sets them as visually hidden `role="group"` nodes, labelled, beside the
  canvas each time a frame settles (`aria`).

## Structure link

`link-structure.ts` (`@molstar/seqx-viewer/link-structure`): `linkStructure(viewer, side, options?)` wires a viewer to a
structure viewer its host translates for, both ways, and returns the disposer. A library over the viewer's public models
and the core's translation (`@molstar/seqx-core/structure`, Structures): no behavior, no port, nothing in the
composition, so a host that never calls it carries none of it. Why a helper: what it adds over the viewer's events
(echoes, throttling, clearing) is per link, and a viewer may have none or two; a host with other needs writes its own
over the same translation.

- **The side** (`StructureSide`, plain data, the host's): `highlight(target | null)`, `select(target | null)`,
  `focus(target, { durationMs })`, and optionally `onHover`, `onSelect`, `onFocus` (each returning its unsubscribe) and
  a `signal`. A target is a Location on units (structured-cloneable).
- **Pending** (EC-169): a hover or a selection whose `toUnits` says `pending` and reaches nothing yet is held, never
  sent as nothing (the last sent stays); the store's nodes arriving (`nodeVersion`) sends it again, then as it resolves.
  Why: a host's per-chain maps arrive after the view, and a cleared highlight there would be a lie.
- **Out**: on each hover (`hover.version`), the hover's Location through `toUnits` under the hover's limits, sent when
  its image changes (an unchanged image sends nothing; an empty one, an unobserved residue, is `null`, EC-118); the
  selection's `shown`, one message a frame (a frame callback: a drag moves it every pointer event), under the
  selection's limits; a jump naming a Location (`view.lastJump`, read as `jumps` moves: a focus, a reveal, and back and
  forward when the step kept one and `followJumps`) as `focus` with the jump's duration (the camera and the sections on
  one clock). A fit names none: nothing sent.
- **In**: the side's hover as a hover of a Location (`HoverModel.setLocation`) owned by the link: drawn projected, no
  card; the pointer taking the hover takes it, and the side's `null` clears only the link's own (a stale `null` does
  nothing, EC-141); a side hover carrying its time (`at`) older than the pointer's last own hover is dropped; the side's
  selection through `fromUnits` (a unit's atoms kept where a section or block shows its atoms form, collapsed onto its
  residues elsewhere, unit by unit) into `selection.set`, then revealed (`revealOnSelect`: `revealLocation`, the zoom
  kept; `focus`; `none`); a `null` clears it; its focus as `focusLocation` (one jump).
- **Echoes** never go back: the link's own hover is not sent (it only ends the pointer's hover sent before it, with
  `null`, sent from inside the side's `onHover` call); a selection equal (`locationsEqual`) to the one last set from the
  side is not sent; a jump made while a message from the side is applied (its selection's reveal, its focus) is not
  sent. So a host keeps Seq\*'s highlight apart from its own hover: a `highlight(null)` arriving while it reports its
  own hover clears what Seq\* lit earlier, never the side's hover. Of the jumps, only those naming a Location go out
  (`focus`, `reveal`, and back and forward over them); `fit`, `zoom` and `move` name none.
- **The model** translated against is the spec's static model with the store's current overlay
  (`StoreModel.currentOverlay`; an empty one before the first frame reads as the model alone, core's Structures).
- **Clearing**: the disposer, the side's `signal` or the viewer retiring (`retiredSignal`) ends the link once: its
  subscriptions and a waiting selection message dropped, `highlight(null)` and `select(null)` sent for the directions
  sent out, the link's hover cleared. A signal aborted before the call ends the link at once.

| Option           | Default                                 | Why the default; the alternatives                                    |
| ---------------- | --------------------------------------- | -------------------------------------------------------------------- |
| `complexes`      | every complex, no `via`                 | the spec names no host's binding; `{ complex, via? }[]` (Paths, Via) |
| `hover`          | `both`                                  | the precedents link all three both ways; `out`, `in`, `none`         |
| `select`         | `both`                                  | as `hover`                                                           |
| `focus`          | `both`                                  | as `hover`                                                           |
| `revealOnSelect` | `reveal`                                | a selection from the structure keeps the zoom; `focus`, `none`       |
| `followJumps`    | `true`                                  | back and forward bring the camera; `false`                           |
| `limits`         | `DefaultSyncLimits` (`hover`, `select`) | a hover a frame is smaller than a selection; any, `over: skip`       |
| `ligand`         | `whole`                                 | a ligand's position is its residue; `atoms`                          |

## Configuration

`config.ts`: `DefaultViewerConfig`, one object of the shape the models read (`viewerConfig(over)` for one group's items
over the defaults; `view.keys` merges by command id). Each item is a choice with its recommended default; all of it is
read when the viewer is made. As items in a spec: `viewerItem('group.item')` (typed by key), read by
`viewerConfigOf(pairs)`; the core's items (`SeqCoreConfig`) are the same pairs. `ParamDefinition`s come with the
settings form.

| Item                                                             | Default                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 | Why the default; the alternatives                                                                                                                                                                                                                                                                                          |
| ---------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `foundations.preset`                                             | `comfortable`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           | the UX foundations' spacing; `compact`                                                                                                                                                                                                                                                                                     |
| `density.letterMinPx`                                            | 7                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       | the least CSS px a position needs for its letter (a number: no reason recorded)                                                                                                                                                                                                                                            |
| `density.dotMinPx`                                               | 2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       | the least px a position needs for a dot (no reason recorded)                                                                                                                                                                                                                                                               |
| `bindings.preset`                                                | `default`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               | the interaction grammar (Bindings); `genome-browser`: a drag on data pans                                                                                                                                                                                                                                                  |
| `store.cacheBytes`                                               | 64 MB                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   | the loaded bytes the store keeps, least recently used evicted, never what is in view (the size: no reason recorded)                                                                                                                                                                                                        |
| `store.inFlight`                                                 | 6                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       | requests in flight at most (no reason recorded)                                                                                                                                                                                                                                                                            |
| `store.tileBins`                                                 | 1,024                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   | bins per tile (no reason recorded)                                                                                                                                                                                                                                                                                         |
| `hover.projectLinked`                                            | true                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    | the hover shows on every linked section; false: only sections of the hovered system                                                                                                                                                                                                                                        |
| `hover.card`                                                     | true                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    | the hover card (the `hoverCard` behavior) shows what the pointer rests on; false: none (the live region still says it)                                                                                                                                                                                                     |
| `hover.cardItems`                                                | 8                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       | the items the card lists at most (a place's features, a group's members): past it the first ones and a row `+k more`; the data keeps every one (`HoverCard.more`), the live region says them all                                                                                                                           |
| `hover.tooltipLines`                                             | 3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       | the lines a card row's value wraps to within the card's capped width, then an ellipsis (a tooltip's text, a long INFO value)                                                                                                                                                                                               |
| `hover.closeCodons`                                              | true                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    | a nucleotide lights its whole codon and its residue (EC-30, EC-144); false: the position alone                                                                                                                                                                                                                             |
| `headers.submenuRestMs` (a behavior's param, in `SeqViewerSpec`) | 300 ms                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  | a submenu's flyout opens once the pointer rests on its row: long enough that a pointer passing over it down the menu opens nothing, short enough to read as hover; a click, → and Enter open it at once; the same rest on another row closes it; 0 turns resting off                                                       |
| `selection.holdNodes`                                            | true                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    | a selection's dynamic nodes stay while selected, so a selected transcript survives a scroll; false: only while shown                                                                                                                                                                                                       |
| `selection.extendScope`                                          | `range`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 | an extend changes the range the anchor is in, the others kept; `selection`: the whole selection becomes that range                                                                                                                                                                                                         |
| `selection.acrossBreaks`                                         | `shown`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 | what a drag or an extend across a collapsed break selects: the positions shown (what is seen); `span`: the hidden stretches too (a gene's whole span to export)                                                                                                                                                            |
| `highlight.color`                                                | `token:cat-5`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           | a channel's color when `highlight.colors` names none: a hex literal or `token:<name>`                                                                                                                                                                                                                                      |
| `highlight.colors`                                               | `{ 'link-hover': 'token:hover' }`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       | colors by channel: a link's hover in the hover's ink                                                                                                                                                                                                                                                                       |
| `cursor.keys`                                                    | `ArrowLeft`, `ArrowRight`, `ArrowUp`, `ArrowDown`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       | back, on, the section above, below; Shift with back or on extends from it                                                                                                                                                                                                                                                  |
| `reveal.cursor`                                                  | true                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    | a cursor moved out of view brings its section to it (the layout's least move); false: it moves unseen                                                                                                                                                                                                                      |
| `reveal.linked`                                                  | `text`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  | a selection made elsewhere scrolls a text section showing none of it to its first part (a text hides most of itself; a window shows the projection where it is); `all`, `none`                                                                                                                                             |
| `viewport.maxPxPerPosition`                                      | 32                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      | the zoom's upper bound (no reason recorded)                                                                                                                                                                                                                                                                                |
| `viewport.groupEqualAxes`                                        | true                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    | sections on equal axes share one viewport (EC-113); false: one each                                                                                                                                                                                                                                                        |
| `viewport.onResize`                                              | `keep-zoom`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             | a section hidden for lack of width comes back as it was (EC-112); `keep-range`                                                                                                                                                                                                                                             |
| `navigator.form`                                                 | `rail`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  | a navigator on every linear section at no row's cost (the author's choice); `none` draws no rail                                                                                                                                                                                                                           |
| `navigator.rowsRail`                                             | `always`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                | a rows band's rail beside a declared overview too: the overview a picture, the rail the precise control (the author's choice); `unlessOverview` (replaced by a navigator over its rows), `none`                                                                                                                            |
| `navigator.positionsRail`                                        | `always`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                | a linear section's positions rail beside a declared navigator too: the navigator a picture, the rail the precise control (the author's choice, 2026-11-21); `unlessNavigator` (replaced by a navigator over its positions), `none`                                                                                         |
| `navigator.railPlace`                                            | `ruler`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 | where the positions rail sits: `ruler` (below the ruler's border), `bottom` (a row under the last track), `navigator` (by a declared horizontal navigator of the section, else `ruler`), `pinned` (a row at the frame's foot, in reach at any scroll; none in an export); or a list of them, a rail each on the one window |
| `navigator.dragOutside`                                          | `sweep`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 | a drag from outside the window sweeps a region and zooms to it (IGV, UCSC, JBrowse); `pan` (centre there, then pan)                                                                                                                                                                                                        |
| `navigator.dragStep`                                             | `far`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   | back undoes a trip, not a nudge (a default to be confirmed by use); `never`, `always`                                                                                                                                                                                                                                      |
| `navigator.doubleClick`, `navigator.zoomPx`                      | `toggle-fit`, 48                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        | the whole, then back in, scale-free; `fit`, `none`; any px                                                                                                                                                                                                                                                                 |
| `navigator.direction`                                            | `display`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               | a drag right moves the view right (decided); `genomic` (positions up whatever the section's direction)                                                                                                                                                                                                                     |
| `navigator.zoomStrip`                                            | `{ minPx: 48, heightPx: 16, overBlocks: 'always' }`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     | a one-system span under `minPx` gets a zoom strip `heightPx` more (a mark above it); over blocks `overBlocks`: `always` the chromosome is the control row, `narrow` only for a window under `minPx`, `never` one row; any px                                                                                               |
| `follow.followerMoves`                                           | `lead`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  | a move asked of a follower moves its leader through the follow, so every section stays in step (the author); `ignore`; `detach` (its own viewport until attached again)                                                                                                                                                    |
| `follow.cyclicArc`                                               | `wrap`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  | a follower's range over a cyclic origin kept as the wrap, two pieces; `hull`                                                                                                                                                                                                                                               |
| `view.historyLimit`                                              | 50                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      | steps back the history of jumps keeps, the oldest dropped first (no reason recorded)                                                                                                                                                                                                                                       |
| `view.backOnFocused`                                             | true                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    | a double-click on what is already focused goes back; false: nothing                                                                                                                                                                                                                                                        |
| `view.jumpMs`                                                    | 240                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     | a jump's animation (focus, fit, reveal, back, forward, a step zoom; the length: no reason recorded); 0: a cut                                                                                                                                                                                                              |
| `view.motion`                                                    | `platform`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              | jumps and stories alike: animated unless the platform asks for reduced motion; `cut`: always a cut; `animate`: always animated                                                                                                                                                                                             |
| `story.transition`                                               | `{ duration: 240, easing: ease-in-out }`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                | a step's transition unless the step says: the base duration of a jump (no reason recorded for its length)                                                                                                                                                                                                                  |
| `story.linger`                                                   | 2,000 ms                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                | how long a played step is shown before the next (no reason recorded)                                                                                                                                                                                                                                                       |
| `story.view`                                                     | `{ path: auto }`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        | pan when the windows overlap or touch, else out, across and in; `pan`, `out-across-in`                                                                                                                                                                                                                                     |
| `story.color`                                                    | `{ space: oklab }`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      | perceptual: a blend reads evenly; `srgb`                                                                                                                                                                                                                                                                                   |
| `story.cut`                                                      | `{ at: 0.5 }`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           | a cut switches mid-window, where the change is least seen (no reason recorded)                                                                                                                                                                                                                                             |
| `story.stepStart`                                                | `spec`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  | a patch step applies over the spec: any step reached directly; `from` per step                                                                                                                                                                                                                                             |
| `story.store`                                                    | `patch`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 | the encoder stores a patch over the spec (it follows edits); `full`, `smaller`                                                                                                                                                                                                                                             |
| `story.maxChain`                                                 | 8                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       | deltas at most between a step and its keyframe: reaching a step costs at most that many applies                                                                                                                                                                                                                            |
| `view.keys`                                                      | `view.back`: Alt+ArrowLeft, Backspace; `view.forward`: Alt+ArrowRight; `view.fit`: 0; `view.zoom-in`: = and +; `view.zoom-out`: -; `view.pan-left`, `view.pan-right`, `view.pan-up`, `view.pan-down`: Shift+Alt+Arrow; `selection.copy-letters`: Mod+c; `selection.bookmark`: Mod+d; `section.reverse`: Shift+R; `section.collapse`: Shift+C; `rows.zoom-in`: ]; `rows.zoom-out`: [; `cursor.page-down`: PageDown; `cursor.page-up`: PageUp; `cursor.scroll-start`: Home; `cursor.scroll-end`: End; `cursor.first-row`: Mod+ArrowUp; `cursor.last-row`: Mod+ArrowDown; `rows.order-next`: Shift+O; `rows.tree`: Shift+T; `help.keys`: ? | the view's commands' keys by command id, the first the shortcut, and the one shown in the toolbar's tooltip (Cmd/Ctrl+D bookmarks the selection: the browser's own bookmark key, taken only while the canvas has focus); an override merges by id (`[]`: none)                                                             |
| `ui.canvas` (in `SeqViewerSpec`)                                 | `true` in the default spec                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              | the canvas stack mounted; a headless viewer (no `ui`) draws through its backend alone                                                                                                                                                                                                                                      |
| `ui.history` (in `SeqViewerSpec`)                                | `uiHistory` in the default spec                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         | a mounted viewer with no `history` given gets a `@molstar/ui/history` stack from it (the standalone viewer's Undo); a headless viewer never does; not listed, not in the bundle                                                                                                                                            |
| `ui.toolbar` (in `SeqViewerSpec`)                                | none                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    | a function mounting the bar above the canvas (`mountToolbar`: Back and Forward, their keys and disabled state); the host places it, as it would a legend: the bar is `@molstar/ui`'s and heavy                                                                                                                             |
| `ui.navigators` (in `SeqViewerSpec`)                             | `navigatorsMount`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       | the host places parts; `mountNavigator` by hand                                                                                                                                                                                                                                                                            |
| `ui.railToggle` (in `SeqViewerSpec`)                             | none (the standalone viewer); `railToggleMount` in an embedded layout                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   | a cell's room is the reader's                                                                                                                                                                                                                                                                                              |
| `ui.legends` (in `SeqViewerSpec`)                                | `legendsMount`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          | a function mounting the Legend part (`mountLegend(element, viewer, { of, omit, form, orient, maxItems })`: the legends of the view or of some tracks, at most 12 entries each, its rows options with accessible names, roving focus, a click selecting a value's items, a hover highlighting them)                         |
| `navigate.wheelZoomStep`                                         | 1.25                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    | the zoom of one wheel notch (100 px of delta) with the zoom modifier, and of a key zoom (no reason recorded)                                                                                                                                                                                                               |
| `navigate.doubleClickZoomStep`                                   | 2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       | Shift+double-click's zoom out about the pointer                                                                                                                                                                                                                                                                            |
| `navigate.keyPanFraction`                                        | 0.1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     | the part of the viewport a key pans (no reason recorded)                                                                                                                                                                                                                                                                   |
| `navigate.coalesceWheel`                                         | true                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    | wheel deltas within a frame applied once, at the frame; false: at each event                                                                                                                                                                                                                                               |
| `navigate.acrossBlocks`                                          | `cross`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 | on a row of blocks (a genome's chromosomes) a pan carries across a gap into the next block (JBrowse's); `stop`: a pan, a zoom and a fit kept in the block the window is in, a jump moving to another (IGV's and UCSC's one chromosome at a time)                                                                           |
| `theme.source`                                                   | `root`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  | a mounted viewer reads its root's CSS (`--msui-*`, `--seq-*`), redrawn when `data-msui-theme` changes; `config`: the core's `theme` item (a headless viewer always)                                                                                                                                                        |
| `theme.cvdSafe`                                                  | false                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   | the reader's color-blind substitution: a scheme not safe under the three deficiencies is drawn as its safe substitute (`theme.cvd-safe` toggles it; a reader's setting, never the spec's)                                                                                                                                  |
| `canvas.height`                                                  | `content`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               | the layout's content, tracks below the page not laid out; a number of CSS px                                                                                                                                                                                                                                               |
| `canvas.maxHeight`                                               | 2,048                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   | the tallest the canvas grows with `content`, CSS px (no reason recorded)                                                                                                                                                                                                                                                   |
| `canvas.dpr`                                                     | `device`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                | the device's pixel ratio; a number for an export or a test                                                                                                                                                                                                                                                                 |
| `canvas.settleMs`                                                | 150                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     | after the last move the view counts as settled: open domains computed, lanes packed (no reason recorded)                                                                                                                                                                                                                   |
| `pointer.rehit`                                                  | true                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    | a resting pointer is hit again once a frame in which the layout was rebuilt, so the hover follows the data (EC-143); false: with the pointer only                                                                                                                                                                          |
| `pointer.pinchFromCtrlWheel`                                     | true                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    | a trackpad's pinch arrives as a Ctrl+wheel with small deltas: the `pinch` gesture; false: a wheel with Ctrl                                                                                                                                                                                                                |
| `pointer.pinchSpeed`                                             | 0.01                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    | a pinch's zoom is `exp(-deltaY · pinchSpeed)`                                                                                                                                                                                                                                                                              |
| `pointer.pinchMaxDelta`                                          | 50                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      | a Ctrl+wheel with a larger delta (px) is a mouse wheel, not a pinch (not on macOS)                                                                                                                                                                                                                                         |
| `pointer.lineHeight`                                             | 16                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      | px a line of a wheel in line mode scrolls; a page is the viewport's own axis                                                                                                                                                                                                                                               |
| `address.limit`                                                  | 10,000                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  | the most nodes or elements an address names: more is an issue, never a truncation                                                                                                                                                                                                                                          |
| `address.maxDepth`                                               | 64                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      | how deep below the root resolving walks                                                                                                                                                                                                                                                                                    |
| `address.refDepth`                                               | 8                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       | how many refs a chain passes through                                                                                                                                                                                                                                                                                       |
| `links.prefix`                                                   | `seq-`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  | this viewer's keys in a link; other keys are another viewer's                                                                                                                                                                                                                                                              |
| `links.maxLength`                                                | 8,192                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   | the longest link read, in characters                                                                                                                                                                                                                                                                                       |
| `links.maxActions`                                               | 16                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      | the most actions one link runs                                                                                                                                                                                                                                                                                             |
| `links.highlightChannel`                                         | `link`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  | the channel a link without `channel` highlights on                                                                                                                                                                                                                                                                         |
| `links.hoverChannel`                                             | `link-hover`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            | a link's hover has its own channel, never the pointer's (EC-182)                                                                                                                                                                                                                                                           |
| `copyLink.prefer`                                                | `['ref', 'ids']`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        | the forms an address is copied in: a ref naming exactly it, then ids and forms                                                                                                                                                                                                                                             |
| `copyLink.numbering`                                             | `ordinal`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               | labels are copied in it; one that does not read back falls back to it                                                                                                                                                                                                                                                      |
| `copyLink.action`                                                | `focus`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 | the action a copied link carries; `reveal`, `select`, `highlight`                                                                                                                                                                                                                                                          |
| `revisions.keepRecent`                                           | 64                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      | the recent roots kept by revision, for a rebase of a stale patch                                                                                                                                                                                                                                                           |
| `edit.rebaseStale`                                               | false                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   | a stale patch is refused (EC-81): a merge the writer did not see is a decision; true: rebased through `merge3`                                                                                                                                                                                                             |
| `preview.dropped`                                                | `issue`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 | a preview that no longer applies is dropped with a warning in the spec's issues; `silent`                                                                                                                                                                                                                                  |
| `bounds.maxBytes`                                                | 64 MB                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   | the spec's size (Language, Bounds in core's `_spec`)                                                                                                                                                                                                                                                                       |
| `bounds.maxDepth`                                                | 256                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     | levels of nodes and nested values (EC-76)                                                                                                                                                                                                                                                                                  |
| `bounds.maxQueryDepth`                                           | 64                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      | nested forms in an address string                                                                                                                                                                                                                                                                                          |
| `bounds.maxColumnLength`                                         | 2^31 − 1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                | a column's length                                                                                                                                                                                                                                                                                                          |
| `bounds.maxNodeIssues`                                           | 100                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     | issues kept a node, the rest counted (EC-86)                                                                                                                                                                                                                                                                               |
| `bounds.maxPlaces`                                               | 10                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      | places kept a code, the rest counted (EC-86)                                                                                                                                                                                                                                                                               |
| `bounds.themeDepth`                                              | 4                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       | a theme form's nesting (`seq.lang.theme.depth`): rules over rules read, a deeper tree is a mistake                                                                                                                                                                                                                         |
| `bounds.tooltipText`                                             | 1,000                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   | a tooltip's template, characters (`seq.lang.tooltip.text`); a formatted text is cut there too                                                                                                                                                                                                                              |
| `bounds.tooltipRows`                                             | 1,000,000                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               | the rows of a tooltip table read whole (addresses, MVS's schemas)                                                                                                                                                                                                                                                          |
| `bounds.maxBlockValues`                                          | 2^24                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    | the values (rows times columns) a table of an answer may hold, held by the store                                                                                                                                                                                                                                           |

The spec's bounds are held by the viewer's patches, previews, story steps and whole roots (`SpecModel.apply` refuses
with `'out-of-bounds'` and `boundsIssues`; `viewer.load` rejects with the issues); the table bound by its store. The
nested sources' depth is the resolve config's (`maxNestedDepth`, 4).
