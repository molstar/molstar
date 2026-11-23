# Stories: states, transitions and an animation the application drives

## Intent

The author, 2026-11-18: "i always looked at "Stories" just at sequence of snapshots, compared to molviewspec, we added
the "diff patches" to make them more efficient. so my goal here is being able to plan an animation, and make everything
looks nice, smooth, and consistent"; then "ultimately it would be good if these "animations" were also controllable by
ports. basically the animation loop will be owned by the "application" and it will decide what to apply and when. we
"just" need to provide a consistent library."

A story is a sequence of view states over one spec, like slides. A step is stored as a full state or, where that is
smaller or should follow the spec, as a patch: the form is storage, not what a step is. Seq\* computes each state, plans the transition between two states from their diff, and samples it at any time;
the application's clock decides what is shown and when. Seq\* runs no animation loop of its own.

## Precedents

- MolViewSpec: a story is a list of full snapshots, each with a transition duration and a linger time; the viewer
  plays them. Here a step is stored as a patch over the spec, and the player is the application.
- [changes.md](changes.md): the history port, by which the host owns the undo stack. The animation port is its twin:
  the host owns the clock.
- [ux-foundations.md](ux-foundations.md) section 7: jumps animate, the hand is followed directly, data does not animate,
  reduced motion cuts.

## Decision

### A step is a full state, stored as a patch where that makes sense

The author, 2026-11-18: "not every step needs to be stored as a patch. some steps might be full states too. the point
is that you can store it as a patch when it makes sense."

- **The `snapshots` container** ([composition.md](composition.md)) holds the story's steps in order. A step has an id
  (`/snapshots/<id>`), a label, its content, its view, and its timing.
- **A step's content is one of two forms**, the writer's choice per step:
  - **a patch** (`patchLanguage(SeqLanguage)`, no `base`): small when a step changes little, and it follows the spec:
    an edit of the spec reaches every step stored this way;
  - **a full state**: a spec tree of its own, for a step that differs much from the spec or must not move when the spec
    is edited (a figure frozen as published). It does not follow edits of the spec; its refs resolve inside its own
    tree.
  A writer converts one form to the other at any time (`diff` and `applyPatch`). Recommendation: a patch over the spec
  by default.
- **Keyframes and deltas** (the author, 2026-11-18, after Braid's rewind and video codecs): a full state and a patch
  over the spec are keyframes, reached in one step; a patch `from` another step is a delta. A step is at most
  `story.maxChain` deltas (8) from a keyframe, so reaching any step costs at most that many applies (walk back to the
  keyframe, apply forward). The storing helper is the encoder: `story.store: smaller` picks per step, within the bound,
  a delta from the previous step, a patch over the spec when the chain reaches the bound, or a full state when that is
  smaller; a step the writer marked frozen stays a full state, one marked to follow the spec stays a patch over it.
- **States in memory are already shared**: revisions share every unchanged subtree, so many states cost what differs.
  The stored form matters for the saved file; sampled animation frames are never stored.
- **A step's state** is `stateOf(spec, step)`: the step's full state, or the spec with the step's patch applied. It is
  computed on demand and cached by the spec's revision and the step. A patch applies over the spec itself by default,
  so any step is reached directly and editing one step never changes another. A step may name another step as its
  start (`from: /snapshots/<id>`, either form), so a writer can store a chain of small diffs; a cycle or a missing start
  is an issue on the step.
- **A step's view** is part of its state: what each section shows (a focus address, or a range per section) and the
  view shown. The viewport is session state elsewhere ([changes.md](changes.md)); a story names it on purpose, as data
  of the step, so a state is a whole picture.
- **Refs in a step stored as a patch** are typed by kind and key and resolved over the step's dry run
  ([patches.md](patches.md), A stored patch). They follow renames and moves like any other ref. A full state's refs are
  its own tree's, checked like a spec's.
- **A patch step that no longer applies** (an op whose target is gone, a sibling renamed, a move now refused) is
  reported on the step, and its state is the ops that still apply, flagged. An edit of the spec is never refused
  because of what it does to a step's dry run; a step's ref left dangling is refused as for any ref (`allowDangling`).
- **Entering a step is navigation**: no history step, no change to the spec. Edits while a patch step is shown patch
  the spec (the base), and the step is applied again over it; while a full-state step is shown, edits are refused with
  a hint to edit the step (the author, 2026-11-18). Editing a step's own content is its own mode and command ("edit this step"), each change one
  history step.

### A transition is planned from the diff of two states

- `plan(from: State, to: State, options?) → TransitionPlan`: plain data, a list of **channels**, each one property that
  changes, its kind, its start and end values and its time window within the transition. An application may inspect a
  plan, retime it, drop or add a channel.
- **Channel kinds and their default rules** (each a config item, [Configuration](#configuration)):

  | Kind | What changes | Default |
  |---|---|---|
  | `view` | a section's range, the view shown | an interpolated path: pan if the ranges overlap, else out, across and in; a change of view crossfades |
  | `number` | a numeric param (a height, a domain, an opacity, a lane gap) | tweened |
  | `color` | a color or a theme's resolved colors | blended in a perceptual space |
  | `enter` / `leave` | a track, layer or section appears or goes | fades and grows (or shrinks) its height; the layout moves the rest |
  | `order` | children reordered | each moves to its new place |
  | `cut` | anything else (a layer type, a data source, a text) | switched at the channel's time (default the midpoint), named in the plan |

- **Data does not animate**: data arriving during a transition draws in place. A channel whose values need data not
  yet held waits on it, and says so in the frame.
- **Reduced motion** makes every channel a cut at the end.

### Sampling is pure; a frame is shown, never committed

- `sample(plan, t) → Frame` for `t` in [0, 1] (eased per channel): deterministic, allocation-free once warm, the same
  function for preview, scrubbing and export.
- **A frame is shown as a transient layer** over the shown state, in the slot a preview uses: shown =
  frame(step(applyPatch(base, preview))). It makes no patch, no revision, no history step; the spec's hash never moves
  during an animation. `sample(plan, 1)` shows the `to` state exactly; leaving the animation leaves the step's state.

### The application owns the clock: the animation port

- **`AnimationPort`** (`seqx-viewer/ports`), its two faces as the history port's:
  - *the application drives Seq\**: `show(frame)` or `progress(plan, t)` when it decides; `clear()` returns to the
    state. Seq\* reports per drawn frame (`drawn(frame)`), and a plan's natural duration (from the steps' timing), so
    the application keeps several views and a structure viewer in step on one clock.
  - *Seq\* asks the application*: Seq\*'s own motions (focus and reveal jumps, a snapshot's view, a zoom step, ux-foundations section 7) are
    plans handed to the port (`request(plan, { duration })`), so a host that animates a structure's camera runs both
    together, or cuts both.
- **Adapters** in each composition root, as for history: the default composition's runs `request`ed plans on
  `requestAnimationFrame`, so a plain viewer animates with no host code; the headless one runs a fake clock frame by
  frame (export); a Mol\* plugin's would use the plugin's animation loop.
- **Timing is data on the step**: the transition into it (`duration`, `easing`) and how long it lingers (`linger`),
  with defaults. It feeds `plan`; the application may override it.

### Export

- A story renders to frames headlessly through the same port: the headless adapter's clock steps `t`, the SVG or
  raster backend draws each frame. What is previewed is what is exported. The render CLI's video and frame-sequence
  outputs are this.

## Configuration

`story.transition.duration` (the base duration of ux-foundations section 7), `story.transition.easing` (ease in and out), `story.linger`,
`story.view.path` (`auto`: pan when the ranges overlap, else out-across-in), `story.color.space` (`oklab`),
`story.cut.at` (0.5), `story.stepStart` (`spec`; `from` per step),
`story.store` (`patch`; `full`, or `smaller` for the encoder), `story.maxChain` (8). Whether motion plays at all is
the viewer's one item, `view.motion` (`platform`: animate unless the platform asks for reduced motion; `cut`;
`animate`), read by jumps and stories alike (decided, 2026-11-19: one item, `story.reducedMotion` dropped).

## In the specs

| Decision | Where |
|---|---|
| The `snapshots` container; a step's id, label, content, view and timing | core › Story steps |
| A step is a patch or a full state, the writer's choice; converted by `diff` and `applyPatch` | core › Story steps |
| Keyframes and deltas; `story.maxChain`; the encoder (`story.store`) | core › The encoder |
| States in memory share subtrees; sampled frames never stored | core › Sampling |
| A step's state (`stateOf`), cached by revision and step; `from` chains, a cycle an issue | core › A step's state |
| A step's view is part of its state | core › Story steps |
| Refs in a patch step follow renames and moves | core › Its addresses are refs of the spec |
| A patch step that no longer applies: the ops that still do, flagged | core › A patch that no longer applies whole |
| Entering a step is navigation; edits patch the base; a full-state step refuses edits | viewer › Story steps |
| A mode and command for editing a step's own content | design-only: not built |
| `plan(from, to)`: plain data, channels and their kinds and rules | core › A plan |
| Data does not animate; a channel waiting on data says so | core › A plan |
| Reduced motion makes every channel a cut at the end (`view.motion`) | core › A plan |
| `sample(plan, t)`: pure, deterministic, allocation-free warm | core › Sampling |
| A frame is shown in the frame slot, never committed | viewer › The frame slot |
| The animation port: the application drives, or Seq\* asks | viewer › Animation |
| Adapters: the default on the page's frames, a headless fake clock | viewer › Adapters |
| Timing is data on the step | core › Story steps |
| Playing a step in the viewer (`playStep`) | viewer › Stories |
| Export: a story rendered to frames headlessly, the render CLI's video | design-only: not built (the headless clock is there; no CLI output) |
| The configuration items | viewer › Configuration |

## Tensions

- **Full states vs small storage**: a patch over the spec repeats what an earlier step changed, and a full state
  repeats the spec. Accepted; the writer picks the form per step, and `from` chains when wanted, at the cost that a
  broken step breaks those after it.
- **A full-state step does not follow the spec**: renaming a track leaves a full state's copy as it was. That is its
  point (a frozen figure); a writer who wants a step to follow stores it as a patch.
- **A step's view is data, the viewport is session**: a story pins views on purpose; editing outside a step never
  writes a view into the spec.
- **Interpolation needs both values**: a channel between a value and a missing one is an `enter`/`leave`, or a cut;
  never invented in between.

## Consequences

- `seqx-core`: `stateOf`, `plan`, `sample` and the channel kinds, DOM-free, tested in Node with property tests
  (`sample(plan, 0)` is `from`, `sample(plan, 1)` is `to`, monotone channels, determinism) and allocation checks.
- `seqx-viewer`: the frame slot over the shown state; the `AnimationPort` with its default and headless adapters; ux-foundations section 7's
  motions through it.
- Stories in Storybook: the p53 story played by the default adapter, scrubbed by a slider (the application's clock),
  and the same plan driving a stand-in structure camera, with headless twins on a fake clock.
- Built after lane N's editing (the base the frame slot sits over); lane K3 rebuilds stored-patch typing by kind and
  key in the second slice. The library and the port are placed at the next gate.
