# Plan: predictive streaming (wanting what a view is about to see)

## Status

- **Current**: proposed, nothing built. Decisions open (see Open).
- Wants that run ahead of the camera, and release that follows what a view can reach, for every kind that streams
  by view: the residency's pages first (their wants the GPU cull's), the tiles and SOG after through the shared
  traversal.
- Builds on: [particle-residency.md](particle-residency.md) (pages by the cull's wants, the release after 10 s, the
  split of a copy lacking pages: step 2d), [streaming.md](streaming.md) (the shared traversal's "prefetch along the
  camera's motion" and hysteresis, the scheduler's priority, budget and cooldown: Design 3 to 5),
  [gl/mesoscale-lod.md](gl/mesoscale-lod.md) (the metric: error in pixels; time fades for a page becoming resident and
  settling: Design 1 and 4).

## Question

The wants are of the view a cull sees. What a view is about to see is guessed only by a page's padding (half its
radius, `PagedLeafPagePadding`); what it can see again by nothing (a page let go 10 s after its last want). Two cases
in the nested yeast 2×2, the first shown by the residency's split fix:

| Case                                     | What happens now                                                              |
| ---------------------------------------- | ----------------------------------------------------------------------------- |
| **A first focus** (a click on a molecule) | the transition flies from the whole view into pages never written; each frame's cull wants only what it reaches, the leaves near arriving behind the camera, near holes (2d's split) or the sample meanwhile |
| **Turning about a focused molecule**     | its slab's other pages let go after 10 s (about 1,000), back in view as it turns (28): read and written again, near holes for the frames they take |
| A wheel zoom, a dolly                    | as a first focus, frame by frame; the end not known                           |
| A pan                                    | the padding alone ahead                                                       |

## What exists

- **The residency's wants** (`paged-cull-wants`, `readWants`): a leaf page of a copy at level 0 in the frustum wanted
  (its copy lacking it: split near and far), within half its radius of it wanted alone; read back at most every
  250 ms where none are pending; written within 32 MB a cull; let go 10 s after the last read that wanted it, or for
  room (least recently wanted first).
- **The camera's intents**: a transition (`Camera.setState` with a duration: focus, reset, a snapshot, an animation)
  holds its end snapshot from its start; the trackball's controls their target, its rotation about it.
- **Planned elsewhere, not built**: the shared traversal's prefetch along the camera's motion and its hysteresis
  ([streaming.md](streaming.md) 3), the scheduler's priority and cooldown (4), the fade on arrival and the settled state
  ([gl/mesoscale-lod.md](gl/mesoscale-lod.md) 4).

## Goals and how they are measured

- **A first focus whole when it lands**: the frames from a transition's end till its view is drawn whole (no near
  hole, no sample where leaves draw), and the near holes' frames during it.
- **Turning about a target costs no reads**: the pages read and written again a turn about a focused molecule; the
  near holes' frames.
- **Within the budget, not at the current view's cost**: the bytes read ahead and never drawn (overfetch), the value
  blocks held; a predicted want never lets go of a page the current view wants.
- Measured on the yeast nested 2×2 (a first focus, a turn) and the isg nested 16³ (a dolly), in Chrome and headless.

## Design

1. **Predicted views, wanted alone.** A predictor gives each frame a few views besides the camera's: frustums and
   levels (`LodView`) a camera will see. The residency's wants kernel runs over them too (the same pass, a view a
   dispatch or a view index a thread): their pages wanted alone (no copy lacking, no split: nothing drawn changes), at
   a priority below the current view's. Why the GPU's wants, not a CPU guess: the residency's lesson (2c against 2d:
   the GPU knows what it draws). Kinds without a GPU want (the tiles, SOG) take the same views through the shared
   traversal ([streaming.md](streaming.md) 3).
2. **A transition's end view** (the first focus, a reset, a snapshot): known at its start (its end snapshot), wanted
   from its first frame; its path's middle as one or two views where it is long (a zoom from the whole view passes
   the copies' coarse levels: wanted at the levels each view draws, by the metric, not leaves everywhere).
3. **The orbit's reach** (turning about a target): the pages a turn about the camera's target can bring into view,
   kept wanted while the target and distance hold: a sphere about the target of `min(radius, distance)` (a focus's
   radius: its slab; a wheel zoom keeps the scene's radius, the distance bounds it), its pages at the levels a view
   at that distance draws. Wanted alone, never let go while in it; not written ahead where never in view (only kept).
4. **Motion ahead** (a pan, a dolly, a turn without a target): the camera's velocity (from its last frames)
   extrapolated by the time a want takes to be drawn (the read back, the read, the write: measured in step 1), one
   view ahead. Hysteresis as the shared traversal's.
5. **Release by reachability**: a page let go where no current or predicted view wanted it for the cooldown, not by
   age alone; for room the least likely first (not in the orbit's reach, then least recently wanted).
6. **Arrival**: what a prediction brings arrives before its view; what arrives late fades in
   ([gl/mesoscale-lod.md](gl/mesoscale-lod.md) 4), the split's near hole filled by it.

## Steps

1. **Baseline**: the cases measured as they are (the goals' metrics), headless with a scripted focus and turn, and in
   Chrome (the browser harness: a click, the transition, a minute's turns).
2. **Predicted views in the wants kernel**: views a dispatch, wanted alone, their priority below the current view's
   (written after it within the upload budget); tests: a page only a predicted view reaches wanted and written, its
   copy not lacking, nothing drawn changed.
3. **A transition's end view**: the camera's transition end handed to the renderer's wants (Canvas3D: the transition's
   end snapshot to a `LodView` and frustum); measured: a first focus's frames till whole.
4. **The orbit's reach**; measured: a turn's pages written again and near holes.
5. **Release by reachability**; measured: the value blocks held, a turn's reads.
6. **Motion ahead**; measured: a pan's and a dolly's near holes against the overfetch.
7. **The shared traversal's predicted views** for the tiles and SOG ([streaming.md](streaming.md) 5).

## Open

- **How far ahead**: the time a want takes to be drawn varies (a remote pack's reads take seconds): a fixed lead or
  one measured a source.
- **The prediction's share of the budget**: a fraction of the value blocks and of the upload a cull, or whatever the
  current view leaves.
- **A transition's path**: its end alone, or samples along it (a long zoom crosses several levels).
- **The particle budget** (`maxParticles`): its coarse bias applied to the predicted views too, or their levels as
  unbudgeted.
- **Where the predictor lives**: the canvas (it knows the camera and its transitions) handing views to the renderer,
  or the shared traversal's (streaming), the residency taking them from there.
