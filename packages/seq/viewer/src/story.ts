/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createStoryCache, frameRoot, frameView, plan, sample, stateOf, stepPatchOf } from '@molstar/seqx-core/story';

// The viewer's story part (_spec/README.md "Composition", Stories and Motion): the core's story steps and transition
// plans, listed by the spec (`story`), so a viewer that lists none bundles none of them: its jumps are cuts and it
// enters no story step (`seq.viewer.step.none`).

/** A story step's state over a root (the spec model's): the core's step cache, state and stored patch. */
export interface SeqStorySteps {
  readonly createStoryCache: typeof createStoryCache;
  readonly stateOf: typeof stateOf;
  readonly stepPatchOf: typeof stepPatchOf;
}

/** How the viewer plans and samples a transition (a jump's, a story step's), and reads story steps. */
export interface SeqStoryPart extends SeqStorySteps {
  readonly plan: typeof plan;
  readonly sample: typeof sample;
  readonly frameRoot: typeof frameRoot;
  readonly frameView: typeof frameView;
}

/** The core's stories and transitions: what the default spec lists. */
export const StoryPart: SeqStoryPart = { plan, sample, frameRoot, frameView, createStoryCache, stateOf, stepPatchOf };
