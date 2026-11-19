/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Stories (_spec/README.md "Stories"): a step's state, a transition planned from the diff of two states, sampled at
// any time; the story's encoder. DOM-free: the application's clock decides what is shown and when.

export {
  DefaultStoryConfig,
  type Easing,
  Easings,
  easingOf,
  type StepTiming,
  type StoryConfig,
} from './config.js';
export { StoryIssues } from './issues.js';
export {
  type Channel,
  type ChannelKind,
  type ColorChannel,
  type CutChannel,
  colorFrom,
  type EnterChannel,
  type LeaveChannel,
  type NumberChannel,
  type OrderChannel,
  type PlanOptions,
  plan,
  type StoryState,
  stripStory,
  type TransitionPlan,
  type ViewChannel,
  type ViewState,
} from './plan.js';
export { ChannelPhase, createFrame, type Frame, frameDone, frameRoot, frameView, sample } from './sample.js';
export {
  applyLoose,
  createStoryCache,
  readStepView,
  type StateOfOptions,
  StepApplyOptions,
  type StepState,
  type StepView,
  type StepViewport,
  type StoryCache,
  stateOf,
  stepPatchOf,
  storySteps,
  timingOf,
} from './state.js';
export { type StepInput, type StoreOptions, storeSteps, withSteps } from './store.js';
