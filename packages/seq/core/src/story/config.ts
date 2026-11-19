/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The stories' configuration items (_spec/README.md "Stories", Configuration) and the easings a channel's time runs
// through. Plain defaults: none is inferred from what a step shows.

/** How a channel's time runs: `ease-in-out` (cubic, the default), `ease-in`, `ease-out`, `linear`. */
export type Easing = 'linear' | 'ease-in' | 'ease-out' | 'ease-in-out';

export const Easings: readonly Easing[] = ['linear', 'ease-in', 'ease-out', 'ease-in-out'];

/** A step's timing: the transition into it and how long it lingers, ms. */
export interface StepTiming {
  readonly duration: number;
  readonly easing: Easing;
  readonly linger: number;
}

export interface StoryConfig {
  readonly transition: {
    /** The transition into a step, ms (240: the base duration of a jump). */
    readonly duration: number;
    readonly easing: Easing;
  };
  /** How long a step is shown before the next, ms (2,000). */
  readonly linger: number;
  /** A view channel's path: `auto` pans when the windows overlap, else out, across and in; `pan`; `out-across-in`. */
  readonly view: { readonly path: 'auto' | 'pan' | 'out-across-in' };
  /** Where colors blend: `oklab` (perceptual, the default) or `srgb`. */
  readonly color: { readonly space: 'oklab' | 'srgb' };
  /** When a cut switches, a fraction of its channel's window (0.5). */
  readonly cut: { readonly at: number };
  /** What a patch step applies over when it names no start: the spec (`spec`). */
  readonly stepStart: 'spec';
  /** How the encoder stores a step: a patch over the spec (`patch`), a full state (`full`), or the `smaller`. */
  readonly store: 'patch' | 'full' | 'smaller';
  /** The most deltas between a step and its keyframe (8). */
  readonly maxChain: number;
}

export const DefaultStoryConfig: StoryConfig = {
  transition: { duration: 240, easing: 'ease-in-out' },
  linger: 2000,
  view: { path: 'auto' },
  color: { space: 'oklab' },
  cut: { at: 0.5 },
  stepStart: 'spec',
  store: 'patch',
  maxChain: 8,
};

/** An easing by name; the default for anything else. */
export function easingOf(x: unknown, fallback: Easing = DefaultStoryConfig.transition.easing): Easing {
  return typeof x === 'string' && (Easings as readonly string[]).includes(x) ? (x as Easing) : fallback;
}
