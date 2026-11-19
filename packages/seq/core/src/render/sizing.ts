/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Track heights (_spec/README.md "Render", Track heights): a track with no height is as tall as its content, in steps
// of a lane (`content`, the default), every content-sized track the section's tallest (`uniform`), or its layer type's
// fixed default (`fixed`); an extension registers more (`ext:` names, `track-sizing` registry entries). A track's own
// `height` overrides any; a host's `trackSizing` port, given the tracks and the named sizing's answer, has the last word.

import type { ResolvedSection, ResolvedTrack } from '@molstar/seqx-core/view';
import type { RenderConfig } from './config.js';

/** What a layer's items need over the whole data (never the view): their packed lanes, and a lane for labels beside. */
export interface LayerContent {
  readonly lanes: number;
  /** Some item's label does not fit inside its mark at the fit view: a label lane above the marks. */
  readonly labels: boolean;
}

/** A section's track as a sizing reads it. */
export interface SizedTrack {
  readonly track: ResolvedTrack;
  /**
   * Its kind: the types of its layers in order (`ranges`, `letters+ranges`): tracks of one kind share the tallest
   * content height under `content`.
   */
  readonly kind: string;
  /** Its content's height, px: its lanes and their padding, at least its other layers' fixed heights; undefined: none. */
  readonly content?: number;
  /** Its layer types' fixed height (the tallest of its layers' types'). */
  readonly fixed: number;
  /** Its own height, when the spec gives one (it overrides any sizing). */
  readonly height?: number;
}

/** A named sizing: a section's tracks' heights from what they hold (one a track, in order). */
export interface TrackSizing {
  readonly name: string;
  heights(tracks: readonly SizedTrack[], config: RenderConfig): readonly number[];
}

/**
 * The host's last word (the composition's `trackSizing` port): a section's tracks and the named sizing's answer (its
 * own heights over it already), returning the heights; an answer of another length or not a list of px is ignored.
 */
export type TrackSizingHook = (
  section: ResolvedSection,
  tracks: readonly SizedTrack[],
  named: readonly number[],
) => readonly number[];

/** The height of `lanes` lanes and their padding (_spec/README.md "Render", Track heights: 22, 36, 50, 64, 78 px). */
export function lanesHeight(lanes: number, config: RenderConfig): number {
  const n = Math.max(1, Math.min(config.maxLanes, lanes));
  return 2 * config.lanePadding + n * config.laneHeight;
}

/** Tracks fit their content; tracks of one kind share the tallest (neighbours do not vary lane by lane). */
export const ContentSizing: TrackSizing = {
  name: 'content',
  heights(tracks) {
    const tallest = new Map<string, number>();
    for (const t of tracks)
      if (t.content !== undefined) tallest.set(t.kind, Math.max(tallest.get(t.kind) ?? 0, t.content));
    return tracks.map((t) => (t.content === undefined ? t.fixed : (tallest.get(t.kind) as number)));
  },
};

/** Every content-sized track the tallest content height in the section. */
export const UniformSizing: TrackSizing = {
  name: 'uniform',
  heights(tracks) {
    let h = 0;
    for (const t of tracks) if (t.content !== undefined) h = Math.max(h, t.content);
    return tracks.map((t) => (t.content === undefined ? t.fixed : h));
  },
};

/** Each track its layer types' fixed default (the second slice's look). */
export const FixedSizing: TrackSizing = {
  name: 'fixed',
  heights: (tracks) => tracks.map((t) => t.fixed),
};

/** The built-in sizings by name: what a section's `heights` names without an extension. */
export const BuiltInSizings: Readonly<Record<string, TrackSizing>> = {
  content: ContentSizing,
  uniform: UniformSizing,
  fixed: FixedSizing,
};

/** Whether a sizing's or a port's answer is a height a track (a list of px of the tracks' length). */
export function heightsFit(h: readonly number[] | undefined, n: number): h is readonly number[] {
  if (!Array.isArray(h) || h.length !== n) return false;
  for (const x of h) if (typeof x !== 'number' || !Number.isFinite(x) || x < 0) return false;
  return true;
}
