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
  /** The lanes are those a source's data has packed so far, not its whole (they grow as tiles settle). */
  readonly partial?: boolean;
}

const Knowns = new WeakMap<object, (() => number)[]>();

/** The sum of the registered types' `known` counts: moves when a source's content grew (the types' list kept). */
export function knownVersion(types: Readonly<Record<string, { known?(): number }>> | undefined): number {
  if (!types) return 0;
  let fns = Knowns.get(types);
  if (!fns) {
    fns = [];
    for (const k in types) {
      const f = types[k].known;
      if (f) fns.push(f);
    }
    Knowns.set(types, fns);
  }
  let v = 0;
  for (let i = 0; i < fns.length; i++) v += fns[i]();
  return v;
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
  /** What its own label lane adds to `content`, px (0 or absent: none): tracks of one kind share lanes, not label lanes. */
  readonly label?: number;
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

/**
 * The height of `lanes` lanes and their padding (_spec/README.md "Render", Track heights: 22, 36, 50, 64, 78 px).
 * @internal
 */
export function lanesHeight(lanes: number, config: RenderConfig): number {
  const n = Math.max(1, lanes);
  return 2 * config.lanePadding + n * config.laneHeight;
}

/**
 * Tracks fit their content; tracks of one kind share the tallest lanes (neighbours do not vary lane by lane), each with
 * its own label lane (a track whose labels all fit inside reserves none).
 */
export const ContentSizing: TrackSizing = {
  name: 'content',
  heights(tracks) {
    const tallest = new Map<string, number>();
    for (const t of tracks)
      if (t.content !== undefined) tallest.set(t.kind, Math.max(tallest.get(t.kind) ?? 0, t.content - (t.label ?? 0)));
    return tracks.map((t) => (t.content === undefined ? t.fixed : (tallest.get(t.kind) as number) + (t.label ?? 0)));
  },
};

/**
 * Every content-sized track the tallest content height in the section.
 * @internal
 */
export const UniformSizing: TrackSizing = {
  name: 'uniform',
  heights(tracks) {
    let h = 0;
    for (const t of tracks) if (t.content !== undefined) h = Math.max(h, t.content);
    return tracks.map((t) => (t.content === undefined ? t.fixed : h));
  },
};

/**
 * Each track its layer types' fixed default (the second slice's look).
 * @internal
 */
export const FixedSizing: TrackSizing = {
  name: 'fixed',
  heights: (tracks) => tracks.map((t) => t.fixed),
};

/**
 * The built-in sizings by name: what a section's `heights` names without an extension.
 * @internal
 */
export const BuiltInSizings: Readonly<Record<string, TrackSizing>> = {
  content: ContentSizing,
  uniform: UniformSizing,
  fixed: FixedSizing,
};

/**
 * A track's label lane (_spec/README.md "Render", Labels inside if they fit, else beside): the spans its layers' labels
 * took, shared by the layers of one band in their order (the layer's order is the first drop priority), cleared before
 * the band is laid out. Kept by the layout, its arrays grown: no allocation once warm.
 * @internal
 */
export class LabelLane {
  n = 0;
  x0 = new Float64Array(16);
  x1 = new Float64Array(16);
  /** The marks the band's layers drew (`mk` of them: x0, x1, lane): what a leader up to the label lane must not cross. */
  mk = 0;
  mx0 = new Float64Array(16);
  mx1 = new Float64Array(16);
  mlane = new Int32Array(16);

  reset() {
    this.n = 0;
    this.mk = 0;
  }

  /** A mark of a layer's lane `lane`, spanning `x0` to `x1`. */
  addMark(x0: number, x1: number, lane: number) {
    const q = this.mk++;
    if (q >= this.mx0.length) {
      const m = 2 * this.mx0.length;
      const a = new Float64Array(m);
      const b = new Float64Array(m);
      const c = new Int32Array(m);
      a.set(this.mx0);
      b.set(this.mx1);
      c.set(this.mlane);
      this.mx0 = a;
      this.mx1 = b;
      this.mlane = c;
    }
    this.mx0[q] = x0;
    this.mx1[q] = x1;
    this.mlane[q] = lane;
  }

  /**
   * Whether a leader from `x` up from the top of lane `lane` and right along the label line to `to` crosses a label
   * placed or a mark in a lane above `lane`.
   */
  crosses(x: number, lane: number, to: number): boolean {
    for (let q = 0; q < this.n; q++) if (this.x0[q] < to && this.x1[q] > x) return true;
    for (let q = 0; q < this.mk; q++) if (this.mlane[q] < lane && this.mx0[q] <= x && x < this.mx1[q]) return true;
    return false;
  }

  /**
   * The first x in `x0` to `x1` with no mark in a lane above `lane` over it: where a leader from an item in `lane` can
   * rise to the label lane (its start, for the first lane); NaN when its whole extent is covered.
   */
  clearX(x0: number, x1: number, lane: number): number {
    let x = x0;
    for (let moved = true; moved && x < x1; ) {
      moved = false;
      for (let q = 0; q < this.mk; q++)
        if (this.mlane[q] < lane && this.mx0[q] <= x && x < this.mx1[q]) {
          x = this.mx1[q];
          moved = true;
        }
    }
    return x < x1 ? x : Number.NaN;
  }

  /** The first x from `x` rightward where a label `w` wide stands `gap` clear of every label placed. */
  free(x: number, w: number, gap: number): number {
    for (let moved = true; moved; ) {
      moved = false;
      for (let q = 0; q < this.n; q++)
        if (x < this.x1[q] + gap && x + w + gap > this.x0[q]) {
          x = this.x1[q] + gap;
          moved = true;
        }
    }
    return x;
  }

  add(x: number, w: number) {
    const q = this.n++;
    if (q >= this.x0.length) {
      const a = new Float64Array(2 * this.x0.length);
      const b = new Float64Array(2 * this.x0.length);
      a.set(this.x0);
      b.set(this.x1);
      this.x0 = a;
      this.x1 = b;
    }
    this.x0[q] = x;
    this.x1[q] = x + w;
  }
}

/** Whether a sizing's or a port's answer is a height a track (a list of px of the tracks' length). */
export function heightsFit(h: readonly number[] | undefined, n: number): h is readonly number[] {
  if (!Array.isArray(h) || h.length !== n) return false;
  for (const x of h) if (typeof x !== 'number' || !Number.isFinite(x) || x < 0) return false;
  return true;
}
