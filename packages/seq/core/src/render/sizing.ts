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

// Scratch of the content sizing: the content-sized tracks' indices and their lanes' px, kept between settles (grown, never
// shrunk): a sizing allocates nothing once warm beyond its answer.
const order: number[] = [];
const laneOf: number[] = [];
const EmptyTracks: readonly SizedTrack[] = [];
let sizedTracks: readonly SizedTrack[] = EmptyTracks;

function byKindThenLanes(a: number, b: number): number {
  const ka = sizedTracks[a].kind;
  const kb = sizedTracks[b].kind;
  return ka === kb ? laneOf[a] - laneOf[b] : ka < kb ? -1 : 1;
}

/**
 * Whether lanes `px` tall (the lanes' px, no padding) are within 2 lanes or 2x of `base`, compared in px (a lane
 * `laneHeight` px, not a float lane count: a fractional `laneHeight` must not flip a boundary), with a tolerance of
 * rounding.
 */
function close(base: number, px: number, laneHeight: number): boolean {
  const eps = 1e-9 * Math.max(1, px);
  return px - base <= 2 * laneHeight + eps || px <= 2 * base + eps;
}

/**
 * Tracks fit their content; tracks of one kind whose lanes (their content less the label lane) are within 2 lanes or 2x
 * of the cluster's smallest share the tallest of them (neighbours do not vary lane by lane: two CDR tracks are equal), a
 * track far from the rest fits its own (one 40-lane track does not draw thirty one-feature tracks at the cap). The
 * clusters are by sorted lane height, anchored on the smallest: a track joins when it is within 2 lanes or 2x of the
 * cluster's smallest, so the result does not depend on the tracks' order (1, 3 and 5 lanes: {1, 3} and {5}). Each track
 * keeps its own label lane (a track whose labels all fit inside reserves none).
 */
export const ContentSizing: TrackSizing = {
  name: 'content',
  heights(tracks, config) {
    const n = tracks.length;
    const out = new Array<number>(n);
    let m = 0;
    for (let i = 0; i < n; i++) {
      const t = tracks[i];
      if (t.content === undefined) {
        out[i] = t.fixed;
        continue;
      }
      order[m++] = i;
      laneOf[i] = t.content - (t.label ?? 0) - 2 * config.lanePadding;
    }
    order.length = m;
    sizedTracks = tracks;
    try {
      order.sort(byKindThenLanes);
      for (let a = 0; a < m; ) {
        const kind = tracks[order[a]].kind;
        const base = laneOf[order[a]];
        let b = a;
        let tallest = 0;
        while (
          b < m &&
          tracks[order[b]].kind === kind &&
          (b === a || close(base, laneOf[order[b]], config.laneHeight))
        ) {
          const t = tracks[order[b]];
          tallest = Math.max(tallest, (t.content as number) - (t.label ?? 0));
          b++;
        }
        for (let k = a; k < b; k++) out[order[k]] = tallest + (tracks[order[k]].label ?? 0);
        a = b;
      }
    } finally {
      sizedTracks = EmptyTracks;
    }
    return out;
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

/** What offers labels to `TakenLabels`: told once they are settled, whether each is kept. */
export interface TakenOwner {
  placed(id: number, kept: boolean): void;
}

/**
 * The labels a section's layers offer at points (a ring's stub labels), settled together once the section's layers are
 * laid out (_spec/README.md "Render", Section layouts): the strongest first (`weight`, ties by offer order), each kept
 * only where its box meets no label kept, no wedge (`block`: a triangle no label may cover, the layers' stubs) and its
 * tip is the clearance from every tip kept; the others dropped. Kept by the layout, its arrays grown: no allocation once warm.
 * @internal
 */
export class TakenLabels {
  n = 0;
  private b = new Float64Array(9 * 16);
  private owners: (TakenOwner | undefined)[] = [];
  private ids = new Int32Array(16);
  private nt = 0;
  private t = new Float64Array(6 * 16);
  private order: number[] = [];

  reset() {
    this.n = 0;
    this.nt = 0;
  }

  /** A wedge no label may cover (a triangle: three points). */
  block(ax: number, ay: number, bx: number, by: number, cx: number, cy: number) {
    const o = 6 * this.nt++;
    if (this.t.length < o + 6) {
      const g = new Float64Array(2 * this.t.length);
      g.set(this.t);
      this.t = g;
    }
    const t = this.t;
    t[o] = ax;
    t[o + 1] = ay;
    t[o + 2] = bx;
    t[o + 3] = by;
    t[o + 4] = cx;
    t[o + 5] = cy;
  }

  /** A label offered: its box, its tip and the clearance it keeps from the tips kept, its weight (the stronger is kept). */
  offer(
    owner: TakenOwner,
    id: number,
    x0: number,
    y0: number,
    x1: number,
    y1: number,
    x: number,
    y: number,
    weight: number,
    tip: number,
  ) {
    const q = this.n++;
    const o = 9 * q;
    if (this.b.length < o + 9) {
      const g = new Float64Array(2 * this.b.length);
      g.set(this.b);
      this.b = g;
      const ids = new Int32Array(2 * this.ids.length);
      ids.set(this.ids);
      this.ids = ids;
    }
    const b = this.b;
    b[o] = x0;
    b[o + 1] = y0;
    b[o + 2] = x1;
    b[o + 3] = y1;
    b[o + 4] = x;
    b[o + 5] = y;
    b[o + 6] = weight;
    b[o + 8] = tip;
    this.owners[q] = owner;
    this.ids[q] = id;
  }

  /** Whether a box meets triangle `k` (separating axes: the box's two and the triangle's three edge normals). */
  private hits(x0: number, y0: number, x1: number, y1: number, k: number): boolean {
    const t = this.t;
    const o = 6 * k;
    for (let e = 0; e < 5; e++) {
      const i = o + 2 * (e - 2);
      const j = e < 4 ? i + 2 : o;
      const nx = e < 2 ? 1 - e : t[j + 1] - t[i + 1];
      const ny = e < 2 ? e : t[i] - t[j];
      let lo = Number.POSITIVE_INFINITY;
      let hi = Number.NEGATIVE_INFINITY;
      let blo = lo;
      let bhi = hi;
      for (let q = 0; q < 4; q++) {
        const d = nx * (q & 1 ? x1 : x0) + ny * (q & 2 ? y1 : y0);
        blo = Math.min(blo, d);
        bhi = Math.max(bhi, d);
        if (q < 3) {
          const p = nx * t[o + 2 * q] + ny * t[o + 2 * q + 1];
          lo = Math.min(lo, p);
          hi = Math.max(hi, p);
        }
      }
      if (bhi <= lo || blo >= hi) return false;
    }
    return true;
  }

  /** Settles the offered labels (see the class), tells each owner, and clears the lists. */
  settle() {
    const n = this.n;
    if (n === 0) {
      this.nt = 0;
      return;
    }
    const b = this.b;
    const order = this.order;
    order.length = n;
    for (let i = 0; i < n; i++) order[i] = i;
    order.sort((p, q) => b[9 * q + 6] - b[9 * p + 6] || p - q);
    const kept: number[] = [];
    const out = new Uint8Array(n);
    for (const i of order) {
      const o = 9 * i;
      let by = -1;
      for (const j of kept) {
        const k = 9 * j;
        if (
          (b[o] < b[k + 2] && b[o + 2] > b[k] && b[o + 1] < b[k + 3] && b[o + 3] > b[k + 1]) ||
          Math.hypot(b[k + 4] - b[o + 4], b[k + 5] - b[o + 5]) < b[o + 8]
        ) {
          by = j;
          break;
        }
      }
      if (by >= 0) continue;
      // the strongest label is always kept (a hub's wedges cover every place by it)
      let covered = false;
      for (let w = 0; w < this.nt && !covered && kept.length > 0; w++)
        covered = this.hits(b[o], b[o + 1], b[o + 2], b[o + 3], w);
      if (covered) continue;
      kept.push(i);
      out[i] = 1;
    }
    for (let i = 0; i < n; i++) this.owners[i]?.placed(this.ids[i], out[i] === 1);
    this.n = 0;
    this.nt = 0;
  }
}

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
