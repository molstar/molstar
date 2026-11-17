/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Paths between systems (_spec/README.md "Resolve", Paths): maps composed without expanding, at most one codon step,
// the shortest unambiguous one, and what a path answers: the image of ranges, of a boundary, an element's status.
// A codon step is a scale: a residue `k` is the triplet `3k, 3k + 1, 3k + 2` of fine places, so every map on the
// protein's side is scaled by 3 and the whole path is one linear map of fine places.

import {
  composeLinear,
  ImplicitLinear,
  isEmptyLinear,
  type Linear,
  linearDir,
  linearsEqual,
  scaleLinear,
} from './linear.js';

/** `1`: a codon step from residues to nucleotides; `-1`: back; `0`: none (positions to positions). */
export type Step = -1 | 0 | 1;

/** A directed step between two systems: a map read forward or backward, or a unit's place on its molecule. */
export interface Edge {
  /** The map's address (a unit's own address for its place on its molecule). */
  readonly map: string;
  readonly from: string;
  readonly to: string;
  /** Whether the map is read backward. */
  readonly backward: boolean;
  readonly step: Step;
  /**
   * The step's linear part: positions to positions, or for a codon step fine places (`3k + j`) of the residue side to
   * the nucleotides (forward), and nucleotides to fine places (backward).
   */
  readonly linear: Linear;
  /** A unit's step to its molecule (`up`) or from it (`down`): free, never counted, under the copy rule. */
  readonly unit?: 'up' | 'down';
}

/** A composed path: one linear map between the systems' positions, or fine places on the side of a codon step. */
export interface Path {
  readonly from: string;
  readonly to: string;
  readonly step: Step;
  readonly linear: Linear;
  /** The steps it composes, in order (empty from a system to itself). */
  readonly edges: readonly Edge[];
}

/** The path from a system to itself. */
export function identityPath(system: string): Path {
  return { from: system, to: system, step: 0, linear: new ImplicitLinear(0, 1, -Infinity, Infinity), edges: [] };
}

/** The path read the other way. */
export function inversePath(p: Path): Path {
  const edges = p.edges
    .slice()
    .reverse()
    .map((e) => ({
      ...e,
      from: e.to,
      to: e.from,
      backward: !e.backward,
      step: -e.step as Step,
      linear: e.linear.inverse(),
    }));
  return { from: p.to, to: p.from, step: -p.step as Step, linear: p.linear.inverse(), edges };
}

/** Composes edges into one path; undefined when it takes more than one codon step (EC-19). */
export function composeEdges(from: string, edges: readonly Edge[]): Path | undefined {
  let step: Step = 0;
  let linear: Linear = new ImplicitLinear(0, 1, -Infinity, Infinity);
  for (const e of edges) {
    if (e.step !== 0) {
      if (step !== 0) return undefined;
      step = e.step;
      // a step up from residues: what came before is on the residue side, scaled to fine places
      if (e.step === 1) linear = scaleLinear(linear, 3);
      linear = composeLinear(linear, e.linear);
    } else linear = composeLinear(linear, step === -1 ? scaleLinear(e.linear, 3) : e.linear);
  }
  return { from, to: edges.length ? edges[edges.length - 1].to : from, step, linear, edges };
}

// ---- what a path answers ----

/**
 * The image of `[start, end, …]` pairs as `[start, end, …]` pairs on the target (unsorted, possibly overlapping). With
 * `out`, it is cleared and filled, so a warm buffer is reused with no allocation (the hover loop's); without, a new
 * array.
 */
export function projectRanges(p: Path, ranges: ArrayLike<number>, out: number[] = []): number[] {
  out.length = 0;
  for (let i = 0; i + 1 < ranges.length; i += 2) {
    const s = ranges[i];
    const e = ranges[i + 1];
    const from = out.length;
    if (p.step === 1) p.linear.targets(3 * s, 3 * e, out);
    else p.linear.targets(s, e, out);
    // nucleotides lift to the residues whose triplets they are in
    if (p.step === -1)
      for (let j = from; j < out.length; j += 2) {
        out[j] = Math.floor(out[j] / 3);
        out[j + 1] = Math.ceil(out[j + 1] / 3);
      }
  }
  return out;
}

/**
 * The continuous image of a position `x` (a viewport's fractional edge): the edge where element `floor(x)` starts on
 * the target, plus the fraction times the scale in the run's direction; undefined where the element maps nowhere. A
 * follower moves smoothly with its leader within a run.
 */
export function projectPoint(p: Path, x: number): number | undefined {
  const k = Math.floor(x);
  const frac = x - k;
  const r = p.linear.forward(p.step === 1 ? 3 * k : k, (p.step === 1 ? 3 * k : k) + 1)[0];
  if (!r) return undefined;
  const edge = r.dir === 1 ? r.to : r.to + 1;
  if (p.step === 1) return edge + frac * 3 * r.dir;
  if (p.step === -1) return (edge + frac * r.dir) / 3;
  return edge + frac * r.dir;
}

/** The direction of a path: −1 where it reverses (a minus-strand gene), else 1 (its first run's for a path that turns). */
export function pathDir(p: Path): 1 | -1 {
  return linearDir(p.linear);
}

/** The fine image of one place, with the direction of its run, or undefined. */
function imageOf(l: Linear, k: number): { to: number; dir: 1 | -1 } | undefined {
  const r = l.forward(k, k + 1)[0];
  return r ? { to: r.to, dir: r.dir } : undefined;
}

/**
 * Where a boundary lands (EC-28, EC-29): the boundary before element `k` maps before its element's first image; if
 * that is unmapped, after element `k − 1`'s image; else nowhere (undefined). From nucleotides to residues, a boundary
 * inside a codon lands after its residue.
 */
export function projectBoundary(p: Path, k: number): number | undefined {
  const fine = p.step === 1 ? 3 * k : k;
  let b: number | undefined;
  const at = imageOf(p.linear, fine);
  if (at) b = at.dir === 1 ? at.to : at.to + 1;
  else {
    const before = imageOf(p.linear, fine - 1);
    if (before) b = before.dir === 1 ? before.to + 1 : before.to;
  }
  if (b === undefined || p.step !== -1) return b;
  return b % 3 === 0 ? b / 3 : Math.floor(b / 3) + 1;
}

export type ElementStatus = 'paired' | 'conflict' | 'unpaired';

/** An element's status along a path: paired when all of it has an image (a residue: its whole triplet). */
export function pathStatus(p: Path, k: number): ElementStatus {
  if (p.step === 1) {
    let n = 0;
    for (const r of p.linear.forward(3 * k, 3 * k + 3)) n += r.length;
    return n === 3 ? 'paired' : 'unpaired';
  }
  return p.linear.forward(k, k + 1).length ? 'paired' : 'unpaired';
}

// ---- search ----

export interface PathGraph {
  /** The steps out of a system: the systems whose steps the search asks are what its answer depends on. */
  edgesFrom(system: string): readonly Edge[];
}

export interface PathBounds {
  /** Counted steps at most (unit steps are free): 4. */
  readonly maxSteps?: number;
  /** Candidates of the shortest length at most, before ambiguity is judged: 16. */
  readonly maxCandidates?: number;
}

export const DefaultPathBounds = { maxSteps: 4, maxCandidates: 16 } as const;

export type PathSearch =
  | { readonly ok: true; readonly path: Path }
  | { readonly ok: false; readonly why: 'none' | 'ambiguous' | 'limit'; readonly candidates: number };

/**
 * The shortest unambiguous path (EC-19 to EC-25): candidates by counted length, unit steps free; at most one codon
 * step; a candidate that pairs nothing is dropped first; equal compositions count once; two that differ are
 * ambiguous. The copy rule is part of the search: no path climbs from a unit to its molecule and down to another
 * unit of it, so a longer legal path beats a shorter illegal one.
 */
export function findPath(graph: PathGraph, from: string, to: string, bounds: PathBounds = {}): PathSearch {
  if (from === to) return { ok: true, path: identityPath(from) };
  const maxSteps = bounds.maxSteps ?? DefaultPathBounds.maxSteps;
  const maxCandidates = bounds.maxCandidates ?? DefaultPathBounds.maxCandidates;
  for (let length = 0; length <= maxSteps; length++) {
    const found: Path[] = [];
    let over = false;
    const edges: Edge[] = [];
    const visited = new Set<string>([from]);
    const climbed = new Set<string>();
    const go = (at: string, counted: number, codon: boolean) => {
      if (over) return;
      if (at === to) {
        if (counted !== length) return;
        const p = composeEdges(from, edges.slice());
        if (p && !isEmptyLinear(p.linear)) {
          if (found.length >= maxCandidates) over = true;
          else found.push(p);
        }
        return;
      }
      for (const e of graph.edgesFrom(at)) {
        if (visited.has(e.to)) continue;
        const cost = e.unit ? 0 : 1;
        if (counted + cost > length) continue;
        if (e.step !== 0 && codon) continue;
        // the copy rule: down from a molecule only if the path never climbed out of it
        if (e.unit === 'down' && climbed.has(e.from)) continue;
        visited.add(e.to);
        edges.push(e);
        if (e.unit === 'up') climbed.add(e.to);
        go(e.to, counted + cost, codon || e.step !== 0);
        if (e.unit === 'up') climbed.delete(e.to);
        edges.pop();
        visited.delete(e.to);
      }
    };
    go(from, 0, false);
    if (over) return { ok: false, why: 'limit', candidates: found.length };
    if (!found.length) continue;
    const distinct: Path[] = [];
    for (const p of found)
      if (!distinct.some((q) => q.step === p.step && linearsEqual(q.linear, p.linear))) distinct.push(p);
    if (distinct.length === 1) return { ok: true, path: distinct[0] };
    return { ok: false, why: 'ambiguous', candidates: distinct.length };
  }
  return { ok: false, why: 'none', candidates: 0 };
}
