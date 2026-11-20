/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `segments` and `structure` maps (_spec/README.md "Resolve", Maps; maps.md, Structure maps): runs of positions,
// `from + i ↦ to + dir·i`, in their cheapest encoding (one run is arithmetic); a structure map adds each run's
// provenance and the residues' status on `to` (conflict, unobserved; what no run holds is unpaired). Unobserved is a
// unit's own fact: where a structure map onto a unit says otherwise, it is reported once and the unit wins.

import { HalfOpenRanges } from '@molstar/core/data/int/half-open-ranges';
import { isColumn, type Value } from '@molstar/core/repr/tree';
import { composeLinear, ImplicitLinear, type Linear, type Run, RunsLinear, runTargets, wrapRun } from './linear.js';

/** A structure map's own facts beside its runs. */
export interface StructureMapData {
  /** Each run's provenance on `from` (`sifts`, `alignment`, `author`), in `from` order; null where none is said. */
  readonly provenance: readonly { readonly start: number; readonly end: number; readonly word: string | null }[];
  /** The residues of `to` it says are paired but differ (`[start, end, …]`, sorted, merged). */
  readonly conflict: HalfOpenRanges;
  /** The residues of `to` it says are unobserved (sorted, merged); a unit's `observed` wins. */
  readonly unobserved: HalfOpenRanges;
}

export type StructureStatus = 'paired' | 'conflict' | 'unpaired' | 'unobserved';

interface Ends {
  readonly length?: number;
  readonly cyclic: boolean;
}

type Strings = { get(i: number): string | null };
const strings = (v: Value | undefined): Strings | undefined =>
  v && typeof (v as Strings).get === 'function' ? (v as Strings) : undefined;

/** A `segments` table's runs, in the order written. */
export function segmentRuns(table: Record<string, Value> | undefined): Run[] {
  if (!table || !isColumn(table.from) || !isColumn(table.to) || !isColumn(table.length)) return [];
  const from = table.from as ArrayLike<number>;
  const to = table.to as ArrayLike<number>;
  const length = table.length as ArrayLike<number>;
  const strand = strings(table.strand);
  const out: Run[] = [];
  for (let i = 0; i < Math.min(from.length, to.length, length.length); i++) {
    if (length[i] < 1) continue;
    out.push({ from: from[i], to: to[i], length: length[i], dir: strand?.get(i) === '-' ? -1 : 1 });
  }
  return out;
}

/**
 * The runs as a map between two systems: what lies past either maps nowhere (EC-21); onto a cyclic target of known
 * length a run wraps, split at the origin (EC-11). One run within its target is arithmetic, else runs.
 */
export function segmentsLinear(runs: readonly Run[], from: Ends, to: Ends): Linear {
  const lf = from.length ?? Infinity;
  const clipped: Run[] = [];
  for (const r of runs) {
    // the source side first: positions past `from` are none of its
    const s = Math.max(r.from, 0);
    const e = Math.min(r.from + r.length, lf);
    if (e <= s) continue;
    const run: Run = { from: s, to: r.to + r.dir * (s - r.from), length: e - s, dir: r.dir };
    if (to.cyclic && to.length !== undefined) clipped.push(...wrapRun(run, to.length));
    else clipped.push(run);
  }
  if (clipped.length === 1 && !(to.cyclic && to.length !== undefined)) {
    const r = clipped[0];
    const [lo, hi] = runTargets(r);
    const tl = to.length ?? Infinity;
    if (lo >= 0 && hi <= tl) return new ImplicitLinear(r.to - r.dir * r.from, r.dir, r.from, r.from + r.length);
  }
  const linear = new RunsLinear(clipped);
  if (to.cyclic || to.length === undefined) return linear;
  return composeLinear(linear, new ImplicitLinear(0, 1, 0, to.length));
}

/** A structure map's provenance and status rows. */
export function structureData(
  segments: Record<string, Value> | undefined,
  status: Value | undefined,
): StructureMapData {
  const runs = segmentRuns(segments);
  const words = strings(segments?.provenance);
  const provenance = runs
    .map((r, i) => ({ start: r.from, end: r.from + r.length, word: words?.get(i) ?? null }))
    .sort((a, b) => a.start - b.start);
  const conflict: number[] = [];
  const unobserved: number[] = [];
  const t = status as Record<string, Value> | undefined;
  if (t && isColumn(t.start) && isColumn(t.end)) {
    const start = t.start as ArrayLike<number>;
    const end = t.end as ArrayLike<number>;
    const word = strings(t.status);
    for (let i = 0; i < Math.min(start.length, end.length); i++) {
      const w = word?.get(i);
      if (w === 'conflict') conflict.push(start[i], end[i]);
      else if (w === 'unobserved') unobserved.push(start[i], end[i]);
    }
  }
  return {
    provenance,
    conflict: HalfOpenRanges.normalize(conflict),
    unobserved: HalfOpenRanges.normalize(unobserved),
  };
}

/** The image of a map on its target, sorted and merged. */
function imageOf(linear: Linear): HalfOpenRanges {
  const out: number[] = [];
  for (const r of linear.runs() ?? []) out.push(...runTargets(r));
  return HalfOpenRanges.normalize(out);
}

/**
 * Whether a structure map onto a unit says other than the unit of what is unobserved, on the residues it pairs: the
 * unit's `observed` (absent: all observed) against the map's `unobserved`.
 */
export function observedDisagrees(
  data: StructureMapData,
  linear: Linear,
  unit: { readonly length?: number; readonly observed?: ArrayLike<number> },
): boolean {
  const image = imageOf(linear);
  const own = HalfOpenRanges.intersect(data.unobserved, image);
  if (unit.length === undefined) return false;
  const all = HalfOpenRanges.normalize([0, unit.length]);
  const seen = unit.observed ? HalfOpenRanges.normalize(unit.observed) : all;
  const theirs = HalfOpenRanges.intersect(HalfOpenRanges.subtract(all, seen), image);
  return HalfOpenRanges.compare(own, theirs) !== 0;
}

/**
 * A residue's status on a structure map's target (`k` a position of `to`): `unobserved` by the unit's `observed` when
 * given (the unit wins), else the map's; then `conflict`; then `paired` where a run holds it, else `unpaired`.
 */
export function structureStatus(
  map: { readonly linear: Linear; readonly structure?: StructureMapData },
  k: number,
  observed?: ArrayLike<number>,
): StructureStatus {
  const unobserved = observed
    ? !HalfOpenRanges.overlaps(HalfOpenRanges.normalize(observed), k, k + 1)
    : !!map.structure && HalfOpenRanges.overlaps(map.structure.unobserved, k, k + 1);
  if (unobserved) return 'unobserved';
  if (map.structure && HalfOpenRanges.overlaps(map.structure.conflict, k, k + 1)) return 'conflict';
  return map.linear.backward(k, k + 1).length ? 'paired' : 'unpaired';
}

/** A position of `from`'s provenance on a structure map (the run holding it), or undefined. */
export function structureProvenance(data: StructureMapData, k: number): string | null | undefined {
  for (const p of data.provenance) if (k >= p.start && k < p.end) return p.word;
  return undefined;
}
