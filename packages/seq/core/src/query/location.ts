/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Location (_spec/README.md "Queries", Location): positions on systems. Per system one part, its ranges sorted,
// merged, half-open index pairs in one Int32Array (`[start0, end0, start1, end1, …]`); rects pair two parts (two
// axes). A Location is a value: built normalized by `location`, never mutated after, compared by `locationsEqual`,
// hashed by `hashLocation`. A cyclic system's wrapped range is two ranges, as any two ranges are. Equality is
// structural: rects are kept as given (sorted, each once), never merged into fewer rects covering the same cells.

import { HalfOpenRanges } from '@molstar/core/data/int/half-open-ranges';
import { hash2, hashFnv32a, hashString } from '@molstar/core/data/util/hash-functions';
import { compareStrings } from './compare.js';

/** Whether a number is a position a Location holds: an integer within int32 (an Int32Array would wrap others). */
export function isPosition(n: number): boolean {
  return HalfOpenRanges.isInt32(n);
}

/** The positions of one system: sorted, merged, half-open `[start, end)` pairs, none empty. */
export interface LocationPart {
  /** The system's address (`/data/p53`), printed canonically. */
  readonly system: string;
  readonly ranges: Int32Array;
}

/** A rectangle over two axes: positions × positions, rows × columns, states × positions. */
export interface LocationRect {
  readonly a: LocationPart;
  readonly b: LocationPart;
}

export interface Location {
  /** One part per system, by system address; none empty. */
  readonly parts: readonly LocationPart[];
  /** Sorted (by `compareRects`) and each once; neither side empty. */
  readonly rects: readonly LocationRect[];
}

export const EmptyLocation: Location = Object.freeze({ parts: Object.freeze([]), rects: Object.freeze([]) });

/** Ranges of one system as flat `[start, end, …]` pairs, in any order, overlapping or empty. */
export interface PartInput {
  readonly system: string;
  readonly ranges: ArrayLike<number>;
}

export interface RectInput {
  readonly a: PartInput;
  readonly b: PartInput;
}

/**
 * Sorts and merges `[start, end, …]` pairs (adjacent ones too), dropping empty ones; a new array. An end that is not a
 * position (`isPosition`) is a programmer error and throws, never wraps: callers given data check it and report.
 */
export function normalizeRanges(ranges: ArrayLike<number>): Int32Array {
  return HalfOpenRanges.normalize(ranges);
}

function part(system: string, ranges: ArrayLike<number>): LocationPart {
  return Object.freeze({ system, ranges: normalizeRanges(ranges) });
}

/** Parts merged by system, sorted by system, empty ones dropped. */
function normalizeParts(parts: readonly PartInput[]) {
  const by = new Map<string, number[]>();
  for (const p of parts) {
    let acc = by.get(p.system);
    if (!acc) by.set(p.system, (acc = []));
    for (let i = 0; i < p.ranges.length; i++) acc.push(p.ranges[i]);
  }
  const out: LocationPart[] = [];
  for (const system of [...by.keys()].sort(compareStrings)) {
    const p = part(system, by.get(system)!);
    if (p.ranges.length) out.push(p);
  }
  return out;
}

/** A normalized, frozen Location. */
export function location(input: {
  readonly parts?: readonly PartInput[];
  readonly rects?: readonly RectInput[];
}): Location {
  const parts = normalizeParts(input.parts ?? []);
  const rects: LocationRect[] = [];
  for (const r of input.rects ?? []) {
    const a = part(r.a.system, r.a.ranges);
    const b = part(r.b.system, r.b.ranges);
    if (a.ranges.length && b.ranges.length) rects.push(Object.freeze({ a, b }));
  }
  return freeze(parts, sortedUnique(rects));
}

/** One system's ranges. */
export function rangesLocation(system: string, ranges: ArrayLike<number>): Location {
  return location({ parts: [{ system, ranges }] });
}

function freeze(parts: readonly LocationPart[], rects: readonly LocationRect[]): Location {
  if (!parts.length && !rects.length) return EmptyLocation;
  return Object.freeze({ parts: Object.freeze(parts.slice()), rects: Object.freeze(rects.slice()) });
}

function sortedUnique(rects: LocationRect[]): LocationRect[] {
  rects.sort(compareRects);
  const out: LocationRect[] = [];
  for (const r of rects) if (!out.length || compareRects(out[out.length - 1], r) !== 0) out.push(r);
  return out;
}

export function isEmptyLocation(l: Location): boolean {
  return l.parts.length === 0 && l.rects.length === 0;
}

/** The part of a system, or undefined. */
export function partOf(l: Location, system: string): LocationPart | undefined {
  for (const p of l.parts) if (p.system === system) return p;
  return undefined;
}

/** How many positions the parts hold (rects not counted). */
export function positionCount(l: Location): number {
  let n = 0;
  for (const p of l.parts) n += HalfOpenRanges.size(p.ranges);
  return n;
}

// ---- identity ----

function compareParts(a: LocationPart, b: LocationPart): number {
  return compareStrings(a.system, b.system) || HalfOpenRanges.compare(a.ranges, b.ranges);
}

/** A total order of rects: by side a, then side b. */
export function compareRects(a: LocationRect, b: LocationRect): number {
  return compareParts(a.a, b.a) || compareParts(a.b, b.b);
}

export function locationsEqual(a: Location, b: Location): boolean {
  if (a === b) return true;
  if (a.parts.length !== b.parts.length || a.rects.length !== b.rects.length) return false;
  for (let i = 0; i < a.parts.length; i++) if (compareParts(a.parts[i], b.parts[i]) !== 0) return false;
  for (let i = 0; i < a.rects.length; i++) if (compareRects(a.rects[i], b.rects[i]) !== 0) return false;
  return true;
}

const hashes = new WeakMap<Location, number>();

function hashPart(h: number, p: LocationPart): number {
  return hash2(hash2(h, hashString(p.system)), hashFnv32a(p.ranges));
}

/** A 32-bit hash, equal for equal Locations; cached per Location (they are values). */
export function hashLocation(l: Location): number {
  let h = hashes.get(l);
  if (h !== undefined) return h;
  h = hash2(l.parts.length, l.rects.length);
  for (const p of l.parts) h = hashPart(h, p);
  for (const r of l.rects) h = hashPart(hashPart(h, r.a), r.b);
  hashes.set(l, h);
  return h;
}

// ---- set operations on parts (per system); rects by the rule each names ----

function combineParts(a: Location, b: Location, keep: (inA: boolean, inB: boolean) => boolean): LocationPart[] {
  const systems = new Set<string>();
  for (const p of a.parts) systems.add(p.system);
  for (const p of b.parts) systems.add(p.system);
  const out: LocationPart[] = [];
  for (const system of [...systems].sort(compareStrings)) {
    const ranges = HalfOpenRanges.combine(
      partOf(a, system)?.ranges ?? HalfOpenRanges.Empty,
      partOf(b, system)?.ranges ?? HalfOpenRanges.Empty,
      keep,
    );
    if (ranges.length) out.push(Object.freeze({ system, ranges }));
  }
  return out;
}

/** Every position and rect of either. */
export function unionLocations(a: Location, b: Location): Location {
  return freeze(
    combineParts(a, b, (x, y) => x || y),
    sortedUnique([...a.rects, ...b.rects]),
  );
}

/** Positions in both; rects as the pairwise intersections of rects over the same two systems. */
export function intersectLocations(a: Location, b: Location): Location {
  const rects: LocationRect[] = [];
  for (const r of a.rects)
    for (const q of b.rects) {
      if (r.a.system !== q.a.system || r.b.system !== q.b.system) continue;
      const ra = HalfOpenRanges.intersect(r.a.ranges, q.a.ranges);
      const rb = HalfOpenRanges.intersect(r.b.ranges, q.b.ranges);
      if (ra.length && rb.length)
        rects.push(
          Object.freeze({
            a: Object.freeze({ system: r.a.system, ranges: ra }),
            b: Object.freeze({ system: r.b.system, ranges: rb }),
          }),
        );
    }
  return freeze(
    combineParts(a, b, (x, y) => x && y),
    sortedUnique(rects),
  );
}

/** The positions of `a` not in `b`. Rects are not subtracted: a caller refuses them (a rect minus a rect is no rect). */
export function subtractLocations(a: Location, b: Location): Location {
  return freeze(
    combineParts(a, b, (x, y) => x && !y),
    a.rects,
  );
}

/** Whether a half-open item `[start, end)` overlaps the ranges (`HalfOpenRanges.overlaps`: boundaries, the system's end). */
export function overlapsRanges(ranges: Int32Array, start: number, end: number, length?: number): boolean {
  return HalfOpenRanges.overlaps(ranges, start, end, length);
}
