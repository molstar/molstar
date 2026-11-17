/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Numberings and labels at the edge (_spec/README.md "Queries", Numberings and labels): the built-ins `index`
// (0-based) and `ordinal` (1-based, the default), the lookup of a declared numbering (a unit's own before its
// molecule's), and labels resolved to positions: a single label exactly, a range's missing end snapped inward. Labels
// need not rise with the index (IMGT): the forms snap first, then apply one rule in index order.

import { compareLabels, parseLabel } from './labels.js';
import type { Numbering, QuerySystem } from './model.js';

/** The numbering forms use without `:numbering`: 1-based. */
export const DefaultNumbering = 'ordinal';
export const BuiltInNumberings: readonly string[] = ['index', 'ordinal'];

const IntegerLabel = /^(0|-?[1-9]\d*)$/;

/**
 * Positions numbered from `first` (`index`: 0, `ordinal`: 1). Its labels rise with the index, so `snap` and `nearest`
 * answer in O(1): a label's integer places it (a suffix puts it just after its integer, a name after every number).
 */
export function sequentialNumbering(id: string, length: number, first: number): Numbering {
  const labelAt = (i: number) => (Number.isInteger(i) && i >= 0 && i < length ? `${i + first}` : undefined);
  const snap = (label: string, side: 'start' | 'end') => {
    const p = parseLabel(label);
    if (p.int === undefined) return side === 'start' ? -1 : length - 1;
    if (side === 'start') {
      const i = Math.max(0, p.int - first + (p.kind === 0 ? 0 : 1));
      return i < length ? i : -1;
    }
    const i = Math.min(length - 1, p.int - first);
    return i >= 0 ? i : -1;
  };
  return {
    id,
    length,
    labelAt,
    indexOf: (label) => {
      if (!IntegerLabel.test(label)) return -1;
      const i = Number(label) - first;
      return i >= 0 && i < length ? i : -1;
    },
    snap,
    nearest: (label, count) => {
      const s = snap(label, 'start');
      const at = s < 0 ? length : s;
      const from = Math.max(0, Math.min(at - (count >> 1), length - count));
      const out: string[] = [];
      for (let i = from; i < Math.min(length, from + count); i++) out.push(`${i + first}`);
      return out;
    },
  };
}

/** A numbering and how its indices shift onto the system's (a molecule's numbering of a unit). */
export interface NumberingOf {
  readonly numbering: Numbering;
  /** The system's index of the numbering's index 0: a unit's position i is its molecule's `start + i`. */
  readonly shift: number;
}

/**
 * The numbering of a system by id: a built-in, the system's own, else its molecule's, shifted to the unit; undefined
 * when none has the id. A unit's own numbering beats its molecule's of the same id.
 */
export function findNumbering(system: QuerySystem, id: string): NumberingOf | undefined {
  if (id === 'index') return { numbering: sequentialNumbering('index', system.length, 0), shift: 0 };
  if (id === 'ordinal') return { numbering: sequentialNumbering('ordinal', system.length, 1), shift: 0 };
  let shift = 0;
  for (let s: QuerySystem | undefined = system; s; ) {
    const n = s.numbering(id);
    if (n) return { numbering: n, shift };
    if (!s.parent) break;
    shift -= s.parent.start;
    s = s.parent.system;
  }
  return undefined;
}

/** The numbering ids a system can use, for an issue's hint. */
export function numberingIds(system: QuerySystem): string[] {
  const ids = new Set<string>(BuiltInNumberings);
  for (let s: QuerySystem | undefined = system; s; s = s.parent?.system) for (const id of s.numberingIds()) ids.add(id);
  return [...ids];
}

/** The system's position of a label, or -1 (absent, or on the molecule outside the unit). */
export function labelIndex(of: NumberingOf, system: QuerySystem, label: string): number {
  const i = of.numbering.indexOf(label);
  if (i < 0) return -1;
  const j = i + of.shift;
  return j >= 0 && j < system.length ? j : -1;
}

/** Counts a visited position; the forms pass one that reads their signal and deadline. */
export type Tick = () => void;

const noTick: Tick = () => {};

/** Whether the numbering's own `snap` and `nearest` answer for the system: its indices are the system's. */
const direct = (of: NumberingOf, system: QuerySystem) => of.shift === 0 && of.numbering.length === system.length;

/**
 * Snaps a label the numbering lacks inward, by the order of labels: a range's start to the first position whose label
 * is at or after it, its end to the last whose label is at or before it; -1 when none is. The numbering's own `snap`
 * answers when it has one, else the unit's positions are scanned, each a tick.
 */
export function snapLabel(
  of: NumberingOf,
  system: QuerySystem,
  label: string,
  side: 'start' | 'end',
  tick: Tick = noTick,
): number {
  const { numbering } = of;
  if (numbering.snap && direct(of, system)) return numbering.snap(label, side);
  if (side === 'start') {
    for (let j = 0; j < system.length; j++) {
      tick();
      const l = numbering.labelAt(j - of.shift);
      if (l !== undefined && compareLabels(l, label) >= 0) return j;
    }
    return -1;
  }
  for (let j = system.length - 1; j >= 0; j--) {
    tick();
    const l = numbering.labelAt(j - of.shift);
    if (l !== undefined && compareLabels(l, label) <= 0) return j;
  }
  return -1;
}

/** Whether any position in `[from, to)` of the system has a label (each position a tick). */
export function anyLabelIn(of: NumberingOf, from: number, to: number, tick: Tick = noTick): boolean {
  for (let j = from; j < to; j++) {
    tick();
    if (of.numbering.labelAt(j - of.shift) !== undefined) return true;
  }
  return false;
}

/** How many positions a hint scans at most: past it, a hint names no labels (it is help, not an answer). */
export const HintScanLimit = 1 << 16;

/**
 * Up to `count` labels of the system nearest a label in the order of labels, for an issue's hint: the numbering's own
 * `nearest`, else a scan keeping only the `count` nearest on each side (no list of every label); none past
 * HintScanLimit positions.
 */
export function nearestLabels(
  of: NumberingOf,
  system: QuerySystem,
  label: string,
  count = 6,
  tick: Tick = noTick,
): string[] {
  const { numbering } = of;
  if (numbering.nearest && direct(of, system)) return numbering.nearest(label, count);
  if (system.length > HintScanLimit) return [];
  // the greatest labels before it (ascending) and the least at or after it (ascending), each at most `count`
  const below: string[] = [];
  const above: string[] = [];
  for (let j = 0; j < system.length; j++) {
    tick();
    const l = numbering.labelAt(j - of.shift);
    if (l === undefined) continue;
    if (compareLabels(l, label) < 0) {
      if (below.length === count && compareLabels(l, below[0]) <= 0) continue;
      insertSorted(below, l);
      if (below.length > count) below.shift();
    } else {
      if (above.length === count && compareLabels(l, above[count - 1]) >= 0) continue;
      insertSorted(above, l);
      if (above.length > count) above.pop();
    }
  }
  const take = Math.min(below.length, Math.max(count >> 1, count - above.length));
  return [...below.slice(below.length - take), ...above].slice(0, count);
}

function insertSorted(list: string[], l: string) {
  let i = list.length;
  while (i > 0 && compareLabels(list[i - 1], l) > 0) i--;
  list.splice(i, 0, l);
}
