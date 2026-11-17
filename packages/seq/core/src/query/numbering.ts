/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Numberings and labels at the edge (_spec/README.md "Queries", Numberings and labels): the built-ins `index`
// (0-based) and `ordinal` (1-based, the default), the lookup of a declared numbering (a unit's own before its
// molecule's), and labels resolved to positions: a single label exactly, a range's missing end snapped inward, a
// backwards range refused unless the system is cyclic, where it wraps as two ranges short of a whole lap.

import { compareLabels } from './labels.js';
import type { Numbering, QuerySystem } from './model.js';

/** The numbering forms use without `:numbering`: 1-based. */
export const DefaultNumbering = 'ordinal';
export const BuiltInNumberings: readonly string[] = ['index', 'ordinal'];

const IntegerLabel = /^(0|-?[1-9]\d*)$/;

/** Positions numbered from `first` (`index`: 0, `ordinal`: 1). */
export function sequentialNumbering(id: string, length: number, first: number): Numbering {
  return {
    id,
    length,
    labelAt: (i) => (Number.isInteger(i) && i >= 0 && i < length ? `${i + first}` : undefined),
    indexOf: (label) => {
      if (!IntegerLabel.test(label)) return -1;
      const i = Number(label) - first;
      return i >= 0 && i < length ? i : -1;
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

/**
 * Snaps a label the numbering lacks inward, by the order of labels: a range's start to the first position whose label
 * is at or after it, its end to the last whose label is at or before it; -1 when none is. Only the unit's positions
 * are scanned.
 */
export function snapLabel(of: NumberingOf, system: QuerySystem, label: string, side: 'start' | 'end'): number {
  if (side === 'start') {
    for (let j = 0; j < system.length; j++) {
      const l = of.numbering.labelAt(j - of.shift);
      if (l !== undefined && compareLabels(l, label) >= 0) return j;
    }
    return -1;
  }
  for (let j = system.length - 1; j >= 0; j--) {
    const l = of.numbering.labelAt(j - of.shift);
    if (l !== undefined && compareLabels(l, label) <= 0) return j;
  }
  return -1;
}

/** Up to `count` labels of the system nearest a label in the order of labels, for an issue's hint. */
export function nearestLabels(of: NumberingOf, system: QuerySystem, label: string, count = 6): string[] {
  const all: string[] = [];
  for (let j = 0; j < system.length; j++) {
    const l = of.numbering.labelAt(j - of.shift);
    if (l !== undefined) all.push(l);
  }
  all.sort(compareLabels);
  let at = all.findIndex((l) => compareLabels(l, label) >= 0);
  if (at < 0) at = all.length;
  const from = Math.max(0, Math.min(at - (count >> 1), all.length - count));
  return all.slice(from, from + count);
}
