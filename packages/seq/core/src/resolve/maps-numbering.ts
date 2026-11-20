/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A numbering's `offset` and `segments` methods (_spec/README.md "Resolve", Numberings): a position plus `by` (an
// author numbering that starts at 94), or runs of consecutive numbers with an insertion code on a run of one (a
// structure's author numbering with gaps and `52A`). Both are arithmetic: nothing per position is held.

import { isColumn, type Value } from '@molstar/core/repr/tree';
import { type Numbering, sequentialNumbering } from '@molstar/seqx-core/query';
import { strings } from './maps-segments.js';

/** A numbering of `length` positions, position `i` labelled `i + by`. */
export function offsetNumbering(id: string, length: number, by: number): Numbering {
  return sequentialNumbering(id, length, by);
}

/**
 * A numbering of runs: positions `[start, end)` labelled `first`, `first + 1`, …, each with its run's insertion code
 * (`52A`); a position in no run has no label. Labels need no order, so it has no `snap`: the forms scan.
 */
export function segmentsNumbering(id: string, length: number, table: Record<string, Value> | undefined): Numbering {
  const runs: { start: number; end: number; first: number; code: string }[] = [];
  if (table && isColumn(table.start) && isColumn(table.end) && isColumn(table.first)) {
    const start = table.start as ArrayLike<number>;
    const end = table.end as ArrayLike<number>;
    const first = table.first as ArrayLike<number>;
    const codes = strings(table.insertion);
    for (let i = 0; i < Math.min(start.length, end.length, first.length); i++)
      runs.push({ start: start[i], end: end[i], first: first[i], code: codes?.get(i) ?? '' });
  }
  runs.sort((a, b) => a.start - b.start);
  const at = (i: number) => {
    // the last run starting at or before `i`
    let lo = 0;
    let hi = runs.length;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (runs[mid].start <= i) lo = mid + 1;
      else hi = mid;
    }
    const r = runs[lo - 1];
    return r && i < r.end ? r : undefined;
  };
  const labelAt = (i: number) => {
    if (!Number.isInteger(i) || i < 0 || i >= length) return undefined;
    const r = at(i);
    return r ? `${r.first + i - r.start}${r.code}` : undefined;
  };
  return {
    id,
    length,
    labelAt,
    indexOf: (label) => {
      const m = /^(-?\d+)([A-Za-z]*)$/.exec(label);
      if (!m) return -1;
      const n = Number(m[1]);
      const code = m[2];
      // a repeated label names its first position; an insertion code is found ignoring case, as labels are
      let found = -1;
      for (const r of runs) {
        const i = r.start + n - r.first;
        if (i < r.start || i >= r.end || i >= length || r.code.toLowerCase() !== code.toLowerCase()) continue;
        if (found < 0 || i < found) found = i;
      }
      return found;
    },
  };
}
