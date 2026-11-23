/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The checks of a map's and a numbering's methods (_spec/README.md, Language, Maps and Numberings): which params a
// method takes, and the runs of `segments` and `structure` maps and of `segments` numberings.

import { CategoryColumn, isColumn, type Node, type Value } from '@molstar/core/repr/tree';
import type { Code } from './check.js';

type Say = (code: Code, message: string, ...tail: (string | number)[]) => void;

const p = (n: Node, key: string): Value | undefined => n.params?.[key];

/**
 * A column of strings as the normalized spec holds it: a category column, or an array where it has nulls (the one
 * reader of `segments`' and `status`' word columns, the checks' and resolve's).
 */
export function strings(v: Value | undefined): { get(i: number): string | null } | undefined {
  if (v instanceof CategoryColumn) return v;
  if (!Array.isArray(v)) return undefined;
  return { get: (i: number) => (typeof v[i] === 'string' ? (v[i] as string) : null) };
}

/** The params each map method takes besides `from`, `to` and `priority`. */
const MapParams: Readonly<Record<string, readonly string[]>> = {
  offset: ['by'],
  codon: ['exons', 'strand', 'table'],
  segments: ['segments', 'source'],
  structure: ['segments', 'status'],
};
/** The params each method needs. */
const MapNeeds: Readonly<Record<string, readonly string[]>> = {
  offset: ['by'],
  codon: ['exons', 'strand'],
  structure: ['segments'],
};
const MethodParams = new Set(Object.values(MapParams).flat());

/** Whether a map's params fit its method: what it needs given, nothing another method takes. */
export function mapParamsFit(n: Node): boolean {
  const method = String(p(n, 'method'));
  const takes = MapParams[method] ?? [];
  for (const k of MapNeeds[method] ?? []) if (p(n, k) === undefined) return false;
  for (const k of MethodParams) if (!takes.includes(k) && p(n, k) !== undefined) return false;
  // a segments map's runs inline or a pairs source's, one of the two
  if (method === 'segments' && (p(n, 'segments') === undefined) === (p(n, 'source') === undefined)) return false;
  // a segment's provenance is a structure map's
  const segments = p(n, 'segments') as Record<string, Value> | undefined;
  return method === 'structure' || segments?.provenance === undefined;
}

/**
 * A `segments` or `structure` map's runs: each at least one position, none overlapping another on `from`, a strand
 * `+` or `-`; a structure map's none overlapping on `to` either (one residue of the structure is one of the sequence;
 * a segments map may be many-to-one, a fusion onto its part, EC-20: _spec/README.md "Resolve", Maps); its status rows
 * each `[start, end)` and `conflict` or `unobserved`.
 */
export function checkMapSegments(n: Node, say: Say): void {
  const t = p(n, 'segments') as Record<string, Value> | undefined;
  if (t && p(n, 'method') === 'structure' && isColumn(t.to) && isColumn(t.length)) {
    const to = t.to as ArrayLike<number>;
    const length = t.length as ArrayLike<number>;
    const strand = strings(t.strand);
    // each run's target `[lo, hi)`: `to` is the image of its first position, downward on `-`
    const lo = (i: number) => (strand?.get(i) === '-' ? to[i] - length[i] + 1 : to[i]);
    const order = Array.from({ length: Math.min(to.length, length.length) }, (_, i) => i)
      .filter((i) => length[i] >= 1)
      .sort((a, b) => lo(a) - lo(b));
    for (let j = 1; j < order.length; j++)
      if (lo(order[j - 1]) + length[order[j - 1]] > lo(order[j])) {
        say('seq.lang.map.segments', 'Two runs of a structure map overlap on `to`', 'segments', order[j]);
        return;
      }
  }
  if (t && isColumn(t.from) && isColumn(t.length)) {
    const from = t.from as ArrayLike<number>;
    const length = t.length as ArrayLike<number>;
    const strand = strings(t.strand);
    const order = Array.from({ length: Math.min(from.length, length.length) }, (_, i) => i).sort(
      (a, b) => from[a] - from[b],
    );
    for (let j = 0; j < order.length; j++) {
      const i = order[j];
      const s = strand?.get(i);
      const bad = length[i] < 1 || (s != null && s !== '+' && s !== '-');
      const overlaps = j > 0 && from[order[j - 1]] + length[order[j - 1]] > from[i];
      if (bad || overlaps) {
        say('seq.lang.map.segments', 'A run holds no position, overlaps another or has no strand', 'segments', i);
        break;
      }
    }
  }
  const st = p(n, 'status') as Record<string, Value> | undefined;
  if (st && isColumn(st.start) && isColumn(st.end)) {
    const start = st.start as ArrayLike<number>;
    const end = st.end as ArrayLike<number>;
    const word = strings(st.status);
    for (let i = 0; i < Math.min(start.length, end.length); i++) {
      const w = word?.get(i);
      if (end[i] < start[i] || (w !== 'conflict' && w !== 'unobserved')) {
        say('seq.lang.map.segments', 'A status row ends before it starts, or says no status', 'status', i);
        break;
      }
    }
  }
}

/** The params each numbering method takes. */
const NumberingParams: Readonly<Record<string, readonly string[]>> = {
  labels: ['labels', 'source', 'field'],
  offset: ['by'],
  segments: ['segments'],
};
const NumberingMethodParams = new Set(Object.values(NumberingParams).flat());

/**
 * Whether a numbering's params fit its method: `labels` one of `labels` and `source` (`field` with a source only);
 * `offset` its `by`; `segments` its runs; nothing another method takes.
 */
export function numberingParamsFit(n: Node): boolean {
  const method = String(p(n, 'method') ?? 'labels');
  const takes = NumberingParams[method] ?? [];
  for (const k of NumberingMethodParams) if (!takes.includes(k) && p(n, k) !== undefined) return false;
  if (method === 'offset') return p(n, 'by') !== undefined;
  if (method === 'segments') return p(n, 'segments') !== undefined;
  const labels = p(n, 'labels');
  const source = p(n, 'source');
  return (labels === undefined) !== (source === undefined) && (p(n, 'field') === undefined || source !== undefined);
}

/** A `segments` numbering's runs: each `[start, end)`, apart from the others; an insertion code on one position. */
export function checkNumberingSegments(n: Node, say: Say): void {
  const t = p(n, 'segments') as Record<string, Value> | undefined;
  if (!t || !isColumn(t.start) || !isColumn(t.end)) return;
  const start = t.start as ArrayLike<number>;
  const end = t.end as ArrayLike<number>;
  const insertion = strings(t.insertion);
  const order = Array.from({ length: Math.min(start.length, end.length) }, (_, i) => i).sort(
    (a, b) => start[a] - start[b],
  );
  for (let j = 0; j < order.length; j++) {
    const i = order[j];
    const coded = insertion?.get(i) != null && end[i] - start[i] !== 1;
    if (end[i] < start[i] || coded || (j > 0 && end[order[j - 1]] > start[i])) {
      say('seq.lang.numbering.segments', 'A run ends before it starts, overlaps another, or codes many', 'segments', i);
      break;
    }
  }
}
