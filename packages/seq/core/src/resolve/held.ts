/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Annotations held whole in the static model (_spec/README.md "Resolve", Pairs and grids): a `pairs` node's ends and a
// `grid` node's values, read once from the spec into typed arrays (typed arrays by reference) and kept with the entry,
// never through the store: they are inline, and read by region only with the Genomes milestone.

import { DefaultMaxBlockValues } from '@molstar/seqx-core/data';
import { type Issue, issueOf, isColumn, type Node, type Value } from '@molstar/core/repr/tree';
import { ResolveIssues } from './issues.js';

/**
 * A `pairs` node's rows as arrays: each pair's two ends, `[aStart, aEnd)` on `a`'s system and `[bStart, bEnd)` on
 * `b`'s; its kind (the `type` column) as an index into `kinds` (−1 none), its key (its `id`, else its row), its label, each end's strand.
 */
export interface PairsData {
  readonly n: number;
  readonly aStart: Int32Array;
  readonly aEnd: Int32Array;
  readonly bStart: Int32Array;
  readonly bEnd: Int32Array;
  readonly kind: Int16Array;
  /** The kinds in order of first appearance. */
  readonly kinds: readonly string[];
  readonly keys: readonly (string | number)[];
  readonly labels?: readonly (string | null)[];
  /** −1 where an end is on its `-` strand, else 1 (EC-206: a reversed side is valid). */
  readonly aStrand: Int8Array;
  readonly bStrand: Int8Array;
  /** Any other category column of the rows, by name, a value a kept pair (a flag a filter reads: `tertiary`). */
  readonly extra?: Readonly<Record<string, readonly (string | null)[]>>;
}

/** A `grid` node's values: `nx × ny`, x-major (`values[i × ny + j]` at x i, y j), NaN missing. */
export interface GridData {
  readonly nx: number;
  readonly ny: number;
  readonly values: Float32Array;
}

const p = (n: Node, key: string): Value | undefined => n.params?.[key];

type Category = { get(i: number): string | null };

/** A category column as normalized, or its values as a list (a tree read without normalizing). */
function category(v: Value | undefined): Category | undefined {
  if (v === undefined || v === null) return undefined;
  if (Array.isArray(v)) return { get: (i) => (typeof v[i] === 'string' ? (v[i] as string) : null) };
  return typeof (v as unknown as Category).get === 'function' ? (v as unknown as Category) : undefined;
}

const PairColumns = new Set(['a_start', 'a_end', 'b_start', 'b_end', 'id', 'type', 'label', 'a_strand', 'b_strand']);

/** A pairs node's rows read into arrays; rows whose ends are not intervals are left out (the checks said so). */
export function readPairs(node: Node): PairsData {
  const rows = (p(node, 'rows') ?? {}) as Record<string, Value>;
  const col = (k: string): ArrayLike<number> => (isColumn(rows[k]) ? (rows[k] as unknown as ArrayLike<number>) : []);
  const as = col('a_start');
  const ae = col('a_end');
  const bs = col('b_start');
  const be = col('b_end');
  const n0 = Math.min(as.length, ae.length, bs.length, be.length);
  const kindCol = category(rows.type);
  const idCol = category(rows.id);
  const aStr = category(rows.a_strand);
  const bStr = category(rows.b_strand);
  const label = Array.isArray(rows.label) ? (rows.label as (string | null)[]) : undefined;
  const aStart = new Int32Array(n0);
  const aEnd = new Int32Array(n0);
  const bStart = new Int32Array(n0);
  const bEnd = new Int32Array(n0);
  const kind = new Int16Array(n0);
  const aStrand = new Int8Array(n0);
  const bStrand = new Int8Array(n0);
  const kinds: string[] = [];
  const kindIndex = new Map<string, number>();
  const keys: (string | number)[] = [];
  const labels: (string | null)[] = [];
  const extraCols = Object.entries(rows)
    .filter(([k]) => !PairColumns.has(k))
    .flatMap(([k, v]) => {
      const c = category(v);
      return c ? [[k, c, [] as (string | null)[]] as const] : [];
    });
  let n = 0;
  for (let i = 0; i < n0; i++) {
    if (!(as[i] >= 0 && ae[i] > as[i] && bs[i] >= 0 && be[i] > bs[i])) continue;
    aStart[n] = as[i];
    aEnd[n] = ae[i];
    bStart[n] = bs[i];
    bEnd[n] = be[i];
    const k = kindCol?.get(i);
    if (k === null || k === undefined) kind[n] = -1;
    else {
      let ki = kindIndex.get(k);
      if (ki === undefined) {
        ki = kinds.length;
        kinds.push(k);
        kindIndex.set(k, ki);
      }
      kind[n] = ki;
    }
    aStrand[n] = aStr?.get(i) === '-' ? -1 : 1;
    bStrand[n] = bStr?.get(i) === '-' ? -1 : 1;
    keys.push(idCol?.get(i) ?? i);
    if (label) labels.push(label[i] ?? null);
    for (const [, c, out] of extraCols) out.push(c.get(i));
    n++;
  }
  return {
    n,
    aStart: aStart.subarray(0, n),
    aEnd: aEnd.subarray(0, n),
    bStart: bStart.subarray(0, n),
    bEnd: bEnd.subarray(0, n),
    kind: kind.subarray(0, n),
    kinds,
    keys,
    ...(label ? { labels } : {}),
    aStrand: aStrand.subarray(0, n),
    bStrand: bStrand.subarray(0, n),
    ...(extraCols.length ? { extra: Object.fromEntries(extraCols.map(([k, , out]) => [k, out])) } : {}),
  };
}

/**
 * A grid node's values on its two systems' lengths: a wrong length or more than `DefaultMaxBlockValues` values (the
 * bound the store checks a block) is dropped with an issue saying the count it expected.
 */
export function readGrid(
  address: string,
  node: Node,
  nx: number | undefined,
  ny: number | undefined,
): { grid?: GridData; issues: Issue[] } {
  const values = p(node, 'values');
  if (!isColumn(values) || nx === undefined || ny === undefined) return { issues: [] };
  const v = values as unknown as Float32Array;
  if (nx * ny > DefaultMaxBlockValues)
    return {
      issues: [
        issueOf(
          ResolveIssues,
          'seq.resolve.grid.size',
          `A grid of ${nx * ny} values is past ${DefaultMaxBlockValues}`,
          {
            address,
            hint: ResolveIssues.codes['seq.resolve.grid.size'].hint,
          },
        ),
      ],
    };
  if (v.length !== nx * ny)
    return {
      issues: [
        issueOf(ResolveIssues, 'seq.resolve.grid.size', `The grid has ${v.length} values, not ${nx * ny}`, {
          address,
          hint: ResolveIssues.codes['seq.resolve.grid.size'].hint,
        }),
      ],
    };
  return { grid: { nx, ny, values: v instanceof Float32Array ? v : Float32Array.from(v) }, issues: [] };
}
