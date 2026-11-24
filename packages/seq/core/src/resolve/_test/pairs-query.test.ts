/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Pairs answer `items` (_spec/README.md "Queries", Pairs answer `items`): a pairs annotation is read by the query forms
// as its items, a pair one item with a row per end on its own system, so `(items :overlapping <location>)` gives the
// pairs with an end in the location, on either system; `extent` and a result's location cover both ends.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { evaluateQuery, type Location, locationOfResult } from '@molstar/seqx-core/query';
import { describe, expect, it } from 'vitest';
import { resolveSpec } from '../model.js';

/** Two chains and contacts between them (a on `x`, b on `y`), and pairs within `x` (a hairpin's). */
function spec() {
  const b = createSeqSpec();
  const x = b.data.sequence('x', { alphabet: 'protein', letters: 'MKTAYIAKQRQISFVKSHFSRQ' });
  const y = b.data.sequence('y', { alphabet: 'protein', letters: 'GSHMLEVLFQGPLGSMENFQKV' });
  b.data.pairs('contacts', {
    a: x,
    b: y,
    rows: {
      a_start: [1, 5, 10, 18],
      a_end: [2, 6, 11, 19],
      b_start: [3, 15, 7, 20],
      b_end: [4, 16, 8, 21],
      id: ['c1', 'c2', 'c3', 'c4'],
      type: ['contact', 'contact', 'hbond', 'contact'],
    },
  });
  b.data.ranges('yf', { target: y, rows: { start: [2], end: [5], id: ['site'] } });
  b.data.pairs('none', { a: x, b: y, rows: { a_start: [], a_end: [], b_start: [], b_end: [] } });
  b.data.pairs('stem', {
    a: x,
    rows: { a_start: [0, 1, 2], a_end: [1, 2, 3], b_start: [20, 19, 18], b_end: [21, 20, 19] },
  });
  return b.build();
}

const Spec = spec();
const model = resolveSpec(Spec);
const q = (a: string) => evaluateQuery(Spec, model, a, { limit: 100 });
const keys = (a: string) => {
  const r = q(a);
  expect(r.issues.filter((i) => i.severity === 'error')).toEqual([]);
  return r.result?.type === 'keys' ? r.result.keys.keys.map((k) => k.key) : undefined;
};
const ranges = (l: Location) => l.parts.map((p) => [p.system, [...p.ranges]]);

describe('pairs answer items', () => {
  it('a pair is one item, a row per end on its own system', () => {
    const a = model.annotation('/data/contacts')!;
    expect(a.target).toBe('/data/x');
    const rows = a.items().items;
    expect(rows.map((r) => [r.key, r.system, r.start, r.end])).toEqual([
      ['c1', '/data/x', 1, 2],
      ['c1', '/data/y', 3, 4],
      ['c2', '/data/x', 5, 6],
      ['c2', '/data/y', 15, 16],
      ['c3', '/data/x', 10, 11],
      ['c3', '/data/y', 7, 8],
      ['c4', '/data/x', 18, 19],
      ['c4', '/data/y', 20, 21],
    ]);
    expect(rows[0].fields?.type).toBe('contact');
  });

  it('items :overlapping gives the pairs with an end in the location, on either system', () => {
    // on a's system
    expect(keys('/data/contacts/(items :overlapping /data/x/(range :from 1 :to 7))')).toEqual(['c1', 'c2']);
    // on b's: the pairs whose other end is there (not those at the same numbers on a)
    expect(keys('/data/contacts/(items :overlapping /data/y/(range :from 1 :to 8))')).toEqual(['c1', 'c3']);
    // both at once, each pair once
    expect(
      keys(
        '/data/contacts/(items :overlapping (merge :of [/data/x/(range :from 2 :to 2) /data/y/(range :from 4 :to 4)]))',
      ),
    ).toEqual(['c1']);
    // a system the pairs are on neither of is unreachable, as for ranges
    const other = q('/data/stem/(items :overlapping /data/y/(all))');
    expect(other.issues.map((i) => i.code)).toEqual(['seq.query.location.unreachable']);
    // pairs on one system: either end
    expect(keys('/data/stem/(items :overlapping /data/x/(range :from 20 :to 21))')).toEqual([0, 1]);
  });

  it('reach is each annotation’s own: a pair’s far end is read beside another annotation on that system', () => {
    // the prefix names the pairs (on x and y) and a ranges annotation on y: the pairs are read on y too
    const both = '/data/(merge :of [/data/contacts /data/yf])';
    expect(keys(`${both}/(items :overlapping /data/y/(range :from 1 :to 8))`)).toEqual(['c1', 'c3', 'site']);
    // a location on both systems: the pairs read on both
    expect(
      keys(
        `${both}/(items :overlapping (merge :of [/data/x/(range :from 19 :to 19) /data/y/(range :from 16 :to 16)]))`,
      ),
    ).toEqual(['c2', 'c4']);
  });

  it('an empty pairs annotation on two systems reaches both: no unreachable issue', () => {
    const r = q('/data/none/(items :overlapping /data/y/(all))');
    expect(r.issues).toEqual([]);
    expect(r.result?.type === 'keys' && r.result.keys.keys).toEqual([]);
  });

  it('item reads a pair’s fields; extent and a result’s location cover both ends', () => {
    expect(keys('/data/contacts/(item :type hbond)')).toEqual(['c3']);
    const e = q('/data/contacts/(extent :of /data/contacts/(item :id c2))');
    expect(e.result?.type === 'location' && ranges(e.result.location)).toEqual([
      ['/data/x', [5, 6]],
      ['/data/y', [15, 16]],
    ]);
    const r = q('/data/contacts/(item :id c3)');
    expect(ranges(locationOfResult(Spec, model, r.result!).location)).toEqual([
      ['/data/x', [10, 11]],
      ['/data/y', [7, 8]],
    ]);
  });
});
