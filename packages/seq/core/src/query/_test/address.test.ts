/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { describe, expect, it } from 'vitest';
import { addressFor, locationOfResult, refResults, resolveAddressIn } from '../address.js';
import { keySet } from '../keys.js';
import { location, rangesLocation } from '../location.js';
import type { QueryResult } from '../result.js';
import { Model, Spec } from './fixture.js';

// Addressing over the hand-written model: what an address names once its refs are followed, what waits and what is
// missing, and the way back, from what is named to its shortest stable address, which reads back to the same thing.

/** The fixture's spec with more refs: a chain, a set, a ref on a dynamic node, a prefix ref, a gone one. */
const Root: Node = {
  ...Spec,
  children: Spec.children!.map((c) =>
    c.kind !== 'refs'
      ? c
      : {
          ...c,
          children: [
            ...c.children!,
            { kind: 'ref', id: 'core', params: { to: '/refs/dbd', label: 'core' } },
            { kind: 'ref', id: 'protein', params: { to: '/data/p53' } },
            { kind: 'ref', id: 'tad', params: { to: '/data/p53/domains/(item :id TAD)' } },
            { kind: 'ref', id: 'sites', params: { to: ['/refs/tad', '/data/p53/domains/(item :id DBD)'] } },
            { kind: 'ref', id: 'mixed', params: { to: ['/refs/tad', '/data/p53'] } },
            { kind: 'ref', id: 't3', params: { to: '/data/genes/(node :key T3)' } },
            { kind: 'ref', id: 'gone', params: { to: '/data/p53/gone' } },
            { kind: 'ref', id: 'loop', params: { to: '/refs/loop2' } },
            { kind: 'ref', id: 'loop2', params: { to: '/refs/loop' } },
          ],
        },
  ),
};

const resolve = (a: string, o = {}) => resolveAddressIn(Root, Model, a, o);
const loc = (r: ReturnType<typeof resolve>) => (r.result?.type === 'location' ? r.result.location : undefined);

describe('resolving an address', () => {
  it('names a location through a ref, a chain of refs and a ref as a prefix', () => {
    const want = rangesLocation('/data/p53', [93, 292]);
    expect(resolve('/refs/dbd')).toMatchObject({ status: 'resolved', via: ['/refs/dbd'], issues: [] });
    expect(loc(resolve('/refs/dbd'))).toEqual(want);
    expect(resolve('/refs/core').via).toEqual(['/refs/core', '/refs/dbd']);
    expect(loc(resolve('/refs/core'))).toEqual(want);
    expect(loc(resolve('/refs/protein/(range :from 94 :to 292)'))).toEqual(want);
  });

  it('names the same positions in two numberings', () => {
    const a = loc(resolve('/data/p53/(range :from 94 :to 292)'));
    const b = loc(resolve('/data/p53/(range :from 93 :numbering index :to 291)'));
    expect(a).toEqual(b);
  });

  it('follows refs nested in a form’s operands', () => {
    const r = resolve('/data/p53/(intersect :of [/refs/dbd (range :from 1 :to 100)])');
    expect(r.status).toBe('resolved');
    expect(loc(r)).toEqual(rangesLocation('/data/p53', [93, 100]));
    expect(r.via).toEqual(['/refs/dbd']);
    const set = resolve('/(extent :of /refs/sites)');
    expect(loc(set)).toEqual(rangesLocation('/data/p53', [0, 61, 93, 292]));
  });

  it('follows refs only where a form takes an operand: a quoted string elsewhere is a value', () => {
    const r = resolve('/data/p53/domains/(item :id "/refs/dbd")');
    expect(r.status).toBe('missing');
    expect(r.via).toEqual([]);
    expect(resolve('/data/p53/domains/(item :label "/refs/tad")').via).toEqual([]);
  });

  it('names a set as the union of its members, one type', () => {
    const r = resolve('/refs/sites');
    expect(r.result).toEqual({
      type: 'keys',
      keys: keySet([
        { owner: '/data/p53/domains', key: 'DBD' },
        { owner: '/data/p53/domains', key: 'TAD' },
      ]),
    });
    expect(resolve('/refs/mixed')).toMatchObject({ status: 'failed', issues: [{ code: 'seq.address.type' }] });
  });

  it('waits for a node data makes that is not held, and for data not loaded', () => {
    expect(resolve('/data/genes/(node :key T1)')).toMatchObject({ status: 'resolved' });
    const t3 = resolve('/refs/t3');
    expect(t3.status).toBe('waiting');
    expect(t3.issues.map((i) => [i.code, i.severity, i.address])).toEqual([
      ['seq.address.waiting', 'info', '/refs/t3'],
    ]);
    expect(t3.issues[0].places).toEqual([{ address: '/refs/t3' }, { address: '/data/genes/(node :key T3)' }]);
    const peaks = resolve('/data/p53/peaks/(item :id c)');
    expect(peaks.status).toBe('waiting');
    expect(peaks.complete).toBe(false);
  });

  it('reports what no longer resolves with what it named, never re-bound', () => {
    const r = resolve('/refs/gone');
    expect(r.status).toBe('missing');
    expect(r.issues).toEqual([
      expect.objectContaining({
        code: 'seq.address.missing',
        address: '/refs/gone',
        places: [{ address: '/refs/gone' }, { address: '/data/p53/gone' }],
      }),
    ]);
    expect(resolve('/data/p53/domains/(item :id GONE)').status).toBe('missing');
  });

  it('refuses a cycle, a chain past the bound, a set followed by segments and a malformed address', () => {
    expect(resolve('/refs/loop').issues.map((i) => i.code)).toEqual(['seq.address.ref.cycle']);
    expect(resolve('/refs/core', { refDepth: 1 }).issues.map((i) => i.code)).toEqual(['seq.address.ref.depth']);
    expect(resolve('/refs/sites/(range :from 1 :to 2)').issues.map((i) => i.code)).toEqual(['seq.address.ref.set']);
    expect(resolve('/data/p53/(range :from 1 :to 2)/(x').status).toBe('failed');
    expect(resolve('/data/(intersect :of [/refs/loop])').issues.map((i) => i.code)).toEqual(['seq.address.ref.cycle']);
  });

  it('is bounded: more matches than the limit is an issue, never a truncation', () => {
    const r = resolve('/refs/sites', { limit: 1 });
    expect(r.status).toBe('failed');
    expect(r.issues.map((i) => i.code)).toEqual(['repr.tree.address.limit_exceeded']);
    expect(r.result).toBeUndefined();
    expect(resolve('/data/*', { limit: 3 }).issues.map((i) => i.code)).toEqual(['repr.tree.address.limit_exceeded']);
  });
});

describe('what is named, as positions', () => {
  it('reads elements by their extents, a system whole, an annotation by its items', () => {
    const at = (a: string) => locationOfResult(Root, Model, resolve(a).result!).location;
    expect(at('/refs/tad')).toEqual(rangesLocation('/data/p53', [0, 61]));
    expect(at('/data/p53')).toEqual(rangesLocation('/data/p53', [0, 393]));
    expect(at('/data/plasmid/features')).toEqual(rangesLocation('/data/plasmid', [0, 5, 20, 40, 90, 100]));
  });
});

describe('the address of what is named', () => {
  const back = (target: QueryResult, options = {}) => addressFor(Root, Model, target, options);
  const at = (l: ReturnType<typeof location>) => back({ type: 'location', location: l });

  it('is a ref when one names exactly it, the first in the spec', () => {
    expect(at(rangesLocation('/data/p53', [93, 292]))).toBe('/refs/dbd');
    expect(back(resolve('/data/p53/domains/(item :id TAD)').result!)).toBe('/refs/tad');
    expect(at(rangesLocation('/data/p53', [93, 292])), 'ids preferred').toBe('/refs/dbd');
    expect(back({ type: 'location', location: rangesLocation('/data/p53', [93, 292]) }, { prefer: ['ids'] })).toBe(
      '/data/p53/(range :from 94 :to 292)',
    );
  });

  it('is ids and forms otherwise, in the numbering asked, reading back to the same thing', () => {
    expect(at(rangesLocation('/data/p53', [9, 10]))).toBe('/data/p53/(at :label 10)');
    expect(back({ type: 'location', location: rangesLocation('/data/p53', [9, 20]) }, { numbering: 'uniprot' })).toBe(
      '/data/p53/(range :from 10 :numbering uniprot :to 20)',
    );
    // an insertion code is a label of its own
    expect(back({ type: 'location', location: rangesLocation('/data/ab', [3, 5]) }, { numbering: 'kabat' })).toBe(
      '/data/ab/(range :from 82A :numbering kabat :to 82B)',
    );
    // labels out of index order (IMGT) read back
    expect(back({ type: 'location', location: rangesLocation('/data/ab', [4, 7]) }, { numbering: 'imgt' })).toBe(
      '/data/ab/(range :from 111A :numbering imgt :to 112)',
    );
    // a repeated label does not read back: ordinal instead
    expect(back({ type: 'location', location: rangesLocation('/data/p53', [2, 3]) }, { numbering: 'dup' })).toBe(
      '/data/p53/(at :label 3)',
    );
    // a hole in the numbering: ordinal instead
    expect(back({ type: 'location', location: rangesLocation('/data/ab', [7, 8]) }, { numbering: 'kabat' })).toBe(
      '/data/ab/(at :label 8)',
    );
  });

  it('wraps a cyclic system’s range over the origin, and joins several under `or`', () => {
    expect(at(rangesLocation('/data/plasmid', [0, 5, 90, 100]))).toBe('/data/plasmid/(range :from 91 :to 5)');
    expect(at(rangesLocation('/data/plasmid', [0, 100]))).toBe('/data/plasmid/(all)');
    expect(at(rangesLocation('/data/p53', [0, 5, 10, 11]))).toBe(
      '/data/p53/(merge :of [(range :from 1 :to 5) (at :label 11)])',
    );
    expect(
      at(
        location({
          parts: [
            { system: '/data/p53', ranges: [0, 5] },
            { system: '/data/ab', ranges: [0, 1] },
          ],
        }),
      ),
    ).toBe('/(merge :of [/data/ab/(at :label 1) /data/p53/(range :from 1 :to 5)])');
  });

  it('names a rect by its two sides', () => {
    const rect = location({
      rects: [{ a: { system: '/data/p53', ranges: [9, 25] }, b: { system: '/data/ab', ranges: [0, 4] } }],
    });
    const a = at(rect)!;
    expect(a).toBe('/(rect :a /data/p53/(range :from 10 :to 25) :b /data/ab/(range :from 1 :to 4))');
    expect(loc(resolve(a))).toEqual(rect);
  });

  it('names elements by their items and nodes by their addresses', () => {
    expect(
      back({
        type: 'keys',
        keys: keySet([
          { owner: '/data/p53/domains', key: 'TET' },
          { owner: '/data/p53/domains', key: 1 },
        ]),
      }),
    ).toBe('/data/p53/domains/(item :id [1 TET])');
    expect(back(resolve('/data/genes/(node :key T1)').result!)).toBe('/data/genes/(node :key T1)');
    expect(back(resolve('/data/1tup/(chain :id A)').result!)).toBe('/data/1tup/(chain :id A)');
    expect(back({ type: 'location', location: location({}) })).toBeUndefined();
  });

  it('resolves each ref once per root and model', () => {
    expect(refResults(Root, Model)).toBe(refResults(Root, Model));
    expect([...refResults(Root, Model).keys()]).toContain('/refs/t3');
  });
});
