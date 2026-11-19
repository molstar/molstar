/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A source's numbering read from the store (_spec/README.md "Render", Numberings from a source), over a hand-held
// store: bounded (past `maxNumberingPositions` nothing is asked for), held labels kept when a tile is evicted (never
// asked again), a failed tile and a source without the sequence's system reported, rows past the end and another
// count of positions reported once held, a field other than `label` and numeric labels; and the query model: pending
// as `waiting`, a unit's molecule's source numbering found through its parent.

import { categoryColumn, type Column } from '@molstar/core/repr/tree';
import type { Coverage, DataEntry, Want } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { findNumbering, type QueryModel, type QuerySystem, resolveAddressIn } from '@molstar/seqx-core/query';
import { resolveSpec } from '@molstar/seqx-core/resolve';
import { describe, expect, it } from 'vitest';
import type { StoreReader } from '../layer-type.js';
import { readSourceNumbering, withSourceNumberings } from '../numberings.js';

const Letters = 'MEEPQSDP';

function model(field?: string) {
  const b = createSeqSpec();
  const src = b.data.source('t', { format: 'json', system: 'h', value: [] });
  b.data.sequence('h', { alphabet: 'protein', letters: Letters }, (s) =>
    s.numbering('kabat', { method: 'labels', source: src, ...(field ? { field } : {}) }),
  );
  const m = resolveSpec(b.build());
  return { m, entry: m.numberingOf('/data/h', 'kabat')! };
}

function loaded(positions: number[], fields: Record<string, Column>, version = 1): DataEntry {
  return {
    status: 'loaded',
    version,
    region: { system: 'h', start: 0, end: Letters.length, level: 1 } as never,
    tables: {
      ranges: {
        form: 'list',
        start: Int32Array.from(positions),
        end: Int32Array.from(positions, (p) => p + 1),
        fields: Object.entries(fields).map(([name, column]) => ({ name, column })),
      },
    },
    nodes: [],
    bytes: 0,
  };
}

/** A store of one tile: its coverage and entry are the test's to change. */
function store(entry: DataEntry | undefined, system = 'h', length = Letters.length) {
  const wants: Want[] = [];
  const s = {
    coverage: { systems: [{ name: system, length }], shapes: [], levels: [1] } as Coverage | undefined,
    entry,
    wants,
    reader: undefined as unknown as StoreReader,
  };
  s.reader = { coverage: () => s.coverage, keysOf: () => ['t|h|1|0'], get: () => s.entry };
  return s;
}

const labels = (...l: (string | null)[]) => categoryColumn(l);
const all = [0, 1, 2, 3, 4, 5, 6, 7];

describe('a source’s numbering', () => {
  it('past the bound: reported, never asked for', () => {
    const { m, entry } = model();
    const s = store(undefined);
    const r = readSourceNumbering(m, entry, s.reader, { maxPositions: 4, want: (w) => s.wants.push(w) });
    expect(r.status).toBe('blank');
    expect(r.status === 'blank' && r.issue.code).toBe('seq.render.numbering.limit');
    expect(s.wants).toEqual([]);
  });

  it('held labels are kept when a tile is evicted: not asked for again', () => {
    const { m, entry } = model();
    const s = store(loaded(all, { label: labels('1', '2', '3', '4', '5', '6', '7', '8') }));
    const want = (w: Want) => s.wants.push(w);
    const held = readSourceNumbering(m, entry, s.reader, { want });
    expect(held.status === 'held' && held.numbering.labelAt(7)).toBe('8');
    expect(s.wants).toHaveLength(1);
    s.entry = { status: 'evicted', version: 2, region: s.entry!.region };
    for (let i = 0; i < 5; i++) expect(readSourceNumbering(m, entry, s.reader, { want })).toBe(held);
    expect(s.wants).toHaveLength(1);
  });

  it('a failed tile, and a source without the sequence’s system, are reported', () => {
    const { m, entry } = model();
    const failed = store({
      status: 'failed',
      version: 1,
      region: {} as never,
      kind: 'failed',
      issue: { code: 'x', severity: 'error', message: 'gone' },
    } as DataEntry);
    const r = readSourceNumbering(m, entry, failed.reader);
    expect(r.status === 'blank' && [r.issue.code, r.issue.message]).toEqual([
      'seq.render.numbering.source',
      'A tile of the source failed: gone',
    ]);
    const other = store(undefined, 'chr1');
    const o = readSourceNumbering(m, entry, other.reader);
    expect(o.status === 'blank' && [o.issue.code, o.issue.message]).toEqual([
      'seq.render.numbering.source',
      'The source has no system h',
    ]);
  });

  it('rows past the end and rows for another count of positions are reported once held', () => {
    const { m, entry } = model();
    const s = store(loaded([0, 1, 2, 9, 12], { label: labels('1', '2', '3', 'x', 'y') }));
    const r = readSourceNumbering(m, entry, s.reader);
    expect(r.status === 'held' && r.issues.map((i) => [i.code, i.message])).toEqual([
      ['seq.render.numbering.rows', '2 rows of the source are past the sequence’s 8 positions'],
      ['seq.render.numbering.rows', 'The source has rows for 3 of 8 positions'],
    ]);
    // the held issues are the same objects on the next read: said once, not made again
    expect(readSourceNumbering(m, entry, s.reader)).toBe(r);
  });

  it('a field other than `label`, and numeric labels', () => {
    const { m, entry } = model('kabat');
    const s = store(
      loaded(all, { label: labels(...all.map(() => 'no')), kabat: Float64Array.of(1, 2, 3, 52, 53, 82, 83, 100) }),
    );
    const r = readSourceNumbering(m, entry, s.reader);
    expect(r.status === 'held' && Array.from({ length: 8 }, (_, i) => r.numbering.labelAt(i))).toEqual([
      '1',
      '2',
      '3',
      '52',
      '53',
      '82',
      '83',
      '100',
    ]);
    expect(r.status === 'held' && r.numbering.indexOf('82')).toBe(5);
  });
});

describe('a source’s numbering in the query model', () => {
  it('pending: a form waits, and its tiles are asked for; held: it resolves', () => {
    const { m } = model();
    const s = store({ status: 'pending', version: 1, region: {} as never });
    const q = withSourceNumberings(m, m, s.reader, { want: (w) => s.wants.push(w) });
    const r = resolveAddressIn(m.spec, q, '/data/h/(at :label 2 :numbering kabat)');
    expect(r.status).toBe('waiting');
    expect(r.issues.map((i) => i.code)).toEqual(['seq.query.numbering.pending', 'seq.address.waiting']);
    expect(s.wants.length).toBeGreaterThan(0);
    for (const w of s.wants) expect(w).toEqual({ source: '/data/t', system: 'h', start: 0, end: 8, level: 1 });
    s.entry = loaded(all, { label: labels('1', '2', '3', '4', '5', '6', '7', '8') });
    const held = resolveAddressIn(
      m.spec,
      withSourceNumberings(m, m, s.reader),
      '/data/h/(at :label 2 :numbering kabat)',
    );
    expect(held.status).toBe('resolved');
  });

  it('a blank one is left out of the ids an unknown numbering lists', () => {
    const { m } = model();
    const s = store(undefined, 'chr1');
    const r = resolveAddressIn(m.spec, withSourceNumberings(m, m, s.reader), '/data/h/(at :label 2 :numbering kabat)');
    expect(r.issues.map((i) => [i.code, i.hint])).toEqual([
      ['seq.query.numbering.unknown', 'numberings: index, ordinal'],
    ]);
  });

  it('wraps a unit’s molecule: the molecule’s source numbering numbers the unit', () => {
    const { m } = model();
    const s = store(loaded(all, { label: labels('1', '2', '3', '4', '5', '6', '7', '8') }));
    // a unit of positions 2 to 6 of the molecule, as a complex makes it
    const unit: QuerySystem = {
      address: '/data/c/u',
      length: 4,
      numbering: () => undefined,
      numberingIds: () => [],
      parent: { system: m.system('/data/h')!, start: 2 },
    };
    const base: QueryModel = { system: (a) => (a === unit.address ? unit : m.system(a)), annotation: () => undefined };
    const q = withSourceNumberings(base, m, s.reader);
    const n = findNumbering(q.system(unit.address)!, 'kabat');
    expect(n && [n.shift, n.numbering.labelAt(2 - n.shift)]).toEqual([-2, '5']);
    expect(findNumbering(base.system(unit.address)!, 'kabat')).toBeUndefined();
  });
});
