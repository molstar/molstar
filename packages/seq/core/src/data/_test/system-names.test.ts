/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A source's systems by a requested name (_spec/README.md "Data", Systems by name): an alignment's columns are named
// `<id>/columns`, and an old source whose records say `columns` is still served under that name (the alias).

import { describe, expect, it } from 'vitest';
import { inlineSeq } from '../seq-format.js';
import { openSeq } from '../seq-runtime.js';
import { DataStore } from '../store.js';
import { coverageSystemOf, sourceNameOf } from '../system-names.js';

const signal = new AbortController().signal;
const ranges = { ranges: { form: 'list', start: Int32Array.of(2, 10), end: Int32Array.of(5, 14) } } as const;

describe('an alignment’s columns from an old source naming them `columns`', () => {
  it('are served by the `seq` reader under `<alignment>/columns`', async () => {
    const old = openSeq(inlineSeq({ name: 'columns', length: 40 }, ranges as never));
    const a = await old.read({ system: 'seed/columns', start: 0, end: 40, level: 1, shapes: ['ranges'] }, signal);
    expect(a.ok && a.tables.ranges?.form === 'list' && Array.from(a.tables.ranges.start)).toEqual([2, 10]);
    // a name of no alignment's columns is not served by the alias
    const other = await old.read({ system: 'seed', start: 0, end: 40, level: 1, shapes: ['ranges'] }, signal);
    expect(other.ok).toBe(false);
  });

  it('are held by the store under `<alignment>/columns`', async () => {
    const store = new DataStore();
    store.attach('/data/old', openSeq(inlineSeq({ name: 'columns', length: 40 }, ranges as never)));
    const want = { source: '/data/old', system: 'seed/columns', start: 0, end: 40, level: 1 };
    store.want([want]);
    await store.whenSettled();
    const keys = store.keysOf(want);
    expect(keys.length).toBeGreaterThan(0);
    expect(keys.every((k) => store.get(k)?.status === 'loaded')).toBe(true);
  });

  it('match a coverage and a names map by the alias last', () => {
    const systems = [
      { name: 'columns', length: 40 },
      { name: 'seed/columns', length: 30 },
    ];
    expect(coverageSystemOf(systems, 'seed/columns')?.length).toBe(30);
    expect(coverageSystemOf(systems.slice(0, 1), 'full/columns')?.length).toBe(40);
    expect(sourceNameOf(new Map([['columns', 'columns']]), 'full/columns')).toBe('columns');
    expect(sourceNameOf(new Map([['columns', 'columns']]), 'full')).toBeUndefined();
  });
});
