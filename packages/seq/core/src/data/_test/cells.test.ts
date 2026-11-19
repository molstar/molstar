/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { resolveSpec, rowsSystemOf } from '@molstar/seqx-core/resolve';
import { resolveView } from '@molstar/seqx-core/view';
import { AlignmentFileRuntime } from '../alignment-readers.js';
import { AlignmentRuntime, openAlignment } from '../alignment-source.js';
import type { CellsCodes } from '../protocol.js';
import { answerBytes, checkAnswer, joinTables, sliceTable } from '../tables.js';

const rows = { label: ['a/3-6', 'b', 'c'], letters: ['AC-D', 'ACCD', '-CCE'] };
const signal = new AbortController().signal;

describe('the cells shape (an alignment’s rows on row tiles, columns whole)', () => {
  it('an alignment source answers a tile’s rows × every column with its rows table; checked as an answer', async () => {
    const s = openAlignment(rows);
    expect(s.coverage.systems).toEqual([
      { name: 'rows', length: 3 },
      { name: 'columns', length: 4 },
    ]);
    const request = { system: 'rows', start: 1, end: 3, level: 1, shapes: ['cells' as const] };
    const a = await s.read(request, signal);
    expect(checkAnswer(a, request)).toBeUndefined();
    if (!a.ok) throw new Error('no answer');
    const c = a.tables.cells as CellsCodes;
    expect(String.fromCharCode(...c.codes)).toBe('ACCD-CCE');
    expect(c.label!.get(0)).toBe('b');
    expect(Array.from(c.first!)).toEqual([1, 1]);
    expect(answerBytes(a.tables, undefined)).toBeGreaterThan(8);
    const whole = await s.read({ ...request, start: 0 }, signal);
    expect(whole.ok && (whole.tables.cells as CellsCodes).first![0]).toBe(3);
    // cut and joined by rows
    const all = (whole.ok && whole.tables.cells) as CellsCodes;
    const a0 = sliceTable(all, 1, { start: 0, end: 1 });
    const a1 = sliceTable(all, 1, { start: 1, end: 3 });
    const joined = joinTables([a0, a1], 1) as CellsCodes;
    expect(Array.from(joined.codes)).toEqual(Array.from(all.codes));
    // past the request, refused
    expect(checkAnswer(whole, request)).toBe('the cells table: outside the region');
  });

  it('rows of another width are refused when opened; a file reader opens the same source', async () => {
    await expect(
      AlignmentRuntime.open(
        { address: '/data/a', format: 'alignment', value: { label: ['a'], letters: ['AC', 'A'] } as never },
        {
          assets: { readText: () => Promise.reject(new Error('no')) },
          signal,
        },
      ),
    ).rejects.toThrow();
    const s = await AlignmentFileRuntime.open(
      { address: '/data/a', format: 'aligned-fasta', value: '>a/3-6\nAC-D\n>b\nACCD\n' },
      { assets: { readText: () => Promise.reject(new Error('no')) }, signal },
    );
    expect(s.coverage.systems[0].length).toBe(2);
  });

  it('an alignment resolves to two systems, its cells an annotation on its columns; a rows layer reads it', () => {
    const b = createSeqSpec();
    const aln = b.data.alignment('aln', { alphabet: 'protein', rows });
    b.views.view('v', (v) =>
      v.section('s', { of: aln }, (s) => s.track('t', {}, (t) => t.layer('r', { type: 'rows', data: aln }))),
    );
    const spec = b.build();
    expect(b.check()).toEqual([]);
    const model = resolveSpec(spec);
    expect(model.systems.get('/data/aln')?.length).toBe(4);
    expect(model.systems.get(rowsSystemOf('/data/aln'))?.length).toBe(3);
    expect(model.annotations.get('/data/aln')).toMatchObject({
      shape: 'cells',
      target: '/data/aln',
      rows: '/data/aln/(rows)',
    });
    const view = resolveView(model, '/views/v')!;
    const layer = view.sections[0].tracks[0].layers[0];
    expect(layer).toMatchObject({ placed: true, data: { shape: 'cells', system: '/data/aln' } });
  });
});
