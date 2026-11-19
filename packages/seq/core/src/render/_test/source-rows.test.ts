/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `(row :key k)` on an alignment read from a source (_spec/README.md "Overlay", The query model): a row names its
// position once the store holds the source's rows (the tiles a rows layer read), as an inline alignment's labels do;
// before, and without the store, the key is unknown.

import { AlignmentRuntime } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { resolveAddressIn } from '@molstar/seqx-core/query';
import { EmptyOverlay, overlayQueryModel } from '@molstar/seqx-core/resolve';
import { describe, expect, it } from 'vitest';
import { settle } from './helpers.js';

function spec() {
  const b = createSeqSpec();
  const file = b.data.source('file', {
    format: 'alignment',
    value: { label: ['a', 'b', 'c', 'b'], letters: ['ACDE', 'AC-E', 'A-DE', 'ACDF'] },
  });
  const aln = b.data.alignment('aln', { alphabet: 'protein', source: file });
  b.views.view('v', (v) => v.section('s', { of: aln }).track('rows').layer('rows', { type: 'rows', data: aln }));
  return b.build();
}

describe('a row form over an alignment read from a source', () => {
  it('names the first row of its key once the rows are held; unknown before and without the store', async () => {
    const s = await settle(spec(), [AlignmentRuntime]);
    const rows = '/data/aln/(rows)';
    expect(overlayQueryModel(s.model, EmptyOverlay).row?.(rows, 'b')).toBeUndefined();
    const q = overlayQueryModel(s.model, EmptyOverlay, s.store);
    // the rows are not read before a frame wants them
    expect(q.row?.(rows, 'b')).toBeUndefined();
    await s.frame({ size: { width: 400, height: 200 }, settled: true });
    expect(q.row?.(rows, 'b')).toBe(1);
    expect(q.row?.(rows, 'c')).toBe(2);
    expect(q.row?.(rows, 'x')).toBeUndefined();
    const r = resolveAddressIn(s.model.spec, q, '/data/aln/(row :key [c a])');
    expect(r.result?.type === 'location' && [...r.result.location.parts[0].ranges]).toEqual([0, 1, 2, 3]);
  });
});
