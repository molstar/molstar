/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `(row :key k)` on an alignment read from a source (_spec/README.md "Overlay", The query model): a row names its
// position once the store holds the source's rows, as an inline alignment's labels do; before, the form waits
// (`seq.query.row.pending`) and asks for the rows through whoever resolves (`RowsSource.want`), so no layer has to
// read them first; a key no row has is unknown.

import { AlignmentRuntime, type Want } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { resolveAddressIn } from '@molstar/seqx-core/query';
import { EmptyOverlay, overlayQueryModel } from '@molstar/seqx-core/resolve';
import { describe, expect, it } from 'vitest';
import { settle } from './helpers.js';

function spec(withLayer = true) {
  const b = createSeqSpec();
  const file = b.data.source('file', {
    format: 'alignment',
    value: { label: ['a', 'b', 'c', 'b'], letters: ['ACDE', 'AC-E', 'A-DE', 'ACDF'] },
  });
  const aln = b.data.alignment('aln', { alphabet: 'protein', source: file });
  b.views.view('v', (v) => {
    const s = v.section('s', { of: aln });
    if (withLayer) s.track('rows').layer('rows', { type: 'rows', data: aln });
  });
  return b.build();
}

const Rows = '/data/aln/(rows)';

describe('a row form over an alignment read from a source', () => {
  it('names the first row of its key once the rows are held; nothing without the store', async () => {
    const s = await settle(spec(), [AlignmentRuntime]);
    expect(overlayQueryModel(s.model, EmptyOverlay).row?.(Rows, 'b')).toBeUndefined();
    const src = { data: s.store };
    const q = overlayQueryModel(s.model, EmptyOverlay, src);
    expect(overlayQueryModel(s.model, EmptyOverlay, src)).toBe(q);
    await s.frame({ size: { width: 400, height: 200 }, settled: true });
    expect(q.rowsPending?.(Rows)).toBe(false);
    expect(q.row?.(Rows, 'b')).toBe(1);
    expect(q.row?.(Rows, 'c')).toBe(2);
    expect(q.row?.(Rows, 'x')).toBeUndefined();
    const r = resolveAddressIn(s.model.spec, q, '/data/aln/(row :key [c a])');
    expect(r.result?.type === 'location' && [...r.result.location.parts[0].ranges]).toEqual([0, 1, 2, 3]);
  });

  it('waits while the rows have not arrived (pending, never unknown), and asks for them', async () => {
    const s = await settle(spec(false), [AlignmentRuntime]);
    const wants: Want[] = [];
    const q = overlayQueryModel(s.model, EmptyOverlay, { data: s.store, want: (w) => wants.push(w) });
    const waiting = resolveAddressIn(s.model.spec, q, '/data/aln/(row :key b)');
    expect(waiting.status).toBe('waiting');
    expect(waiting.issues.map((i) => i.code)).toEqual(['seq.query.row.pending', 'seq.address.waiting']);
    // the whole rows system, asked for by the form (no layer reads the rows here)
    expect(wants.length).toBeGreaterThan(0);
    expect(wants.every((w) => w.source === '/data/file' && w.start === 0 && w.end === 4)).toBe(true);
    s.store.want(wants);
    await s.store.whenSettled();
    const r = resolveAddressIn(s.model.spec, q, '/data/aln/(row :key b)');
    expect(r.status).toBe('resolved');
    expect(r.result?.type === 'location' && [...r.result.location.parts[0].ranges]).toEqual([1, 2]);
  });

  it('says a key no row has is unknown once the rows are held', async () => {
    const s = await settle(spec(false), [AlignmentRuntime]);
    const q = overlayQueryModel(s.model, EmptyOverlay, { data: s.store, want: (w) => s.store.want([w]) });
    resolveAddressIn(s.model.spec, q, '/data/aln/(row :key x)');
    await s.store.whenSettled();
    const r = resolveAddressIn(s.model.spec, q, '/data/aln/(row :key x)');
    expect(r.status).toBe('failed');
    expect(r.issues.map((i) => i.code)).toEqual(['seq.query.row.unknown']);
  });
});
