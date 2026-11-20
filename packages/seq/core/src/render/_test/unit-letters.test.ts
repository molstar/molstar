/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A unit on a source's sequence, starting past 0 (EC-169): its letters are its sequence's over its range, read on the
// sequence's positions through the unit's step (pending while unread, never an empty plan); a form reading its
// standard components waits.

import type { Node } from '@molstar/core/repr/tree';
import { SeqRuntime } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { evaluateQuery, rangesLocation } from '@molstar/seqx-core/query';
import { projectRanges, resolveSpec } from '@molstar/seqx-core/resolve';
import { resolveView } from '@molstar/seqx-core/view';
import { describe, expect, it } from 'vitest';
import { layerSource } from '../layer-data.js';

function spec(): Node {
  const b = createSeqSpec();
  const src = b.data.source('src', { format: 'seq', uri: 's.seq' });
  const s = b.data.sequence('s', { alphabet: 'protein', source: src });
  let unit = '';
  b.data.complex('c', {}, (c) => {
    unit = c.unit('A', { of: s, range: [2, 8] });
  });
  b.views.view('v', (v) =>
    v
      .section('u', { of: unit as never })
      .track('t')
      .layer('l', { type: 'letters' }),
  );
  return b.build();
}

describe('a unit past 0 on a source’s sequence (EC-169)', () => {
  const root = spec();
  const model = resolveSpec(root, { runtimes: [SeqRuntime] });

  it('its letters are read on its sequence’s positions over its range, through the unit’s step', () => {
    const layer = resolveView(model, '/views/v')!.sections[0].tracks[0].layers[0];
    const src = layerSource(model, layer)!;
    expect([src.source, src.system, src.dataSystem]).toEqual(['/data/src', 's', '/data/s']);
    // the unit's positions 0–6 are the sequence's 2–8
    expect(projectRanges(src.back!, [0, 6])).toEqual([2, 8]);
  });

  it('a form reading its standard components waits for the letters; the modified rows need none', () => {
    const r = evaluateQuery(root, model, '/data/c/A/(component :id ALA)', { limit: 10 });
    expect(r.complete).toBe(false);
    expect(r.pending).toEqual(rangesLocation('/data/c/A', [0, 6]));
    expect(evaluateQuery(root, model, '/data/c/A/(component :standard false)', { limit: 10 }).complete).toBe(true);
  });
});
