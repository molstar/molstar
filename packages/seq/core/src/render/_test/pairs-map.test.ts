/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A segments map as a pairs layer's data (_spec/README.md "Resolved view", Bridges; "Render", The pairs layer): its
// runs are pairs anywhere, in tracks an arc a run between its `from` and `to` ranges, keyed by its row.

import { describe, expect, it } from 'vitest';
import { createSeqSpec } from '../../lang/builder.js';
import { resolveSpec, segmentPairs } from '../../resolve/index.js';
import type { PairsLay } from '../index.js';
import { settle } from './helpers.js';

function spec() {
  const b = createSeqSpec();
  const x = b.data.sequence('x', { alphabet: 'dna', letters: 'ACGT'.repeat(30) });
  // a duplication within one sequence: [0, 20) to [60, 80) forward, [30, 40) to [109, 99] reversed
  const m = b.data.map('dup', {
    from: x,
    to: x,
    method: 'segments',
    segments: { from: [0, 30, 50], to: [60, 109, 0], length: [20, 10, 0], strand: ['+', '-', '+'] },
  });
  b.views.view('v', (v) => {
    const s = v.section('s', { of: x });
    s.track('letters').layer('letters', { type: 'letters' });
    s.track('arcs', { height: 60 }).layer('arcs', { type: 'pairs', data: m });
  });
  return b.build();
}

describe('a segments map as pairs', () => {
  it('reads a run a pair: `a` its from range, `b` its to range and strand, keyed by its row (an empty run left out)', () => {
    const p = segmentPairs({
      from: Int32Array.of(0, 30, 50),
      to: Int32Array.of(60, 109, 0),
      length: Int32Array.of(20, 10, 0),
      strand: ['+', '-', '+'],
    } as never);
    expect(p.n).toBe(2);
    expect([...p.aStart]).toEqual([0, 30]);
    expect([...p.aEnd]).toEqual([20, 40]);
    expect([...p.bStart]).toEqual([60, 100]);
    expect([...p.bEnd]).toEqual([80, 110]);
    expect([...p.bStrand]).toEqual([1, -1]);
    expect(p.keys).toEqual([0, 1]);
  });

  it('draws its runs as arcs in a track, each hit as its run', async () => {
    const model = resolveSpec(spec());
    expect(model.issues.filter((i) => i.severity === 'error')).toEqual([]);
    const s = await settle(spec());
    expect(s.view.issues).toEqual([]);
    const f = await s.renderer.frame({ size: { width: 800, height: 300 }, settled: true });
    const lay = f.layout.sections[0].sectionLays[0].lay as unknown as PairsLay;
    expect(lay.n).toBe(2);
    expect(lay.prep!.owner).toBe('/data/dup');
    expect([...lay.pair.subarray(0, 2)].sort()).toEqual([0, 1]);
  });
});
