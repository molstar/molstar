/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A layer's data filter (_spec/README.md "Themes", Filters): the query after its annotation leaves rows out as the
// shared item reader reads them (ranges, secondary, variants alike), the track shows how many in its "n filtered" chip
// and a hit on the chip names the layer; a query form not read row by row filters nothing and says so.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { hitView, type RangesLay, type TrackLayout, type VariantsLay } from '../index.js';
import { P53 } from './fixtures.js';
import { settle } from './helpers.js';

const Size = { width: 800, height: 120 };

/** A track's "n filtered" chip: its `seq:filtered` state. */
const filtered = (t: TrackLayout) => t.header?.states.find((s) => s.name === 'seq:filtered');
const rows = {
  start: [0, 63, 93, 318, 323],
  end: [61, 94, 292, 360, 355],
  id: ['a', 'b', 'c', 'd', 'e'],
  category: ['activation', 'linker', 'binding', 'linker', 'oligomer'],
};

function spec(type: 'ranges' | 'variants', query: string): Node {
  const b = createSeqSpec();
  const p = b.data.sequence('p', { alphabet: 'protein', letters: P53 });
  const d = b.data.ranges('d', { target: p, rows });
  b.views.view('main', (v) =>
    v
      .section('p', { of: p })
      .track('t', { label: 'Domains' })
      .layer('l', { type, data: `${d}${query}` }),
  );
  return b.build();
}

async function layout(type: 'ranges' | 'variants', query: string) {
  const s = await settle(spec(type, query));
  const f = await s.frame({ size: Size, settled: true });
  return f;
}

describe('a layer’s data filter', () => {
  it('leaves the rows out as the items are read, and the track counts them', async () => {
    const f = await layout('ranges', '/(except :of (item :category linker))');
    const t = f.layout.sections[0].tracks[0];
    const lay = t.layers[0].lay as RangesLay;
    expect(lay.keys.slice(0, lay.items).sort()).toEqual(['a', 'c', 'e']);
    const chip = filtered(t)!;
    expect(chip).toMatchObject({ count: 2, text: '2 filtered', address: '/views/main/p/t/l', mark: false });
    // on a line of its own under the label (the track is two lines high)
    expect(chip.y).toBeGreaterThan(t.labelY);
    const hit = hitView(f.layout, chip.x + 1, chip.box.y + chip.box.height / 2);
    expect(hit).toMatchObject({ target: 'header-state', state: 'seq:filtered', layer: '/views/main/p/t', count: 2 });
  });

  it('reads `item` by a list, `or` and `and`, in the other item layers too', async () => {
    const a = await layout('ranges', '/(item :category [linker oligomer])');
    expect((a.layout.sections[0].tracks[0].layers[0].lay as RangesLay).items).toBe(3);
    const b = await layout('variants', '/(merge :of [(item :id a) (item :id e)])');
    const v = b.layout.sections[0].tracks[0].layers[0].lay as VariantsLay;
    expect([...v.rows.keys.slice(0, v.rows.n)].sort()).toEqual(['a', 'e']);
    // a lone `or` of items is the track's own rows (the base), not what the user hid: not counted
    expect(filtered(b.layout.sections[0].tracks[0])).toBeUndefined();
  });

  it('no filter: no chip; a form not read row by row: every row, and said', async () => {
    const plain = await layout('ranges', '');
    expect(filtered(plain.layout.sections[0].tracks[0])).toBeUndefined();
    const range = await layout('ranges', '/(range :from 1 :to 100)');
    expect((range.layout.sections[0].tracks[0].layers[0].lay as RangesLay).items).toBe(5);
    expect(range.layout.issues.map((i) => i.code)).toContain('seq.render.filter');
  });
});
