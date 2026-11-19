/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { expect } from 'vitest';
import { itWithin } from '@molstar/test-harness/timing';
import { conservation, P53, SecondaryRows, VariantRows } from './fixtures.js';
import { settle } from './helpers.js';

// A frame of each small layer at the fit view under the frame target (16 ms at 1440 CSS px and dpr 2; _spec/README.md
// "Render", Layer types): in the quick lane the frame runs unbounded and its marks are counted, in the heavy lane it is
// timed (the best of three runs). The draw and hover loops' allocations are the opt-in tier's (`frame.quick.ts`).

const FRAME_MS = 16;
const Size = { width: 1440, height: 200 } as const;

type Builder = ReturnType<typeof createSeqSpec>;

function spec(layer: (b: Builder, p53: Parameters<Builder['data']['ranges']>[1]['target']) => Record<string, unknown>) {
  const b = createSeqSpec({ title: 'frame' });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
  const params = layer(b, p53);
  b.views.view('main', (v) => {
    v.section('protein', { of: p53 })
      .track('t')
      .layer('l', params as never);
  });
  return b.build();
}

const cases: readonly (readonly [string, () => ReturnType<typeof spec>])[] = [
  ['ruler', () => spec(() => ({ type: 'ruler' }))],
  [
    'secondary',
    () => spec((b, p) => ({ type: 'secondary', data: b.data.ranges('ss', { target: p, rows: SecondaryRows }) })),
  ],
  [
    'variants',
    () => spec((b, p) => ({ type: 'variants', data: b.data.ranges('vs', { target: p, rows: VariantRows }) })),
  ],
  [
    'values, a line',
    () =>
      spec((b, p) => ({
        type: 'values',
        encoding: 'line',
        data: b.data.values('c', { target: p, values: conservation() }),
      })),
  ],
  [
    'values, an area',
    () =>
      spec((b, p) => ({
        type: 'values',
        encoding: 'area',
        data: b.data.values('c', { target: p, values: conservation() }),
      })),
  ],
  [
    'ranges, dashed with arrows',
    () =>
      spec((b, p) => ({
        type: 'ranges',
        line: 'dashed',
        arrow: 'both',
        data: b.data.ranges('ss', { target: p, rows: SecondaryRows }),
      })),
  ],
];

for (const [name, build] of cases)
  itWithin(`${name}: a frame at the fit view`, async (within) => {
    const s = await settle(build());
    await s.frame({ size: Size, dpr: 2 });
    for (const window of [
      { start: 0.5, end: 392.5 },
      { start: 3, end: 390 },
    ]) {
      const f = within(FRAME_MS, () => s.renderer.frame({ size: Size, dpr: 2, viewport: window }));
      expect(f.layout.issues).toEqual([]);
      expect(f.lists.passes.data.size).toBeGreaterThan(0);
    }
  });
