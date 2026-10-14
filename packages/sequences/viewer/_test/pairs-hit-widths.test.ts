/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// VL-24, LC-01 (seqspec:interaction#hover): under a width theme, of two arcs at one place in one color the thinner is
// on top, in the hit as in the draw, by the layer's own widths: another pairs track drawn after it (its widths
// reversed) changes nothing.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { expect, it } from 'vitest';
import { Hit } from '../src/draw/api.js';
import { drawn, layerOf, record } from './draw-fixtures.js';

const scale = (field: string, range: [number, number]) => ({ theme: 'field_scale', params: { field, range } });
/** Two arcs from 10 to 40, counts 1 and 8 (track A), and with `second` track B with the counts reversed. */
function doc(second: boolean) {
  const b = createSeqViewSpec();
  const a = b.data.sequence({ ref: 'a', value: 'MKVLATGHEK'.repeat(6) });
  const pairs = (ref: string, count: number[]) =>
    b.data.annotation({
      ref,
      geometry: 'pairs',
      target: a,
      value: { start: [10, 10], b_target: ['a', 'a'], b_start: [40, 40], count },
    } as never);
  const [pa, pb] = [pairs('pa', [1, 8]), pairs('pb', [8, 1])];
  const v = b.view({ ref: 'v', label: 'V' }).block({ sequence: a });
  const width = scale('count', [1, 16]);
  v.track({ label: 'A', height: 4 }).pairs({ data: pa, display: 'arcs', width } as never);
  if (second) v.track({ label: 'B', height: 4 }).pairs({ data: pb, display: 'arcs', width } as never);
  return b.getDocument();
}

it('the hit on a width-themed pairs track reads its own widths, whatever pairs track drew last', async () => {
  const items: string[] = [];
  for (const second of [false, true]) {
    const t = await drawn(doc(second), { width: 600, height: 400 });
    const { section } = layerOf(t.ctx, 'A');
    const { frame } = record(t.ctx, 'A');
    const [top, h] = [frame.v.top, frame.v.height];
    if (second) record(t.ctx, 'B');
    const ax = t.ctx.layout.axis(section);
    if (!ax) throw new Error('no axis');
    const x = (ax.x(0, ax.shownOf(0, 10)) + ax.x(0, ax.shownOf(0, 41))) / 2;
    const hit = Hit.create();
    const at: string[] = [];
    for (let y = Math.floor(top); y < top + h; y++)
      if (t.ctx.layout.hitTest(section, x, y, hit) && hit.item >= 0) at.push(`${y}:${hit.item}`);
    items.push(at.join(' '));
    t.ctx.dispose();
  }
  // The thin arc (item 0) is hit where it is drawn, over the wide one.
  expect(items[0]).toMatch(/:0/);
  expect(items[1]).toBe(items[0]);
});
