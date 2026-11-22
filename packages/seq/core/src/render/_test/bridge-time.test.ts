/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A bridge's time targets (_spec/README.md "Render", Bridges): a frame of 4,096 ribbons under 8 ms (a pan of one side:
// that side and the bridge laid out again), a hover over them under 2 ms. The bounds are the heavy lane's (`itWithin`,
// best of three); the quick lane runs each unbounded. A pan of one side lays out that side and the bridge only
// (counted, quick).

import { itWithin } from '@molstar/test-harness/timing';
import { describe, expect, it, vi } from 'vitest';
import { type FrameInput, hitView, TracksLayout } from '../index.js';
import { RibbonsType, type RibbonsLay } from '../pairs-ribbons.js';
import { bridgeSpec, randomBridgePairs } from './bridge-fixtures.js';
import { settle } from './helpers.js';

const Size = { width: 1440, height: 900 };
const N = 20_000;

const at = (k: number, settled = true): FrameInput => ({
  size: Size,
  dpr: 2,
  settled,
  viewports: { '/views/v/top': { start: 0, end: N }, '/views/v/bottom': { start: k, end: N - 3 + k } },
});

describe('4,096 ribbons', () => {
  const pairs = randomBridgePairs(4096, 23);

  itWithin('a pan frame of one side (that side and the bridge laid out again)', async (within) => {
    const s = await settle(bridgeSpec({ pairs, lengths: [N, N] }));
    const inputs = [0, 1, 2, 3].map((k) => at(k));
    for (const i of inputs) await s.frame(i);
    const lay = (await s.frame(inputs[0])).layout.bridges!.list[0].lays[0].lay as RibbonsLay;
    expect(lay.n).toBe(4096);
    let k = 0;
    const pan = () => {
      const f = s.renderer.frame(inputs[++k % inputs.length]);
      if (!f.rebuilt) throw new Error('a pan frame lays out again');
    };
    for (let i = 0; i < 8; i++) pan();
    within(8, pan);
  });

  itWithin('a hover over them', async (within) => {
    const s = await settle(bridgeSpec({ pairs, lengths: [N, N] }));
    const f = await s.frame(at(0));
    const lay = f.layout.bridges!.list[0].lays[0].lay as RibbonsLay;
    const points: [number, number][] = [];
    for (let i = 0; i < lay.n; i += 41) {
      const o = 6 * i;
      points.push([
        (lay.rib[o] + lay.rib[o + 1] + lay.rib[o + 3] + lay.rib[o + 4]) / 4,
        (lay.rib[o + 2] + lay.rib[o + 5]) / 2,
      ]);
    }
    for (const p of points) hitView(f.layout, p[0], p[1]);
    for (const p of points) within(2, () => hitView(f.layout, p[0], p[1]));
  });
});

describe('a pan of one side', () => {
  it('lays out that side and the bridge only; a frame with nothing moved, neither', async () => {
    const s = await settle(bridgeSpec({ pairs: randomBridgePairs(200, 5), lengths: [N, N] }));
    await s.frame(at(0));
    await s.frame(at(0));
    const place = vi.spyOn(TracksLayout, 'place');
    const ribbons = vi.spyOn(RibbonsType, 'layout');
    s.renderer.frame(at(5));
    expect(place).toHaveBeenCalledTimes(1);
    expect(place.mock.calls[0][0].address).toBe('/views/v/bottom');
    expect(ribbons).toHaveBeenCalledTimes(1);
    place.mockClear();
    ribbons.mockClear();
    s.renderer.frame({ ...at(5), emphasis: { hover: { layer: '', section: '', system: '' } } });
    expect(place).not.toHaveBeenCalled();
    expect(ribbons).not.toHaveBeenCalled();
    place.mockRestore();
    ribbons.mockRestore();
  });
});
