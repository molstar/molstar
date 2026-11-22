/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A bridge laid out and drawn (_spec/README.md "Render", Bridges): its band between its sides, a ribbon a pair with its
// ends at the data's positions, crossed by GA's rule (EC-206), a stub, an unplaced pair said, the hit closed form.

import { describe, expect, it } from 'vitest';
import { bridgeAt, hitView } from '../index.js';
import { renderSvg } from '../render-svg.js';
import type { RibbonsLay } from '../pairs-ribbons.js';
import { bridgeSpec, type BridgeSpecOptions } from './bridge-fixtures.js';
import { settle } from './helpers.js';

const Size = { width: 900, height: 600 };

async function frame(o: BridgeSpecOptions, input: Record<string, unknown> = {}) {
  const s = await settle(bridgeSpec(o));
  const f = await s.renderer.frame({ size: Size, settled: true, ...input });
  return { s, f, lay: f.layout.bridges?.list[0]?.lays[0]?.lay as RibbonsLay | undefined };
}

/** The x of the boundary before position `p` on a side (its first section, `top` or `bottom`), no blocks. */
function xOf(f: Awaited<ReturnType<typeof frame>>['f'], side: 0 | 1, p: number): number {
  return f.layout.sections[side].axis!.x(p);
}

describe('a bridge between two sections', () => {
  it('takes its band between its sides, the sections’ gap above and below, its height by default 96', async () => {
    const { f, lay } = await frame({ pairs: [{ a: [10, 50], b: [20, 60] }] });
    const [top, bottom] = f.layout.sections;
    const br = f.layout.bridges!.list[0];
    expect(br.box.y).toBe(top.y + top.height);
    expect(br.box.y + br.box.height).toBe(bottom.y);
    const gap = f.layout.config.sectionGap;
    expect(bottom.y - (top.y + top.height)).toBe(96 + 2 * gap);
    expect(lay!.n).toBe(1);
  });

  it('draws a ribbon a pair with its ends at the data’s positions on each side', async () => {
    const { f, lay } = await frame({ pairs: [{ a: [10, 50], b: [20, 60] }] });
    const r = lay!.rib;
    expect(r[0]).toBeCloseTo(xOf(f, 0, 10), 6);
    expect(r[1]).toBeCloseTo(xOf(f, 0, 50), 6);
    expect(r[2]).toBe(lay!.edgeA);
    expect(r[3]).toBeCloseTo(xOf(f, 1, 20), 6);
    expect(r[4]).toBeCloseTo(xOf(f, 1, 60), 6);
    expect(r[5]).toBe(lay!.edgeB);
  });

  it('is crossed exactly where GA’s rule says, by geometry: a reversed run, or a reversed side, untwisted by both (EC-206)', async () => {
    const crossed = async (o: Partial<BridgeSpecOptions>) => {
      const { lay } = await frame({ pairs: [{ a: [10, 50], b: [20, 60], strand: o.pairs?.[0]?.strand }], ...o });
      const r = lay!.rib;
      return Math.sign(r[1] - r[0]) !== Math.sign(r[4] - r[3]);
    };
    expect(await crossed({})).toBe(false);
    expect(await crossed({ pairs: [{ a: [10, 50], b: [20, 60], strand: '-' }] })).toBe(true);
    expect(await crossed({ bottomReversed: true })).toBe(true);
    expect(await crossed({ pairs: [{ a: [10, 50], b: [20, 60], strand: '-' }], bottomReversed: true })).toBe(false);
    expect(await crossed({ topReversed: true, bottomReversed: true })).toBe(false);
  });

  it('places the data’s `a` end on the side that shows its system, whichever side is the bridge’s `a`', async () => {
    const { f, lay } = await frame({ pairs: [{ a: [10, 50], b: [20, 60] }], swapped: true });
    // the bridge's `a` is the bottom (y): the pair's `a` end (on x) is turned onto the top
    const r = lay!.rib;
    expect(lay!.orient[0]).toBe(-1);
    expect(r[0]).toBeCloseTo(xOf(f, 1, 20), 6);
    expect(r[3]).toBeCloseTo(xOf(f, 0, 10), 6);
  });

  it('says an end on a system neither side shows, once, and draws the rest', async () => {
    const { f, lay } = await frame({ pairs: [{ a: [10, 50], b: [20, 60], on: ['x', 'z'] }] });
    expect(lay!.n).toBe(0);
    expect(lay!.unplaced).toBe(1);
    expect(f.layout.issues.map((i) => i.code)).toContain('seq.render.bridge.unplaced');
  });

  it('draws a stub where an end’s system is shown on no block (a block’s range leaves it out)', async () => {
    const { lay } = await frame({ pairs: [{ a: [10, 50], b: [250, 280] }], bottom: [{ range: [0, 200] }] });
    expect(lay!.n).toBe(1);
    expect(lay!.flags[0] & 2).toBe(2);
    const r = lay!.rib;
    // a quarter of the band toward the far side, tapering to the placed end's middle
    expect(r[5] - r[2]).toBeCloseTo(0.25 * (lay!.edgeB - lay!.edgeA), 6);
    expect(r[3]).toBe(r[4]);
  });

  it('widens an end under a device px about its middle; both under it, one curve', async () => {
    const { lay } = await frame({ pairs: [{ a: [10, 11], b: [20, 21] }], lengths: [40_000, 30_000] });
    expect(lay!.flags[0] & 1).toBe(1);
    expect(Math.abs(lay!.rib[1] - lay!.rib[0])).toBeCloseTo(1, 6);
  });

  it('clamps an end past the window one data box past the side’s edge, continuous on a pan (EC-129)', async () => {
    const o = { pairs: [{ a: [10, 50], b: [250, 290] }] } as BridgeSpecOptions;
    const { f, lay } = await frame(o, { viewports: { '/views/v/bottom': { start: 0, end: 50 } } });
    const box = f.layout.sections[1].axis!;
    // both its corners clamped: the end widened to a device px about the clamp
    expect((lay!.rib[3] + lay!.rib[4]) / 2).toBe(box.x0 + 2 * box.width);
    const near = await frame(o, { viewports: { '/views/v/bottom': { start: 0, end: 200 } } });
    expect(near.lay!.rib[3]).toBeCloseTo(near.f.layout.sections[1].axis!.x(250), 6);
  });

  it('hits a ribbon at the point’s height in closed form, its pair both ends; the band is the bridge’s', async () => {
    const { f, lay } = await frame({ pairs: [{ a: [10, 50], b: [20, 60] }] });
    const r = lay!.rib;
    const y = (r[2] + r[5]) / 2;
    // at half height S is a half: the sides' middles
    const x = (r[0] + r[3]) / 4 + (r[1] + r[4]) / 4;
    const h = hitView(f.layout, x, y);
    expect(h?.bridge).toBe(0);
    expect(h?.section).toBe('br');
    expect(h?.element?.key).toBe('p0');
    expect(h?.extent?.parts.map((p) => p.system)).toEqual(['/data/x', '/data/y']);
    expect(bridgeAt(f.layout, x, y)).toBe(0);
    expect(hitView(f.layout, f.layout.bridges!.list[0].box.x + 1, y)).toBeUndefined();
    expect(bridgeAt(f.layout, x, 2)).toBe(-1);
  });

  it('reads a segments map as its pairs: a run a ribbon of equal length, a reversed run crossed', async () => {
    const { f, lay } = await frame({
      pairs: [
        { a: [10, 50], b: [20, 60] },
        { a: [100, 130], b: [200, 230], strand: '-' },
      ],
      map: true,
    });
    expect(f.layout.issues.filter((i) => i.severity !== 'info')).toEqual([]);
    expect(lay!.n).toBe(2);
    const byPair = (k: number) => 6 * [0, 1].find((i) => lay!.pair[i] === k)!;
    const o = byPair(0);
    expect(lay!.rib[o + 3]).toBeCloseTo(xOf(f, 1, 20), 6);
    const r = byPair(1);
    expect(Math.sign(lay!.rib[r + 1] - lay!.rib[r])).not.toBe(Math.sign(lay!.rib[r + 4] - lay!.rib[r + 3]));
    expect(
      hitView(f.layout, (lay!.rib[o] + lay!.rib[o + 3]) / 2 + 2, (lay!.rib[o + 2] + lay!.rib[o + 5]) / 2)?.element,
    ).toEqual({ owner: '/data/runs', key: 0 });
  });

  it('one genome against itself: a pair a ribbon from its `a` end above; with `mirror` its swap too', async () => {
    const pairs = [{ a: [10, 50], b: [200, 240], on: ['x', 'x'] }] as BridgeSpecOptions['pairs'];
    const plain = await frame({ pairs, self: true });
    expect(plain.lay!.n).toBe(1);
    expect(plain.lay!.rib[0]).toBeCloseTo(xOf(plain.f, 0, 10), 6);
    const mirrored = await frame({ pairs, self: true, layer: { mirror: true } });
    expect(mirrored.lay!.n).toBe(2);
    expect([...mirrored.lay!.orient.subarray(0, 2)].sort()).toEqual([-1, 1]);
  });

  it('exports a bridge as a group with its sides, its ribbons as cubic paths', async () => {
    const { svg } = await renderSvg(bridgeSpec({ pairs: [{ a: [10, 50], b: [20, 60] }] }), { size: Size });
    expect(svg).toContain('data-bridge="/views/v/br" data-a="/views/v/top" data-b="/views/v/bottom"');
    expect(svg).toContain('data-count="1"');
    expect(svg).toMatch(/<path d="M[^"]*C[^"]*L[^"]*C[^"]*Z"/);
  });

  it('leaves a hidden bridge out: no band, the sections stacked with their gap', async () => {
    const { f } = await frame({ pairs: [{ a: [10, 50], b: [20, 60] }], bridge: { hidden: true } });
    expect(f.layout.bridges).toBeUndefined();
    const [top, bottom] = f.layout.sections;
    expect(bottom.y - (top.y + top.height)).toBe(f.layout.config.sectionGap);
  });
});
