/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `renderPng` over resvg: the PNG is the SVG's picture (probe pixels of each track's rects in both themes, at a scale),
// deterministic, with the bundled fonts only (text draws ink, the same on every machine).

import { type DrawLists, renderPng, renderSvg, walkPass } from '@molstar/seqx-core/render';
import { describe, expect, it } from 'vitest';
import { resvgRasterizer } from '../rasterize.js';
import { decodePng, p53Tracks } from './fixture.js';

interface Probe {
  /** The track (clip block) the rect is in. */
  readonly block: number;
  readonly x: number;
  readonly y: number;
  readonly rgb: readonly [number, number, number];
}

/**
 * Up to `per` probes in each track of the data pass: an opaque rect's inner corner, in no later rect of the pass and
 * inside the figure (a letter's ink is at its middle, a clip may cut a rect's far side).
 */
function probesOf(lists: DrawLists, per = 4): Probe[] {
  const rects: { block: number; x: number; y: number; w: number; h: number; fill: number }[] = [];
  let block = -1;
  walkPass(lists.passes.data, {
    rect: (l, i) => rects.push({ block, x: l.x[i], y: l.y[i], w: l.w[i], h: l.h[i], fill: l.fill[i] }),
    path: () => {},
    glyphs: () => {},
    text: () => {},
    image: () => {},
    clip: () => block++,
    unclip: () => {},
  });
  const out: Probe[] = [];
  const taken = new Map<number, number>();
  rects.forEach((r, k) => {
    const x = r.x + 1.5;
    const y = r.y + 1.5;
    if ((r.fill & 0xff) !== 0xff || r.w < 4 || r.h < 4) return;
    if (x >= lists.width || y >= lists.height) return;
    if (rects.slice(k + 1).some((o) => x >= o.x && x < o.x + o.w && y >= o.y && y < o.y + o.h)) return;
    if ((taken.get(r.block) ?? 0) >= per) return;
    taken.set(r.block, (taken.get(r.block) ?? 0) + 1);
    out.push({ block: r.block, x, y, rgb: [(r.fill >>> 24) & 255, (r.fill >>> 16) & 255, (r.fill >>> 8) & 255] });
  });
  return out;
}

describe('renderPng', () => {
  for (const theme of ['light', 'dark'] as const)
    for (const scale of [1, 2]) {
      it(`${theme} at ${scale}x: the probe pixels of every track are the SVG's fills`, async () => {
        const r = await renderPng(p53Tracks(), {
          size: { width: 640, height: 'auto' },
          theme,
          scale,
          rasterize: resvgRasterizer,
        });
        expect(r.issues).toEqual([]);
        const img = decodePng(r.png);
        expect([img.width, img.height]).toEqual([r.pixels.width, r.pixels.height]);
        expect(img.width).toBe(Math.round(640 * scale));
        const probes = probesOf(r.lists);
        // letters, domains and conservation: a track each
        expect(new Set(probes.map((p) => p.block)).size).toBe(3);
        for (const p of probes) {
          const px = img.at(Math.floor(p.x * scale), Math.floor(p.y * scale));
          for (let c = 0; c < 3; c++)
            expect(
              Math.abs(px[c] - p.rgb[c]),
              `${theme} ${scale}x track ${p.block} at ${p.x},${p.y}`,
            ).toBeLessThanOrEqual(2);
          expect(px[3]).toBe(255);
        }
      });
    }

  it('draws the theme background: the SVG paints it, so the PNG has no clear pixels', async () => {
    const base = { size: { width: 200, height: 40 }, rasterize: resvgRasterizer };
    const dark = decodePng((await renderPng(p53Tracks(), { ...base, theme: 'dark' })).png);
    expect(dark.at(0, 0)).toEqual([0x12, 0x14, 0x19, 255]);
    const light = decodePng((await renderPng(p53Tracks(), base)).png);
    expect(light.at(0, 0)).toEqual([0xf1, 0xef, 0xea, 255]);
  });

  it('is deterministic: the same spec gives the same bytes', async () => {
    const o = { size: { width: 480, height: 'auto' as const }, rasterize: resvgRasterizer };
    const a = await renderPng(p53Tracks(), o);
    const b = await renderPng(p53Tracks(), o);
    expect(a.svg).toBe(b.svg);
    expect(Buffer.from(a.png).equals(Buffer.from(b.png))).toBe(true);
  });

  it('is the SVG export picture: the same SVG is rasterized', async () => {
    const o = { size: { width: 480, height: 'auto' as const } };
    const seen: string[] = [];
    const png = await renderPng(p53Tracks(), {
      ...o,
      rasterize: (svg, request) => {
        seen.push(svg);
        return resvgRasterizer(svg, request);
      },
    });
    expect(seen).toEqual([(await renderSvg(p53Tracks(), o)).svg]);
    expect(png.svg).toBe(seen[0]);
  });

  it('draws the text with the bundled fonts: a letter leaves ink in its cell', async () => {
    const r = await renderPng(p53Tracks(), {
      size: { width: 640, height: 'auto' },
      viewport: { start: 0, end: 30 },
      rasterize: resvgRasterizer,
    });
    const img = decodePng(r.png);
    expect(r.lists.passes.data.glyphs.n).toBeGreaterThan(0);
    // the first run of letters: around its first glyph's baseline the pixels are more than the cell's fill
    const g = r.lists.passes.data.glyphs;
    const x0 = Math.floor(g.x[g.start[0]]);
    const y0 = Math.floor(g.y[0]);
    const seen = new Set<string>();
    for (let y = y0 - 9; y <= y0; y++) for (let x = x0 - 4; x <= x0 + 4; x++) seen.add(img.at(x, y).join(','));
    expect(seen.size).toBeGreaterThan(2);
  });
});
