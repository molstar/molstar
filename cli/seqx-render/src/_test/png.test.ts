/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `renderPng` over resvg: the PNG is the SVG's picture (probe pixels of each track's rects in both themes, at a scale),
// deterministic, with the bundled fonts only (text draws ink, the same on every machine).

import {
  bundledWidth,
  DefaultFonts,
  type DrawLists,
  DrawLists as Lists,
  drawListsToSvg,
  escapeXml,
  walkPass,
} from '@molstar/seqx-core/render';
import { renderPng } from '@molstar/seqx-core/render/export-png';
import { renderSvg } from '@molstar/seqx-core/render/render-svg';
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

/** The dark pixels' horizontal extent in a decoded picture (the ink of black text on white). */
function inkWidth(img: ReturnType<typeof decodePng>): number {
  let lo = img.width;
  let hi = -1;
  for (let y = 0; y < img.height; y++)
    for (let x = 0; x < img.width; x++)
      if (img.at(x, y)[0] < 128) {
        lo = Math.min(lo, x);
        hi = Math.max(hi, x);
      }
  return hi - lo + 1;
}

describe('fonts', () => {
  // the token table's own stacks, as the SVG writes them: the roles must reach the bundled families
  const ink = async (stack: string, role: 'sans' | 'mono') => {
    const svg =
      '<svg xmlns="http://www.w3.org/2000/svg" width="400" height="80" viewBox="0 0 400 80">' +
      '<rect width="400" height="80" fill="#ffffff"/>' +
      `<text x="20" y="60" font-family="${escapeXml(stack)}" font-size="40" fill="#000000">iiiiMMMM</text></svg>`;
    const png = await resvgRasterizer(svg, { width: 400, height: 80, scale: 1 });
    return { width: inkWidth(decodePng(png)), expected: bundledWidth(`40px ${role}`, 'iiiiMMMM') };
  };

  it('the mono role draws in the mono face (a fixed advance) and the sans role in the sans', async () => {
    const mono = await ink(DefaultFonts.mono, 'mono');
    const sans = await ink(DefaultFonts.sans, 'sans');
    // the ink is the advances less the side bearings of the first and last letters (under 0.2 em)
    for (const f of [mono, sans]) expect(f.expected - f.width, JSON.stringify(f)).toBeGreaterThanOrEqual(0);
    expect(mono.expected - mono.width).toBeLessThan(8);
    expect(sans.expected - sans.width).toBeLessThan(8);
    // the two roles differ (the wrong font gives one width for both)
    expect(mono.width - sans.width).toBeGreaterThan(15);
  });
});

describe('images', () => {
  it('stay crisp: a 2 × 2 image drawn at 40 × 40 keeps exact cell colors, not a blend', async () => {
    const lists = new Lists();
    lists.width = 40;
    lists.height = 40;
    // red, green / blue, white
    const data = Uint8Array.from([255, 0, 0, 255, 0, 255, 0, 255, 0, 0, 255, 255, 255, 255, 255, 255]);
    lists.writer('data').image(0, 0, 40, 40, { width: 2, height: 2, data });
    const svg = drawListsToSvg(lists);
    expect(svg).toContain('image-rendering="optimizeSpeed"');
    const img = decodePng(await resvgRasterizer(svg, { width: 40, height: 40, scale: 1 }));
    // away from the cells' centers, where a smoothed image would blend its neighbors
    expect(img.at(15, 15)).toEqual([255, 0, 0, 255]);
    expect(img.at(24, 15)).toEqual([0, 255, 0, 255]);
    expect(img.at(15, 24)).toEqual([0, 0, 255, 255]);
    expect(img.at(24, 24)).toEqual([255, 255, 255, 255]);
  });
});
