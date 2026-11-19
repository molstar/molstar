/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `renderPng`: the SVG export handed to a rasterizer port at its own size and the scale × dpr asked, its result and
// size returned; the render's issues kept. (The pixels against the SVG's picture: the render CLI's tests, over resvg.)

import { describe, expect, it } from 'vitest';
import { type RasterizeRequest, renderPng, svgSize } from '../export-png.js';
import { renderSvg } from '../render-svg.js';
import { p53TracksSpec } from './fixtures.js';

describe('renderPng', () => {
  it('rasterizes the SVG export at the scale and dpr asked and reports the pixel size', async () => {
    const o = { size: { width: 320, height: 100 }, viewport: { start: 0, end: 30 } };
    const seen: { svg: string; request: RasterizeRequest }[] = [];
    const bytes = Uint8Array.of(1, 2, 3);
    const r = await renderPng(p53TracksSpec(), {
      ...o,
      scale: 2,
      dpr: 1.5,
      rasterize: (svg, request) => {
        seen.push({ svg, request });
        return bytes;
      },
    });
    const svg = (await renderSvg(p53TracksSpec(), { ...o, dpr: 1.5 })).svg;
    expect(seen).toEqual([{ svg, request: { width: 320, height: 100, scale: 3 } }]);
    expect(r.svg).toBe(svg);
    expect(r.png).toBe(bytes);
    expect(r.pixels).toEqual({ width: 960, height: 300 });
    expect(r.settled).toBe(true);
    expect(r.issues).toEqual([]);
  });

  it('awaits an asynchronous rasterizer, and takes a height auto from the SVG', async () => {
    const r = await renderPng(p53TracksSpec(), {
      size: { width: 320, height: 'auto' },
      rasterize: async () => Uint8Array.of(9),
    });
    expect([...r.png]).toEqual([9]);
    expect(r.pixels.width).toBe(320);
    expect(r.pixels.height).toBe(Math.round(svgSize(r.svg).height));
    expect(r.pixels.height).toBeGreaterThan(0);
  });
});

describe('svgSize', () => {
  it('reads the root’s width and height, 0 for none', () => {
    expect(
      svgSize('<svg xmlns="http://www.w3.org/2000/svg" width="12.5" height="7" viewBox="0 0 12.5 7"></svg>'),
    ).toEqual({
      width: 12.5,
      height: 7,
    });
    expect(svgSize('<g/>')).toEqual({ width: 0, height: 0 });
  });
});
