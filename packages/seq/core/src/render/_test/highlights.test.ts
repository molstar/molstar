/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Overlapping highlight channels stay apart (_spec/README.md "Render", Emphasis): where several light a position the
// band splits into equal stripes in channel order, at most `emphasis.highlightStripes`, the rest in the last; one alone
// takes the whole band; the styles `tint`, `band` and `outline`; tracks, text lines and rings, canvas lists and SVG,
// both themes.

import { rangesLocation } from '@molstar/seqx-core/query';
import { describe, expect, it } from 'vitest';
import { type EmphasisHighlight, recordDrawLists, renderSvg } from '../index.js';
import { p53TextSpec, p53TracksSpec } from './fixtures.js';
import { ringSpec } from './rings-fixtures.js';

const Protein = '/data/p53';
const Red = 0xff0000ff;
const Green = 0x00ff00ff;
const Blue = 0x0000ffff;
const Gold = 0xffcc00ff;
const Pink = 0xff66ccff;

const on = (system: string, range: [number, number], color: number, style?: EmphasisHighlight['style']) =>
  ({ location: rangesLocation(system, range), color, ...(style ? { style } : {}) }) satisfies EmphasisHighlight;

interface R {
  readonly x: number;
  readonly y: number;
  readonly w: number;
  readonly h: number;
  readonly paint: string;
}

/** The highlight pass's rects, in drawing order: `rect x y w h fill stroke [width]`. */
function listRects(lists: Parameters<typeof recordDrawLists>[0]): R[] {
  return recordDrawLists(lists, ['highlight'])
    .split('\n')
    .filter((l) => l.startsWith('rect '))
    .map((l) => {
      const [, x, y, w, h, ...paint] = l.split(' ');
      return { x: Number(x), y: Number(y), w: Number(w), h: Number(h), paint: paint.join(' ') };
    });
}

/** The same rects from the SVG: the fill and its opacity. */
function svgRects(svg: string): R[] {
  const out: R[] = [];
  for (const m of svg.matchAll(/<rect ([^>]*)\/>/g)) {
    const a = Object.fromEntries([...m[1].matchAll(/([\w-]+)="([^"]*)"/g)].map((x) => [x[1], x[2]]));
    if (a['fill-opacity'] !== '0.3') continue;
    out.push({
      x: Number(a.x),
      y: Number(a.y),
      w: Number(a.width),
      h: Number(a.height),
      paint: `${a.fill} ${a['fill-opacity']}`,
    });
  }
  return out;
}

const Tracks = { size: { width: 800, height: 240 }, viewport: { start: 0, end: 30 } };

async function tracks(highlights: EmphasisHighlight[], theme: 'light' | 'dark' = 'light', config = {}) {
  const r = await renderSvg(p53TracksSpec(), { ...Tracks, theme, config, emphasis: { highlights } });
  return { rects: listRects(r.lists), svg: svgRects(r.svg) };
}

describe.each(['light', 'dark'] as const)('overlapping channels in tracks (%s)', (theme) => {
  it('two channels: the whole band alone, equal stripes together, in channel order', async () => {
    const { rects, svg } = await tracks([on(Protein, [4, 10], Red), on(Protein, [6, 12], Green)], theme);
    const top = rects[0].y;
    const H = rects[0].h;
    // red alone, red over green, green under red, green alone
    expect(rects.map((r) => r.h)).toEqual([H, H / 2, H / 2, H]);
    expect(rects.map((r) => r.y)).toEqual([top, top, top + H / 2, top]);
    expect(rects.map((r) => r.paint.split(' ')[0])).toEqual(['#ff00004d', '#ff00004d', '#00ff004d', '#00ff004d']);
    // the stripes cover the same stretch; the lone parts the rest
    expect(rects[1].x).toBe(rects[2].x);
    expect(rects[1].w).toBe(rects[2].w);
    expect(rects[0].x + rects[0].w).toBe(rects[1].x);
    expect(rects[1].x + rects[1].w).toBe(rects[3].x);
    // the SVG draws the same geometry
    expect(svg.map((r) => [r.x, r.y, r.w, r.h])).toEqual(rects.map((r) => [r.x, r.y, r.w, r.h]));
    expect(svg.map((r) => r.paint)).toEqual(['#ff0000 0.3', '#ff0000 0.3', '#00ff00 0.3', '#00ff00 0.3']);
  });

  it('five channels on one stretch: four stripes, the fifth in the last', async () => {
    const hs = [Red, Green, Blue, Gold, Pink].map((c) => on(Protein, [4, 10], c));
    const { rects, svg } = await tracks(hs, theme);
    const top = rects[0].y;
    const H = rects[0].h * 4;
    expect(rects.length).toBe(5);
    expect(rects.map((r) => r.h)).toEqual([H / 4, H / 4, H / 4, H / 4, H / 4]);
    expect(rects.map((r) => r.y)).toEqual([0, 1, 2, 3, 3].map((k) => top + (k * H) / 4));
    expect(svg.map((r) => r.y)).toEqual(rects.map((r) => r.y));
  });

  it('the stripe count is a config value', async () => {
    const hs = [Red, Green, Blue].map((c) => on(Protein, [4, 10], c));
    const { rects } = await tracks(hs, theme, { emphasis: { highlightStripes: 2 } });
    expect(rects.map((r) => r.h)).toEqual([rects[0].h, rects[0].h, rects[0].h]);
    expect(rects[0].h * 2).toBeGreaterThan(rects[0].h);
    expect(new Set(rects.map((r) => r.y)).size).toBe(2);
  });

  it('channels that do not overlap each take the whole band', async () => {
    const { rects } = await tracks(
      [on(Protein, [4, 8], Red), on(Protein, [8, 12], Green), on(Protein, [20, 22], Blue)],
      theme,
    );
    expect(rects.length).toBe(3);
    expect(new Set(rects.map((r) => r.h)).size).toBe(1);
    expect(new Set(rects.map((r) => r.y)).size).toBe(1);
  });
});

describe('the styles', () => {
  it('a band is a bar along the foot; an outline a stroke without fill; the tint the default', async () => {
    const { rects } = await tracks([
      on(Protein, [4, 8], Red),
      on(Protein, [10, 14], Green, 'band'),
      on(Protein, [16, 20], Blue, 'outline'),
    ]);
    expect(rects.length).toBe(3);
    const [tint, band, outline] = rects;
    expect(tint.paint).toBe('#ff00004d');
    expect(band.paint).toBe('#00ff00ff');
    expect(band.h).toBe(3);
    expect(band.y + band.h).toBe(tint.y + tint.h);
    expect(outline.paint).toBe('none stroke #0000ffff 1');
  });

  it('overlapping channels of different styles each draw in their own stripe', async () => {
    const { rects } = await tracks([on(Protein, [4, 10], Red), on(Protein, [4, 10], Green, 'band')]);
    expect(rects.length).toBe(2);
    const [tint, band] = rects;
    // the tint fills the upper half; the band is the foot of the lower half
    expect(band.h).toBe(3);
    expect(band.y + band.h).toBe(tint.y + 2 * tint.h);
  });
});

describe('text lines and rings', () => {
  it('a text line splits its band into two stripes where the channels overlap', async () => {
    const r = await renderSvg(p53TextSpec(), {
      size: { width: 800, height: 'auto' },
      emphasis: { highlights: [on(Protein, [4, 10], Red), on(Protein, [6, 12], Green)] },
    });
    const rects = listRects(r.lists);
    expect(rects.map((x) => x.paint.split(' ')[0])).toEqual(['#ff00004d', '#ff00004d', '#00ff004d', '#00ff004d']);
    expect(rects[1].h).toBe(rects[0].h / 2);
    expect(rects[2].y).toBe(rects[0].y + rects[0].h / 2);
    const svg = svgRects(r.svg);
    expect(svg.map((x) => [x.x, x.y, x.w, x.h])).toEqual(rects.map((x) => [x.x, x.y, x.w, x.h]));
  });

  it('a ring splits its band radially: a path a stripe, in order, in the SVG too', async () => {
    const hs = [on('/data/construct', [10, 100], Red), on('/data/construct', [50, 150], Green)];
    const r = await renderSvg(ringSpec(), { size: { width: 600, height: 600 }, emphasis: { highlights: hs } });
    const lines = recordDrawLists(r.lists, ['highlight'])
      .split('\n')
      .filter((l) => l.startsWith('path '));
    expect(lines.map((l) => l.slice(l.lastIndexOf(' ') + 1))).toEqual([
      '#ff00004d',
      '#ff00004d',
      '#00ff004d',
      '#00ff004d',
    ]);
    const fills = [...r.svg.matchAll(/<path [^>]*fill="(#[0-9a-f]{6})" fill-opacity="0.3"\/>/g)].map((m) => m[1]);
    expect(fills).toEqual(['#ff0000', '#ff0000', '#00ff00', '#00ff00']);
  });
});
