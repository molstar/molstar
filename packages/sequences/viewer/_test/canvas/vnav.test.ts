/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The vertical overview and the scrollbar as pixels (CS2), in Chromium: the two passes' canvases lie in the box the
// rows' band takes at the section's right edge; the overview's ground and its map (the rows layer's `overview`, drawn
// into an image kept while its key holds); its window washed outside and outlined in two tones; the scrollbar's thumb
// at its share, darker while dragged; both pairs of grips as whole device pixels, always the
// columns' on the overview and the rows' on its window and on the thumb with `verticalZoom`, the one whose edge zone
// holds the pointer lit, both over the body, none outside (the old stack tested the horizontal rail's only). Skin
// colors made distinct to read them exactly.

import type { SeqSet } from '@molstar/seq-core/model';
import { afterEach, describe, expect, it } from 'vitest';
import type { AxisView } from '../../src/axis.js';
import { SCROLLBAR } from '../../src/hit.js';
import { Inputs, PointerType, type SectionLayout } from '../../src/index.js';
import { WINDOW_FRAME, WindowKind, windowFrame } from '../../src/windows.js';
import { syntheticFamily } from '../families.js';
import { rowLoci } from '../overlay-docs.js';
import { canvasesOf, close, type Opened, opened, pixelsOf, settle } from './open.js';

let t: Opened | null = null;
afterEach(() => {
  if (t) close(t);
  t = null;
});

const SKIN =
  '--seq-surface:rgb(0 255 0);--seq-surface-alt:rgb(255 255 0);--seq-ink-faint:rgb(200 200 200);' +
  '--seq-ink-quiet:rgb(50 50 50);--seq-grip:rgb(255 0 0);--seq-grip-hover:rgb(0 0 255);' +
  '--seq-selection-border:rgb(255 0 255);--seq-highlight:rgb(0 255 255)';
const [GROUND, FAINT, QUIET, REST, LIT] = ['0,255,0', '200,200,200', '50,50,50', '255,0,0', '0,0,255'];

async function msa(dpr: number, o: { zoom?: boolean; rows?: number; height?: number; skin?: string } = {}) {
  const x = await opened(syntheticFamily(o.rows ?? 120, 8), {
    width: 801,
    height: o.height ?? 300,
    dpr,
    render: { verticalZoom: o.zoom ?? false },
    style: o.skin ?? SKIN,
  });
  settle(x);
  return x;
}

/** One pass's canvas as a reader: "r,g,b" opaque, "r,g,b/a" translucent, "" clear. */
function reader(c: HTMLCanvasElement) {
  const px = pixelsOf(c);
  return (x: number, y: number) => {
    const i = (y * c.width + x) * 4;
    return px[i + 3] === 0
      ? ''
      : px[i + 3] === 255
        ? `${px[i]},${px[i + 1]},${px[i + 2]}`
        : `${px[i]},${px[i + 1]},${px[i + 2]}/${px[i + 3]}`;
  };
}
const dprOf = (d: number) => d;
const VNAV = 4;
const OVER = 5;

/** The geometry in device px of the box's local frame: the overview's and the scrollbar's windows, spans. */
function geometry(o: Opened, dpr: number) {
  const l = o.ctx.layout.get(o.section)?.current as SectionLayout;
  const v = o.ctx.layout.axis(o.section) as AxisView;
  const ox = Math.round((l.width - l.vnavWidth - SCROLLBAR) * dpr);
  const oy = Math.round(l.scrolledTop * dpr);
  const f = new Float64Array(WINDOW_FRAME);
  const s = new Float64Array(WINDOW_FRAME);
  windowFrame(v, l, WindowKind.overview, f);
  windowFrame(v, l, WindowKind.scrollbar, s);
  const px = (cssPx: number, origin: number) => Math.round(cssPx * dpr) - origin;
  return {
    l,
    f,
    s,
    ox,
    oy,
    // The overview window as drawn, local device px.
    x0: px(f[8], ox),
    x1: px(f[9], ox),
    y0: px(f[10], oy),
    y1: px(f[11], oy),
    // The span the grips stay in.
    ovx: Math.round((l.width - l.vnavWidth) * dpr) - ox,
    right: Math.ceil(l.width * dpr) - ox,
    bottom: Math.ceil(l.height * dpr) - oy,
    lw: Math.max(1, Math.round(dpr)),
    t3: Math.max(2, Math.round(3 * dpr)),
  };
}

describe('the box and the overview', () => {
  for (const dpr of [1, 2])
    it(`the passes lie in the rows' band at the right edge: the scrollbar and the overview (DPR ${dpr})`, async () => {
      t = await msa(dpr);
      const o = t;
      const c = canvasesOf(o);
      const g = geometry(o, dpr);
      for (const k of [VNAV, OVER]) {
        expect([c[k].width, c[k].height]).toEqual([g.right, g.bottom]);
        expect([c[k].style.left, c[k].style.top]).toEqual([`${g.ox / dpr}px`, `${g.oy / dpr}px`]);
        expect([Number.parseFloat(c[k].style.width) * dpr, Number.parseFloat(c[k].style.height) * dpr]).toEqual([
          g.right,
          g.bottom,
        ]);
      }
      // The pinned rows above and the navigator's canvases are the full width; nothing else covers the rows' band.
      expect(c[0].width).toBe(Math.ceil(801 * dpr));
    });

  it('no box while the rows fit: both canvases empty', async () => {
    t = await msa(1, { rows: 6 });
    const c = canvasesOf(t);
    expect([c[VNAV].width, c[VNAV].height, c[OVER].width, c[OVER].height]).toEqual([0, 0, 0, 0]);
  });

  for (const dpr of [1, 2])
    it(`the ground under the map, the map in the rows' colors, the scrollbar's strip left clear (DPR ${dpr})`, async () => {
      t = await msa(dpr);
      const g = geometry(t, dpr);
      const px = reader(canvasesOf(t)[VNAV]);
      expect(px(2, 4)).toBe('');
      expect(px(g.ovx - 1, 4)).toBe('');
      const colors = new Set<string>();
      for (let y = 0; y < g.bottom; y += 3) for (let x = g.ovx; x < g.right; x += 2) colors.add(px(x, y));
      // The map has the rows' cells (several classes), not just the ground.
      expect(colors.size).toBeGreaterThan(4);
      // Drawn again on a scroll without rebuilding the image: the same pixels.
      const before = pixelsOf(canvasesOf(t)[VNAV]);
      t.ctx.viewports.get(t.section)?.scrollTo(0, 400, 'api');
      settle(t);
      expect(pixelsOf(canvasesOf(t)[VNAV])).toEqual(before);
    });

  it('the window outside washed, in two tones inside (a light line outside a dark one), the rows lit by a selection marked', async () => {
    t = await msa(2);
    const o = t;
    const g = geometry(o, 2);
    const px = reader(canvasesOf(o)[OVER]);
    // Below the left grip (centered on the window's height).
    const mid = g.y0 + 3 * dprOf(2);
    // The window's left edge: light, then dark, each `lw` wide; just inside, nothing drawn (the map shows through).
    expect([px(g.x0, mid), px(g.x0 + g.lw, mid), px(g.x0 + 2 * g.lw + 1, mid)]).toEqual([
      expect.stringMatching(/^255,255,255\/24\d$/),
      expect.stringMatching(/^1[67],19,2[34]\/217$/),
      '',
    ]);
    // Below the window the map is washed with the ground (translucent), inside it is not.
    expect(px(g.ovx + 30, g.y1 + 20)).toMatch(/^255,255,0\/158$/);
    expect(px(g.ovx + 30, (g.y0 + g.y1) >> 1)).toBe('');
  });

  it('the rows a selection or a highlight lights are marked across the overview, at least two device px each', async () => {
    t = await msa(2);
    const o = t;
    const g = geometry(o, 2);
    const set = o.ctx.graph.assembled.value?.set as SeqSet;
    o.ctx.selection.set(rowLoci(set, 60, 0, 10), 'api');
    o.ctx.highlight.set('host', rowLoci(set, 30, 0, 10));
    settle(o);
    const px = reader(canvasesOf(o)[OVER]);
    const at = (row: number) => Math.floor((row * g.bottom) / 120);
    const x = g.ovx + 20;
    expect([px(x, at(60)), px(x, at(60) + 1)]).toEqual(['255,0,255', '255,0,255']);
    expect([px(x, at(30)), px(x, at(30) + 1)]).toEqual(['0,255,255', '0,255,255']);
    // Elsewhere the map shows through (washed, or the window's own: nothing drawn here).
    expect(px(x, at(90))).toMatch(/^255,255,0\/158$/);
    // Cleared, the marks go.
    o.ctx.selection.clear('api');
    o.ctx.highlight.clear('host');
    settle(o);
    expect(reader(canvasesOf(o)[OVER])(x, at(60))).toMatch(/^255,255,0\/158$/);
  });
});

describe('the scrollbar', () => {
  for (const dpr of [1, 2])
    it(`the thumb at its share of the span, in the faint tone, darker while it is dragged (DPR ${dpr})`, async () => {
      t = await msa(dpr);
      const o = t;
      o.ctx.viewports.get(o.section)?.scrollTo(0, 900, 'api');
      settle(o);
      const g = geometry(o, dpr);
      const px = reader(canvasesOf(o)[OVER]);
      const bx = Math.round((g.l.width - g.l.vnavWidth - SCROLLBAR) * dpr) - g.ox + Math.round(SCROLLBAR * dpr) / 2;
      const inset = g.s[11] - g.s[10] >= 4 ? 1 : 0;
      const top = Math.round((g.s[10] + inset) * dpr) - g.oy;
      const bottom = Math.round((g.s[11] - inset) * dpr) - g.oy;
      const x = Math.floor(bx);
      expect([px(x, top), px(x, (top + bottom) >> 1), px(x, bottom - 1)]).toEqual([FAINT, FAINT, FAINT]);
      expect(px(x, top - 2)).toBe('');
      // The strip's edges: 2 CSS px of margin each side.
      expect([px(1, (top + bottom) >> 1), px(Math.round((SCROLLBAR - 1) * dpr), (top + bottom) >> 1)]).toEqual([
        '',
        '',
      ]);
      // A press on the thumb and a move: the scrollbar owns the drag, the thumb darkens.
      const p = Inputs.pointer();
      p.section = o.ctx.layout.sectionIndex(o.section);
      p.type = PointerType.down;
      p.x = g.l.width - g.l.vnavWidth - SCROLLBAR / 2;
      p.y = g.l.scrolledTop + (top + bottom) / 2 / dpr;
      p.button = 0;
      p.buttons = 1;
      o.ctx.input.pointer(p);
      p.type = PointerType.move;
      p.y += 20;
      o.ctx.input.pointer(p);
      settle(o);
      expect(o.ctx.ui.drag.value?.owner).toBe('scrollbar');
      const px2 = reader(canvasesOf(o)[OVER]);
      const moved = Math.round((geometry(o, dpr).s[10] + 1) * dpr) - g.oy;
      expect(px2(x, moved + 3)).toBe(QUIET);
    });
});

describe('the end grips', () => {
  /** The pill's four pixel columns or rows through its center, and the halo outside them. */
  for (const dpr of [1, 2])
    it(`the overview's column grips sit at the window's ends and its row grips with verticalZoom, on whole device pixels (DPR ${dpr})`, async () => {
      t = await msa(dpr, { zoom: true });
      const o = t;
      const g = geometry(o, dpr);
      const px = reader(canvasesOf(o)[OVER]);
      const cy = (g.y0 + g.y1) >> 1;
      // Left and right grips (columns): at the window's ends, on the span's far ends at the fit.
      const a = Math.max(g.ovx, g.x0 - (g.t3 >> 1));
      const e = Math.min(g.right - g.t3, g.x1 - (g.t3 >> 1));
      expect([a, e + g.t3]).toEqual([g.ovx, g.right]);
      for (const x of [a, e]) expect([px(x, cy), px(x + g.t3 - 1, cy)], `column grip at ${x}`).toEqual([REST, REST]);
      // Row grips: a pill GRIP_ACROSS wide across the window's top and bottom, centered on x.
      const rx = (g.x0 + g.x1) >> 1;
      const ty = Math.max(0, g.y0 - (g.t3 >> 1));
      const by = Math.min(g.bottom - g.t3, g.y1 - (g.t3 >> 1));
      for (const y of [ty, by]) expect([px(rx, y), px(rx, y + g.t3 - 1)], `row grip at ${y}`).toEqual([REST, REST]);
      // Crisp: the ground right outside a pill's end along its axis, across the halo.
      expect(px(rx, ty + g.t3 + g.lw - 1)).toBe(GROUND);
      expect([px(rx - Math.round(10 * dpr), ty), px(rx + Math.round(10 * dpr), ty)]).not.toContain(REST);
    });

  it('without verticalZoom the window has column grips only, and the thumb none', async () => {
    t = await msa(2);
    const g = geometry(t, 2);
    const px = reader(canvasesOf(t)[OVER]);
    const rx = (g.x0 + g.x1) >> 1;
    const ty = Math.max(0, g.y0 - (g.t3 >> 1));
    expect([px(rx, ty), px(rx, ty + g.t3 - 1)]).not.toContain(REST);
    const cy = (g.y0 + g.y1) >> 1;
    expect(px(g.ovx, cy)).toBe(REST);
    // No grip on the thumb (its strip holds only the thumb's tone).
    for (let y = 0; y < g.bottom; y++) expect(px(Math.round(SCROLLBAR), y)).not.toBe(REST);
  });

  it("the thumb's row grips (verticalZoom) at its ends, centered across the strip", async () => {
    t = await msa(2, { zoom: true });
    const g = geometry(t, 2);
    const px = reader(canvasesOf(t)[OVER]);
    const inset = g.s[11] - g.s[10] >= 4 ? 1 : 0;
    const top = Math.round((g.s[10] + inset) * 2) - g.oy;
    const x = Math.round(SCROLLBAR * 2) >> 1;
    // Rest: the first grip's rows at the thumb's top are the grip tone (the pill is as wide as the strip allows).
    const rest = [px(x, Math.max(0, top - (g.t3 >> 1))), px(x, Math.max(0, top - (g.t3 >> 1)) + g.t3 - 1)];
    expect(rest).toEqual([REST, REST]);
  });

  /** A pointer moved over the section, the lit window part set with the cursor, the overlay redrawn. */
  const hover = (o: Opened, x: number, y: number) => {
    const p = Inputs.pointer();
    p.type = PointerType.move;
    p.section = o.ctx.layout.sectionIndex(o.section);
    p.x = x;
    p.y = y;
    o.ctx.input.pointer(p);
    settle(o);
  };
  it('the grip whose edge zone holds the pointer is lit; both over the body; none outside', async () => {
    t = await msa(1, { zoom: true });
    const o = t;
    // Scrolled into the rows so the window floats, not at the span's ends.
    o.ctx.viewports.get(o.section)?.scrollTo(0, 800, 'api');
    settle(o);
    const g = geometry(o, 1);
    const cx = (g.x0 + g.x1) >> 1;
    const cy = (g.y0 + g.y1) >> 1;
    const top = g.y0 - (g.t3 >> 1);
    const bot = g.y1 - (g.t3 >> 1);
    const colorAt = () => {
      const px = reader(canvasesOf(o)[OVER]);
      return [px(cx, Math.max(0, top) + 1), px(cx, bot + 1)];
    };
    // Top edge zone (the window's y0 in surface px), bottom edge zone, the body, outside the window.
    hover(o, g.ox + cx, g.oy + g.y0 + 1);
    expect(colorAt()).toEqual([LIT, REST]);
    hover(o, g.ox + cx, g.oy + g.y1 - 1);
    expect(colorAt()).toEqual([REST, LIT]);
    hover(o, g.ox + cx, g.oy + cy);
    expect(colorAt()).toEqual([LIT, LIT]);
    hover(o, 100, g.oy + 20);
    expect(colorAt()).toEqual([REST, REST]);
  });
});
