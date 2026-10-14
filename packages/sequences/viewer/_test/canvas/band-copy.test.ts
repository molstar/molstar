/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Copies on pan for the rows band (seqspec:viewer#the-render-backend; CS2 letters) as pixels in Chromium: a lazy
// alignment's band (one entry, its tiles from the in-process tile server) and an in-memory alignment's dense band are
// copied, and after random pans and scrolls the copied frame equals a full redraw, pixel for pixel: lettered rows at
// DPR 1 and 2 (through the glyph atlas and through fillText), rows of colored cells, thin dense rows, and at fractional
// zooms (a tile drawn at its whole scale, in parts: two full draws at different scrolls agree). A frame drawn while
// what a move uncovered is pending copies too (a pending tile's coarse level drawn in its place alone), and each frame
// redrawn as tiles and labels arrive equals a full redraw.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { tilesStream } from '@molstar/seq-core/stream/tiles';
import { tileServer } from '@molstar/seq-core/stream/tiles/server';
import { describe, expect, it } from 'vitest';
import { Canvas2DBackend, Canvas2DShared } from '../../src/canvas/index.js';
import { type DrawResources, SkinToken } from '../../src/draw/backend.js';
import { Dirty, SectionComposer } from '../../src/draw/index.js';
import {
  createSeqViewer,
  DefaultSeqViewerSpec,
  DefaultSkin,
  manualScheduler,
  type SeqViewerContext,
  type Skin,
} from '../../src/index.js';
import type { Streams } from '../../src/streams.js';
import { fixtureMeasurer } from '../../src/testing/index.js';
import { holdBack } from '../arrival-fixtures.js';

const URI = 'https://example.org/lazy';
const [W, H] = [700, 300];

function lazyDoc() {
  const b = createSeqViewSpec();
  const source = { uri: URI, protocol: 'seqstar-tiles/1' };
  const aln = b.data.alignment({ ref: 'lazy', row_count: 3000, column_count: 4000, source } as never);
  const v = b.view({ ref: 'g', label: 'G', height: 'box' });
  v.block({ alignment: aln });
  v.track({ label: 'Rows' }).rows();
  return b.getDocument();
}
/** An in-memory alignment of 400 rows: its band dense below `denseRowHeight`. */
function memoryDoc() {
  const b = createSeqViewSpec();
  const L = 'ACDEFGHIKLMNPQRSTVWY-';
  const value = Array.from({ length: 400 }, (_, r) =>
    Array.from(
      { length: 900 },
      (_, c) => L[((Math.imul(r + 7, 2654435761) ^ Math.imul(c + 3, 40503)) >>> 0) % L.length],
    ).join(''),
  );
  const aln = b.data.alignment({ ref: 'aln', value });
  const v = b.view({ ref: 'g', label: 'G', height: 'box' });
  v.block({ alignment: aln });
  v.track({ label: 'Rows' }).rows();
  return b.getDocument();
}

/** The in-process tile server behind the lazy alignment. */
function lazyLoaders() {
  const server = tileServer({ rowCount: 3000, columnCount: 4000 });
  return {
    load: () => Promise.reject(new Error('no fixture')),
    stream: tilesStream((u) => (u === URI ? server.transport : null)),
  };
}
async function open(doc: ReturnType<typeof lazyDoc>, dpr: number, lazy: boolean) {
  const ctx = createSeqViewer(DefaultSeqViewerSpec(), {
    scheduler: manualScheduler(),
    text: fixtureMeasurer(),
    ...(lazy ? { loaders: lazyLoaders() } : {}),
  });
  const done = ctx.setDocument(doc);
  for (const ref of ctx.viewports.sections.value) ctx.viewports.get(ref)?.setSize(W, H, dpr);
  const r = await done;
  if (!r.ok) throw new Error(JSON.stringify(r.issues));
  await ctx.graph.update();
  return ctx;
}
function resources(ctx: SeqViewerContext): DrawResources {
  const tokens = new Map<number, string>();
  for (const [name, id] of Object.entries(SkinToken))
    tokens.set(id, (DefaultSkin as Skin)[name as keyof Skin] as string);
  return {
    font: (id) => ctx.text.fontString(id),
    fontPx: (id) => Number.parseFloat(/([0-9.]+)px/.exec(ctx.text.fontString(id))?.[1] ?? '10'),
    token: (t) => tokens.get(t) ?? '#ff00ff',
    glyph: () => null,
  };
}
/** A composer drawing into a visible canvas through a back target (the main pass's, as the surface's). */
function surface(ctx: SeqViewerContext, section: string, w: number, h: number, atlas: boolean) {
  const c = new OffscreenCanvas(w, h);
  const g = c.getContext('2d', { willReadFrequently: true }) as OffscreenCanvasRenderingContext2D;
  const be = new Canvas2DBackend([g], new Canvas2DShared(resources(ctx), { atlas }));
  const composer = new SectionComposer(ctx, section);
  return { composer, frame: () => composer.frame(be, Dirty.tracks), px: () => g.getImageData(0, 0, w, h).data };
}
function random(seed: number) {
  let s = seed >>> 0;
  return () => {
    s = (s * 1664525 + 1013904223) >>> 0;
    return s / 2 ** 32;
  };
}
/** The stream manager's passes and the tile server's answers, until nothing is in flight; `each` after each pass. */
async function arrive(ctx: SeqViewerContext, each?: (k: number) => void) {
  const st = ctx.streams as Streams;
  for (let k = 0; k < 8; k++) {
    st.pass();
    await new Promise((ok) => setTimeout(ok, 0));
    each?.(k);
    if (st.stats.current.inFlight === 0 && k > 1) break;
  }
}
/** Device px of `a` and `b` (RGBA) that differ, and the first one's (x, y). */
function differ(a: Uint8ClampedArray, b: Uint8ClampedArray, w: number) {
  let n = 0;
  let at = -1;
  for (let i = 0; i < a.length; i++)
    if (a[i] !== b[i]) {
      n++;
      if (at < 0) at = i >> 2;
    }
  return { n, x: at % w, y: Math.floor(at / w) };
}

interface Run {
  lazy: boolean;
  dpr: number;
  cell: number;
  rowScale: number;
  atlas: boolean;
  seed: number;
  count: number;
}
/** Random pans and scrolls, each frame compared with a full redraw: the copies made, and those while pending. */
async function moves(o: Run): Promise<number> {
  const ctx = await open(o.lazy ? lazyDoc() : memoryDoc(), o.dpr, o.lazy);
  const section = ctx.viewports.sections.value[0];
  const vp = ctx.viewports.get(section);
  if (!vp) throw new Error('no viewport');
  // From the middle, so no move is clamped at an edge.
  vp.setView(o.cell, 0, 'api');
  vp.setRowScale(o.rowScale, 0);
  vp.panBy(1500 * o.cell, 600, 'api');
  await ctx.graph.update();
  const copied = surface(ctx, section, W * o.dpr, H * o.dpr, o.atlas);
  const full = surface(ctx, section, W * o.dpr, H * o.dpr, o.atlas);
  copied.frame();
  await arrive(ctx);
  copied.composer.redraw();
  copied.frame();
  const rnd = random(o.seed);
  let pendingCopies = 0;
  const compare = (where: string) => {
    full.composer.redraw();
    full.frame();
    const d = differ(copied.px(), full.px(), W * o.dpr);
    expect(d.n, `${where}, first at (${d.x}, ${d.y})`).toBe(0);
  };
  for (let k = 0; k < o.count; k++) {
    const far = rnd() < 0.1;
    const dx = rnd() < 0.7 ? (rnd() - 0.5) * (far ? 2000 : 160) : 0;
    const dy = dx === 0 || rnd() < 0.4 ? (rnd() - 0.5) * (far ? 800 : 120) : 0;
    vp.panBy(dx, dy, 'api');
    if (rnd() < 0.25) {
      // A frame while what the move uncovered is pending (the band drawn whole), then a frame per arrival: each redraws
      // the regions that changed (their columns; a labels chunk, the band whole) and equals a full redraw.
      const before = copied.composer.stats.copies;
      copied.frame();
      pendingCopies += copied.composer.stats.copies - before;
      compare(`move ${k}, pending`);
      await arrive(ctx, (p) => {
        copied.frame();
        compare(`move ${k}, arrival ${p}`);
      });
      continue;
    }
    await arrive(ctx);
    copied.frame();
    compare(`move ${k} (${dx.toFixed(1)}, ${dy.toFixed(1)})`);
  }
  const n = copied.composer.stats.copies;
  ctx.dispose();
  // Some moves' frames copied while what they uncovered was pending (a far pan's does not fit).
  expect(pendingCopies).toBeGreaterThan(0);
  return n;
}

describe('a rows band copies on pan once its regions arrived, equal to a full redraw', () => {
  it.each([1, 2])('lazy lettered rows (18 px, cell 14) at DPR %i, letters from the atlas', async (dpr) => {
    const n = await moves({ lazy: true, dpr, cell: 14, rowScale: 1, atlas: true, seed: dpr, count: 40 });
    expect(n).toBeGreaterThan(20);
  });
  it('lazy lettered rows at DPR 2, letters by fillText', async () => {
    const n = await moves({ lazy: true, dpr: 2, cell: 14, rowScale: 1, atlas: false, seed: 5, count: 30 });
    expect(n).toBeGreaterThan(15);
  });
  it.each([1, 2])('lazy rows of cells (cell 4) at DPR %i', async (dpr) => {
    const n = await moves({ lazy: true, dpr, cell: 4, rowScale: 1, atlas: true, seed: 10 + dpr, count: 30 });
    expect(n).toBeGreaterThan(15);
  });
  it('lazy dense rows (4 px, 1 px cells), coarse tiles', async () => {
    const n = await moves({ lazy: true, dpr: 1, cell: 1, rowScale: 4 / 18, atlas: true, seed: 13, count: 30 });
    expect(n).toBeGreaterThan(15);
  });
  it.each([1, 2])("an in-memory alignment's dense band (rows 4 px) at DPR %i", async (dpr) => {
    const n = await moves({ lazy: false, dpr, cell: 3, rowScale: 4 / 18, atlas: true, seed: 20 + dpr, count: 30 });
    expect(n).toBeGreaterThan(15);
  });
  // A sample's width or a dense row's height off whole device px: each tile at its whole scale, so a copy equals a redraw.
  it.each([
    [1, 4.3, 1],
    [2, 15.12, 1],
    [1, 15.12, 1],
    [1, 3, 0.2],
    [2, 2.7, 0.15],
    [2, 13.7, 1.1],
  ])('at DPR %i, cell %f, row scale %f the band copies too', async (dpr, cell, rowScale) => {
    const n = await moves({ lazy: true, dpr, cell, rowScale, atlas: true, seed: 30 + dpr, count: 30 });
    expect(n).toBeGreaterThan(15);
  });
  // Reduced fine tiles (translucent by their share) over the coarse level: a tile that was pending leaves nothing
  // under its neighbors once it arrives.
  it.each([
    [1, 0.5, 0.5 / 18],
    [2, 0.3, 0.6 / 18],
  ])('reduced tiles at DPR %i, cell %f, row scale %f', async (dpr, cell, rowScale) => {
    const n = await moves({ lazy: true, dpr, cell, rowScale, atlas: true, seed: 60 + dpr, count: 30 });
    expect(n).toBeGreaterThan(15);
  });
  it('an in-memory dense band at a fractional zoom (cell 2.6, rows 3.6 px) copies', async () => {
    const n = await moves({ lazy: false, dpr: 1, cell: 2.6, rowScale: 0.2, atlas: true, seed: 40, count: 30 });
    expect(n).toBeGreaterThan(15);
  });
});

describe("a tile's pixels do not depend on where the view starts", () => {
  // Two full draws at different scrolls (whole device px apart) agree in their overlap, at fractional zooms: a tile is
  // drawn at its whole scale and clipped, never resampled from the first column or row in view.
  it.each([
    [1, 4.3, 1, true],
    [2, 15.12, 1, true],
    [1, 15.12, 1, true],
    [1, 3, 0.2, true],
    [2, 2.6, 0.2, false],
  ])('at DPR %i, cell %f, row scale %f (lazy %s)', async (dpr, cell, rowScale, lazy) => {
    const ctx = await open(lazy ? lazyDoc() : memoryDoc(), dpr, lazy);
    const section = ctx.viewports.sections.value[0];
    const vp = ctx.viewports.get(section);
    if (!vp) throw new Error('no viewport');
    vp.setView(cell, 0, 'api');
    vp.setRowScale(rowScale, 0);
    vp.panBy(1500 * cell, 600, 'api');
    await ctx.graph.update();
    const [w, h] = [W * dpr, H * dpr];
    const a = surface(ctx, section, w, h, true);
    await arrive(ctx);
    a.frame();
    const l = ctx.layout.get(section)?.current;
    if (!l) throw new Error('no layout');
    const [L, R, T] = [Math.round(l.dataLeft * dpr), Math.round(l.dataRight * dpr), Math.round(l.scrolledTop * dpr)];
    const s = vp.live.current;
    let [x0, y0, pa] = [s.scrollX, s.scrollY, a.px()];
    const rnd = random(50 + dpr);
    for (let k = 0; k < 12; k++) {
      // Whole device px in both directions, less than a third of the data area's.
      vp.panBy(Math.round((rnd() - 0.5) * 0.6 * (R - L)) / dpr, Math.round((rnd() - 0.5) * 0.6 * (h - T)) / dpr, 'api');
      await arrive(ctx);
      a.composer.redraw();
      a.frame();
      const pb = a.px();
      // The previous draw's device px (x, y) shows at (x + ox, y + oy) in this one.
      const ox = Math.round((x0 - s.scrollX) * dpr);
      const oy = Math.round((y0 - s.scrollY) * dpr);
      let n = 0;
      let first = '';
      for (let y = Math.max(T, T - oy); y < Math.min(h, h - oy); y++)
        for (let x = Math.max(L, L - ox); x < Math.min(R, R - ox); x++) {
          const [i, j] = [4 * (y * w + x), 4 * ((y + oy) * w + x + ox)];
          if (pa[i] === pb[j] && pa[i + 1] === pb[j + 1] && pa[i + 2] === pb[j + 2] && pa[i + 3] === pb[j + 3])
            continue;
          if (n++ === 0) first = `(${x}, ${y})`;
        }
      expect(n, `pan ${k} by (${ox}, ${oy}), first at ${first}`).toBe(0);
      [x0, y0, pa] = [s.scrollX, s.scrollY, pb];
    }
    ctx.dispose();
  });
});

describe('a frame redrawn for arrivals alone equals a full redraw', () => {
  /**
   * The lazy alignment at a zoom, every region in view arrived (`hold`: the stream wrapped first, before anything is
   * drawn); a surface drawing it, and a full redraw to compare with.
   */
  async function at(dpr: number, cell: number, rowScale: number, hold: boolean) {
    const ctx = await open(lazyDoc(), dpr, true);
    const section = ctx.viewports.sections.value[0];
    const vp = ctx.viewports.get(section);
    if (!vp) throw new Error('no viewport');
    vp.setView(cell, 0, 'api');
    vp.setRowScale(rowScale, 0);
    vp.panBy(1500 * cell, 600, 'api');
    await ctx.graph.update();
    await arrive(ctx);
    const held = hold ? holdBack(ctx, section) : null;
    const drawn = surface(ctx, section, W * dpr, H * dpr, true);
    const full = surface(ctx, section, W * dpr, H * dpr, true);
    const compare = (where: string) => {
      full.composer.redraw();
      full.frame();
      const d = differ(drawn.px(), full.px(), W * dpr);
      expect(d.n, `${where}, first at (${d.x}, ${d.y})`).toBe(0);
    };
    return { ctx, section, vp, drawn, compare, held };
  }

  // Fine tiles held back over three columns draw the coarse level in their place alone; released one column a frame,
  // each frame redraws that column's, and nothing of the coarse level stays under the others (reduced tiles,
  // translucent by their share).
  it.each([
    [1, 0.5, 0.5 / 18],
    [2, 0.3, 0.6 / 18],
  ])('fine tiles held back, then arriving, at DPR %i, cell %f, row scale %f', async (dpr, cell, rowScale) => {
    const t = await at(dpr, cell, rowScale, true);
    const h = t.held as NonNullable<typeof t.held>;
    const lv = Array.from(h.levels);
    expect(lv[2] !== lv[0] || lv[3] !== lv[1], `levels ${lv}`).toBe(true);
    const c = h.center();
    for (const tx of [c - 1, c, c + 1]) h.hold(tx);
    t.drawn.frame();
    t.compare('held');
    for (const tx of [c, c + 1, c - 1]) {
      h.release(tx);
      t.drawn.frame();
      t.compare(`column ${tx} arrived`);
    }
    t.ctx.dispose();
  });

  // A jump to rows whose labels and tiles have not arrived: each pass's arrivals redraw what changed, the row labels
  // with a labels chunk.
  it.each([1, 2])("a labels chunk's arrival redraws the row labels, at DPR %i", async (dpr) => {
    const t = await at(dpr, 14, 1, false);
    t.drawn.frame();
    for (const row of [1300, 2600]) {
      t.vp.scrollTo(t.vp.live.current.scrollX, row * 18, 'api');
      t.drawn.frame();
      t.compare(`rows from ${row}, pending`);
      let passes = 0;
      await arrive(t.ctx, (p) => {
        passes++;
        t.drawn.frame();
        t.compare(`rows from ${row}, pass ${p}`);
      });
      expect(passes).toBeGreaterThan(0);
    }
    t.ctx.dispose();
  });
});
