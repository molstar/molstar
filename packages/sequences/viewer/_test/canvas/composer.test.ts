/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Copies on pan as pixels in Chromium (CS1a; seqspec:viewer#the-render-backend): after random moves on both axes, the
// composer's copied frame on the Canvas2D backend (a back target, presented whole) equals a full redraw, pixel for
// pixel, with RN1's renderers at each level of detail and at DPR 1 and 2; stacked scrolls copy; a wrapped section with
// its gutter's line numbers (`none`: drawn once per frame) and a wide gutter (an msa's 40 px) equal a full redraw too,
// and so do range labels sticking at the view's edges (seqspec:sources#levels-and-the-zoom-gate). The overlay as pixels
// (CS1b): a stacked block's selection stays in its own lane.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { describe, expect, it } from 'vitest';
import { Canvas2DBackend, Canvas2DShared } from '../../src/canvas/index.js';
import { type DrawResources, SkinToken } from '../../src/draw/backend.js';
import { Dirty, SectionComposer } from '../../src/draw/index.js';
import {
  createSeqViewer,
  DefaultSeqViewerSpec,
  DefaultSkin,
  manualScheduler,
  type RenderProps,
  type SeqViewerContext,
  type Skin,
} from '../../src/index.js';
import { fixtureMeasurer } from '../../src/testing/index.js';
import { seqLoci, twoBlocks } from '../overlay-docs.js';
import { allTypes } from '../prepare-fixtures.js';

async function open(doc: ReturnType<typeof allTypes>, w: number, h: number, dpr: number, render = {}) {
  const ctx = createSeqViewer(
    { ...DefaultSeqViewerSpec(), render },
    { scheduler: manualScheduler(), text: fixtureMeasurer() },
  );
  const done = ctx.setDocument(doc);
  for (const ref of ctx.viewports.sections.value) ctx.viewports.get(ref)?.setSize(w, h, dpr);
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
function surface(ctx: SeqViewerContext, section: string, w: number, h: number) {
  const c = new OffscreenCanvas(w, h);
  const g = c.getContext('2d', { willReadFrequently: true }) as OffscreenCanvasRenderingContext2D;
  const be = new Canvas2DBackend([g], new Canvas2DShared(resources(ctx)));
  const composer = new SectionComposer(ctx, section);
  return { composer, frame: () => composer.frame(be, Dirty.tracks), px: () => g.getImageData(0, 0, w, h).data };
}
/** A seeded generator, so a failure reproduces. */
function random(seed: number) {
  let s = seed >>> 0;
  return () => {
    s = (s * 1664525 + 1013904223) >>> 0;
    return s / 2 ** 32;
  };
}
interface Run {
  doc: ReturnType<typeof allTypes>;
  dpr: number;
  cells: number[];
  count: number;
  seed: number;
  /** Vertical moves only (wrapped). */
  vertical?: boolean;
  render?: Partial<RenderProps>;
}
/** Random pans and scrolls at each cell size, each compared with a full redraw; the copies made. */
async function moves(o: Run): Promise<number> {
  const [W, H] = [700, 260];
  const ctx = await open(o.doc, W, H, o.dpr, o.render);
  const section = ctx.viewports.sections.value[0];
  const vp = ctx.viewports.get(section);
  if (!vp) throw new Error('no viewport');
  const copied = surface(ctx, section, W * o.dpr, H * o.dpr);
  const full = surface(ctx, section, W * o.dpr, H * o.dpr);
  const rnd = random(o.seed);
  for (const cell of o.cells) {
    if (!o.vertical) vp.setView(cell, 0, 'api');
    await ctx.graph.update();
    copied.frame();
    for (let k = 0; k < o.count; k++) {
      const far = rnd() < 0.1;
      const dx = o.vertical ? 0 : (rnd() - 0.5) * (far ? 2000 : 160);
      const dy = o.vertical || rnd() < 0.6 ? (rnd() - 0.5) * (far ? 800 : 120) : 0;
      vp.panBy(dx, dy, 'api');
      copied.frame();
      full.composer.redraw();
      full.frame();
      const [a, b] = [copied.px(), full.px()];
      let diff = 0;
      let at = -1;
      for (let i = 0; i < a.length; i++)
        if (a[i] !== b[i]) {
          diff++;
          if (at < 0) at = i >> 2;
        }
      const where = `cell ${cell}, move ${k} (${dx.toFixed(1)}, ${dy.toFixed(1)}), first at (${at % (W * o.dpr)}, ${Math.floor(at / (W * o.dpr))})`;
      expect(diff, where).toBe(0);
    }
  }
  const n = copied.composer.stats.copies;
  ctx.dispose();
  return n;
}

/** Letters, emphasis, values and color over a long protein, wrapped: its lines numbered in the gutter. */
function wrapped() {
  const b = createSeqViewSpec();
  const a = b.data.sequence({ ref: 'a', value: 'MKVLATGHEKWRSPQIDNCY'.repeat(40) });
  const d = b.data.annotation({
    ref: 'd',
    geometry: 'ranges',
    target: a,
    value: { start: [12, 300], stop: [70, 420] },
  });
  const v = b.view({ ref: 'v', label: 'V', layout: 'text', height: 'box' });
  v.block({ sequence: a });
  const t = v.track({ label: 'Letters' });
  t.letters();
  t.emphasis({ data: d });
  return b.getDocument();
}

/** A protein filtered to 30 sites: a break hairline between each two. */
function filtered() {
  const b = createSeqViewSpec();
  const a = b.data.sequence({ ref: 'a', value: 'MKVLATGHEKWRSPQIDNCY'.repeat(40) });
  const start = Array.from({ length: 30 }, (_, k) => k * 26);
  const site = b.data.annotation({
    ref: 'site',
    geometry: 'ranges',
    target: a,
    value: { start, stop: start.map((x) => x + 12) },
  });
  const v = b.view({ ref: 'v', label: 'V', positions: site });
  v.block({ sequence: a });
  v.track({ label: 'Letters' }).letters();
  v.track({ label: 'Values' }).values({
    data: b.data.annotation({
      ref: 'vals',
      geometry: 'values',
      target: a,
      value: { value: start.map((_, k) => k / 30) },
    }),
  });
  return b.getDocument();
}

/** Long labeled ranges over a protein, in two lanes: each label cut by the view's edges as it pans, sticking there. */
function labeled() {
  const b = createSeqViewSpec();
  const a = b.data.sequence({ ref: 'a', value: 'MKVLATGHEKWRSPQIDNCY'.repeat(30) });
  const r = b.data.annotation({
    ref: 'r',
    geometry: 'ranges',
    target: a,
    value: {
      start: [20, 250, 300, 330, 420],
      stop: [260, 290, 380, 336, 590],
      label: ['Kinase domain', 'Linker', 'SH2 domain', 'Site', 'C-terminal tail'],
    },
  });
  const v = b.view({ ref: 'v', label: 'V' });
  v.block({ sequence: a });
  v.track({ label: 'Letters' }).letters();
  v.track({ label: 'Features' }).ranges({ data: r, labels: true });
  v.track({ label: 'Brackets' }).ranges({ data: r, labels: true, display: 'brackets' } as never);
  return b.getDocument();
}

describe('copies equal a full redraw, pixel for pixel', () => {
  it.each([1, 2])("range labels sticking at the view's edges pan at DPR %i", async (dpr) => {
    const n = await moves({ doc: labeled(), dpr, cells: [14, 6, 2.5], count: 40, seed: 13 + dpr });
    expect(n).toBeGreaterThan(80);
  });

  it.each([1, 2])('every RN1 type side by side at DPR %i; unstable layers redrawn whole', async (dpr) => {
    const n = await moves({ doc: allTypes({ length: 600 }), dpr, cells: [14, 4, 1.5, 0.6], count: 30, seed: dpr });
    expect(n).toBeGreaterThan(60);
  });
  it('stacked scrolls and pans copy', async () => {
    const doc = structuredClone(allTypes({ length: 600 }));
    type Node = { kind: string; params?: Record<string, unknown>; children?: Node[] };
    const walk = (n: Node): void => {
      if (n.kind === 'section') n.params = { ...n.params, arrangement: 'stacked' };
      n.children?.forEach(walk);
    };
    walk(doc.root as Node);
    const n = await moves({ doc, dpr: 2, cells: [12, 3], count: 30, seed: 3 });
    expect(n).toBeGreaterThan(30);
  });
  it.each(['none', 'aboveEnd'] as const)(
    "wrapped scrolls with %s group labels (the gutter's numbers once)",
    async (labels) => {
      const n = await moves({
        doc: wrapped(),
        dpr: 2,
        cells: [0],
        count: 40,
        seed: 5,
        vertical: true,
        render: { wrappedGroupLabels: labels },
      });
      expect(n).toBeGreaterThan(20);
    },
  );
  it("a filtered section's break hairlines", async () => {
    const n = await moves({ doc: filtered(), dpr: 2, cells: [14, 3], count: 30, seed: 11 });
    expect(n).toBeGreaterThan(20);
  });
  it("a wide gutter (an msa's 40 px) pans", async () => {
    const n = await moves({
      doc: allTypes({ length: 600 }),
      dpr: 1,
      cells: [10, 2],
      count: 30,
      seed: 7,
      render: { gutterLeft: 40 },
    });
    expect(n).toBeGreaterThan(30);
  });
});

describe('the overlay as pixels', () => {
  it("stacked: a block's selection and hover stay in its own lane after scrolls, at DPR 2", async () => {
    const [W, H, dpr] = [600, 300, 2];
    const ctx = await open(twoBlocks({ stacked: true }), W, H, dpr);
    const section = ctx.viewports.sections.value[0];
    const set = ctx.graph.assembled.value?.set;
    if (!set) throw new Error('no set');
    ctx.selection.set(seqLoci(set, 'a', 10, 13), 'api');
    const canvas = () => new OffscreenCanvas(W * dpr, H * dpr).getContext('2d', { willReadFrequently: true });
    const targets = [0, 1, 2, 3].map(canvas) as OffscreenCanvasRenderingContext2D[];
    const g = targets[1];
    const be = new Canvas2DBackend(targets, new Canvas2DShared(resources(ctx)));
    const composer = new SectionComposer(ctx, section);
    const vp = ctx.viewports.get(section);
    for (const dy of [0, 7, 13]) {
      vp?.panBy(0, dy, 'api');
      composer.frame(be, Dirty.tracks | Dirty.overlay);
      const l = ctx.layout.get(section)?.current;
      const ax = ctx.layout.axis(section);
      if (!l || !ax) throw new Error('no layout');
      const px = g.getImageData(0, 0, W * dpr, H * dpr).data;
      const x = Math.round(((ax.x(0, 10 as never) + ax.x(0, 13 as never)) / 2) * dpr);
      const alpha = (y: number) => px[(Math.round(y * dpr) * W * dpr + x) * 4 + 3];
      const mid = (b: number) => {
        const [e0, e1] = [l.blockEntries[2 * b], l.blockEntries[2 * b + 1] - 1];
        const top = l.scrolledTop + l.entryTop[e0] - (vp?.live.current.scrollY ?? 0);
        return (top + l.scrolledTop + l.entryTop[e1] + l.entryHeight[e1] - (vp?.live.current.scrollY ?? 0)) / 2;
      };
      expect(alpha(mid(0)), `scrolled ${dy}`).toBeGreaterThan(0);
      if (mid(1) < l.height) expect(alpha(mid(1)), `scrolled ${dy}`).toBe(0);
    }
    ctx.dispose();
  });
});
