/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST3b: progressive loading (seqspec:sources#levels-and-the-zoom-gate; seqspec:sources#the-viewers-side;
// seqspec:viewer#labels). On the overlay, over the main pass's cells: a pending region with nothing drawn in its place
// is a skeleton (the pending shade, a shimmer while the composer is live), a failed one an error band (tinted, hatched,
// "failed to load"), one just arrived the shade fading out over 150 ms, regions landing together alike; `seqv:retry`
// asks failed regions again; values keep their coarser level under a pending one (no skeleton); the `stream` label says
// "loading…" and "failed" and names a streamed item and its fields; field themes reach streamed items through the
// manifest's vocabulary; an items layer's lanes are
// held by its regions; a miniature naming a rows track in a section without an overview column loads its overview.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { describe, expect, it } from 'vitest';
import { GENE_CATEGORIES, GRCH38, geneAt, hg38Streams, streamedDoc } from '../../core/_test/format/stream-provider.js';
import { FADE_MS } from '../src/draw/progress.js';
import { KernelRenderers, type LayerRenderer, SeqCommands } from '../src/index.js';
import type { Streams } from '../src/streams.js';
import { driven } from './input-fixtures.js';
import { lazyLoaders, opened } from './layout-fixtures.js';
import { CacheBytes, GENOME, MaxInFlight } from './stream-fixtures.js';
import { Dirty, progressGenome, rectsOf, SkinToken } from './stream-progress-fixtures.js';

const flush = () => new Promise((ok) => setTimeout(ok, 0));
const AT = 100_000_000;
/** The shimmer's strips: the surface under its full alpha. */
const sheen = (a: number[]) => a.filter((x) => x < 255).length;

describe('skeletons, fade-in and errors on the overlay', () => {
  it('a pending region is a skeleton: the shade, shimmering only while live and in motion', async () => {
    const g = await progressGenome(AT, AT + 100_000, { items: true, values: true, live: true });
    g.ctx.config.set(MaxInFlight, 64);
    g.m.pass();
    // A third into a sweep: the band over the data area.
    g.clock[0] = 1300;
    g.frame(Dirty.tracks);
    const be = g.frame(Dirty.overlay);
    expect(rectsOf(be, SkinToken.pending).n).toBeGreaterThan(2);
    const at = rectsOf(be, SkinToken.surface);
    expect(sheen(at.alphas)).toBeGreaterThan(0);
    expect(g.composer.animating).toBe(true);
    // The sweep moves with the clock (one sweep a second).
    const x0 = Array.from(be.rectItems.x.subarray(0, be.rectItems.n));
    g.clock[0] = 1600;
    const x1 = Array.from(g.frame(Dirty.overlay).rectItems.x.subarray(0, be.rectItems.n));
    expect(x1).not.toEqual(x0);
    // Reduced motion, or an export's composer: the shade alone, no frame asked.
    g.ctx.config.set((await import('../src/index.js')).SeqViewerConfig.Viewport.ReducedMotion, true);
    expect(sheen(rectsOf(g.frame(Dirty.overlay), SkinToken.surface).alphas)).toBe(0);
    expect(rectsOf(g.backend, SkinToken.pending).n).toBeGreaterThan(2);
    expect(g.composer.animating).toBe(false);
    g.ctx.dispose();
  });

  it('an arrival fades in over 150 ms, regions landing together alike; then nothing is asked of the overlay', async () => {
    const g = await progressGenome(AT, AT + 100_000, { items: true, values: true, live: true });
    g.ctx.config.set(MaxInFlight, 64);
    g.clock[0] = 1000;
    await g.release();
    g.frame(Dirty.tracks);
    // The first frame to draw them: the shade whole over each.
    expect(new Set(rectsOf(g.frame(Dirty.overlay), SkinToken.pending).alphas)).toEqual(new Set([255]));
    g.clock[0] = 1050;
    const fading = rectsOf(g.frame(Dirty.overlay), SkinToken.pending);
    expect(fading.n).toBeGreaterThan(2);
    expect(new Set(fading.alphas)).toEqual(new Set([Math.round(255 * (1 - 50 / FADE_MS))]));
    expect(g.composer.animating).toBe(true);
    g.clock[0] = 1000 + FADE_MS + 1;
    expect(rectsOf(g.frame(Dirty.overlay), SkinToken.pending).n).toBe(0);
    expect(g.composer.animating).toBe(false);
    g.ctx.dispose();
  });

  it('a failed region is an error band with its hint; `seqv:retry` asks it again', async () => {
    const g = await progressGenome(AT, AT + 100_000, { items: true, values: true, live: true });
    g.ctx.config.set(MaxInFlight, 64);
    g.failing.on = true;
    g.m.pass();
    await flush();
    const ref = (k: string) => {
      const m = g.ctx.graph.assembled.value?.view;
      for (let i = 0; m && i < m.resolved.layers.count; i++)
        if (m.set.strings.get(m.resolved.layers.type[i]) === k) return m.set.strings.get(m.resolved.layers.ref[i]);
      return '';
    };
    expect(g.m.hasFailed(ref('ranges'))).toBe(true);
    g.frame(Dirty.tracks);
    const be = g.frame(Dirty.overlay);
    const tint = rectsOf(be, SkinToken.conflict);
    expect(tint.n).toBeGreaterThan(2);
    expect(new Set(tint.alphas).size).toBe(1);
    expect(be.pathItems.n).toBeGreaterThan(10);
    expect(be.textItems.text.slice(0, be.textItems.n)).toContain('failed to load');
    // Failed is no skeleton: nothing shimmers, nothing is asked of the overlay.
    expect(sheen(rectsOf(be, SkinToken.surface).alphas)).toBe(0);
    expect(g.composer.animating).toBe(false);
    // Retry: the cooldown ends and the next pass asks them again; answered, the bands are gone.
    g.failing.on = false;
    const asked = g.asked();
    expect(await g.ctx.commands.run({ command: SeqCommands.Retry.id, params: {} }, 'api')).toBe(true);
    await g.release();
    expect(g.asked()).toBeGreaterThan(asked);
    expect(g.m.hasFailed(ref('ranges'))).toBe(false);
    g.frame(Dirty.tracks);
    expect(rectsOf(g.frame(Dirty.overlay), SkinToken.conflict).n).toBe(0);
    g.ctx.dispose();
  });

  it('values keep a coarser level under a pending one: no skeleton there (coarse under fine)', async () => {
    const g = await progressGenome(AT, AT + 100_000, { values: true });
    g.ctx.config.set(MaxInFlight, 64);
    g.m.pass();
    // The coarsest level is asked first: answered alone, it is drawn while the view's level loads.
    const top = g.s.phylop.requests.filter((r) => r.level === 5).length;
    expect(top).toBeGreaterThan(0);
    g.s.phylop.release(top);
    await flush();
    g.frame(Dirty.tracks);
    expect(g.backend.rectItems.n).toBeGreaterThan(100);
    expect(g.s.phylop.held).toBeGreaterThan(0);
    expect(rectsOf(g.frame(Dirty.overlay), SkinToken.pending).n).toBe(0);
    g.ctx.dispose();
  });

  it('lazy rows: a tile with nothing in its place is a skeleton until its tiles arrive', async () => {
    const uri = 'https://example.org/lazy';
    const b = createSeqViewSpec({ title: 'lazy' });
    const aln = b.data.alignment({
      ref: 'lazy',
      row_count: 3000,
      column_count: 400,
      source: { uri, protocol: 'seqstar-tiles/1' },
    });
    b.view({ ref: 'v', label: 'V', style: 'msa', height: 'box' })
      .block({ alignment: aln })
      .track({ label: 'Rows' })
      .rows();
    const { composed, recording } = await import('./composer-fixtures.js');
    const t = await composed(b.getDocument(), { loaders: lazyLoaders(uri, 3000, 400) });
    t.composer.live = true;
    t.vp?.setView(1.5, 0, 'api');
    await t.ctx.graph.update();
    const be = recording();
    t.frame(be, Dirty.tracks);
    expect(rectsOf(t.frame(be, Dirty.overlay), SkinToken.pending).n).toBeGreaterThan(0);
    expect(t.composer.animating).toBe(true);
    for (let k = 0; k < 6; k++) {
      t.scheduler.advance(16);
      await flush();
    }
    t.scheduler.advance(FADE_MS + 16);
    t.frame(be, Dirty.tracks);
    t.frame(be, Dirty.overlay);
    t.scheduler.advance(FADE_MS + 16);
    expect(rectsOf(t.frame(be, Dirty.overlay), SkinToken.pending).n).toBe(0);
    t.ctx.dispose();
  });
});

describe('streamed items: hover, labels, field themes, lanes', () => {
  async function browser(color?: unknown) {
    const s = hg38Streams({ chromosomes: GRCH38.slice(0, 3), letters: true, gated: true });
    const doc = streamedDoc({ items: true, values: false });
    if (color) {
      const view = (doc.root.children ?? []).find((c) => c.kind === 'view');
      const genes = JSON.stringify(view).includes('Genes') ? view : null;
      const walk = (n: { kind: string; params?: Record<string, unknown>; children?: unknown[] }) => {
        if (n.kind === 'layer' && n.params?.type === 'ranges') n.params.color = color;
        for (const c of n.children ?? []) walk(c as never);
      };
      if (genes) walk(genes as never);
    }
    const t = await driven(doc, {
      renderers: [...KernelRenderers] as LayerRenderer[],
      loaders: { load: () => Promise.reject(new Error('none')), stream: s.stream },
      size: { [GENOME]: [1000, 300] },
    });
    t.ctx.config.set(MaxInFlight, 64);
    t.vp(GENOME).zoomTo(1, AT, AT + 100_000);
    await t.ctx.graph.update();
    const answer = async () => {
      for (let k = 0; k < 3; k++) {
        (t.ctx.streams as Streams).pass();
        for (const p of [s.genome, s.genes, s.variants, s.phylop]) p.release();
        await flush();
      }
    };
    return { t, s, answer };
  }
  /** Down the data area to the track labeled `name`, at x; the hover. */
  function point(t: Awaited<ReturnType<typeof browser>>['t'], name: string, x: number) {
    const l = t.layout(GENOME);
    const m = t.ctx.graph.assembled.value?.view;
    const h = t.ctx.hover.live.current;
    const label = (k: number) => m?.set.strings.get(m.resolved.tracks.label[k]);
    for (let y = 1; y < 300; y += 1) {
      t.pointer('move', GENOME, l.dataLeft + x, y);
      if (h.track >= 0 && label(h.track) === name && h.part >= 0) break;
    }
    return h;
  }

  it("the stream label: loading… over a pending region, made again when it lands; a gene's label and fields", async () => {
    const { t, answer } = await browser();
    const l = t.layout(GENOME);
    const m = t.ctx.graph.assembled.value?.view;
    const label = (k: number) => m?.set.strings.get(m.resolved.tracks.label[k]);
    const h = t.ctx.hover.live.current;
    (t.ctx.streams as Streams).pass();
    for (let y = 1; y < 300 && (h.track < 0 || label(h.track) !== 'Genes'); y += 2)
      t.pointer('move', GENOME, l.dataLeft + 300, y);
    const lines = () => (t.ctx.labels.current.value?.lines ?? []).map((x) => `${x.key}: ${x.value}`);
    expect(lines()).toContain('Genes: loading…');
    await answer();
    expect(lines()).not.toContain('Genes: loading…');
    // A gene under the pointer: its number in the layer as the hover's part, its label and category as lines.
    let found = -1;
    for (let x = 10; x < 900 && found < 0; x += 7) found = point(t, 'Genes', x).part;
    expect(found).toBeGreaterThanOrEqual(0);
    const got = lines();
    expect(got.some((x) => /^Genes: G\d+$/.test(x))).toBe(true);
    expect(got.some((x) => GENE_CATEGORIES.some((c) => x === `category: ${c}`))).toBe(true);
    t.ctx.dispose();
  });

  it('a hovered gene is outlined over its exons, in its lane', async () => {
    const { t, answer } = await browser();
    await answer();
    let part = -1;
    for (let x = 10; x < 900 && part < 0; x += 7) part = point(t, 'Genes', x).part;
    expect(part).toBeGreaterThanOrEqual(0);
    const { recording } = await import('./composer-fixtures.js');
    const { SectionComposer } = await import('../src/draw/composer.js');
    const c = new SectionComposer(t.ctx, GENOME);
    const be = recording();
    c.frame(be, Dirty.tracks);
    be.reset();
    c.frame(be, Dirty.overlay);
    // The outline: stroked rects (HoverOutline's boxes).
    let stroked = 0;
    for (let i = 0; i < be.rectItems.n; i++) if (be.rectItems.stroke[i] > 0) stroked++;
    expect(stroked).toBeGreaterThan(0);
    c.dispose();
    t.ctx.dispose();
  });

  it("a field theme reaches streamed items: genes by category in the manifest's vocabulary order", async () => {
    const palette = { kind: 'categorical', colors: ['#ff0000', '#00ff00', '#0000ff'] };
    const { t, answer } = await browser({ theme: 'field_palette', params: { field: 'category', palette } });
    await answer();
    const { recording } = await import('./composer-fixtures.js');
    const { SectionComposer } = await import('../src/draw/composer.js');
    const c = new SectionComposer(t.ctx, GENOME);
    const be = recording();
    c.frame(be, Dirty.tracks);
    const colors = new Set(Array.from(be.rectItems.color.subarray(0, be.rectItems.n)));
    // The genes in view (chr2, the fixture's), each in its category's color, the palette in the vocabulary's order.
    const RGB = [0xff0000, 0x00ff00, 0x0000ff];
    const want = new Set<number>();
    for (let k = Math.floor(AT / 37_000) - 1; k * 37_000 < AT + 100_000; k++) {
      const x = geneAt(1, k);
      if (x.stop > AT && x.start < AT + 100_000) want.add(RGB[GENE_CATEGORIES.indexOf(x.category)]);
    }
    expect(want.size).toBeGreaterThan(1);
    for (const c of RGB) expect(colors.has(c), `#${c.toString(16)}`).toBe(want.has(c));
    c.dispose();
    t.ctx.dispose();
  });

  it("an items layer's lanes are held by its regions: a tour under a small cache keeps the map bounded", async () => {
    const g = await progressGenome(AT, AT + 100_000, { items: true });
    g.ctx.config.set(MaxInFlight, 64);
    g.ctx.config.set(CacheBytes, 256 * 1024);
    const m = g.ctx.graph.assembled.value?.view;
    let ref = '';
    for (let k = 0; m && k < m.resolved.layers.count; k++)
      if (m.set.strings.get(m.resolved.layers.type[k]) === 'ranges') ref = m.set.strings.get(m.resolved.layers.ref[k]);
    const layer = g.m.get(ref) as unknown as { laneOf: Map<string, unknown>; all: { keys: string[] | null }[] };
    let most = 0;
    for (let k = 0; k < 40; k++) {
      g.vp.zoomTo(1, AT + k * 3_000_000, AT + k * 3_000_000 + 100_000);
      g.layout(GENOME);
      await g.release();
      most = Math.max(most, layer.laneOf.size);
    }
    // Every key held is a held region's; no more than the regions held carry.
    const held = new Set(layer.all.flatMap((r) => r.keys ?? []));
    expect(layer.laneOf.size).toBe(held.size);
    expect(most).toBeLessThan(400);
    g.ctx.dispose();
  });
});

describe('a miniature naming a rows track in a section without an overview column', () => {
  it("loads the overview's tiles for it (MAP_WIDTH's levels), and none when nothing names it", async () => {
    const uri = 'https://example.org/lazy';
    const doc = (named: boolean) => {
      const b = createSeqViewSpec({ title: 'lazy' });
      const aln = b.data.alignment({
        ref: 'lazy',
        row_count: 3000,
        column_count: 400,
        source: { uri, protocol: 'seqstar-tiles/1' },
      });
      const v = b.view({ ref: 'v', label: 'V', style: 'tracks', height: 'box' }).block({ alignment: aln });
      v.track({ ref: 't.rows', label: 'Rows' }).rows();
      v.track({ label: 'Map' }).layer({
        type: 'seqv:map',
        ...(named ? { options: { tracks: ['t.rows'] } } : {}),
        fallback: 'omit',
      } as never);
      return b.getDocument();
    };
    const overviewTiles = async (named: boolean) => {
      const t = await opened(doc(named), {
        renderers: [...KernelRenderers] as LayerRenderer[],
        loaders: lazyLoaders(uri, 3000, 400),
        render: { verticalOverview: 0 },
      });
      t.vp('v.main').setView(1.5, 0, 'api');
      await t.ctx.graph.update();
      expect(t.ctx.layout.axis('v.main')?.vnavWidth).toBe(0);
      (t.ctx.streams as Streams).pass();
      const n = t.ctx.streams.stats.current.requests;
      t.ctx.dispose();
      return n;
    };
    expect(await overviewTiles(true)).toBeGreaterThan(await overviewTiles(false));
  });
});
