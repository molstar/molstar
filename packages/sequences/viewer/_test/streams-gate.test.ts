/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST4b: the zoom gate in the request path (seqspec:sources#streamed-sources; seqspec:sources#the-viewers-side), and
// what ST3b's questions asked of it: an items data request answered with a level's bins is held as that level's
// density, and where a device px holds at most 64 positions its part in view is asked again in regions of 4,096, whose
// items are drawn; a streamed item is hovered as its range (lit across the section) and a click selects that range; a
// click on an error band retries;
// glyph themes over a field reach streamed items; `Stream.ShimmerFps` paces the shimmer (0: a still band); a failed
// region with a coarser level drawn in its place is marked, not banded.

import type { DataPatch, StreamRequest } from '@molstar/seq-core/format';
import { describe, expect, it } from 'vitest';
import {
  GRCH38,
  hg38Streams,
  itemsManifest,
  sampledItems,
  scripted,
  streamedDoc,
  VARIANT_LEVELS,
  VOCABULARY,
} from '../../core/_test/format/stream-provider.js';
import { SkinToken } from '../src/draw/backend.js';
import { Dirty } from '../src/draw/composer.js';
import { KernelRenderers, type LayerRenderer, SeqViewerConfig } from '../src/index.js';
import type { Streams } from '../src/streams.js';
import { composed, recording } from './composer-fixtures.js';
import { driven } from './input-fixtures.js';
import { GENOME } from './stream-fixtures.js';
import { progressGenome, rectsOf } from './stream-progress-fixtures.js';

const flush = () => new Promise((ok) => setTimeout(ok, 0));
const THREE = GRCH38.slice(0, 3);
const KEYS = THREE.map((c) => c[0]);
const AT = 100_000_000;

/** Variants whose provider gates every data request over 4,096 positions: level 0's `count` bins instead. */
function gatedVariants() {
  return scripted(
    itemsManifest(THREE, VARIANT_LEVELS, VOCABULARY.variants),
    (r: StreamRequest): DataPatch => {
      if (r.level >= 0 || r.force || r.stop - r.start <= 4096) return sampledItems(r, KEYS, 'variants', VARIANT_LEVELS);
      const binned = sampledItems({ ...r, level: 0 }, KEYS, 'variants', VARIANT_LEVELS);
      return { ...binned, level: 0 };
    },
    { gated: true },
  );
}

async function gatedGenome(start: number, stop: number) {
  const s = hg38Streams({ chromosomes: THREE, letters: true, gated: true });
  const variants = gatedVariants();
  const loaders = {
    load: () => Promise.reject(new Error('none')),
    stream: (src: { stream: string | null }) => (src.stream === 'variants' ? variants.provider : s.stream(src)),
  };
  const t = await composed(streamedDoc({ items: true, values: false }), { width: 1000, height: 300, dpr: 2, loaders });
  t.ctx.config.set(SeqViewerConfig.Stream.MaxInFlight, 64);
  t.vp?.zoomTo(1, start, stop);
  t.layout(GENOME);
  const answer = async () => {
    for (let k = 0; k < 3; k++) {
      (t.ctx.streams as Streams).pass();
      for (const p of [s.genome, s.genes, variants]) p.release();
      await flush();
    }
  };
  return { t, variants, answer };
}

/** The variants layer's ref. */
function variantsRef(t: Awaited<ReturnType<typeof composed>>): string {
  const m = t.ctx.graph.assembled.value?.view;
  let ref = '';
  for (let k = 0; m && k < m.resolved.layers.count; k++)
    if (m.set.strings.get(m.resolved.layers.type[k]) === 'variants') ref = m.set.strings.get(m.resolved.layers.ref[k]);
  return ref;
}

describe('the zoom gate', () => {
  it('a data request answered with bins is density; nothing finer is asked while a device px holds over 64', async () => {
    const { t, variants, answer } = await gatedGenome(AT, AT + 400_000);
    await answer();
    expect(variants.requests.every((r) => r.level >= 0 || r.stop - r.start > 4096)).toBe(true);
    const be = t.frame(recording(), Dirty.tracks);
    expect(be.glyphItems.n).toBe(0);
    expect(be.rectItems.n).toBeGreaterThan(50);
    t.ctx.dispose();
  });

  it('a gated region says to zoom in or click; forced, its part in view is asked past the gate and its items drawn', async () => {
    const { t, variants, answer } = await gatedGenome(AT, AT + 400_000);
    await answer();
    const texts = (be: ReturnType<typeof recording>) => be.textItems.text.slice(0, be.textItems.n);
    const before = t.frame(recording(), Dirty.tracks);
    expect(texts(before)).toContain('zoom in to see variants, or click to load them');
    const ref = variantsRef(t);
    expect(t.ctx.streams.force(ref, 1, AT + 50_000)).toBe(true);
    expect(t.ctx.streams.force(ref, 1, AT + 50_000)).toBe(false);
    await answer();
    // The region of 65,536 positions is asked in regions of 4,096, each forced; none is asked whole.
    const forced = variants.requests.filter((r) => r.force === true);
    expect(forced.length).toBe(16);
    expect(forced.every((r) => r.level < 0 && r.stop - r.start === 4096)).toBe(true);
    for (const r of forced) expect(Math.floor(r.start / 65_536)).toBe(Math.floor((AT + 50_000) / 65_536));
    expect(variants.requests.every((r) => r.level >= 0 || r.stop - r.start <= 4096 || r.force !== true)).toBe(true);
    // Every request for data carries the region bound (`Stream.RegionBytes`).
    expect(variants.requests.filter((r) => r.level < 0).every((r) => r.max_bytes === 1 << 20)).toBe(true);
    // Its items are drawn in the density's place (the other regions in view stay gated, hint and all).
    const be = t.frame(recording(), Dirty.tracks);
    expect(be.rectItems.n).not.toBe(before.rectItems.n);
    expect(texts(be)).toContain('zoom in to see variants, or click to load them');
    t.ctx.dispose();
  });

  it("a click on a region cut by the view's edge loads the part in view, clipped to the region", async () => {
    const { t, variants, answer } = await gatedGenome(AT, AT + 400_000);
    await answer();
    // Region 1525 ends at 100,007,936: of it only [AT, 100,007,936) is in view (two regions of 4,096).
    expect(t.ctx.streams.force(variantsRef(t), 1, AT + 10)).toBe(true);
    await answer();
    const forced = variants.requests.filter((r) => r.force === true);
    expect(forced.length).toBeGreaterThan(0);
    expect(forced.length).toBeLessThanOrEqual(2);
    for (const r of forced) {
      expect(r.stop - r.start).toBe(4096);
      expect(r.stop).toBeGreaterThan(AT);
      expect(r.stop).toBeLessThanOrEqual(100_007_936);
    }
    t.ctx.dispose();
  });

  it("zoomed in, a gated region's part in view is asked in regions of 4,096, and its items drawn", async () => {
    const { t, variants, answer } = await gatedGenome(AT + 40_000, AT + 45_000);
    await answer();
    const fine = variants.requests.filter((r) => r.level < 0 && r.stop - r.start === 4096);
    expect(fine.length).toBeGreaterThan(0);
    expect(fine.length).toBeLessThan(6);
    // Every finer region lies inside the view's gated data region.
    for (const r of fine) expect(Math.floor(r.start / 65_536)).toBe(Math.floor((AT + 40_000) / 65_536));
    const be = t.frame(recording(), Dirty.tracks);
    expect(be.glyphItems.n).toBeGreaterThan(0);
    t.ctx.dispose();
  });
});

describe('streamed items hovered, selected and retried', () => {
  async function browser(gated: boolean) {
    const s = hg38Streams({ chromosomes: THREE, letters: true, gated: true });
    const failing = { on: false };
    const t = await driven(streamedDoc({ items: true, values: false }), {
      renderers: [...KernelRenderers] as LayerRenderer[],
      loaders: {
        load: () => Promise.reject(new Error('none')),
        stream: (src) => {
          const p = s.stream(src);
          return (
            p && { fetch: (r, sig) => (failing.on && r.system ? Promise.reject(new Error('503')) : p.fetch(r, sig)) }
          );
        },
      },
      size: { [GENOME]: [1000, 300] },
    });
    t.ctx.config.set(SeqViewerConfig.Stream.MaxInFlight, 64);
    t.vp(GENOME).zoomTo(1, AT, AT + 100_000);
    await t.ctx.graph.update();
    failing.on = gated;
    for (let k = 0; k < 3; k++) {
      (t.ctx.streams as Streams).pass();
      for (const p of [s.genome, s.genes, s.variants]) p.release();
      await flush();
    }
    return { t, s, failing };
  }
  /** Down the data area at x to the Genes track; the hover. */
  const onGenes = (t: Awaited<ReturnType<typeof browser>>['t'], x: number, item = false) => {
    const l = t.layout(GENOME);
    const m = t.ctx.graph.assembled.value?.view;
    const h = t.ctx.hover.live.current;
    let y = 1;
    for (; y < 300; y++) {
      t.pointer('move', GENOME, l.dataLeft + x, y);
      const genes = h.track >= 0 && m?.set.strings.get(m.resolved.tracks.label[h.track]) === 'Genes';
      if (genes && (!item || h.part >= 0)) break;
    }
    return { h, l, y };
  };

  it("a gene's hover is its range on the unit, and a click selects that range", async () => {
    const { t } = await browser(false);
    let found = false;
    for (let x = 10; x < 900 && !found; x += 7) {
      const { h, l, y } = onGenes(t, x, true);
      if (h.part < 0) continue;
      found = true;
      expect(h.span[1] - h.span[0]).toBeGreaterThan(1000);
      expect(h.index >= h.span[0] && h.index < h.span[1]).toBe(true);
      const [a, b] = [h.span[0], h.span[1]];
      t.pointer('down', GENOME, l.dataLeft + x, y);
      t.pointer('up', GENOME, l.dataLeft + x, y);
      // Its double click zooms: the click's selection waits for the double-click time.
      t.advance(500);
      const sel = t.ctx.selection.current.value.seq as { elements: { indices: { min?: number } }[] } | null;
      expect(sel).toBeTruthy();
      const { OrderedSet } = await import('@molstar/seq-core/model');
      const ix = sel?.elements[0].indices as never;
      expect([OrderedSet.min(ix), OrderedSet.max(ix) + 1]).toEqual([a, b]);
    }
    expect(found).toBe(true);
    t.ctx.dispose();
  });

  it("a click on a failed region's band retries it, and selects nothing", async () => {
    const { t, s, failing } = await browser(true);
    const { l, y } = onGenes(t, 300);
    const m = t.ctx.graph.assembled.value?.view;
    let ref = '';
    for (let k = 0; m && k < m.resolved.layers.count; k++)
      if (m.set.strings.get(m.resolved.layers.type[k]) === 'ranges') ref = m.set.strings.get(m.resolved.layers.ref[k]);
    expect(t.ctx.streams.hasFailed(ref)).toBe(true);
    failing.on = false;
    const asked = s.genes.requests.length;
    t.pointer('down', GENOME, l.dataLeft + 300, y);
    t.pointer('up', GENOME, l.dataLeft + 300, y);
    expect(t.ctx.selection.current.value.seq ?? null).toBe(null);
    (t.ctx.streams as Streams).pass();
    expect(s.genes.requests.length).toBeGreaterThan(asked);
    t.ctx.dispose();
  });
});

describe("glyph themes, the shimmer's pace and covered failures", () => {
  it('a field_map over `alt` gives streamed variants their glyphs', async () => {
    const doc = streamedDoc({ items: true, values: false });
    const walk = (n: { kind: string; params?: Record<string, unknown>; children?: unknown[] }) => {
      if (n.kind === 'layer' && n.params?.type === 'variants')
        n.params.glyph = {
          theme: 'field_map',
          params: { field: 'alt', map: { A: 'diamond', C: 'square', G: 'star', T: 'cross' }, default: 'circle' },
        };
      for (const c of n.children ?? []) walk(c as never);
    };
    walk(doc.root as never);
    const s = hg38Streams({ chromosomes: THREE, letters: true });
    const t = await composed(doc, {
      width: 1000,
      height: 300,
      dpr: 2,
      loaders: { load: () => Promise.reject(new Error('none')), stream: s.stream },
    });
    t.vp?.zoomTo(1, AT + 40_000, AT + 45_000);
    t.layout(GENOME);
    for (let k = 0; k < 3; k++) {
      (t.ctx.streams as Streams).pass();
      await flush();
    }
    const be = t.frame(recording(), Dirty.tracks);
    const ids = new Set(Array.from(be.glyphItems.glyph.subarray(0, be.glyphItems.n)));
    expect(be.glyphItems.n).toBeGreaterThan(2);
    expect(ids.size).toBeGreaterThan(1);
    expect(ids.has(0)).toBe(false);
    t.ctx.dispose();
  });

  it('`Stream.ShimmerFps` paces the shimmer; 0 leaves a still band and asks for no frame', async () => {
    const g = await progressGenome(AT, AT + 100_000, { items: true, live: true });
    g.ctx.config.set(SeqViewerConfig.Stream.MaxInFlight, 64);
    g.m.pass();
    g.clock[0] = 1300;
    g.frame(Dirty.tracks);
    g.frame(Dirty.overlay);
    expect(g.composer.animating).toBe(true);
    expect(g.composer.animStep[0]).toBeCloseTo(1000 / 30);
    g.ctx.config.set(SeqViewerConfig.Stream.ShimmerFps, 0);
    const be = g.frame(Dirty.overlay);
    expect(rectsOf(be, SkinToken.pending).n).toBeGreaterThan(2);
    expect(rectsOf(be, SkinToken.surface).alphas.filter((a) => a < 255).length).toBe(0);
    expect(g.composer.animating).toBe(false);
    g.ctx.dispose();
  });

  it("a failed values region with its coarse level drawn is marked along the band's foot, not banded", async () => {
    const g = await progressGenome(AT, AT + 100_000, { values: true, live: true });
    g.ctx.config.set(SeqViewerConfig.Stream.MaxInFlight, 64);
    // The coarsest level answered; every finer request fails.
    g.failing.when = (r) => r.level !== 5;
    g.m.pass();
    g.s.phylop.release();
    await flush();
    expect(g.ctx.streams.stats.current.errors).toBeGreaterThan(0);
    g.frame(Dirty.tracks);
    const be = g.frame(Dirty.overlay);
    const marks = rectsOf(be, SkinToken.conflict);
    expect(marks.n).toBeGreaterThan(0);
    expect(be.pathItems.n).toBe(0);
    expect(marks.alphas.every((a) => a === 255)).toBe(true);
    g.ctx.dispose();
  });
});

describe("navigation over a genome of blocks (David's JBrowse 2 report, ST4b)", () => {
  it('a locus, the whole genome, a whole chromosome and the locus again: the locus draws as it did', async () => {
    const s = hg38Streams({ chromosomes: THREE, letters: true });
    const variants = gatedVariants();
    variants.gated = false;
    const loaders = {
      load: () => Promise.reject(new Error('none')),
      stream: (src: { stream: string | null }) => (src.stream === 'variants' ? variants.provider : s.stream(src)),
    };
    const t = await composed(streamedDoc({ items: true, sequence: true }), {
      width: 800,
      height: 400,
      dpr: 1,
      loaders,
    });
    const m = t.ctx.streams as Streams;
    const settle = async () => {
      t.layout(GENOME);
      for (let k = 0; k < 12; k++) {
        m.pass();
        t.scheduler.advance(16);
        await flush();
      }
    };
    const draw = () => {
      const be = t.frame(recording(), Dirty.tracks);
      return {
        rects: be.rectItems.n,
        glyphs: be.glyphItems.n,
        texts: be.textItems.text.slice(0, be.textItems.n).sort(),
      };
    };
    const at = async (go: () => void) => {
      go();
      await settle();
      return draw();
    };
    const locus = () => t.vp?.zoomTo(1, AT, AT + 26_000);
    const first = await at(locus);
    expect(first.rects).toBeGreaterThan(50);
    const genome = await at(() => t.vp?.fit());
    expect(genome.texts).toContain('zoom in to see letters');
    await at(() => t.vp?.zoomTo(1, 0, THREE[1][1]));
    const again = await at(locus);
    expect(again).toEqual(first);
    expect(m.stats.current.inFlight).toBe(0);
    t.ctx.dispose();
  });
});
