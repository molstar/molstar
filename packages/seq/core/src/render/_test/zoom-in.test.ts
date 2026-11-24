/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A layer whose wants together pass the store's bound on a layer's tiles or read cost (_spec/README.md "Render", Data):
// its band hatched with "Zoom in to read: ⌘/Ctrl+scroll, or +", an info issue said, distinct from pending (which settles by
// itself); zoomed in, the layer reads and draws again. The drawn list is what the canvas and the SVG both draw.

import type { Coverage, FormatRuntime, OpenedSource } from '@molstar/seqx-core/data';
import { SeqRuntime } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import type { Frame } from '../index.js';
import { EmptyNote } from '../layout/tracks.js';
import { ZoomInText } from '../tracks-draw.js';
import { BigRuntime, bigSpec } from './fixtures.js';
import { settle } from './helpers.js';

function texts(f: Frame): string[] {
  const out: string[] = [];
  for (const p of Object.values(f.lists.passes))
    for (let i = 0; i < p.texts.n; i++) out.push(f.lists.string(p.texts.str[i]));
  return out;
}

const Size = { width: 720, height: 300 };

describe('the zoom-in state', () => {
  it('a layer past maxTilesPerLayer draws its band hatched with “Zoom in to read: ⌘/Ctrl+scroll, or +”, said once', async () => {
    // 2 tiles a layer at most: 30 kb at its level (16 positions a bin, 16 kb a tile) with its margins are more
    const s = await settle(bigSpec(), [BigRuntime], {}, {}, { maxTilesPerLayer: 2, summaries: false });
    const f = await s.frame({ size: Size, viewport: { start: 0, end: 30_000 } });
    expect(texts(f)).toContain(ZoomInText);
    expect(f.layout.issues.some((i) => i.code === 'seq.render.zoom_in')).toBe(true);
    expect(s.store.issues.value.some((i) => i.code === 'seq.store.too_many_tiles')).toBe(true);
    // zoomed in, the layers read again: no text, no issue
    const near = await s.frame({ size: Size, viewport: { start: 1000, end: 1400 } });
    expect(texts(near)).not.toContain(ZoomInText);
    expect(near.layout.issues.some((i) => i.code === 'seq.render.zoom_in')).toBe(false);
  });

  it('a pending band says nothing (it settles by itself)', async () => {
    const s = await settle(bigSpec(), [BigRuntime]);
    const f = s.renderer.frame({ size: Size, viewport: { start: 0, end: 3_000_000 } });
    expect(texts(f)).not.toContain(ZoomInText);
  });

  /** A one-level source over a long system (a VCF over chr17's 83 Mb, a matrix of items only): no coarser level. */
  function oneLevel(shape: 'ranges' | 'matrix', length: number, readCost?: () => number | undefined) {
    const coverage: Coverage = {
      systems: [{ name: 'chr', length }],
      shapes: [{ system: 'chr', shape, ...(shape === 'matrix' && { states: ['a', 'b'] }) }],
      levels: [1],
    };
    const opened: OpenedSource = {
      coverage,
      read: async () => ({ ok: true, tables: {} }),
      ...(readCost ? { readCost } : {}),
      dispose() {},
    };
    const runtime: FormatRuntime = { ...SeqRuntime, open: async () => opened };
    const b = createSeqSpec({ title: 'one level' });
    const src = b.data.source('src', { format: 'seq', value: 'x' as never });
    const chr = b.data.sequence('chr', { alphabet: 'dna', source: src, name: 'chr' });
    const node =
      shape === 'ranges'
        ? b.data.ranges('v', { target: chr, source: src })
        : b.data.matrix('v', { target: chr, source: src } as never);
    b.views.view('v', (v) => {
      v.section('chr', { of: chr }).track('v', { height: 40 }).layer('v', { type: shape, data: node });
    });
    return { spec: b.build(), runtime };
  }

  for (const [what, shape, length] of [
    ['a variants track over 83 Mb at one level of 16 kb tiles (about 5,000)', 'ranges', 83_000_000],
    ['a matrix of items only over 300 kb', 'matrix', 300_000],
  ] as const)
    it(`${what}: one want past its own tile bound says zoom in, never an endless wait`, async () => {
      const { spec, runtime } = oneLevel(shape, length);
      const s = await settle(spec, [runtime]);
      const f = await s.frame({ size: Size });
      expect(texts(f)).toContain(ZoomInText);
      expect(f.layout.issues.filter((i) => i.code === 'seq.render.zoom_in').length).toBe(1);
      // zoomed in, it reads
      const near = await s.frame({ size: Size, viewport: { start: 0, end: 20_000 } });
      expect(texts(near)).not.toContain(ZoomInText);
    });

  // the hatch says "the source has data this level cannot draw"; where its index says nothing is stored over the whole
  // span, the track says "no items" in its quiet row instead (never a hatch over data that is simply not there)
  it('over a span an index knows is empty: the “no items” row, never the hatch; over one with data, the hatch', async () => {
    const empty = oneLevel('ranges', 83_000_000, () => 0);
    const s = await settle(empty.spec, [empty.runtime]);
    const f = await s.frame({ size: Size, settled: true });
    expect(texts(f)).toContain(EmptyNote);
    expect(texts(f)).not.toContain(ZoomInText);
    expect(f.layout.issues.some((i) => i.code === 'seq.render.zoom_in')).toBe(false);
    // the same source with bytes stored in the span (or its cost unknown): the hatch, no note
    for (const cost of [() => 4096, () => undefined]) {
      const full = oneLevel('ranges', 83_000_000, cost);
      const r = await settle(full.spec, [full.runtime]);
      const g = await r.frame({ size: Size, settled: true });
      expect(texts(g)).toContain(ZoomInText);
      expect(texts(g)).not.toContain(EmptyNote);
      expect(g.layout.issues.filter((i) => i.code === 'seq.render.zoom_in').length).toBe(1);
    }
  });
});
