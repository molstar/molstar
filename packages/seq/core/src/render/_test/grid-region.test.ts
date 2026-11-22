/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A grid read by region (_spec/README.md "Render", The grid layer: by region): store tiles drawn as image tiles at the
// finest level whose bin is at least a device pixel, the mirror of a symmetric grid's served tile, a pending tile a
// pending rect with nothing held under it and the coarser held tile once one is (EC-169, EC-223); a cell hit with its
// bin (Default 4).

import type { GridDense, SeqData } from '@molstar/seqx-core/data';
import { SeqRuntime, seqToValue } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import type { GridLay } from '../grid.js';
import { gridCellHover, hitView } from '../index.js';
import { settle } from './helpers.js';

const N = 2048;

/** A contact map over chr (2,048 positions): a band along the diagonal, served symmetric at levels 1 and 4. */
function hic(): SeqData {
  const dense = (level: number): GridDense => {
    const n = N / level;
    const v = new Float32Array(n * n);
    for (let i = 0; i < n; i++)
      for (let j = 0; j < n; j++) v[i * n + j] = (level * level * 100) / (1 + Math.abs(i - j));
    return { form: 'dense', x: 0, y: 0, nx: n, ny: n, value: v };
  };
  const request = (level: number) => ({
    system: 'chr',
    start: 0,
    end: N,
    level,
    shapes: ['grid' as const],
    y: { system: 'chr', start: 0, end: N },
  });
  return {
    coverage: {
      systems: [{ name: 'chr', length: N }],
      shapes: [{ system: 'chr', shape: 'grid', symmetric: true, aggregate: 'sum' }],
      levels: [1, 4],
    },
    pairs: [1, 4].map((level) => ({
      request: request(level),
      answer: { ok: true as const, tables: { grid: dense(level) } },
    })),
  };
}

function spec() {
  const b = createSeqSpec({ title: 'Hi-C by region' });
  const src = b.data.source('src', { format: 'seq', value: seqToValue(hic()) as never });
  const chr = b.data.sequence('chr', { alphabet: 'dna', source: src, name: 'chr' });
  const g = b.data.grid('contacts', { x: chr, source: src } as never);
  b.views.view('main', (v) => {
    const top = v.section('chr', { of: chr });
    top.track('ruler').layer('ruler', { type: 'ruler' });
    v.section('plot', { layout: 'grid', of: chr, follow: top.address }).track('plot', {}, (t) => {
      t.layer('grid', { type: 'grid', data: g });
    });
  });
  return b.build();
}

const Size = { width: 600, height: 900 };

function gridLay(f: {
  layout: { sections: readonly { tracks: readonly { layers: readonly { lay: unknown }[] }[] }[] };
}) {
  return f.layout.sections[1].tracks[0].layers[0].lay as GridLay;
}

describe('a grid by region', () => {
  it('draws pending rects before its tiles land, then its store tiles as image tiles, a cell hit with its bin', async () => {
    const s = await settle(spec(), [SeqRuntime]);
    const first = gridLay(s.renderer.frame({ size: Size }));
    // nothing read yet: pending rects, never empty
    expect(first.pend.length).toBeGreaterThan(0);
    expect(first.images.length).toBe(0);
    const f = await s.frame({ size: Size });
    const lay = gridLay(f);
    expect(lay.pend.length).toBe(0);
    // 2,048 positions on about 600 px: 4 positions a bin read; tiles of 256 bins, their mirrors drawn too
    expect(lay.binLevel).toBe(4);
    expect(lay.images.length).toBe(4);
    // tiles meet across their edges: no gap on x
    const xs = new Set<number>();
    for (let k = 0; k < lay.images.length; k++) xs.add(lay.ib[4 * k] + lay.ib[4 * k + 2]);
    expect([...xs].some((x) => lay.ib.some((v, k) => k % 4 === 0 && v === x))).toBe(true);
    const h = hitView(f.layout, lay.box.x + lay.box.width / 3, lay.box.y + lay.box.height / 4);
    expect(h?.target).toBe('cell');
    expect(h?.bin).toBe(4);
    // its hover: its bin on each axis and its value (a mirror's cell as the served one's), its level
    const card = gridCellHover(lay, 10, 9, () => 'chr')!;
    expect(card.title).toBe('chr:9–12 × chr:9–12');
    expect(card.rows).toContainEqual(['Value', '1,600']);
    expect(card.rows).toContainEqual(['Bins', '4 bp bins']);
    // across the diagonal: the mirror of the served tile
    expect(gridCellHover(lay, 1000, 9, () => 'chr')!.rows).toContainEqual(['Value', '6.43']);
  });

  it('zoomed in, reads level 1 and draws the coarser held tile while its own tiles are pending', async () => {
    const s = await settle(spec(), [SeqRuntime]);
    await s.frame({ size: Size });
    const near = { size: Size, viewport: { start: 0, end: 300 } };
    const pending = gridLay(s.renderer.frame(near));
    expect(pending.binLevel).toBe(1);
    // level 4's tile is held: drawn under, no pending rect
    expect(pending.images.length).toBeGreaterThan(0);
    expect(pending.pend.length).toBe(0);
    // a cell of a tile not loaded at the level read: loading
    expect(gridCellHover(pending, 10, 10, () => 'chr')!.rows).toContainEqual(['Value', 'loading']);
    const lay = gridLay(await s.frame(near));
    expect(lay.binLevel).toBe(1);
    expect(lay.images.length).toBeGreaterThan(0);
  });

  it('its wants are bounded together: past `maxTilesPerLayer` none is read and its band says zoom in', async () => {
    const s = await settle(spec(), [SeqRuntime], {}, {}, { maxTilesPerLayer: 2 });
    const f = await s.frame({ size: Size });
    const l = f.layout.sections[1].tracks[0].layers[0] as { lay: GridLay; zoomIn?: boolean };
    expect(l.zoomIn).toBe(true);
    expect(l.lay.images.length).toBe(0);
    expect(f.layout.issues.map((i) => i.code)).toContain('seq.render.zoom_in');
  });
});
