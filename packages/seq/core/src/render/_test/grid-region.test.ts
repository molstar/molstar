/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A grid read by region (_spec/README.md "Render", The grid layer: by region): store tiles drawn as image tiles at the
// finest level whose bin is at least a device pixel, the mirror of a symmetric grid's served tile, a pending tile a
// pending rect with nothing held under it and the coarser held tile once one is (EC-169, EC-223); a cell hit with its
// bin (Default 4).

import { type GridDense, SeqRuntime, type SeqData, seqToValue } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import type { GridLay } from '../grid.js';
import { gridCellHover, hitView } from '../index.js';
import { renderSvg } from '../render-svg.js';
import { settle } from './helpers.js';
import { gridRegionSpec as spec } from './region-render-fixtures.js';

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

  it('past the layer bound at its level, reads coarser before it says zoom in (Decided 16)', async () => {
    const s = await settle(spec(), [SeqRuntime], {}, {}, { maxTilesPerLayer: 2 });
    const near = { size: Size, viewport: { start: 0, end: 300 } };
    const f = await s.frame(near);
    const l = f.layout.sections[1].tracks[0].layers[0] as { lay: GridLay; zoomIn?: boolean };
    // level 1 wants 2 × 2 tiles over the view: level 4 one
    expect(l.zoomIn).toBeUndefined();
    expect(l.lay.binLevel).toBe(4);
    expect(l.lay.images.length).toBeGreaterThan(0);
  });

  it('two region grids on one pair, the upper one map and the lower another, both drawn (EC-59)', async () => {
    const s = await settle(spec(true), [SeqRuntime]);
    const f = await s.frame({ size: Size });
    const layers = f.layout.sections[1].tracks[0].layers;
    expect(layers.length).toBe(2);
    const [upper, lower] = layers.map((l) => l.lay as GridLay);
    expect(upper.images.length).toBeGreaterThan(0);
    expect(lower.images.length).toBeGreaterThan(0);
    // each keeps its own triangle of its own map: a cell above the diagonal is the first's, below the second's
    expect(upper.triangle).toBe('upper');
    expect(lower.triangle).toBe('lower');
    // the first map's value above the diagonal (bins 10 and 2 at level 4), the second's three times it below
    expect(gridCellHover(upper, 40, 8, () => 'chr')!.rows).toContainEqual(['Value', '177.8']);
    expect(gridCellHover(lower, 8, 40, () => 'chr')!.rows).toContainEqual(['Value', '533.3']);
  });

  it('an export waits for its 2D tiles and draws them (EC-125)', async () => {
    const out = await renderSvg(spec(), { size: { width: 600, height: 900 }, runtimes: [SeqRuntime] });
    expect(out.settled).toBe(true);
    expect(out.svg).toContain('<image');
  });

  it('a tile whose last bin runs past its system’s end is drawn on both axes, its cells at their true positions, clipped (EC-223)', async () => {
    // a length no multiple of the level (4): the last bin holds positions past the end, on x and on y
    const lay = await pastTheEnd(1001, 4);
    expect(lay.binLevel).toBe(4);
    expect(lay.images.length).toBe(1);
    const scale = lay.ic[2] / 1001;
    // the image spans its 251 cells at their positions (1,004), its clip the system's 1,001: no cell squeezed
    // (a pixel's snapping either side: past 1,002, never the 1,001 a squeezed tile would be)
    expect(lay.ib[2] / scale).toBeGreaterThan(1002);
    expect(lay.ib[2] / scale).toBeLessThan(1006);
    expect(lay.ic[0] + lay.ic[2]).toBeLessThanOrEqual(lay.box.x + lay.box.width + 1);
    expect(lay.ic[1] + lay.ic[3]).toBeLessThanOrEqual(lay.box.y + lay.box.height + 1);
  });

  it('a section’s range cutting bins at both edges: its edge tiles drawn, clipped to the range, cells never squeezed (EC-223)', async () => {
    // 1,000 to 2,999 at level 16: both edges inside a bin
    const lay = await pastTheEnd(4001, 16, [1000, 2999]);
    expect(lay.binLevel).toBe(16);
    expect(lay.images.length).toBeGreaterThan(0);
    const scale = lay.box.width / 1999;
    // the drawn part is the range (the clip), each image wider than its clip where a bin is cut
    let clipped = 0;
    for (let k = 0; k < lay.images.length; k++) {
      expect(lay.ic[4 * k]).toBeGreaterThanOrEqual(lay.box.x - 1);
      expect(lay.ic[4 * k] + lay.ic[4 * k + 2]).toBeLessThanOrEqual(lay.box.x + lay.box.width + 1);
      clipped += lay.ic[4 * k + 2];
      // a cell is 16 positions wide whatever the clip: the image's width over its cells
      expect(lay.ib[4 * k + 2] / lay.images[k].width / scale).toBeCloseTo(16, 0);
    }
    expect(clipped / scale).toBeCloseTo(1999, -1);
  });
});

/** A grid read by region over one system of length `n` at level `l` (a section's `range` when given), settled. */
async function pastTheEnd(n: number, l: number, range?: [number, number]): Promise<GridLay> {
  const cells = Math.ceil(n / l);
  const grid: GridDense = {
    form: 'dense',
    x: 0,
    y: 0,
    nx: cells,
    ny: cells,
    value: new Float32Array(cells * cells).fill(1),
  };
  const data: SeqData = {
    coverage: {
      systems: [{ name: 'chr', length: n }],
      shapes: [{ system: 'chr', shape: 'grid', y: ['chr'] }],
      levels: [l],
    },
    pairs: [
      {
        request: {
          system: 'chr',
          start: 0,
          end: n,
          level: l,
          shapes: ['grid'],
          y: { system: 'chr', start: 0, end: n },
        },
        answer: { ok: true, tables: { grid } },
      },
    ],
  };
  const b = createSeqSpec({ title: 'Past the end' });
  const src = b.data.source('src', { format: 'seq', value: seqToValue(data) as never });
  const chr = b.data.sequence('chr', { alphabet: 'dna', source: src, name: 'chr' });
  const g = b.data.grid('g', { x: chr, source: src } as never);
  b.views.view('main', (v) => {
    const top = v.section('chr', { of: chr, ...(range ? { range } : {}) });
    top.track('ruler').layer('ruler', { type: 'ruler' });
    v.section('plot', { layout: 'grid', of: chr, follow: top.address }).track('plot', {}, (t) =>
      t.layer('g', { type: 'grid', data: g } as never),
    );
  });
  const s = await settle(b.build(), [SeqRuntime]);
  await s.frame({ size: Size });
  return gridLay(await s.frame({ size: Size }));
}
