/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The matrix layer (_spec/README.md "Render", The matrix layer): states as rows in the band, items at their own
// extents (a long deletion solid, drawn under the short items, EC-222), missing values transparent (EC-220), coarse bins
// as image tiles, pending across the band before its tiles land (EC-169), a cell hit as its item and state.

import { SeqRuntime } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import {
  colorHex,
  createTokens,
  DarkTokens,
  LightTokens,
  type SeqTokenName,
  SeqTokenNames,
} from '@molstar/seqx-core/render';
import { describe, expect, it } from 'vitest';
import { hitView, matrixCellHover } from '../index.js';
import type { MatrixLay } from '../matrix.js';
import { settle } from './helpers.js';
import { GenotypesLength as N, Many, States, sourced } from './region-render-fixtures.js';

const S = States.length;

function inline() {
  const b = createSeqSpec({ title: 'Methylation' });
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGTACGT' });
  // 8 positions × 2 states, a missing value at position 2 of the second state
  const values = Float32Array.from([0, 1, 0.5, 0.5, 0.2, Number.NaN, 1, 0, 0.1, 0.9, 0.3, 0.3, 0.7, 0.2, 0.4, 0.6]);
  const m = b.data.matrix('m', { target: s, values, states: ['tumour', 'normal'] });
  b.views.view('main', (v) => {
    v.section('s', { of: s }).track('m', { height: 40 }, (t) => {
      t.layer('m', { type: 'matrix', data: m });
    });
  });
  return b.build();
}

type Frame = { layout: { sections: readonly { lays: readonly (readonly { lay: unknown }[])[] }[] } };
/** The first strip's matrix lay (a short sequence wraps: its strips' lays, not its tracks'). */
const layOf = (f: Frame) => f.layout.sections[0].lays[0][0].lay as MatrixLay;

describe('the matrix layer', () => {
  it('colours genotype codes by tokens of their own in both themes, the value ramp changing none of them', async () => {
    const size = { width: 800, height: 200 };
    const near = { size, viewport: { start: 0, end: 2000 } };
    const cells = async (tokens: typeof LightTokens) => {
      const s = await settle(sourced(), [SeqRuntime]);
      s.renderer.tokens = tokens;
      const lay = layOf(await s.frame(near));
      const used = new Set<number>();
      for (let k = 0; k < lay.n; k++) used.add(lay.rc[k]);
      return used;
    };
    for (const t of [LightTokens, DarkTokens]) {
      const used = await cells(t);
      // the fixture holds codes 0, 1, 2 and uncalled (255)
      for (const name of ['seq-genotype-0', 'seq-genotype-1', 'seq-genotype-2', 'seq-missing'] as const) {
        expect(used.has(t.colors[name]), `${t.theme} ${name}`).toBe(true);
      }
      expect(used.size).toBe(4);
      // a theme changing the ramp repaints no genotype
      const hex = Object.fromEntries(SeqTokenNames.map((n) => [n, colorHex(t.colors[n])])) as Record<
        SeqTokenName,
        string
      >;
      for (const n of ['seq-ramp-0', 'seq-ramp-1', 'seq-ramp-2', 'seq-ramp-3', 'seq-ramp-4'] as const)
        hex[n] = '#123456';
      expect(await cells(createTokens(t.theme, hex))).toEqual(used);
    }
  });

  it('draws inline values a state a row, a missing one in the missing colour, a cell hit as its position and state (EC-220, EC-221)', async () => {
    const s = await settle(inline(), [SeqRuntime]);
    const f = await s.frame({ size: { width: 400, height: 200 } });
    const lay = layOf(f);
    expect(lay.states).toBe(2);
    expect(lay.items).toBe(true);
    // 16 cells, one missing: drawn in the missing colour, never the ramp's low end (dark over dark) nor 0
    expect(lay.n).toBe(16);
    const missing = LightTokens.colors['seq-missing'];
    expect(Array.from(lay.rc.subarray(0, lay.n)).filter((c) => c === missing)).toHaveLength(1);
    expect(lay.rowScale).toBeCloseTo(lay.box.height / 2);
    const h = hitView(f.layout as never, lay.axis!.x(3.5), lay.box.y + lay.box.height * 0.75);
    expect(h?.target).toBe('cell');
    expect(h?.position).toBe(3);
    expect(h?.yPosition).toBe(1);
    expect(h?.ySystem).toBe(lay.statesSystem);
    expect(h?.bin).toBeUndefined();
    expect(matrixCellHover(lay, 2, 1, () => 's')!.text).toBe('s:3, Sample normal · no data');
    expect(matrixCellHover(lay, 2, 0, () => 's')!.rows).toContainEqual(['Sample', 'tumour · 0.2']);
  });

  it('draws items at their extents, a long deletion solid and under the short ones; an uncalled code drawn missing (EC-222)', async () => {
    const s = await settle(sourced(), [SeqRuntime]);
    const size = { width: 800, height: 200 };
    const near = { size, viewport: { start: 0, end: 2000 } };
    const first = layOf(s.renderer.frame(near));
    // nothing read: pending across the band, never empty
    expect(first.pend.length).toBeGreaterThan(0);
    expect(first.n).toBe(0);
    const lay = layOf(await s.frame(near));
    expect(lay.pend.length).toBe(0);
    expect(lay.items).toBe(true);
    expect(lay.states).toBe(S);
    // the deletion: the first rects, one a state, a thousand positions wide
    const w = lay.axis!.x(1500) - lay.axis!.x(500);
    expect(lay.rx[2]).toBeCloseTo(w);
    expect(lay.rx[4 * 3 + 2]).toBeCloseTo(w);
    expect(lay.n).toBe(21 * S);
    // the uncalled state of the first SNP: the missing token
    const missing = LightTokens.colors['seq-missing'];
    let found = false;
    for (let k = 0; k < lay.n; k++) if (lay.rc[k] === missing) found = true;
    expect(found).toBe(true);
    // a hover: the shortest item there, its state, its genotype
    const name = () => 'chr';
    expect(matrixCellHover(lay, 500, 0, name)!.text).toBe('chr:501, Sample NA12878 · 1/1 (homozygous alternate)');
    expect(matrixCellHover(lay, 550, 0, name)!.text).toBe('chr:501–1,500, Sample NA12878 · 1/1 (homozygous alternate)');
    expect(matrixCellHover(lay, 550, 1, name)!.rows[1][1]).toMatch(/ · 0\/0 \(homozygous reference\)$/);
    expect(matrixCellHover(lay, 0, 1, name)!.rows[1][1]).toMatch(/ · \.\/\. \(uncalled\)$/);
  });

  it('zoomed out, draws the bins as an image tile of bins × states, a hit with its bin', async () => {
    const s = await settle(sourced(), [SeqRuntime]);
    const f = await s.frame({ size: { width: 800, height: 200 } });
    const lay = layOf(f);
    expect(lay.items).toBe(false);
    expect(lay.level).toBe(256);
    expect(lay.n).toBe(0);
    expect(lay.images.length).toBe(1);
    expect(lay.images[0].width).toBe(N / 256);
    expect(lay.images[0].height).toBe(S);
    // each bin's majority at its share's strength (200 of 255 here): the colour's alpha
    expect(lay.images[0].data[3]).toBe(200);
    const h = hitView(f.layout as never, lay.box.x + lay.box.width / 2, lay.box.y + 1);
    expect(h?.bin).toBe(256);
    expect(h?.yPosition).toBe(0);
    // a bin's hover: its range, its state, its majority with its share, its level
    expect(matrixCellHover(lay, 300, 1, () => 'chr')!.text).toBe(
      'chr:257–512, Sample NA12891 · 1/1 (homozygous alternate), 78% of the bin, 256 bp bins',
    );
  });

  it('on a ring, draws its states as thin rings of runs at the outer circumference’s columns; past `rings.imageRuns` its rows reduced, said once (EC-108)', async () => {
    const size = { width: 640, height: 900 };
    const s = await settle(sourced(Many, 'rings'), [SeqRuntime]);
    const f = await s.frame({ size });
    expect(f.layout.issues.map((i) => i.code)).toEqual([]);
    const data = f.lists.passes.data;
    // the image is not on the page: its rows' runs are, a sector each (1,024 bins a state, a class a bin)
    expect(data.images.n).toBe(0);
    expect(data.paths.n).toBeGreaterThanOrEqual(40 * 1000);
    const few = await settle(sourced(Many, 'rings'), [SeqRuntime], { rings: { imageRuns: 4096 } });
    const g = await few.frame({ size });
    expect(g.layout.issues.map((i) => [i.code, i.severity])).toEqual([['seq.render.ring.image', 'info']]);
    expect(g.lists.passes.data.paths.n).toBeLessThanOrEqual(4096 + 64);
  });

  it('past the store’s bounds at every level, its band says zoom in (the layer’s wants bounded together)', async () => {
    const s = await settle(sourced(), [SeqRuntime], {}, {}, { maxTilesPerLayer: 0 });
    const f = await s.frame({ size: { width: 800, height: 200 } });
    const l = f.layout.sections[0].tracks[0].layers[0] as { lay: MatrixLay; zoomIn?: boolean };
    expect(l.zoomIn).toBe(true);
    expect(l.lay.n + l.lay.images.length).toBe(0);
    expect(f.layout.issues.map((i) => i.code)).toContain('seq.render.zoom_in');
  });
});
