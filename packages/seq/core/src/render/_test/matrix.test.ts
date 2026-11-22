/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The matrix layer (_spec/README.md "Render", The matrix layer): states as rows in the band, items at their own
// extents (a long deletion solid, drawn under the short items, EC-222), missing values transparent (EC-220), coarse bins
// as image tiles, pending across the band before its tiles land (EC-169), a cell hit as its item and state.

import type { MatrixBins, MatrixItems, SeqData } from '@molstar/seqx-core/data';
import { SeqRuntime, seqToValue } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { LightTokens } from '@molstar/seqx-core/render';
import { describe, expect, it } from 'vitest';
import { hitView } from '../index.js';
import type { MatrixLay } from '../matrix.js';
import { settle } from './helpers.js';

const N = 262_144;
const States = ['NA12878', 'NA12891', 'NA12892', 'HG00096'];
const S = States.length;
/** A ring's 40 states (a heatmap ring of samples, Circos). */
const Many = Array.from({ length: 40 }, (_, i) => `s${i}`);

/** Genotypes over chr: a SNP every 100 positions in 0..2,000, a deletion over 500..1,500 in the first state. */
function genotypes(names: readonly string[] = States): SeqData {
  const S = names.length;
  const start: number[] = [];
  const end: number[] = [];
  for (let p = 0; p < 2000; p += 100) start.push(p), end.push(p + 1);
  start.push(500), end.push(1500);
  const n = start.length;
  const codes = new Uint8Array(n * S);
  for (let i = 0; i < n; i++) for (let s = 0; s < S; s++) codes[i * S + s] = (i + s) % 3;
  // the deletion: every copy in the first state, none in the others; one state uncalled at the first SNP
  codes[(n - 1) * S] = 2;
  codes[1] = 255;
  const items: MatrixItems = {
    form: 'items',
    start: Int32Array.from(start),
    end: Int32Array.from(end),
    states: S,
    codes,
  };
  const binsN = N / 256;
  const classes = new Uint8Array(binsN * S);
  const shares = new Uint8Array(binsN * S).fill(200);
  for (let b = 0; b < binsN; b++) for (let s = 0; s < S; s++) classes[b * S + s] = (b + s) % 3;
  const bins: MatrixBins = { form: 'bins', start: 0, states: S, classes, shares };
  const request = (level: number) => ({ system: 'chr', start: 0, end: N, level, shapes: ['matrix' as const] });
  return {
    coverage: {
      systems: [{ name: 'chr', length: N }],
      shapes: [{ system: 'chr', shape: 'matrix', states: [...names] }],
      levels: [1, 256],
    },
    pairs: [
      { request: request(1), answer: { ok: true as const, tables: { matrix: items } } },
      { request: request(256), answer: { ok: true as const, tables: { matrix: bins } } },
    ],
  };
}

function sourced(names: readonly string[] = States, layout: 'tracks' | 'rings' = 'tracks') {
  const b = createSeqSpec({ title: 'Genotypes' });
  const src = b.data.source('src', { format: 'seq', value: seqToValue(genotypes(names)) as never });
  const chr = b.data.sequence('chr', { alphabet: 'dna', source: src, name: 'chr' });
  const m = b.data.matrix('gt', { target: chr, source: src } as never);
  b.views.view('main', (v) => {
    v.section('chr', { of: chr, layout }).track('gt', { height: 80 }, (t) => {
      t.layer('gt', { type: 'matrix', data: m });
    });
  });
  return b.build();
}

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
  it('draws inline values a state a row, a missing one transparent, a cell hit as its position and state', async () => {
    const s = await settle(inline(), [SeqRuntime]);
    const f = await s.frame({ size: { width: 400, height: 200 } });
    const lay = layOf(f);
    expect(lay.states).toBe(2);
    expect(lay.items).toBe(true);
    // 16 cells, one missing
    expect(lay.n).toBe(15);
    expect(lay.rowScale).toBeCloseTo(lay.box.height / 2);
    const h = hitView(f.layout as never, lay.axis!.x(3.5), lay.box.y + lay.box.height * 0.75);
    expect(h?.target).toBe('cell');
    expect(h?.position).toBe(3);
    expect(h?.yPosition).toBe(1);
    expect(h?.ySystem).toBe(lay.statesSystem);
    expect(h?.bin).toBeUndefined();
  });

  it('draws items at their extents, a long deletion solid and under the short ones; an uncalled code drawn missing', async () => {
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
    const h = hitView(f.layout as never, lay.box.x + lay.box.width / 2, lay.box.y + 1);
    expect(h?.bin).toBe(256);
    expect(h?.yPosition).toBe(0);
  });

  it('on a ring, draws its states as thin rings of runs; past `rings.imageRuns` its rows reduced, said once', async () => {
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
});
