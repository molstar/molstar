/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A matrix layer whose blocks' annotations have different state counts (IX1, found in the pre-review): `prepareMatrix`
// took the last placement's count for the whole layer, and the tiles read `values[row · states + state]` with it, so a
// block whose annotation had another count read across its rows. The layer's bands are the most states any block has,
// each block's values are read with its own stride and a block with fewer states leaves its lower bands empty; the
// `value` label reads a block's value with its stride too, and none at a band past its states.
// (A resolved layer's matrices all share their states; a host's state, a hand-made `LayerInput`, meets this.)

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { type SeqSet, SeqSets } from '@molstar/seq-core/model';
import { describe, expect, it } from 'vitest';
import { drawMatrixTiles, Hit, HoverKind, type MatrixState, prepareMatrix } from '../src/draw/api.js';
import { HitRegion, KernelRenderers } from '../src/index.js';
import { LayerFrame, RecordingBackend } from '../src/testing/index.js';
import { opened } from './layout-fixtures.js';

const N = 12;
/** The value of state q of position i in block b's annotation (a distinct number per cell; never missing). */
const value = (b: number, i: number, q: number) => 1 + b + i + 10 * q;

function twoMatrices() {
  const states = [5, 3];
  const { set, issues } = SeqSets.build({
    molecules: [
      { key: 'a', residues: 'A'.repeat(N) },
      { key: 'b', residues: 'C'.repeat(N) },
    ],
    units: [
      { key: 'A', molecule: 'a' },
      { key: 'B', molecule: 'b' },
    ],
    annotations: states.map((n, b) => ({
      key: `m${b}`,
      geometry: 'matrix',
      target: 'AB'[b],
      states: Array.from({ length: n }, (_, q) => `s${q}`),
      columns: {
        index: Array.from({ length: N }, (_, i) => i),
        value: Array.from({ length: N }, (_, i) => Array.from({ length: n }, (_, q) => value(b, i, q))),
      },
    })) as never,
  });
  expect(issues.filter((i) => i.severity === 'error')).toEqual([]);
  const anns = [0, 1].map((b) => {
    const a = set.annotations.find((x) => x && set.strings.get(x.key) === `m${b}`);
    if (!a) throw new Error('no annotation');
    return a;
  });
  return { set, anns };
}

describe('prepareMatrix over blocks with different state counts', () => {
  it('the bands are the most states; each block keeps the stride of its own annotation', () => {
    const { set, anns } = twoMatrices();
    const input = {
      set,
      section: { blocks: [0, 1].map(() => ({ length: N, base: 0 })) },
      layer: {
        // Block 1's placement comes last: its 3 states are what the layer used to take.
        placements: anns.map((a, block) => ({ annotation: a.id, lane: -1, block, path: -1, slot: a.slot[0] })),
        paths: [],
      },
      annotations: anns,
      params: {},
    };
    const s = prepareMatrix(input as never);
    expect([s.states, Array.from(s.stride)]).toEqual([5, [5, 3]]);
    expect(s.labels).toEqual(['s0', 's1', 's2', 's3', 's4']);
    expect(s.rows.map((r) => r?.length)).toEqual([N, N]);
  });
});

/** A document with the two matrices' layers in two tracks of one view (blocks A and B); the first's state, joined. */
async function joinedLayer() {
  const b = createSeqViewSpec({ title: 'Strides' });
  const seqs = ['A', 'B'].map((k) => b.data.sequence({ ref: `s${k}`, label: k, value: 'AC'.repeat(N / 2) }));
  const units = ['A', 'B'].map((k, i) => b.data.unit({ ref: `u${k}`, sequence: seqs[i], label: k } as never));
  const ms = [5, 3].map((n, k) =>
    b.data.annotation({
      ref: `m${k}`,
      geometry: 'matrix',
      target: units[k],
      states: Array.from({ length: n }, (_, q) => `s${q}`),
      value: {
        index: Array.from({ length: N }, (_, i) => i),
        value: Array.from({ length: N }, (_, i) => Array.from({ length: n }, (_, q) => value(k, i, q))),
      },
    } as never),
  );
  const v = b.view({ ref: 'v', label: 'V' }).block({ unit: units[0] }).block({ unit: units[1] });
  for (const [k, m] of ms.entries()) v.track({ ref: `v.m${k}`, label: `M${k}` }).matrix({ data: m } as never);
  const t = await opened(b.getDocument(), { renderers: [...KernelRenderers] });
  const m = t.ctx.graph.assembled.value?.view;
  const prepared = [0, 1].map((k) => m?.prepared.findIndex((p, i) => p.renderer === 'matrix' && i === k) ?? -1);
  const [s0, s1] = prepared.map((k) => m?.prepared[k].state as MatrixState);
  const matrix = t.ctx.registry.get('renderers', 'matrix');
  if (!m || !matrix || !s0 || !s1) throw new Error('no matrix layers');
  // Both blocks in one state, as a host's single layer would hold them.
  const joined: MatrixState = {
    ...s0,
    rows: [s0.rows[0], s1.rows[1]],
    ann: [s0.ann[0], s1.ann[1]],
    states: 5,
    stride: Int32Array.of(5, 3),
  };
  return { t, m, matrix, prepared, joined };
}

describe("the tiles read each block's values with its own stride", () => {
  it('block B (3 states) shows its values in its first 3 bands and nothing below, beside block A (5)', async () => {
    const { t, m, matrix, prepared, joined } = await joinedLayer();
    const rec = new RecordingBackend();
    const frame = new LayerFrame(t.ctx, matrix, 'v.main', prepared[0]);
    frame.draw(rec);
    frame.v.backend = rec;
    rec.reset();
    drawMatrixTiles(frame.v, joined);
    const images = rec.images.content.filter((c): c is Uint32Array => c instanceof Uint32Array);
    // One image per block (both fit one tile), a row of N samples per band.
    const mine = images.filter((c) => c.length === 5 * N);
    expect(mine.length).toBe(2);
    const theme = m.prepared[prepared[0]].themes[0];
    const px = (v: number) => {
      const c = theme?.value({ ...frame.v.loc, value: v, state: 0, index: 0 } as never) ?? 0;
      return ((255 << 24) | ((c & 255) << 16) | (c & 0xff00) | ((c >> 16) & 255)) >>> 0;
    };
    const [a, bImg] = mine;
    for (let i = 0; i < N; i++) {
      for (let q = 0; q < 5; q++) expect(a[q * N + i], `A ${i} ${q}`).toBe(px(value(0, i, q)));
      for (let q = 0; q < 3; q++) expect(bImg[q * N + i], `B ${i} ${q}`).toBe(px(value(1, i, q)));
      for (let q = 3; q < 5; q++) expect(bImg[q * N + i], `B ${i} ${q} empty`).toBe(0);
    }
  });
});

describe("the value label reads each block's values with its own stride", () => {
  it("block B (3 states): its value at a state it has, none past its states; block A's at its fifth", async () => {
    const { t, m, prepared, joined } = await joinedLayer();
    const k = prepared[0];
    (m.prepared[k] as { state: unknown }).state = joined;
    const r = m.resolved;
    const set = t.ctx.graph.assembled.value?.set as SeqSet;
    const hit = Hit.create();
    /** The value line of a hover at block `block`, position `j`, band `q` of the joined layer. */
    const line = (block: number, j: number, q: number) => {
      const g = r.sections.blockStart[0] + block;
      Object.assign(hit, {
        section: 0,
        region: HitRegion.data,
        kind: HoverKind.seq,
        layer: k,
        track: r.layers.track[k],
        block,
        index: r.blocks.base[g] + j,
        shown: j,
        state: q,
      });
      Object.assign(hit.seq, { set, unit: r.blocks.unit[g], index: j });
      t.ctx.hover.setFromHit(hit, 'pointer');
      const label = t.ctx.labels.compute(t.ctx.hover.live.current);
      return (label?.lines ?? []).filter((l) => l.key === `s${q}`).map((l) => l.value);
    };
    expect(line(1, 4, 2)).toEqual([`${value(1, 4, 2)}`]);
    expect(line(1, 4, 3)).toEqual([]);
    expect(line(1, 4, 4)).toEqual([]);
    expect(line(0, 4, 4)).toEqual([`${value(0, 4, 4)}`]);
    t.ctx.dispose();
  });
});
