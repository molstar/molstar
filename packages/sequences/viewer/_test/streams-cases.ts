/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The stream manager's region-cache cases, for their heap guard and its quick twin (streams.heap.ts, streams.quick.ts):
// a pass over case 9's cached view (18,000 lazy rows, a remote source) and over an alignment in memory (local tiles) with
// the view moving among held regions, and what a draw reads per frame of a layer (its lanes in view, a block's levels,
// its regions and their states). A step is one frame's reads; every move's view is loaded first, so the measured steps
// find their regions held.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { TILE } from '@molstar/seq-core/format';
import { case9 } from '../../core/_test/format/cases.js';
import type { StreamManager } from '../src/index.js';
import { blockLevels, lanesInView, type Streams } from '../src/streams.js';
import { type Opened, opened, preparedRenderers } from './layout-fixtures.js';

const P = { renderers: preparedRenderers() };
const flush = () => new Promise((ok) => setTimeout(ok, 0));
let sink = 0;

function inMemory(rows: number, columns: number) {
  const AA = 'ACDEFGHIKLMNPQRSTVWY';
  const lines = ['# STOCKHOLM 1.0'];
  for (let k = 0; k < rows; k++) {
    let r = '';
    for (let c = 0; c < columns; c++) r += (k * 7 + c) % 7 === 0 ? '-' : AA[(c * 7 + k) % 20];
    lines.push(`R${k}/1-${columns} ${r}`);
  }
  lines.push('//', '');
  const b = createSeqViewSpec();
  const aln = b.data.alignment({
    ref: 'big',
    source: { value: lines.join('\n'), format: 'stockholm', max_bytes: 1 << 24 },
  });
  const v = b.view({ ref: 'large', label: 'Large', style: 'bare', height: 'box' });
  v.block({ alignment: aln });
  v.track({ label: 'Rows' }).rows();
  return b.getDocument();
}

/** Every move's view loaded first, so the measured passes find their regions held. */
async function settled(t: Opened, section: string, moves: readonly (readonly [number, number])[]) {
  const s = t.ctx.streams as Streams;
  for (let i = 0; i < 2 * moves.length; i++) {
    t.vp(section).panBy(moves[i % moves.length][0], moves[i % moves.length][1], 'pointer');
    t.layout(section);
    for (let k = 0; k < 4; k++) {
      s.pass();
      await flush();
    }
  }
  return s;
}

const MOVES = [
  [0, 7],
  [7, 0],
  [0, -7],
  [-7, 0],
] as const;

/** Case 9's view at row 4,711 with its moves' regions held, and a step: a move among them and a pass. */
export async function case9Pass() {
  const t = await opened(case9(), { ...P, size: { 'msa.main': [1000, 500] } });
  const vp = t.vp('msa.main');
  vp.scrollTo(0, 4711 * 18, 'api');
  const s = await settled(t, 'msa.main', MOVES);
  const step = (i: number) => {
    const m = MOVES[i & 3];
    vp.panBy(m[0], m[1], 'pointer');
    s.pass();
  };
  return { s, step };
}

/** An alignment in memory zoomed out (2,000 rows × 900 columns, local tiles held): a move among held regions and a pass. */
export async function memoryPass() {
  const t = await opened(inMemory(2000, 900), { ...P, size: { 'large.main': [600, 400] } });
  const vp = t.vp('large.main');
  vp.zoomTo(0, 0, 900);
  vp.setRowScale(0.05, 0);
  const s = await settled(t, 'large.main', MOVES);
  const computed = s.stats.current.computed;
  const step = (i: number) => {
    const m = MOVES[i & 3];
    vp.panBy(m[0], m[1], 'pointer');
    s.pass();
  };
  return { s, step, computed };
}

/** What a draw reads per frame of case 9's layer: its lanes in view, a block's levels, its regions and their states. */
export async function drawReads() {
  const t = await opened(case9(), { ...P, size: { 'msa.main': [1000, 500] } });
  const s = await settled(t, 'msa.main', MOVES);
  const m = t.ctx.graph.assembled.value?.view;
  let layer = '';
  for (let k = 0; m && k < m.resolved.layers.count; k++) {
    const ref = m.set.strings.get(m.resolved.layers.ref[k]);
    if ((s as StreamManager).get(ref)) layer = ref;
  }
  const view = s.get(layer);
  const l = t.layout('msa.main');
  const v = t.ctx.layout.axis('msa.main');
  if (!view || !v) throw new Error('no stream');
  const lanes = new Float64Array(6);
  const lv = new Int32Array(4);
  const step = (i: number) => {
    if (lanesInView(l, v, lanes)) blockLevels(v, 0, view.budget, lanes, lv);
    const ty = Math.floor(lanes[0] / (TILE * 2 ** lv[0])) + (i & 1);
    const tile = view.tile(lv[0], lv[1], ty, i & 1);
    sink = (sink + view.state(lv[2], lv[3], i & 3, 0) + (tile ? tile.cls[i & 0xffff] : 0) + view.version) & 0xffff;
  };
  return { step };
}
