/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The cases of the channels' hover guards (IN3; channels.heap.ts) and their quick twin (channels.quick.ts), a step an
// element change: `gene`, a nucleotide or a residue hovered by turns, set through the hover manager, then the channels'
// rebuild (closed, projected into every section's hover set, its primary) and the overlay's reads; `alignment`, a
// million element changes over 10,000 distinct alignment rows and their positions (a row's residue on its own columns
// read without the row's cached map, which allocates on a miss). The heap guard and the quick twin share them.

import { CoordinateSystemKind as K, type SeqSet } from '@molstar/seq-core/model';
import { gene, rows } from '../../core/_test/format/project-fixtures.js';
import { ChannelMarks, Hit, HitRegion, HoverKind, MarkBit } from '../src/index.js';
import { shown } from './projection-fixtures.js';

export const ROWS = 10_000;

let sink = 0;
/** What the steps summed, for the guards' non-vacuity. */
export const sunk = () => sink;

/** A hit on a residue of `unit`, reused. */
function seqHit(set: SeqSet, unit: number) {
  const hit = Hit.create();
  Object.assign(hit, { section: 0, region: HitRegion.data, kind: HoverKind.seq, block: 0 });
  Object.assign(hit.seq, { set, unit });
  Object.assign(hit.system, { kind: K.unit, id: unit, row: -1 });
  return hit;
}

/** A nucleotide or a residue hovered by turns, the channels synced, a section's marks read. */
export async function genes(): Promise<{ step: (i: number) => void; ch: ChannelMarks }> {
  const g = await shown(gene({ chromosome: true }));
  // The channels alone: the hover label's churn per element change is IN2's guard (input.heap.ts), not this one.
  g.ctx.behaviors.get('hoverLabel')?.remove();
  const ch = new ChannelMarks(g.ctx);
  const [gu, pu] = [g.unit('gene/unit'), g.unit('protein/unit')];
  const [nt, aa] = [seqHit(g.set, gu), seqHit(g.set, pu)];
  const step = (i: number) => {
    const h = i & 1 ? aa : nt;
    const x = i & 1 ? i % 20 : i % 120;
    h.index = h.shown = x;
    h.seq.index = x as never;
    g.ctx.hover.setFromHit(h, 'pointer');
    ch.sync();
    const m = ch.marks(i % 3);
    sink += (m?.bitsAt(0, x) ?? 0) + (m?.primary(0) ?? 0) + (m?.any(0, 40, 80, MarkBit.hover) ? 1 : 0);
  };
  return { step, ch };
}

/** Element changes over `ROWS` distinct rows and positions of an alignment. */
export async function alignment(): Promise<{ step: (i: number) => void; ch: ChannelMarks }> {
  const r = await shown(rows(ROWS));
  r.ctx.behaviors.get('hoverLabel')?.remove();
  const ch = new ChannelMarks(r.ctx);
  const hit = Hit.create();
  Object.assign(hit, { section: 0, region: HitRegion.data, kind: HoverKind.aligned, block: 0 });
  Object.assign(hit.aligned, { set: r.set, alignment: 0 });
  Object.assign(hit.system, { kind: K.row, id: 0, row: 0 });
  const step = (i: number) => {
    const row = (i * 7919) % ROWS; // every row, in no order
    const column = i % 5;
    hit.dataRow = hit.system.row = row;
    hit.index = hit.shown = column;
    hit.aligned.row = row as never;
    hit.aligned.column = column as never;
    r.ctx.hover.setFromHit(hit, 'pointer');
    ch.sync();
    sink += ch.marks(0)?.hover.count ?? 0;
  };
  return { step, ch };
}
