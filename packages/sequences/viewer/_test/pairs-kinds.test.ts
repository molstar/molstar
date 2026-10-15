/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Several kinds of contact at one cell (the `pairs` renderer; seqspec:sources#levels-and-the-zoom-gate): when pair ends
// of different `category` share a cell, one glyph per distinct kind is drawn side by side in it, shrunk to fit and past
// what fits the first ones and a "+n"; each glyph is its own hit target (the hit's `part` is its kind), and hovering
// one lights only that kind's pairs at that end, with their arcs. A cell of one kind draws and hits as before.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { describe, expect, it } from 'vitest';
import { Hit, HoverKind, type SeqViewerContext } from '../src/index.js';
import { CELL_MIN_GLYPH, CellIO, cellLayout } from '../src/pairs-index.js';
import type { PairsState } from '../src/prepare.js';
import { drawn, layerOf, record } from './draw-fixtures.js';
import { lit } from './overlay-fixtures.js';

/** A cell state for `cellLayout`: cell 0 holds `kinds` kinds. */
const cellOf = (kinds: number) => ({ cellKinds: Int32Array.of(kinds) }) as unknown as PairsState;
function layout(kinds: number, room: number) {
  const io = new Float64Array(CellIO.SIZE);
  io[CellIO.width] = room;
  io[CellIO.widest] = 7;
  const shown = cellLayout(cellOf(kinds), 0, io);
  return { shown, pitch: io[CellIO.pitch], hidden: io[CellIO.hidden], marker: io[CellIO.marker] };
}

describe("a cell's kinds (`cellLayout`)", () => {
  it('one kind, or no room for two glyphs: not split', () => {
    expect(layout(1, 30).shown).toBe(0);
    expect(layout(3, 2 * CELL_MIN_GLYPH - 1).shown).toBe(0);
  });
  it("every kind fits: one glyph each, at most the lone glyph's width, shrunk to fit", () => {
    expect(layout(2, 30)).toEqual({ shown: 2, pitch: 7, hidden: 0, marker: 0 });
    expect(layout(4, 20)).toEqual({ shown: 4, pitch: 5, hidden: 0, marker: 0 });
    expect(layout(5, 5 * CELL_MIN_GLYPH)).toEqual({ shown: 5, pitch: CELL_MIN_GLYPH, hidden: 0, marker: 0 });
  });
  it("any number of kinds: the first ones that fit and a '+n' marker for the rest", () => {
    for (const [kinds, room] of [
      [10, 30],
      [25, 60],
      [120, 100],
      [3, 12],
    ]) {
      const l = layout(kinds, room);
      expect(l.shown).toBeGreaterThanOrEqual(1);
      expect(l.shown + l.hidden).toBe(kinds);
      expect(l.pitch).toBeGreaterThanOrEqual(CELL_MIN_GLYPH);
      // Glyphs and marker within the cell.
      expect(l.pitch * l.shown + l.marker).toBeLessThanOrEqual(room + 1e-9);
    }
    expect(layout(10, 30)).toMatchObject({ shown: 4, hidden: 6 });
  });
  it('too little room for glyphs and a marker: not split', () => {
    expect(layout(6, 12).shown).toBe(0);
  });
});

const KINDS = ['hbond', 'hydrophobic', 'hbond', 'pi_stacking'];
/** Residue 5 touches four partners (residues 6 on) with three kinds; each partner has one. */
function doc(o: { n?: number; kinds?: string[] } = {}) {
  const kinds = o.kinds ?? KINDS;
  const b = createSeqViewSpec({ title: 'Kinds' });
  const a = b.data.sequence({ ref: 'a', label: 'A', value: 'ACDEFGHIKLMNPQRSTVWY'.repeat((o.n ?? 20) / 20) });
  const partners = kinds.map((_, i) => 6 + i);
  const pairs = b.data.annotation({
    ref: 'pp',
    geometry: 'pairs',
    target: a,
    value: {
      start: partners.map(() => 5),
      b_target: partners.map(() => 'a'),
      b_start: partners,
      category: kinds,
    },
  } as never);
  const v = b.view({ ref: 'v', label: 'V', height: 'box' }).block({ sequence: a });
  v.track({ label: 'Letters' }).letters();
  v.track({ label: 'Contacts' }).pairs({
    data: pairs,
    arcs: 'hover',
    glyph: { theme: 'field_map', params: { map: { hbond: 'diamond', hydrophobic: 'circle', pi_stacking: 'square' } } },
  } as never);
  return b.getDocument();
}
const stateOf = (ctx: SeqViewerContext) => {
  const { k } = layerOf(ctx, 'Contacts');
  return ctx.graph.assembled.value?.view?.prepared[k].state as PairsState;
};

describe('several kinds at one cell (draw and hit)', () => {
  it('indexes the kinds of a cell: residue 5 holds three, in kind order, each with its first end', async () => {
    const t = await drawn(doc(), { width: 600, height: 200 });
    const s = stateOf(t.ctx);
    expect(s.count).toBe(4);
    const c = s.endCell[s.endAt[0]];
    expect(s.cellKinds[c]).toBe(3);
    const kinds = Array.from(s.kindList.subarray(s.cellKindAt[c], s.cellKindAt[c] + 3));
    expect(kinds).toEqual([...kinds].sort((x, y) => x - y));
    // Pairs 0 and 2 are the same kind; the kinds' first ends are pairs 0, 1 and 3.
    expect(s.kind[0]).toBe(s.kind[2]);
    const reps = Array.from(s.kindRep.subarray(s.cellKindAt[c], s.cellKindAt[c] + 3), (q) => s.endOf[q] >> 1);
    expect(reps).toEqual([0, 1, 3]);
  });

  it("draws one glyph per kind side by side in the cell, a lone kind's end as before", async () => {
    const t = await drawn(doc(), { width: 600, height: 200 });
    const g = record(t.ctx, 'Contacts').rec.glyphItems;
    // Three at residue 5's cell and the four partners' ends.
    expect(g.n).toBe(3 + 4);
    const x = Array.from(g.x.subarray(0, g.n));
    const cell = (t.ctx.layout.axis(layerOf(t.ctx, 'Contacts').section)?.cell ?? 0) as number;
    const split = x.slice(0, 3);
    expect(cell).toBeGreaterThanOrEqual(3 * CELL_MIN_GLYPH);
    // Side by side, a pitch apart, within one cell, the kinds' glyph shapes (diamond, circle, square) in order.
    expect(split[1] - split[0]).toBeCloseTo(split[2] - split[1], 6);
    expect(split[1] - split[0]).toBeGreaterThanOrEqual(CELL_MIN_GLYPH - 1e-9);
    expect(split[2] - split[0]).toBeLessThanOrEqual(cell);
    expect(new Set(Array.from(g.glyph.subarray(0, 3))).size).toBe(3);
    // The glyphs are those of the kinds, in order of kind, none repeated for the second hbond.
    expect(Array.from(g.radius.subarray(0, 3)).every((r) => r < 3.5)).toBe(true);
  });

  it("past what fits: the first kinds and a '+n' marker; a narrow cell draws one glyph", async () => {
    const kinds = Array.from({ length: 12 }, (_, i) => `k${i}`);
    const wide = await drawn(doc({ kinds, n: 20 }), { width: 600, height: 200 });
    const s = stateOf(wide.ctx);
    expect(s.cellKinds[s.endCell[s.endAt[0]]]).toBe(12);
    const r = record(wide.ctx, 'Contacts').rec;
    const marker = Array.from({ length: r.textItems.n }, (_, i) => r.textItems.text[i]).filter((x) =>
      x.startsWith('+'),
    );
    expect(marker.length).toBe(1);
    expect(Number(marker[0].slice(1))).toBeGreaterThan(0);
    // The kinds that fit, then each partner's end (residues 6 to 17).
    expect(r.glyphItems.n).toBe(12 - Number(marker[0].slice(1)) + 12);
  });

  it('a cell too narrow for its kinds and a marker draws as before: a glyph per end', async () => {
    const narrow = await drawn(doc({ n: 40 }), { width: 600, height: 200 });
    expect(record(narrow.ctx, 'Contacts').rec.glyphItems.n).toBe(4 + 4);
  });

  it("hits the glyph under the pointer: its kind is the hit's part and its pair the item", async () => {
    const t = await drawn(doc(), { width: 600, height: 200 });
    const { section } = layerOf(t.ctx, 'Contacts');
    const { rec } = record(t.ctx, 'Contacts');
    const s = stateOf(t.ctx);
    const h = Hit.create();
    const xs = Array.from(rec.glyphItems.x.subarray(0, 3));
    const y = s.frame[0];
    const found: number[][] = [];
    for (const x of xs) {
      expect(t.ctx.layout.hitTest(section, x, y, h)).toBe(true);
      expect(h.kind).toBe(HoverKind.item);
      found.push([h.part, h.item]);
    }
    // Three kinds, three parts, each the item of a pair of that kind.
    expect(new Set(found.map((f) => f[0])).size).toBe(3);
    for (const [i, [part]] of found.entries()) expect(s.kind[[0, 1, 3][i]]).toBe(part);
    // A lone end hits as before: no part.
    t.ctx.layout.hitTest(section, rec.glyphItems.x[3], y, h);
    expect([h.kind, h.part]).toEqual([HoverKind.item, -1]);
  });
});

describe('hovering a glyph of a split cell', () => {
  it("lights only that kind's pairs at that end, with their arcs", async () => {
    const t = await lit(doc(), { width: 600, height: 200 });
    const { section } = layerOf(t.ctx, 'Contacts');
    const { rec } = record(t.ctx, 'Contacts');
    const s = stateOf(t.ctx);
    const xs = Array.from(rec.glyphItems.x.subarray(0, 3));
    const y = s.frame[0];
    const base = t.overlay().pathItems.n;
    expect(section).toBe(t.section);
    // The arcs on the overlay per glyph: hbond two pairs, hydrophobic one, pi stacking one.
    const arcs: number[] = [];
    for (const x of xs) {
      t.hover(x, y);
      arcs.push(t.overlay().pathItems.n - base);
    }
    expect(arcs).toEqual([2, 1, 1]);
    // A cell of one kind lights the pair as before.
    t.hover(rec.glyphItems.x[3], y);
    expect(t.overlay().pathItems.n - base).toBe(1);
  });
});
