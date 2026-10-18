/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// RN2's draws on the recording backend (seqspec:viewer#extension-points; seqspec:sources#levels-and-the-zoom-gate):
// every feature, pair type draws at each level of detail and below the row's minimum; features narrower than a pixel
// are drawn; `ranges` in lanes, with labels cut by their room, moving with their range on a pan and sticking 4 px
// inside the edge that cuts it, brackets, boundaries, the hairline across a gap and (wrapped) underlines on the lanes
// of the line; `points` and `variants` as glyphs (stacked, with their alternative's letter) and as density;
// `secondary`'s shapes; `pairs`' overlay (a lit arc is the arc drawn, `arcs: "hover"` draws on the overlay only and
// leaves nothing stale), the matrix hover through the real hit, and `logo`'s stacks equal to a brute force.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { describe, expect, it } from 'vitest';
import { Hit, HoverKind, type MatrixState, NO_COLOR, type PairsState, SkinToken } from '../src/draw/api.js';
import { IndexVisits } from '../src/intervals.js';
import { asText, matrixDoc } from './conformance-fixtures.js';
import { drawn, entryOf, layerOf, record, rects, texts } from './draw-fixtures.js';
import { allTypes, rangesDoc, transcripts } from './prepare-fixtures.js';
import { densityDoc } from './rn2-fixtures.js';
import { featuresDoc, glyphs, lit, logoDoc, MarkBit, overlayOf, paths } from './rn2b-fixtures.js';

type Opened = Awaited<ReturnType<typeof drawn>>;
const issues = (t: Opened) => t.ctx.graph.issues.value.filter((i) => i.code.startsWith('seqv.renderer'));
async function zoom(t: Opened, label: string, cell: number, x = 0) {
  t.ctx.viewports.get(layerOf(t.ctx, label).section)?.setView(cell, x, 'api');
  await t.ctx.graph.update();
}
/** Items a layer drew inside its entry's band (device px) in a color that draws: rects, glyphs, texts and paths. */
function inked(t: Opened, label: string, tiny = false) {
  const { rec } = record(t.ctx, label, undefined, tiny);
  const e = entryOf(t.ctx, label);
  const [y0, y1] = [Math.floor(e.top * e.dpr) - 1, Math.ceil((e.top + e.height) * e.dpr) + 1];
  let n = 0;
  for (const [, y, , h, c] of rects(rec)) if (y >= y0 && y + h <= y1 && c !== NO_COLOR) n++;
  for (const [, y, , , c] of glyphs(rec)) if (y >= y0 && y <= y1 && c !== NO_COLOR) n++;
  for (const x of texts(rec)) if (x.y >= y0 && x.y <= y1 && x.color !== NO_COLOR) n++;
  for (let p = 0; p < rec.pathItems.n; p++) if (rec.pathItems.color[p] !== NO_COLOR) n++;
  return n;
}

describe('every RN2 type draws', () => {
  const LABELS = ['Contacts', 'Features', 'Sites', 'Variants', 'Secondary'];
  it.each([1, 2])("at each level of detail and below the row's minimum, at DPR %i", async (dpr) => {
    const t = await drawn(allTypes({ length: 400 }), { width: 900, height: 420, dpr });
    for (const cell of [14, 4, 1.5, 0.3]) {
      await zoom(t, 'Letters', cell);
      for (const label of LABELS)
        for (const tiny of [false, true])
          expect(inked(t, label, tiny), `${label} at ${cell}${tiny ? ' (tiny)' : ''}`).toBeGreaterThan(0);
      expect(issues(t)).toEqual([]);
    }
  });
});

describe('reductions count every item', () => {
  it('a feature and a point narrower than a pixel are drawn', async () => {
    const t = await drawn(rangesDoc([[5000, 5001]], 20000), { width: 900, height: 200 });
    await zoom(t, 'Letters', 0.03);
    for (const label of ['Features', 'Sites']) expect(inked(t, label), label).toBeGreaterThan(0);
  });
});

describe('ranges', () => {
  it("one lane: a rect over each range's positions (one device px apart from the next), labeled and centered", async () => {
    const t = await drawn(featuresDoc(), { width: 900, height: 420 });
    await zoom(t, 'Ranges', 10);
    const { rec, frame } = record(t.ctx, 'Ranges');
    const e = entryOf(t.ctx, 'Ranges');
    const v = frame.v;
    expect(rects(rec).map(([x, y, w, h]) => [x, y, w, h])).toEqual(
      [
        [10, 30],
        [50, 56],
      ].map(([a, b]) => [v.xd(0, a), e.top + 2 + 1, v.xd(0, b) - v.xd(0, a) - 1, 22 - 2 * 2 - 2]),
    );
    const names = texts(rec);
    expect(names.map((d) => d.text)).toEqual(['Alpha domain', 'Gamma']);
    expect(names[0].x).toBe((v.xd(0, 10) + v.xd(0, 30) - 1) / 2);
  });

  it('labels are cut with an ellipsis where the feature is narrower, and left out under three graphemes', async () => {
    const t = await drawn(featuresDoc({ n: 600 }), { width: 900, height: 420 });
    await zoom(t, 'Ranges', 2);
    const names = texts(record(t.ctx, 'Ranges').rec).map((d) => d.text);
    expect(names).toHaveLength(1);
    expect(names[0].endsWith('…') && names[0].length >= 3).toBe(true);
    expect(names[0]).not.toBe('Alpha domain');
  });

  it.each([1, 2])(
    'a label moves with its range while it fits, and sticks 4 px inside the edge that cuts it, at DPR %i',
    async (dpr) => {
      // R0 is [100, 130): 300 px at a cell of 10, in a data area of 734.
      const t = await drawn(rangesDoc([[100, 130]], 600), { width: 900, height: 420, dpr });
      const draw = async (x: number) => {
        await zoom(t, 'Letters', 10, x);
        const { rec, frame } = record(t.ctx, 'Features');
        const v = frame.v;
        const name = texts(rec)[0];
        const [L, R] = [Math.round(v.left * dpr), Math.round(v.right * dpr)];
        const half = (v.text.label('R0', name.font).width * dpr) / 2;
        // Its bar is a device px short of its right end.
        const [a0, a1] = [v.xd(0, 100), v.xd(0, 130) - 1];
        return { x: name.x, bar: rects(rec)[0][0], L, R, half, a0, a1, c: (a0 + a1) / 2 };
      };
      // In view whole: centered on its bar; a pan moves label and bar alike.
      let prev = await draw(700);
      for (const x of [700, 703, 717, 760, 851, 990]) {
        const d = await draw(x);
        expect(d.a0 >= d.L && d.a1 <= d.R).toBe(true);
        expect(d.x).toBe(d.c);
        expect(d.x - prev.x).toBe(d.bar - prev.bar);
        prev = d;
      }
      // Cut by an edge: still centered on the whole range while that keeps it 4 px inside the visible part.
      for (const x of [430, 450, 500, 1050, 1100, 1135]) {
        const d = await draw(x);
        expect(d.a0 < d.L || d.a1 > d.R).toBe(true);
        expect(d.x).toBe(d.c);
      }
      // Past that, stuck 4 px inside the edge, whatever the pan.
      for (const x of [1150, 1200, 1261]) {
        const d = await draw(x);
        expect(d.c).toBeLessThan(d.L + 4 * dpr + d.half);
        expect(d.x).toBe(d.L + 4 * dpr + d.half);
      }
      for (const x of [290, 350, 421]) {
        const d = await draw(x);
        expect(d.c).toBeGreaterThan(d.R - 4 * dpr - d.half);
        expect(d.x).toBe(d.R - 4 * dpr - d.half);
      }
    },
  );

  it("a cut range's label is cut with an ellipsis to its visible part, inside it", async () => {
    const t = await drawn(featuresDoc({ n: 600 }), { width: 900, height: 420 });
    // Alpha domain [10, 30) scrolled to its last 60 px.
    await zoom(t, 'Ranges', 10, 240);
    const { rec, frame } = record(t.ctx, 'Ranges');
    const v = frame.v;
    const name = texts(rec)[0];
    expect(name.text).toBe('Alpha d…');
    const w = v.text.label('Alpha d', name.font).width + v.text.label('…', name.font).width;
    expect(name.x).toBeCloseTo(v.left + 4 + w / 2, 6);
    expect(name.x + w / 2).toBeLessThanOrEqual(v.xd(0, 30) - 1 - 4);
  });

  it("overlapping ranges take lanes: the second under the first, each a lane's rows with its inset", async () => {
    const t = await drawn(featuresDoc(), { width: 900, height: 420 });
    await zoom(t, 'Ranges', 10);
    const { rec, frame } = record(t.ctx, 'Lanes');
    const e = entryOf(t.ctx, 'Lanes');
    const v = frame.v;
    // 22 px, two lanes: 2 px of padding, lanes of 9, each inset by 1.
    expect(rects(rec).map(([x, y, w, h]) => [x, y, w, h])).toEqual([
      [v.xd(0, 10), e.top + 3, v.xd(0, 30) - v.xd(0, 10) - 1, 7],
      [v.xd(0, 50), e.top + 3, v.xd(0, 56) - v.xd(0, 50) - 1, 7],
      [v.xd(0, 20), e.top + 12, v.xd(0, 44) - v.xd(0, 20) - 1, 7],
    ]);
    // The rows the hit finds are these (draw and hit share them): the lane under the point, its feature.
    const hit = Hit.create();
    for (const [x, y, w, h] of rects(rec)) {
      t.ctx.layout.hitTest(layerOf(t.ctx, 'Lanes').section, x + w / 2, y + h / 2, hit);
      expect(hit.kind).toBe(HoverKind.item);
    }
  });

  it("brackets: a line on the lane's foot and a tick at each end", async () => {
    const t = await drawn(featuresDoc(), { width: 900, height: 420 });
    await zoom(t, 'Ranges', 10);
    const { rec, frame } = record(t.ctx, 'Brackets');
    const v = frame.v;
    const e = entryOf(t.ctx, 'Brackets');
    const x0 = v.xd(0, 10);
    const w = v.xd(0, 30) - x0 - 1;
    const y1 = e.top + 19;
    expect(
      rects(rec)
        .map(([x, y, ww, h]) => [x, y, ww, h])
        .slice(0, 3),
    ).toEqual([
      [x0, y1 - 1, w, 1],
      [x0, y1 - 4, 1, 4],
      [x0 + w - 1, y1 - 4, 1, 4],
    ]);
  });

  it("a boundary is a glyph on the edge before its position (a tick), in its theme's color", async () => {
    const t = await drawn(featuresDoc(), { width: 900, height: 420 });
    await zoom(t, 'Ranges', 14);
    const { rec, frame } = record(t.ctx, 'Cuts');
    const g = glyphs(rec);
    expect(g.map(([x]) => x)).toEqual([frame.v.xd(0, 20), frame.v.xd(0, 48)]);
    expect(g.map(([, , r, id]) => [r, id])).toEqual([
      [5, 7],
      [5, 7],
    ]);
    expect(g[0][4]).not.toBe(g[1][4]);
    expect(rects(rec)).toEqual([]);
  });

  it("a hairline joins an item's fragments across their gap, at the lane's middle", async () => {
    const t = await drawn(transcripts(), { width: 900, height: 300 });
    await zoom(t, 'Genes', 3);
    const { rec, frame } = record(t.ctx, 'Genes');
    const v = frame.v;
    const lines = rects(rec).filter(([, , , h]) => h === 1);
    // t1's exons [10, 30), [60, 80), [120, 140): a hairline across each gap; t2's across [50, 90); t3's across [160, 170).
    const gaps = lines.map(([x, , w]) => [x, x + w]).sort((a, b) => a[0] - b[0]);
    expect(gaps).toEqual(
      [
        [30, 60],
        [50, 90],
        [80, 120],
        [160, 170],
      ].map(([a, b]) => [v.xd(0, a), v.xd(0, b)]),
    );
  });

  it('wrapped: underlines on the lanes of the line, compacted to those with a feature on it, cut at group gaps', async () => {
    const t = await drawn(asText(featuresDoc()), { width: 230, height: 900 });
    const { k, section } = layerOf(t.ctx, 'Lanes');
    const { rec, frame } = record(t.ctx, 'Lanes');
    const l = t.ctx.layout.get(section)?.current;
    const track = t.ctx.graph.assembled.value?.view?.resolved.layers.track[k];
    if (!l) throw new Error('no layout');
    const lines = new Map<number, Set<number>>();
    for (const [x, y, w, h] of rects(rec)) {
      let e = 0;
      while (
        e < l.entryCount &&
        !(l.entryTrack[e] === track && y >= l.entryTop[e] && y < l.entryTop[e] + l.entryHeight[e])
      )
        e++;
      expect(e, `the rect at ${x}, ${y} is in an entry of the layer`).toBeLessThan(l.entryCount);
      expect(h).toBe(3);
      expect(w).toBeGreaterThanOrEqual(2);
      lines.set(e, (lines.get(e) ?? new Set()).add(y));
    }
    expect(lines.size).toBeGreaterThan(2);
    for (const [e, ys] of lines) {
      // The entry's height is its lanes × 5 (`wrappedLaneHeight`) and 2 (`wrappedTrackGap`, `lineHeights`): every
      // lane of it has a feature on the line.
      const n = (l.entryHeight[e] - 2) / 5;
      expect([...ys].sort((a, b) => a - b)).toEqual(
        Array.from({ length: n }, (_, j) => Math.round(l.entryTop[e] + l.entryHeight[e] - 5 - 5 * j)).reverse(),
      );
    }
    // A rect never spans a group gap (it hits back to its item): one per group a feature touches. [10, 30) is in two
    // groups of 10 positions, [20, 44) in three, [50, 56) in one.
    expect(rects(rec)).toHaveLength(6);
    void frame;
  });

  it('density: more items in view than device columns draw as one band of columns, fewer as items', async () => {
    const t = await drawn(densityDoc({ type: 'ranges', items: 3000, len: 200 }), { width: 900, height: 300 });
    await zoom(t, 'D', 0.04);
    const out = rects(record(t.ctx, 'D').rec);
    expect(out.length).toBeGreaterThan(100);
    expect(out.every(([, , w]) => w === 1)).toBe(true);
    await zoom(t, 'D', 20);
    const items = rects(record(t.ctx, 'D').rec);
    expect(items.length).toBeGreaterThan(0);
    expect(items.every(([, , w]) => w > 1)).toBe(true);
  });
});

describe('points and variants', () => {
  it("a glyph at the center of each point's cell, in its theme's color", async () => {
    const t = await drawn(featuresDoc(), { width: 900, height: 420 });
    await zoom(t, 'Points', 14);
    const { rec, frame } = record(t.ctx, 'Points');
    const e = entryOf(t.ctx, 'Points');
    const v = frame.v;
    expect(glyphs(rec).map(([x, y, r, id]) => [x, y, r, id])).toEqual(
      [5, 6, 7, 30].map((i) => [(v.xd(0, i) + v.xe(0, i)) / 2, Math.round(e.top + 11), (e.height - 2) / 3.5, 0]),
    );
    const colors = glyphs(rec).map((g) => g[4]);
    expect(colors[0]).toBe(colors[2]);
    expect(colors[0]).not.toBe(colors[1]);
  });

  it("variants stack by their place in their position, the alternative's letter on each mark", async () => {
    const t = await drawn(featuresDoc(), { width: 900, height: 420 });
    await zoom(t, 'Variants', 14);
    const { rec, frame } = record(t.ctx, 'Variants');
    const e = entryOf(t.ctx, 'Variants');
    const g = glyphs(rec);
    // The band (15.9 px) holds a radius of 6.3 (a variant's mark grows to 8 px with its cell); a row of 13.6 px: [12]
    // twice (a second row down), [40] once.
    expect(g.map(([x, y, r]) => [x, Math.round((y - Math.round(e.top + 2 + 6.3)) * 10) / 10, r])).toEqual([
      [frame.v.xd(0, 12) + 7, 0, 6.3],
      [frame.v.xd(0, 12) + 7, 13.6, 6.3],
      [frame.v.xd(0, 40) + 7, 0, 6.3],
    ]);
    expect(texts(rec).map((d) => d.text)).toEqual(['A', 'C', 'D']);
    // On the alphabetic baseline that centers a capital's ink on its mark (ui 600 at 8 px, DPR 1).
    const cap = t.ctx.text.ink(t.ctx.text.font('ui', 8, 600), 64)[0];
    expect(texts(rec).map((d) => d.y)).toEqual(g.map(([, y]) => Math.round(y + cap / 2)));
  });

  it('past 64 items in view, thinner than two device px a cell: density, one band, no glyphs; zoomed in: glyphs', async () => {
    const t = await drawn(densityDoc({ items: 3000 }), { width: 900, height: 300 });
    await zoom(t, 'D', 0.04);
    const dense = record(t.ctx, 'D').rec;
    expect(glyphs(dense)).toEqual([]);
    expect(dense.rectItems.n).toBeGreaterThan(100);
    await zoom(t, 'D', 14);
    const near = record(t.ctx, 'D').rec;
    expect(near.rectItems.n).toBe(0);
    expect(glyphs(near).length).toBeGreaterThan(0);
  });

  it('the tiny path is density always', async () => {
    const t = await drawn(featuresDoc(), { width: 900, height: 420 });
    await zoom(t, 'Points', 14);
    for (const label of ['Points', 'Variants']) {
      const { rec } = record(t.ctx, label, undefined, true);
      expect(glyphs(rec), label).toEqual([]);
      expect(rects(rec).length, label).toBeGreaterThan(0);
    }
  });
});

describe('secondary', () => {
  it('a coil line, a helix cylinder, a strand arrow, a turn arch and a bar for what is neither, in their rows', async () => {
    const t = await drawn(featuresDoc(), { width: 900, height: 420 });
    await zoom(t, 'Secondary', 14);
    const { rec, frame } = record(t.ctx, 'Secondary');
    const e = entryOf(t.ctx, 'Secondary');
    const v = frame.v;
    const cy = Math.round(e.top + 11);
    const h = Math.round(22 * 0.31);
    const [coil, other] = rects(rec);
    expect(coil).toEqual([
      v.xd(0, v.visLo[0]),
      cy,
      v.xd(0, v.visHi[0]) - v.xd(0, v.visLo[0]),
      1,
      SkinToken.inkFaint,
      255,
    ]);
    expect(other[1]).toBe(cy - Math.round(h * 0.45));
    const [helix, strand, turn] = paths(rec);
    // The helix: a cylinder, its body's top edge between its rounded ends (seqspec:format#layer-types).
    const rx = Math.min(Math.max(1, Math.round(h * 0.6)), (v.xd(0, 12) - v.xd(0, 4)) >> 1);
    expect(helix.slice(4, 8)).toEqual([v.xd(0, 4) + rx, cy - h, v.xd(0, 12) - rx, cy - h]);
    // The strand: seven corners and a close; the arrow's tip at its end, on the coil line.
    expect(strand.slice(4, 6)).toEqual([v.xd(0, 14), cy - Math.round(h * 0.55)]);
    expect(strand.slice(4 + 6, 4 + 8)).toEqual([v.xd(0, 24), cy]);
    // The turn: an arch over the coil line, from its start to its end.
    expect([turn[4], turn[5], turn[turn.length - 2]]).toEqual([v.xd(0, 26), cy, v.xd(0, 26)]);
  });

  it('narrower than three device px, and on the tiny path, a feature is a bar', async () => {
    const t = await drawn(featuresDoc({ n: 6000 }), { width: 900, height: 420 });
    await zoom(t, 'Secondary', 0.12);
    expect(record(t.ctx, 'Secondary').rec.pathItems.n).toBe(0);
    expect(rects(record(t.ctx, 'Secondary').rec).length).toBeGreaterThan(4);
    await zoom(t, 'Secondary', 14);
    expect(record(t.ctx, 'Secondary').rec.pathItems.n).toBe(3);
    expect(record(t.ctx, 'Secondary', undefined, true).rec.pathItems.n).toBe(0);
    expect(rects(record(t.ctx, 'Secondary', undefined, true).rec).length).toBeGreaterThan(4);
  });
});

describe('pairs: the overlay', () => {
  /** The arcs among a recording's paths (strokes; the wedges at the edges are fills). */
  const strokes = (ps: number[][]) => ps.filter((p) => p[0] === 1).length;
  const stateOf = (t: Opened) => {
    const { k } = layerOf(t.ctx, 'Contacts');
    return t.ctx.graph.assembled.value?.view?.prepared[k].state as PairsState;
  };
  /** The first end of block 0 within the first 60 positions: its pair and run. */
  const endOf = (s: PairsState) => {
    let e = 0;
    while (s.block[e] !== 0 || s.start[e] > 60) e++;
    return { k: e >> 1, run: [s.start[e], s.stop[e]] as [number, number] };
  };

  it.each([1, 2])('a lit arc is the arc drawn, stronger (2.5 px, half device px), at DPR %i', async (dpr) => {
    const t = await drawn(allTypes({ length: 160 }), { width: 1200, height: 420, dpr });
    await zoom(t, 'Letters', 6);
    const main = paths(record(t.ctx, 'Contacts').rec).filter((p) => p[0] === 1);
    const { run } = endOf(stateOf(t));
    const over = paths(overlayOf(t.ctx, 'Contacts', lit(MarkBit.selection, [run])));
    expect(over.length).toBeGreaterThan(0);
    for (const arc of over) {
      expect(arc[1]).toBe(5 * dpr);
      expect(
        main.map((p) => p.slice(2).join()),
        'the same curve',
      ).toContain(arc.slice(2).join());
    }
    expect(main.every((p) => p[1] === 2 * dpr)).toBe(true);
  });

  it('every channel lights: the hover, a highlight, the selection, the accent', async () => {
    const t = await drawn(allTypes({ length: 160, arcs: 'hover' }), { width: 1200, height: 420 });
    await zoom(t, 'Letters', 6);
    const { run } = endOf(stateOf(t));
    for (const bit of [MarkBit.hover, MarkBit.highlight, MarkBit.selection, MarkBit.accent])
      expect(paths(overlayOf(t.ctx, 'Contacts', lit(bit, [run]))).length, `bit ${bit}`).toBeGreaterThan(0);
    expect(paths(overlayOf(t.ctx, 'Contacts', lit(MarkBit.hover, []))).length).toBe(0);
  });

  it('`arcs: "hover"`: no arcs on the main pass, the lit ones on the overlay (1.5 px), drawn and hit only while lit', async () => {
    const t = await drawn(allTypes({ length: 160, arcs: 'hover' }), { width: 1200, height: 420 });
    await zoom(t, 'Letters', 6);
    const s = stateOf(t);
    expect(paths(record(t.ctx, 'Contacts').rec).filter((p) => p[0] === 1)).toEqual([]);
    const { k, run } = endOf(s);
    const over = paths(overlayOf(t.ctx, 'Contacts', lit(MarkBit.hover, [run])));
    expect(over.length).toBeGreaterThan(0);
    expect(over.every((p) => p[1] === 3)).toBe(true);
    // Recorded for the hit: exactly the pairs drawn.
    const drawnPairs = [...s.drawn].flatMap((d, i) => (d ? [i] : []));
    expect(drawnPairs).toContain(k);
    expect(drawnPairs).toHaveLength(strokes(over));
    // Another end lit: the first's mark is gone, never stale.
    const other = overlayOf(t.ctx, 'Contacts', lit(MarkBit.hover, [[s.start[2 * 3 + 1], s.stop[2 * 3 + 1]]]));
    expect(s.drawn[3]).toBe(1);
    expect(strokes(paths(other))).toBe([...s.drawn].filter((d) => d).length);
    // Nothing lit: nothing drawn, nothing recorded.
    expect(paths(overlayOf(t.ctx, 'Contacts', lit(MarkBit.hover, []))).length).toBe(0);
    expect(s.drawn.some((d) => d)).toBe(false);
  });

  it('an item hover of the layer lights its pair alone, on top, not the pairs that share its ends', async () => {
    const { pairsDoc } = await import('./rn2-fixtures.js');
    const t = await drawn(pairsDoc({ count: 12_000, long: false }), { width: 1200, height: 300 });
    await zoom(t, 'Letters', 3, 3 * 15_000);
    const s = stateOf(t);
    const { k: layer } = layerOf(t.ctx, 'Contacts');
    // A pair in view with another pair at one of its ends.
    const at = new Map<number, number>();
    for (let e = 0; e < 2 * s.count; e++) at.set(s.start[e], (at.get(s.start[e]) ?? 0) + 1);
    let k = 0;
    while (!(s.start[2 * k] > 15_050 && s.start[2 * k + 1] < 15_300) || (at.get(s.start[2 * k]) ?? 0) < 2) k++;
    const ends: [number, number][] = [0, 1].map((side) => [s.start[2 * k + side], s.stop[2 * k + side]]);
    const marks = lit(
      MarkBit.hover,
      ends.sort((a, b) => a[0] - b[0]),
    ) as ReturnType<typeof lit> & {
      hoverLayer: number;
      hoverItem: number;
    };
    // The hover's runs alone light every pair at either end.
    expect(strokes(paths(overlayOf(t.ctx, 'Contacts', marks)))).toBeGreaterThan(1);
    marks.hoverLayer = layer;
    marks.hoverItem = s.set.annotations[s.annotation[k]]?.item[s.row[k]] ?? -1;
    const one = paths(overlayOf(t.ctx, 'Contacts', marks)).filter((p) => p[0] === 1);
    expect(one).toHaveLength(1);
    expect([...s.drawn].flatMap((d, i) => (d ? [i] : []))).toEqual([k]);
    // Another layer's item: the runs again.
    marks.hoverLayer = layer + 1;
    expect(strokes(paths(overlayOf(t.ctx, 'Contacts', marks)))).toBeGreaterThan(1);
  });

  it("a layer drawn from two annotations: the hovered pair is its annotation's, though the other has an item of that id", async () => {
    const b = createSeqViewSpec({ title: 'Two pair sets' });
    const a = b.data.sequence({ ref: 'a', label: 'A', value: 'ACDEFGHIKLMNPQRSTVWY'.repeat(10) });
    // One pair each: both are item 0 of their annotation.
    const sets = [
      [20, 60],
      [40, 90],
    ].map(([from, to], k) =>
      b.data.annotation({
        ref: `pp${k}`,
        geometry: 'pairs',
        target: a,
        value: { start: [from], b_target: ['a'], b_start: [to] },
      } as never),
    );
    const v = b.view({ ref: 'v', label: 'V', height: 'box' }).block({ sequence: a });
    v.track({ label: 'Letters' }).letters();
    v.track({ label: 'Contacts' }).pairs({ data: sets, arcs: 'always' } as never);
    const t = await drawn(b.getDocument(), { width: 1200, height: 300 });
    await zoom(t, 'Letters', 5);
    const s = stateOf(t);
    const { k: layer } = layerOf(t.ctx, 'Contacts');
    expect(s.count).toBe(2);
    const ids = [0, 1].map((k) => s.set.annotations[s.annotation[k]]?.item[s.row[k]]);
    expect(ids).toEqual([0, 0]);
    for (const k of [0, 1]) {
      const ends: [number, number][] = [0, 1].map((side) => [s.start[2 * k + side], s.stop[2 * k + side]]);
      const marks = lit(MarkBit.hover, ends) as ReturnType<typeof lit> & {
        hoverLayer: number;
        hoverAnnotation: number;
        hoverItem: number;
      };
      Object.assign(marks, { hoverLayer: layer, hoverAnnotation: s.annotation[k], hoverItem: 0 });
      overlayOf(t.ctx, 'Contacts', marks);
      expect(
        [...s.drawn].flatMap((d, i) => (d ? [i] : [])),
        `pair ${k}`,
      ).toEqual([k]);
    }
  });

  it('finds the pairs through the end index: a thousand pairs, a few visits for one lit end', async () => {
    const { pairsDoc } = await import('./rn2-fixtures.js');
    const t = await drawn(pairsDoc({ count: 1000, arcs: 'hover' }), { width: 1200, height: 300 });
    await zoom(t, 'Letters', 0.05);
    const before = IndexVisits.ends;
    const rec = overlayOf(t.ctx, 'Contacts', lit(MarkBit.hover, [[5000, 5001]]));
    expect(IndexVisits.ends - before).toBeLessThan(100);
    expect(paths(rec).length).toBeLessThan(10);
  });
});

describe('pairs: ends', () => {
  it('a theme that gives no color: the ends and their density in the same quiet ink', async () => {
    const { pairsDoc } = await import('./rn2-fixtures.js');
    const t = await drawn(pairsDoc({ count: 2000, colorBy: 'score' }), { width: 900, height: 300 });
    const ink = (cell: number, tiny = false) => {
      t.ctx.viewports.get(layerOf(t.ctx, 'Contacts').section)?.setView(cell, 0, 'api');
      return t.ctx.graph.update().then(() => {
        const { rec } = record(t.ctx, 'Contacts', undefined, tiny);
        return [...new Set([...rects(rec).map((r) => r[4]), ...glyphs(rec).map((g) => g[4])])];
      });
    };
    expect(await ink(10)).toEqual([SkinToken.inkQuiet]);
    expect(await ink(0.05)).toEqual([SkinToken.inkQuiet]);
    expect(await ink(10, true)).toEqual([SkinToken.inkQuiet]);
  });
});

describe('matrix: the hover through the real hit', () => {
  it('the band under the pointer names the state of the image row drawn there', async () => {
    const t = await drawn(matrixDoc(), { width: 900, height: 300 });
    await zoom(t, 'M', 8);
    const { rec, frame } = record(t.ctx, 'M');
    const { k, section } = layerOf(t.ctx, 'M');
    const e = entryOf(t.ctx, 'M');
    const state = t.ctx.graph.assembled.value?.view?.prepared[k].state as MatrixState;
    const theme = t.ctx.graph.assembled.value?.view?.prepared[k].themes[0];
    const hit = Hit.create();
    const i = 37;
    const x = (frame.v.xd(0, i) + frame.v.xe(0, i)) / 2;
    const m = state.ann[0]?.matrix;
    if (!m || !theme) throw new Error('no matrix');
    // The tile holding position i, its pixels as the backend kept them: a row per state over the 100 samples.
    const img = rec.images.content.find((c) => c instanceof Uint32Array && c.length === 300) as Uint32Array | undefined;
    expect(img).toBeDefined();
    for (let band = 0; band < 3; band++) {
      t.ctx.layout.hitTest(section, x, e.top + ((band + 0.5) * e.height) / 3, hit);
      expect([hit.kind, hit.state]).toEqual([HoverKind.seq, band]);
      const loc = {
        ...frame.v.loc,
        value: m.values[(state.rows[0]?.[i] ?? 0) * 3 + hit.state],
        state: hit.state,
        index: i,
      };
      const c = theme.value(loc as never);
      const px = ((255 << 24) | ((c & 255) << 16) | (c & 0xff00) | ((c >> 16) & 255)) >>> 0;
      expect((img as Uint32Array)[hit.state * 100 + i]).toBe(px);
    }
  });
});

describe('logo', () => {
  /** The old walk: a letter's offset is the heights below it in its column's order, found letter by letter. */
  function expected(s: MatrixState, i: number, band: number) {
    const S = s.states;
    const h = s.heights[0] as Float32Array;
    const o = s.order[0] as Uint8Array;
    return Array.from({ length: S }, (_, q) => {
      let below = 0;
      let rank = 0;
      while (rank < S && o[i * S + rank] !== q) below += h[i * S + o[i * S + rank++]];
      return { q, below: below * band, hq: h[i * S + q] * band };
    });
  }

  it("each letter stands on the ones below it, stretched over its height and its cell (the old stack's look)", async () => {
    for (const dpr of [1, 2]) {
      const t = await drawn(logoDoc(), { width: 1200, height: 300, dpr });
      const { k } = layerOf(t.ctx, 'Logo');
      const s = t.ctx.graph.assembled.value?.view?.prepared[k].state as MatrixState;
      const { rec, frame } = record(t.ctx, 'Logo');
      const e = entryOf(t.ctx, 'Logo');
      const bottom = Math.round((e.top + e.height - 2) * dpr);
      const band = (e.height - 4) * dpr;
      expect(rects(rec), 'no bars at letter widths').toEqual([]);
      const b = rec.textItems;
      const font = t.ctx.text.font('mono', 20, 600);
      const adv = t.ctx.text.advances(font);
      let j = 0;
      for (let c = 0; c < 30; c++) {
        const x0 = frame.v.xd(0, c);
        const w = frame.v.xe(0, c) - x0;
        const lw = Math.max(1, w - dpr);
        for (const { q, below, hq } of expected(s, c, band).sort((x, y) => x.below - y.below)) {
          if (hq < dpr) continue;
          const at = `column ${c}, state ${q}, DPR ${dpr}`;
          expect(b.text[j], at).toBe(String.fromCharCode(s.codes[q]));
          expect(b.y[j], at).toBeCloseTo(bottom - below, 6);
          expect(b.x[j], at).toBeCloseTo(x0 + (w - lw) / 2, 6);
          // Stretched: the advance to the cell less a CSS px, the capitals (0.72 of the size) to the letter's height.
          expect(b.sx[j] * adv[s.codes[q]] * dpr, at).toBeCloseTo(lw, 6);
          expect(b.sy[j] * 20 * 0.72 * dpr, at).toBeCloseTo(hq, 6);
          j++;
        }
      }
      expect(j).toBe(b.n);
      expect(j).toBeGreaterThan(50);
    }
  });

  it("thin: one bar per device column, the tallest stack in its top letter's color", async () => {
    const t = await drawn(logoDoc(), { width: 215, height: 300 });
    const { k } = layerOf(t.ctx, 'Logo');
    const s = t.ctx.graph.assembled.value?.view?.prepared[k].state as MatrixState;
    const { rec, frame } = record(t.ctx, 'Logo');
    expect(frame.v.cell).toBeLessThan(2);
    const e = entryOf(t.ctx, 'Logo');
    const out = rects(rec);
    expect(out.length).toBeGreaterThan(10);
    expect(out.length).toBeLessThanOrEqual(30);
    const bottom = Math.round(e.top + e.height - 2);
    const h = s.heights[0] as Float32Array;
    for (const [x, y, , hh] of out) {
      const c = Array.from({ length: 30 }, (_, i) => i).find((i) => frame.v.xd(0, i) === x);
      expect(c).toBeDefined();
      let total = 0;
      for (let q = 0; q < s.states; q++) total += h[(c as number) * s.states + q];
      expect(y + hh).toBe(bottom);
      expect(Math.abs(hh - total * (e.height - 4))).toBeLessThanOrEqual(1);
    }
    expect(texts(rec)).toEqual([]);
  });

  it("thin, several positions a device column: the tallest of each column's positions, from the pyramid or raw", async () => {
    const t = await drawn(logoDoc(12, 8000), { width: 900, height: 300 });
    const { k, section } = layerOf(t.ctx, 'Logo');
    const s = t.ctx.graph.assembled.value?.view?.prepared[k].state as MatrixState;
    const h = s.heights[0] as Float32Array;
    const total = (i: number) => {
      let x = 0;
      for (let q = 0; q < s.states; q++) x += h[i * s.states + q];
      return x;
    };
    // Under 4 positions a column, raw; the fit (about 11 a column) from the pyramid's bins of 8, each in the column
    // of its start.
    const bins: number[] = [];
    for (const cell of [0.6, 0.01]) {
      t.ctx.viewports.get(section)?.setView(cell, 0, 'api');
      await t.ctx.graph.update();
      const { rec, frame } = record(t.ctx, 'Logo');
      const e = entryOf(t.ctx, 'Logo');
      const v = frame.v;
      const level = Math.floor(Math.log2(1 / v.cell));
      const bin = level >= 2 ? 2 ** level : 1;
      bins.push(bin);
      const out = rects(rec);
      expect(out.length).toBeGreaterThan(100);
      for (const [x, y, w, hh] of out) {
        let best = -1;
        for (let i = 0; i < 8000; i++)
          if (v.xd(0, i - (i % bin)) >= x && v.xd(0, i - (i % bin)) < x + w && (best < 0 || total(i) > total(best)))
            best = i;
        expect(Math.abs(hh - total(best) * (e.height - 4)), `column at ${x}, cell ${cell}`).toBeLessThanOrEqual(1);
        expect(y + hh).toBe(Math.round(e.top + e.height - 2));
      }
    }
    expect(bins).toEqual([1, 8]);
  });
});
