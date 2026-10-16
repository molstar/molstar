/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// An alignment's rows as text (seqspec:viewer#rows-in-text): a text row per row on every line, each row's label left
// of it and its last residue number on the line right of it (CLUSTAL W's own counts at 60 a line, in the row's
// numbering), the group labels
// by the reference row, the theme's cells under the letters in the ink that reads on them, a box on one row's range,
// hits on a row's residue; the `espript` theme's two levels.

import { SeqSets } from '@molstar/seq-core/model';
import { describe, expect, it } from 'vitest';
import { Dirty } from '../src/draw/index.js';
import { NO_COLOR } from '../src/draw/backend.js';
import { HoverKind, Hit, type RenderProps } from '../src/index.js';
import { CellInk } from '../src/themes.js';
import { renderSvg } from '../src/svg/index.js';
import { recording } from './composer-fixtures.js';
import { lit } from './overlay-fixtures.js';
import { globinRows, figureDoc } from './text-rows-fixtures.js';
import { at, context } from './theme-fixtures.js';

/** A line per 60 columns, as the fixture's own blocks: one group of 60 a line (no room for two), no group labels. */
async function sixty(): Promise<Partial<RenderProps>> {
  const t = await lit(figureDoc(), { width: 1100, height: 1400 });
  const cell = t.ctx.layout.axis(t.section)?.axis?.fixedCell ?? 0;
  t.ctx.dispose();
  return { wrappedGroup: 60, wrappedMaxWidth: Math.ceil(60 * cell) + 2, wrappedGroupLabels: 'none' };
}
const texts = (r: ReturnType<typeof recording>) =>
  Array.from({ length: r.textItems.n }, (_, i) => ({
    text: r.textItems.text[i],
    x: r.textItems.x[i],
    y: r.textItems.y[i],
    color: r.textItems.color[i],
  }));

async function figure(given?: Partial<RenderProps>) {
  const render = given ?? (await sixty());
  const t = await lit(figureDoc(), { width: 1100, height: 1400, render });
  const l = t.ctx.layout.get(t.section)?.current;
  if (!l) throw new Error('no layout');
  return { ...t, l, tr: l.textRows };
}

describe('the rows as text', () => {
  it('a text row per row on each line; labels left, the last residue number of each row on the line right', async () => {
    const { l, tr, frame } = await figure();
    const rows = globinRows();
    expect(tr.layer).toBeGreaterThanOrEqual(0);
    expect(tr.count).toBe(rows.length);
    expect(l.lines.count).toBe(2);
    expect([l.lines.stop[0] - l.lines.start[0], l.lines.stop[1] - l.lines.start[1]]).toEqual([60, 58]);
    // The rows' band on each line: a text row per row.
    for (let e = 0; e < l.entryCount; e++)
      if (l.entryTrack[e] === tr.track) expect(l.entryHeight[e]).toBe(rows.length * tr.height);
    const drawn = texts(frame(recording(), Dirty.all));
    for (const r of rows) expect(drawn.filter((d) => d.text === r.name).length, r.name).toBe(2);
    // Right of the lines: each row's number, in row order, on each line; CLUSTAL W printed the same.
    const right = drawn.filter((d) => d.x > l.dataLeft && /^\d+$/.test(d.text));
    for (let k = 0; k < 2; k++) {
      const y0 = l.scrolledTop + l.lines.top[k];
      const y1 = l.scrolledTop + l.lines.bottom[k];
      const line = right.filter((d) => d.y >= y0 && d.y < y1).sort((a, b) => a.y - b.y);
      // CLUSTAL W prints a row's residues so far; the row's numbering (`NAME/27-137`) starts at its first.
      const first = (r: (typeof rows)[number]) => Number(/\/(\d+)-/.exec(r.name)?.[1] ?? 1);
      expect(line.map((d) => Number(d.text))).toEqual(rows.map((r) => first(r) + r.ends[k] - 1));
    }
    // Each label left of its line, in the room the layout keeps.
    const label = drawn.find((d) => d.text === rows[0].name);
    expect(label?.x).toBe(Math.round(l.dataLeft - tr.labelRoom));
  });

  it('the group labels number the reference row: its last residue in each group, none where it has none', async () => {
    const { l, frame } = await figure({ wrappedGroup: 10 });
    const ref = globinRows()[0];
    const first = 27;
    const drawn = texts(frame(recording(), Dirty.all)).filter((d) => d.y < l.scrolledTop + l.lines.top[0] + 14);
    const labels = drawn.filter((d) => /^\d+$/.test(d.text)).map((d) => Number(d.text));
    const expected: number[] = [];
    let residues = 0;
    for (let g = 0; g * 10 < l.lines.stop[0]; g++) {
      let last = -1;
      for (let c = g * 10; c < g * 10 + 10; c++)
        if (ref.gapped[c] !== '-' && ref.gapped[c] !== '.') last = first + residues++;
      if (last >= 0) expected.push(last);
    }
    expect(labels).toEqual(expected);
  });

  it('on cells (`text_cells`): an identical column black, its letters light; a gap quiet, on no cell', async () => {
    const { frame, l, tr } = await figure();
    const r = frame(recording(), Dirty.all);
    const top = l.scrolledTop + l.entryTop[l.trackFirst[tr.track]];
    const black: number[][] = [];
    for (let i = 0; i < r.rectItems.n; i++)
      if (r.rectItems.color[i] === 0x000000 && r.rectItems.y[i] >= top && r.rectItems.y[i] < top + tr.count * tr.height)
        black.push([r.rectItems.x[i], r.rectItems.y[i], r.rectItems.w[i], r.rectItems.h[i]]);
    expect(black.length).toBeGreaterThan(8);
    const drawn = texts(r);
    for (const [x, y, w, h] of black) {
      const letter = drawn.find((d) => d.x === x + w / 2 && d.y >= y && d.y <= y + h);
      expect(letter?.color).toBe(CellInk.light);
    }
  });

  it("a box on row 2's range: drawn on row 2's text row, in its color", async () => {
    const { frame, l, tr } = await figure();
    const r = frame(recording(), Dirty.all);
    const top = l.scrolledTop + l.entryTop[l.trackFirst[tr.track]];
    const ys: number[] = [];
    for (let i = 0; i < r.rectItems.n; i++) if (r.rectItems.color[i] === 0xe3120b) ys.push(r.rectItems.y[i]);
    for (let i = 0; i < r.pathItems.n; i++)
      if (r.pathItems.color[i] === 0xe3120b) ys.push(r.pathItems.coords[r.pathItems.cstart[i] + 1]);
    expect(ys.length).toBeGreaterThan(0);
    for (const y of ys) {
      expect(y).toBeGreaterThanOrEqual(top + 2 * tr.height - 1);
      expect(y).toBeLessThanOrEqual(top + 3 * tr.height + 1);
    }
  });

  it("hits: a cell is its row's residue; a gap the row at its column", async () => {
    const { ctx, section, l, tr } = await figure();
    const rows = globinRows();
    const top = l.scrolledTop + l.entryTop[l.trackFirst[tr.track]];
    const h = Hit.create();
    // Row 3 (BRUPA), column 30: a residue; row 1 (9PEZI) column 22: a gap.
    const cellX = (c: number) => {
      const v = ctx.layout.axis(section);
      if (!v) throw new Error('no axis');
      // The first line's x of column c (a mono cell, groups of 60: no gap before 60).
      return l.dataLeft + l.lines.left[0] + (c + 0.5) * (v.axis?.fixedCell ?? 0);
    };
    expect(ctx.layout.hitTest(section, cellX(30), top + 3.5 * tr.height, h)).toBe(true);
    expect([h.kind, h.dataRow, h.index]).toEqual([HoverKind.aligned, 3, 30]);
    const pos = rows[3].gapped.slice(0, 30).replace(/[-.]/g, '').length;
    expect(h.row.index).toBe(pos);
    expect(ctx.layout.hitTest(section, cellX(22), top + 1.5 * tr.height, h)).toBe(true);
    expect(rows[1].gapped[22]).toBe('-');
    expect([h.kind, h.dataRow, h.index]).toEqual([HoverKind.aligned, 1, 22]);
  });

  it('exports alike: SVG holds the labels, numbers and frames, nothing undefined', async () => {
    const xml = await renderSvg(figureDoc(), { width: 1000 });
    for (const s of ['HBA_HUMAN/27-137', 'N-segment', 'α1', '>137<']) expect(xml, s).toContain(s);
    expect(xml).not.toMatch(/NaN|Infinity|undefined/);
  });
});

describe('the espript theme', () => {
  // Four rows: column 0 identical (A); column 1 K, K, R, E: KRH 75% (similar for K and R, none for E); column 2 I, L,
  // V, D: ILMV 75%; column 3 a gap and three W (identical among residues).
  const set = SeqSets.build({
    molecules: [{ key: 'p', residues: 'AKIW' }],
    units: [{ key: 'A', molecule: 'p' }],
    alignments: [{ key: 'e', gapped: ['AKI-', 'AKLW', 'ARVW', 'AEDW'], flavor: 'aligned', alphabet: 'protein' }],
  }).set;
  const code = (ch: string) => ch.charCodeAt(0);
  it('identical black, similar grey, others none', () => {
    const c = context(set, { blocks: [{ alignment: 'e' }] }).color('espript');
    const of = (column: number, ch: string) => c.value(at({ column, letter: code(ch) }));
    expect([of(0, 'A'), of(1, 'K'), of(1, 'R'), of(1, 'E'), of(2, 'I'), of(2, 'D'), of(3, 'W')]).toEqual([
      0x000000,
      0xbdbdbd,
      0xbdbdbd,
      NO_COLOR,
      0xbdbdbd,
      NO_COLOR,
      0x000000,
    ]);
    const strict = context(set, { blocks: [{ alignment: 'e' }] }).color('espript', { similarity: 0.8 });
    expect(strict.value(at({ column: 1, letter: code('K') }))).toBe(NO_COLOR);
  });
});
