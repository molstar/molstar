/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// RN1's draws on the recording backend (seqspec:viewer#extension-points; seqspec:viewer#the-render-backend,
// seqspec:sources#levels-and-the-zoom-gate): every RN1 type draws at each level of detail and below the row's minimum;
// case 14's residues over their codons, one letter or three by the room and `residue_form`, His12 over its larger part;
// reductions count every item; `color` layers with and without values; wrapped group labels above each group's end,
// inline before its start, or none; the ruler over several blocks side by side (names at their starts, no overlap,
// round numbers by jumps equal to those found position by position); values' light min–max band where columns
// aggregate; values colored by their values (seqspec:format#layer-types): a line each position's color over its
// cell, changing halfway between two positions, an area per column, one color with no theme; in SVG, and on rings.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { describe, expect, it } from 'vitest';
import { TextAlign } from '../src/draw/backend.js';
import type { PerPosition } from '../src/prepare.js';
import { renderSvg } from '../src/svg/index.js';
import { columns, drawn, entryOf, layerOf, record, recordLabels, rects, texts } from './draw-fixtures.js';
import { allTypes, BAND_SCORES, colorDoc, gene, linesDoc, rangesDoc, scaleDoc } from './prepare-fixtures.js';

const issues = (t: Awaited<ReturnType<typeof drawn>>) =>
  t.ctx.graph.issues.value.filter((i) => i.code.startsWith('seqv.renderer'));
async function zoom(t: Awaited<ReturnType<typeof drawn>>, section: string, cell: number, x = 0) {
  t.ctx.viewports.get(section)?.setView(cell, x, 'api');
  await t.ctx.graph.update();
}
/** Items a layer drew inside its entry's band (device px). */
function inked(t: Awaited<ReturnType<typeof drawn>>, label: string, tiny = false) {
  const { rec } = record(t.ctx, label, undefined, tiny);
  const e = entryOf(t.ctx, label);
  const [y0, y1] = [Math.floor(e.top * e.dpr), Math.ceil((e.top + e.height) * e.dpr)];
  let n = 0;
  const r = rec.rectItems;
  for (let i = 0; i < r.n; i++) if (r.y[i] >= y0 - 1 && r.y[i] + r.h[i] <= y1 + 1 && r.w[i] > 0) n++;
  const x = rec.textItems;
  for (let i = 0; i < x.n; i++) if (x.y[i] >= y0 && x.y[i] <= y1) n++;
  n += rec.pathItems.n;
  return n;
}

describe('every RN1 type draws', () => {
  const LABELS = ['', 'Letters', 'Emphasis', 'Values', 'Color'];
  it.each([1, 2])("at each level of detail and below the row's minimum, at DPR %i", async (dpr) => {
    const t = await drawn(allTypes({ length: 400 }), { width: 900, height: 420, dpr });
    const { section } = layerOf(t.ctx, 'Letters');
    for (const cell of [14, 4, 1.5, 0.3]) {
      await zoom(t, section, cell);
      for (const label of LABELS)
        for (const tiny of [false, true])
          expect(inked(t, label, tiny), `${label || 'ruler'} at ${cell}${tiny ? ' (tiny)' : ''}`).toBeGreaterThan(0);
      expect(issues(t)).toEqual([]);
    }
  });
});

describe('case 14: residues over their codons', () => {
  it('one letter below a three-letter fit, three above it; `residue_form` forces either; His12 over its larger part', async () => {
    for (const [form, cell, expected] of [
      // Sized as the letters are (the row's one font from its cells): three letters until a codon's three cells no
      // longer hold them at the smallest size, then one.
      ['auto', 4, 'H'],
      ['auto', 16, 'His'],
      ['three', 8, 'His'],
      ['one', 16, 'H'],
    ] as const) {
      const t = await drawn(gene(form), { width: 600, height: 200 });
      const { section } = layerOf(t.ctx, 'Protein');
      // Bases 60 on in view: the intron's end and His12's larger part.
      await zoom(t, section, cell, cell * 60);
      const { rec, frame } = record(t.ctx, 'Protein');
      const drawnTexts = texts(rec).map((d) => d.text);
      expect(drawnTexts, `${form} at ${cell}`).toContain(expected);
      expect(drawnTexts, `${form} at ${cell}`).not.toContain(expected === 'H' ? 'His' : 'H');
      if (expected === 'His') {
        // His12 (bases 45, then 76 and 77): over [76, 78), centered there.
        const his = texts(rec).find((d) => d.text === 'His');
        expect(his?.x).toBeCloseTo((frame.v.xd(0, 76) + frame.v.xd(0, 78)) / 2, 5);
      }
      t.ctx.dispose();
    }
  });
});

describe("case 14: a three-letter code keeps inside its residue's drawn part", () => {
  /** MARGIN CSS px at least each side of the part, whatever the zoom: a codon whole, or split by the intron. */
  it('never reaches the part it is centered over, at every zoom, whole and split codons; else one letter', async () => {
    const MARGIN = 2;
    let threes = 0;
    let ones = 0;
    let his = 0;
    for (const dpr of [1, 2]) {
      const t = await drawn(gene('auto'), { width: 600, height: 200, dpr });
      const { section } = layerOf(t.ctx, 'Protein');
      for (let cell = 4; cell <= 24; cell += 0.25) {
        // Bases 60 on in view: the intron's end and His12's larger part [76, 78), its other base at 45.
        await zoom(t, section, cell, cell * 60);
        const { rec, frame } = record(t.ctx, 'Protein');
        const v = frame.v;
        const drawnTexts = texts(rec).filter((d) => /^[A-Z][a-z]{0,2}$/.test(d.text));
        // The row's one form.
        expect(new Set(drawnTexts.map((d) => d.text.length)).size, `cell ${cell}`).toBeLessThanOrEqual(1);
        for (const d of drawnTexts) {
          if (d.text.length === 1) {
            ones++;
            continue;
          }
          threes++;
          if (d.text === 'His') his++;
          // The text's own box: its width, centered at x; the codon part it is centered on (3 bases, His12's 2).
          const half = (3 * v.text.cell(d.font, null) * dpr) / 2;
          let ok = false;
          for (const [a, b] of [
            [76, 78],
            ...Array.from({ length: 11 }, (_, k) => [12 + 3 * k, 15 + 3 * k] as const),
            ...Array.from({ length: 9 }, (_, k) => [78 + 3 * k, 81 + 3 * k] as const),
          ]) {
            // Clipped to the view, as the part is drawn.
            const x0 = Math.max(v.xd(0, v.visLo[0]), v.xd(0, a));
            const x1 = Math.min(v.xd(0, v.visHi[0]), v.xd(0, b));
            if (Math.abs((x0 + x1) / 2 - d.x) > 0.5) continue;
            ok = true;
            expect(d.x - half, `cell ${cell} dpr ${dpr} ${d.text} left`).toBeGreaterThanOrEqual(x0 + MARGIN * dpr - 1);
            expect(d.x + half, `cell ${cell} dpr ${dpr} ${d.text} right`).toBeLessThanOrEqual(x1 - MARGIN * dpr + 1);
          }
          expect(ok, `cell ${cell}: ${d.text} at ${d.x} is over a known part`).toBe(true);
        }
      }
      t.ctx.dispose();
    }
    // Both forms are met across the zooms.
    expect(threes).toBeGreaterThan(0);
    expect(ones).toBeGreaterThan(0);
    // His12 (the split codon, two cells in its larger part) in three letters at some zoom, inside its part.
    expect(his).toBeGreaterThan(0);
  });
});

describe('reductions count every item', () => {
  it('a value narrower than a pixel is drawn', async () => {
    const t = await drawn(rangesDoc([[5000, 5001]], 20000), { width: 900, height: 200 });
    await zoom(t, layerOf(t.ctx, 'Values').section, 0.03);
    expect(inked(t, 'Values')).toBeGreaterThan(0);
  });
});

describe('color layers', () => {
  it('without values (no per-position values held): the covered cells alone; with values: every cell', async () => {
    const t = await drawn(colorDoc(), { width: 900, height: 200 });
    const { k, section } = layerOf(t.ctx, 'Bands');
    await zoom(t, section, 10);
    const state = t.ctx.graph.assembled.value?.view?.prepared[k].state as PerPosition;
    expect(state.value).toEqual([null]);
    const at = (label: string, i: number) => {
      const { rec, frame } = record(t.ctx, label);
      const x = (frame.v.xd(0, i) + frame.v.xe(0, i)) >> 1;
      return columns(rec, x, x + 1)[0] !== -2;
    };
    expect([15, 35, 25, 5].map((i) => at('Bands', i))).toEqual([true, true, false, false]);
    expect([5, 25, 59].map((i) => at('Shade', i))).toEqual([true, true, true]);
  });
});

describe('wrapped group labels', () => {
  /** 40 residues in an offset numbering (101 on), under a ruler. */
  const doc = (n = 40) => {
    const b = createSeqViewSpec();
    b.data
      .sequence({ ref: 'a', value: 'MKVLATGHEK'.repeat(5).slice(0, n) })
      .numbering({ name: 'auth', method: 'offset', first: 101 });
    const v = b.view({ ref: 'v', label: 'V', layout: 'text', numbering: 'auth' });
    v.track({ label: '' }).ruler();
    v.track({ label: 'Letters' }).letters();
    return b.getDocument();
  };

  it.each([1, 2])(
    "above-end: right-aligned to the end of each group's last letter, smaller, at DPR %i",
    async (dpr) => {
      const t = await drawn(doc(), { width: 900, height: 200, dpr });
      const { section } = layerOf(t.ctx, 'Letters');
      const labels = texts(recordLabels(t.ctx, section));
      const { rec, frame } = record(t.ctx, 'Letters');
      expect(labels.map((d) => d.text)).toEqual(['110', '120', '130', '140']);
      for (const [g, d] of labels.entries()) {
        expect(d.align).toBe(TextAlign.right);
        expect(d.x).toBe(frame.v.xe(0, 10 * g + 9));
      }
      const px = (font: number) => t.ctx.text.fontString(font);
      const letter = texts(rec).find((d) => d.text === 'M');
      expect(Number.parseFloat(/(\d+)px/.exec(px(labels[0].font))?.[1] ?? '')).toBeLessThan(
        Number.parseFloat(/(\d+)px/.exec(px(letter?.font ?? 0))?.[1] ?? '') * 0.7,
      );
      // No inline labels on the letters.
      expect(texts(rec).filter((d) => /^\d+$/.test(d.text))).toEqual([]);
    },
  );

  it('above-end: a short last group keeps its label unless it would reach the previous one', async () => {
    const labels = async (n: number) => {
      // At 15 px a letter, where 141's label reaches 140's.
      const t = await drawn(doc(n), { width: 900, height: 400, render: { wrappedFontPx: 15 } });
      return texts(recordLabels(t.ctx, layerOf(t.ctx, 'Letters').section)).map((d) => d.text);
    };
    expect(await labels(41)).toEqual(['110', '120', '130', '140']);
    expect(await labels(45)).toEqual(['110', '120', '130', '140', '145']);
  });

  it.each([1, 2])(
    'inline-start: before its group, right-aligned 1 px before its first letter, at DPR %i',
    async (dpr) => {
      const t = await drawn(doc(), { width: 900, height: 200, dpr, render: { wrappedGroupLabels: 'inlineStart' } });
      const { section } = layerOf(t.ctx, 'Letters');
      const { rec, frame } = record(t.ctx, 'Letters');
      const sups = texts(rec).filter((d) => /^\d+$/.test(d.text));
      expect(sups.map((d) => d.text)).toEqual(['101', '111', '121', '131']);
      for (const [g, d] of sups.entries()) {
        expect(d.align).toBe(TextAlign.right);
        expect(d.x).toBe(frame.v.xd(0, 10 * g) - Math.round(dpr));
        expect(d.x).toBeGreaterThan(g > 0 ? frame.v.xe(0, 10 * g - 1) : 0);
      }
      expect(texts(recordLabels(t.ctx, section))).toEqual([]);
    },
  );

  it('none: no labels anywhere', async () => {
    const t = await drawn(doc(), { width: 900, height: 200, render: { wrappedGroupLabels: 'none' } });
    const { section } = layerOf(t.ctx, 'Letters');
    expect(texts(recordLabels(t.ctx, section))).toEqual([]);
    expect(texts(record(t.ctx, 'Letters').rec).filter((d) => /^\d+$/.test(d.text))).toEqual([]);
  });
});

describe('the ruler over several blocks side by side', () => {
  const PROTEIN = 'MSEEKLKELAEKAIRLAKEGNPEEALKLLEEALRLAKESGDPELIARVLELLAEVTKQSVHAYRLDDGKWFLQPGGRVSLEMAPTDLQAAH';
  const residues = (n: number, seed: number) =>
    Array.from({ length: n }, (_, i) => PROTEIN[(i * 7 + seed * 13) % PROTEIN.length]).join('');
  /** Blocks of the given lengths side by side (named `<prefix>1`, `<prefix>2`, …), under a ruler. */
  function blocks(lengths: readonly number[], numbering?: { first: number }, prefix = 'b') {
    const b = createSeqViewSpec();
    const seqs = lengths.map((n, k) => {
      const s = b.data.sequence({ ref: `${prefix}${k + 1}`, label: `${prefix}${k + 1}`, value: residues(n, k) });
      if (numbering) s.numbering({ name: 'auth', method: 'offset', first: numbering.first });
      return s;
    });
    const v = b.view({ ref: 'v', label: 'V', ...(numbering ? { numbering: 'auth' } : {}) });
    for (const s of seqs) v.block({ sequence: s });
    v.track({ label: '' }).ruler();
    return b.getDocument();
  }
  type Label = { text: string; x0: number; x1: number };
  /** The ruler's labels as device px extents (names left-aligned, numbers centered). */
  function labels(t: Awaited<ReturnType<typeof drawn>>): Label[] {
    const { rec } = record(t.ctx, '');
    const dpr = entryOf(t.ctx, '').dpr;
    return texts(rec).map((d) => {
      const w = t.ctx.text.label(d.text, d.font).width * dpr;
      const x0 = d.align === TextAlign.center ? d.x - w / 2 : d.x;
      return { text: d.text, x0, x1: x0 + w };
    });
  }
  const overlaps = (ls: readonly Label[]) => {
    const s = [...ls].sort((a, b) => a.x0 - b.x0);
    return s
      .slice(1)
      .filter((l, i) => l.x0 < s[i].x1)
      .map((l) => l.text);
  };
  const isName = (l: Label) => /^(b|chr)\d+$/.test(l.text);

  it('narrow blocks: every block named at its start (over the gap), no numbers, nothing overlapping', async () => {
    const t = await drawn(blocks(Array(16).fill(120_000), undefined, 'chr'), { width: 900, height: 120 });
    const drawnLabels = labels(t);
    const { frame } = record(t.ctx, '');
    const v = frame.v;
    expect(drawnLabels.filter(isName).map((l) => l.text)).toEqual(Array.from({ length: 16 }, (_, k) => `chr${k + 1}`));
    for (const [k, l] of drawnLabels.filter(isName).entries())
      if (k < 15) expect(l.x0).toBe(v.xd(k, 0) + 2);
      else expect(l.x0).toBeGreaterThan(v.xd(14, v.axis?.shown[14] ?? 0));
    expect(drawnLabels.filter((l) => !isName(l))).toEqual([]);
    expect(overlaps(drawnLabels)).toEqual([]);
  });

  it('wide blocks: names and round numbers in each, numbers repeating per block but never overlapping', async () => {
    const t = await drawn(blocks([400, 300]), { width: 1200, height: 120 });
    const drawnLabels = labels(t);
    expect(drawnLabels.filter(isName).map((l) => l.text)).toEqual(['b1', 'b2']);
    const b2 = drawnLabels.find((l) => l.text === 'b2') as Label;
    expect(drawnLabels.filter((l) => !isName(l) && l.x1 <= b2.x0).length).toBeGreaterThan(2);
    expect(drawnLabels.filter((l) => !isName(l) && l.x0 >= b2.x1).length).toBeGreaterThan(1);
    expect(overlaps(drawnLabels)).toEqual([]);
  });

  it("zoomed in across a block boundary: the next block's name is never covered by the last numbers before it", async () => {
    const t = await drawn(blocks([400, 300]), { width: 700, height: 120 });
    const { section } = layerOf(t.ctx, '');
    for (const cell of [3, 5, 9, 14]) {
      // The end of block 1 at the data area's middle.
      await zoom(t, section, cell, Math.max(0, 400 * cell - 300));
      const drawnLabels = labels(t);
      expect(
        drawnLabels.some((l) => l.text === 'b2'),
        `b2 at cell ${cell}`,
      ).toBe(true);
      expect(overlaps(drawnLabels), `cell ${cell}`).toEqual([]);
    }
  });

  it('round numbers found by jumps equal those found position by position (an offset numbering)', async () => {
    const t = await drawn(blocks([5000], { first: -37 }), { width: 1000, height: 120 });
    const { section } = layerOf(t.ctx, '');
    for (const [cell, x] of [
      [0.2, 0],
      [2, 1234],
      [12, 30000],
    ]) {
      await zoom(t, section, cell, x);
      const nums = labels(t).map((l) => Number(l.text));
      expect(nums.length, `cell ${cell}`).toBeGreaterThan(2);
      const step = nums[1] - nums[0];
      expect(
        nums.every((n, i) => n % step === 0 && (i === 0 || n - nums[i - 1] === step)),
        `${nums}`,
      ).toBe(true);
    }
  });
});

describe('values zoomed out over 400,000 positions', () => {
  it('a light min–max band under the mean where columns aggregate, none per cell', async () => {
    const t = await drawn(scaleDoc(), { width: 800, height: 420, dpr: 2 });
    const { section } = layerOf(t.ctx, 'Values');
    const light = () => {
      const r = record(t.ctx, 'Values').rec.rectItems;
      let n = 0;
      for (let i = 0; i < r.n; i++) if (r.alpha[i] < 255) n += r.w[i] * r.h[i];
      return n;
    };
    await zoom(t, section, 0.004);
    expect(light()).toBeGreaterThan(2000);
    await zoom(t, section, 12, 1000);
    expect(light()).toBe(0);
  });
});

describe('values colored by their values', () => {
  const BAND = [0xff7d45, 0xffdb13, 0x65cbf3, 0x0053d6];
  /** The layer's line as subpaths: its color and its points (device px). */
  const lines = (t: Awaited<ReturnType<typeof drawn>>, label: string) => {
    const p = record(t.ctx, label).rec.pathItems;
    return Array.from({ length: p.n }, (_, q) => {
      const [c0, c1] = [p.cstart[q], q + 1 < p.n ? p.cstart[q + 1] : p.nc];
      const pts: [number, number][] = [];
      for (let c = c0; c < c1; c += 2) pts.push([p.coords[c], p.coords[c + 1]]);
      return { color: p.color[q], pts };
    });
  };

  it("a line: each position's band over its own cell, changing halfway between two positions", async () => {
    const t = await drawn(linesDoc(), { width: 1000, height: 400, dpr: 1 });
    await zoom(t, layerOf(t.ctx, 'Line').section, 8);
    const subs = lines(t, 'Line');
    expect(new Set(subs.map((x) => x.color))).toEqual(new Set(BAND));
    // Each band's points at one height (the scores are runs of one value): a point at a band's height is in its color.
    const heightOf = new Map<number, number>();
    for (const { color, pts } of subs)
      for (const [, y] of pts.slice(1, -1)) {
        expect(heightOf.get(color) ?? y, `color ${color.toString(16)}`).toBe(y);
        heightOf.set(color, y);
      }
    expect(heightOf.size).toBe(4);
    const ys = BAND.map((c) => heightOf.get(c) as number);
    expect([...ys].sort((a, b) => b - a)).toEqual(ys); // higher scores higher up
    let splits = 0;
    for (let q = 1; q < subs.length; q++) {
      const [a, b] = [subs[q - 1], subs[q]];
      const [end, start] = [a.pts[a.pts.length - 1], b.pts[0]];
      if (end[0] !== start[0] || end[1] !== start[1]) continue;
      // A change of color: from the last position of one band halfway to the first of the next.
      splits++;
      expect(a.color).not.toBe(b.color);
      const [prev, next] = [a.pts[a.pts.length - 2], b.pts[1]];
      expect([end[0] - prev[0], end[1] * 2]).toEqual([next[0] - end[0], prev[1] + next[1]]);
      expect([prev[1], next[1]]).toEqual([heightOf.get(a.color), heightOf.get(b.color)]);
    }
    expect(splits).toBeGreaterThan(4);
  });

  it('an area per column in its color, a line with no theme in one color, and both where columns aggregate', async () => {
    const t = await drawn(linesDoc(4000), { width: 1000, height: 400, dpr: 1 });
    const { section } = layerOf(t.ctx, 'Line');
    await zoom(t, section, 8);
    const area = (alpha: number) =>
      new Set(rects(record(t.ctx, 'Area').rec).flatMap(([, , , , c, a]) => (a === alpha ? [c] : [])));
    expect(area(77).size).toBe(BAND_SCORES.length);
    expect(new Set(lines(t, 'Plain').map((x) => x.color)).size).toBe(1);
    // Zoomed out: a column's mean picks its color, the min–max band under it in the same.
    await zoom(t, section, 0.2);
    const line = new Set(lines(t, 'Line').map((x) => x.color));
    expect([...line].every((c) => BAND.includes(c)) && line.size > 1).toBe(true);
    const band = rects(record(t.ctx, 'Line').rec).filter(([, , , , , a]) => a === 89);
    expect(band.length).toBeGreaterThan(0);
    expect(band.every(([, , , , c]) => BAND.includes(c))).toBe(true);
    expect(area(77).size).toBeGreaterThan(BAND_SCORES.length);
    expect(new Set(lines(t, 'Plain').map((x) => x.color)).size).toBe(1);
  });

  it("zoomed out, a column's color is its mean's: scores alternating between the end bands draw in the second", async () => {
    const t = await drawn(linesDoc(4000, { alternate: true }), { width: 1000, height: 400, dpr: 1 });
    const { section } = layerOf(t.ctx, 'Line');
    await zoom(t, section, 8);
    expect(new Set(lines(t, 'Line').map((x) => x.color))).toEqual(new Set([BAND[0], BAND[3]]));
    // 10 positions a column: 20.5 and 95.5 in turn average about 58, the second band's.
    await zoom(t, section, 0.1);
    expect(new Set(lines(t, 'Line').map((x) => x.color))).toEqual(new Set([BAND[1]]));
    const band = rects(record(t.ctx, 'Line').rec).filter(([, , , , , a]) => a === 89);
    expect(band.length).toBeGreaterThan(0);
    expect(new Set(band.map(([, , , , c]) => c))).toEqual(new Set([BAND[1]]));
  });

  it.each([
    ['linear', undefined],
    ['on rings', 'circular'],
  ] as const)("in an SVG export, %s: the line's strokes in the four bands", async (_, style) => {
    const xml = await renderSvg(linesDoc(160, { style }), { width: 800 });
    const strokes = new Set([...xml.matchAll(/stroke="(#[0-9a-f]{6})"/gi)].map((m) => m[1].toLowerCase()));
    expect(BAND.map((c) => `#${c.toString(16).padStart(6, '0')}`).filter((c) => !strokes.has(c))).toEqual([]);
    expect(xml).not.toMatch(/NaN|Infinity|undefined/);
  });
});
