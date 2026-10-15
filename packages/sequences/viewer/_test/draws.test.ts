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
// Strip lanes in text (seqspec:viewer#the-text-view): a `values` layer a strip under each line by run of one color
// within a group, nothing at a missing value, under the underline lanes; `color` alone a strip, beside the letters
// their cells. A ligand's pill never vanishes: narrower than its box a filled cell; in text a chip holding its whole
// label, every chip of a row kept through a scroll's redraw of its cut row (a chip's clear is its own), and the hover
// and the selection drawn on the chip.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { describe, expect, it } from 'vitest';
import { Pass, SkinToken, TextAlign } from '../src/draw/backend.js';
import { Dirty } from '../src/draw/index.js';
import type { PerPosition } from '../src/prepare.js';
import { renderSvg } from '../src/svg/index.js';
import { LayerFrame, Op, RecordingBackend } from '../src/testing/index.js';
import { COPY_BACK, recording } from './composer-fixtures.js';
import { lit, passRects } from './overlay-fixtures.js';
import { columns, drawn, entryOf, layerOf, record, recordLabels, rects, texts } from './draw-fixtures.js';
import {
  allTypes,
  BAND_SCORES,
  BANDS,
  chipsDoc,
  colorDoc,
  gene,
  linesDoc,
  rangesDoc,
  scaleDoc,
} from './prepare-fixtures.js';

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

describe("wrapped text in the layer's colors (wrappedColors)", () => {
  const doc = () => {
    const b = createSeqViewSpec();
    b.data.sequence({ ref: 'a', value: 'MKVLATGHEK'.repeat(4) });
    const v = b.view({ ref: 'v', label: 'V', layout: 'text' });
    v.track({ label: 'Letters' }).letters({
      color: { theme: 'letter', params: { colors: 'Hydrophobicity' } },
    } as never);
    return b.getDocument();
  };
  const colors = async (render?: object) => {
    const t = await drawn(doc(), { width: 900, height: 200, render: { wrappedGroupLabels: 'none', ...render } });
    return new Set(texts(record(t.ctx, 'Letters').rec).map((d) => d.color));
  };

  it('are the ink by default, as text is: one color for every letter', async () => {
    expect((await colors()).size).toBe(1);
  });

  it("are each letter's color from the theme with wrappedColors", async () => {
    const tinted = await colors({ wrappedColors: true });
    expect(tinted.size).toBeGreaterThan(3);
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

  it("zoomed in, the line holds its end values to the cells' edges, as the area: a block's ends and around a gap", async () => {
    const t = await drawn(linesDoc(), { width: 1000, height: 400, dpr: 1 });
    const { section } = layerOf(t.ctx, 'Area');
    await zoom(t, section, 8);
    const { rec, frame } = record(t.ctx, 'Area');
    const p = rec.pathItems;
    const xs = new Set<number>();
    for (let c = 0; c < p.nc; c += 2) xs.add(p.coords[c]);
    const area = rects(rec).filter((r) => r[5] === 77);
    const v = frame.v;
    // The block's first cell (in view at the left), and the cells either side of the missing value at 49.
    expect(xs.has(Math.min(...area.map((r) => r[0])))).toBe(true);
    expect(xs.has(v.xd(0, 49))).toBe(true);
    expect(xs.has(v.xd(0, 50))).toBe(true);
    // No point inside the missing cell: the line breaks there.
    expect([...xs].some((x) => x > v.xd(0, 49) && x < v.xd(0, 50))).toBe(false);
    // Zoomed out (a cell under a device px), columns as before: no point at a cell's edge is added.
    await zoom(t, section, 0.1);
    const out = record(t.ctx, 'Area').rec.pathItems;
    expect(out.n).toBeGreaterThan(0);
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

describe('strip lanes in text', () => {
  const BAND = [0xff7d45, 0xffdb13, 0x65cbf3, 0x0053d6];
  /**
   * 100 residues as text: an underlined region, then pLDDT-like values in four bands (5 positions each, so a group of
   * 10 holds two runs), missing at 49 and 99 (`Plddt`), a `color` layer alone (`Shade`) and one beside the letters
   * (`Behind`).
   */
  const doc = () => {
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: 'a', value: 'MKVLATGHEK'.repeat(10) });
    const r = b.data.annotation({ ref: 'r', geometry: 'ranges', target: a, value: { start: [3], stop: [60] } });
    const value = Array.from({ length: 100 }, (_, i) => (i % 50 === 49 ? null : BAND_SCORES[Math.floor(i / 5) % 4]));
    const d = b.data.annotation({ ref: 'd', geometry: 'values', target: a, value: { value } } as never);
    const palette = { theme: 'field_palette', params: { field: 'value', palette: BANDS } };
    const v = b.view({ ref: 'v', label: 'V', layout: 'text' });
    v.track({ label: 'Letters' }).letters();
    v.track({ label: 'Region' }).ranges({ data: r, display: 'underline' } as never);
    v.track({ label: 'Plddt' }).values({ data: d, domain: [0, 100], color: palette } as never);
    v.track({ label: 'Shade' }).color({ data: d, color: palette } as never);
    const behind = v.track({ label: 'Behind' });
    behind.letters();
    behind.color({ data: d, color: palette } as never);
    return b.getDocument();
  };
  /** The rects of a track's layer (its first, or its `type` layer) drawn through the kit's composer. */
  const strips = (t: Awaited<ReturnType<typeof drawn>>, label: string, type?: string) => {
    const first = layerOf(t.ctx, label);
    const m = t.ctx.graph.assembled.value?.view;
    if (!m) throw new Error('no view');
    let k = first.k;
    while (type && m.set.strings.get(m.resolved.layers.type[k]) !== type) k++;
    const entry = t.ctx.registry.get('renderers', m.prepared[k].renderer);
    if (!entry) throw new Error(`no renderer for ${label}`);
    const rec = new RecordingBackend();
    new LayerFrame(t.ctx, entry, first.section, k).draw(rec);
    return rects(rec);
  };

  it.each([1, 2])(
    'values: a strip at the bottom of its band, a rect per run of one color in a group, none at a missing value, DPR %i',
    async (dpr) => {
      const t = await drawn(doc(), { width: 900, height: 600, dpr, render: { wrappedGroupLabels: 'none' } });
      const { section } = layerOf(t.ctx, 'Plddt');
      const l = t.ctx.layout.get(section)?.current;
      const ax = t.ctx.layout.axis(section);
      if (!l || !ax) throw new Error('no layout');
      const { k } = layerOf(t.ctx, 'Plddt');
      const tr = t.ctx.graph.assembled.value?.view?.resolved.layers.track[k];
      const bottoms = new Set<number>();
      for (let e = 0; e < l.entryCount; e++)
        if (l.entryTrack[e] === tr) bottoms.add(Math.round((l.scrolledTop + l.entryTop[e] + l.entryHeight[e]) * dpr));
      const r = strips(t, 'Plddt');
      // Every rect 3 px high at its band's bottom, in a band's color.
      expect(r.every(([, y, , h]) => h === 3 * dpr && bottoms.has(y + h))).toBe(true);
      expect(new Set(r.map((x) => x[4]))).toEqual(new Set(BAND));
      // Runs of 5 within groups of 10: 20 runs on 100 positions, less the two missing ends (49 and 99 shorten theirs).
      expect(r.length).toBe(20);
      // A missing value draws nothing: no rect over position 49's cell on its line (the first).
      const lines = l.lines;
      expect(lines.stop[0]).toBeGreaterThan(49);
      const { frame } = record(t.ctx, 'Plddt');
      frame.v.lineStart = lines.start[0];
      frame.v.lineStop = lines.stop[0];
      frame.v.onlyBlock = lines.block[0];
      const [x0, x1] = [frame.v.xd(0, 49), frame.v.xe(0, 49)];
      const top = Math.min(...r.map((x) => x[1]));
      expect(r.filter(([x, y, w]) => y === top && x < x1 && x + w > x0)).toEqual([]);
      expect(r.filter(([x, y, w]) => y === top && x < x0 && x + w > x0 - 4 * dpr).length).toBe(1);
    },
  );

  it('under the underline lanes, in track order; color alone a strip, beside the letters their cells', async () => {
    const t = await drawn(doc(), { width: 900, height: 600, render: { wrappedGroupLabels: 'none' } });
    const under = strips(t, 'Region');
    const strip = strips(t, 'Plddt');
    const firstLine = (r: number[][]) => Math.min(...r.map((x) => x[1]));
    // On the first line the strip lies below its underline.
    const u = under.filter((x) => x[1] === firstLine(under));
    expect(firstLine(strip)).toBeGreaterThanOrEqual(Math.max(...u.map((x) => x[1] + x[3])));
    // Color alone: the same strip as the values, one band lower.
    const shade = strips(t, 'Shade');
    expect(shade.map((x) => [x[0], x[2], x[3], x[4]])).toEqual(strip.map((x) => [x[0], x[2], x[3], x[4]]));
    expect(firstLine(shade) - firstLine(strip)).toBe(entryOf(t.ctx, 'Shade').top - entryOf(t.ctx, 'Plddt').top);
    // Beside the letters: a cell a position over the letters' band (less its 2 px padding).
    const behind = strips(t, 'Behind', 'color');
    const band = entryOf(t.ctx, 'Behind');
    expect(behind.length).toBe(98);
    expect(behind.every(([, y, , h]) => y === Math.round(band.top + 2) && h === band.height - 4)).toBe(true);
  });

  it('values with no color theme: a strip in the ink beside nothing, nothing beside the letters', async () => {
    const b = createSeqViewSpec();
    const a = b.data.sequence({ ref: 'a', value: 'MKVLATGHEK'.repeat(4) });
    const value = Array.from({ length: 40 }, (_, i) => i / 40);
    const d = b.data.annotation({ ref: 'd', geometry: 'values', target: a, value: { value } });
    const v = b.view({ ref: 'v', label: 'V', layout: 'text' });
    v.track({ label: 'Plain' }).values({ data: d });
    const both = v.track({ label: 'Both' });
    both.letters();
    both.values({ data: d });
    const t = await drawn(b.getDocument(), { width: 900, height: 300, render: { wrappedGroupLabels: 'none' } });
    const plain = strips(t, 'Plain');
    // One run a group: 4 groups on one line.
    expect([plain.length, new Set(plain.map((x) => x[4])).size]).toEqual([4, 1]);
    expect(strips(t, 'Both', 'values')).toEqual([]);
  });
});

describe("a ligand's pill", () => {
  const LIGANDS = ['HEM', 'ERY', 'NI', 'ATP', 'NAG', 'SO4'];
  const tracks = () => {
    const doc = chipsDoc(LIGANDS) as unknown as {
      root: { children: { kind: string; params: Record<string, unknown> }[] };
    };
    const view = doc.root.children.find((c) => c.kind === 'view');
    if (view) view.params.layout = 'tracks';
    return doc as unknown as Parameters<typeof drawn>[0];
  };

  it.each([1, 2])(
    'at the fit of a narrow view, a filled cell at least 2 device px wide per ligand, at DPR %i',
    async (dpr) => {
      // The blocks' gaps (28 px) take most of the narrowest: a ligand's three cells a device px or two there.
      for (const width of [360, 400, 520, 900]) {
        const t = await drawn(tracks(), { width, height: 200, dpr });
        const { rec, frame } = record(t.ctx, 'Letters');
        const e = entryOf(t.ctx, 'Letters');
        const faint = rects(rec).filter(([, , w, , c]) => c === SkinToken.inkFaint && w >= 2 * dpr);
        // One box (or filled cell) per ligand in view (the fit leaves the last out where its gaps do not fit), inside
        // the data area.
        let shown = 0;
        for (let b = 1; b <= LIGANDS.length; b++) if (frame.v.visHi[b] > frame.v.visLo[b]) shown++;
        expect(shown, `width ${width}`).toBeGreaterThanOrEqual(LIGANDS.length - 1);
        expect(faint.length, `width ${width}`).toBe(shown);
        expect(faint.every(([x, , w]) => x >= Math.floor(e.left * dpr) && x + w <= Math.ceil(e.right * dpr))).toBe(
          true,
        );
        t.ctx.dispose();
      }
    },
  );

  it('in text, a chip holding its whole label', async () => {
    const t = await drawn(chipsDoc(LIGANDS), { width: 900, height: 300, render: { wrappedGroupLabels: 'none' } });
    const labels = texts(record(t.ctx, 'Letters').rec)
      .map((d) => d.text)
      .filter((x) => LIGANDS.includes(x));
    expect(labels).toEqual(LIGANDS);
  });
});

describe('chips through redraws', () => {
  const LIGANDS = ['HEM', 'ERY', 'NI', 'ZN', 'ATP', 'NAG', 'SO4', 'GOL'];
  /** The faint rects (pills' borders) of a recording that no later clear covers: what stays on the canvas. */
  function surviving(r: RecordingBackend): number[][] {
    const kept: number[][] = [];
    const b = r.rectItems;
    for (let k = 0; k < r.opCount; k++) {
      const code = r.ops[4 * k];
      if (code === Op.rects) {
        for (let i = r.ops[4 * k + 1]; i < r.ops[4 * k + 1] + r.ops[4 * k + 2]; i++)
          if (b.color[i] === SkinToken.inkFaint) kept.push([b.x[i], b.y[i], b.w[i], b.h[i]]);
      } else if (code === Op.clear) {
        const a = r.ops[4 * k + 3];
        const [x, y, w, h] = [r.args[a], r.args[a + 1], r.args[a + 2], r.args[a + 3]];
        for (let j = kept.length - 1; j >= 0; j--) {
          const q = kept[j];
          if (q[0] < x + w && q[0] + q[2] > x && q[1] < y + h && q[1] + q[3] > y) kept.splice(j, 1);
        }
      }
    }
    return kept;
  }

  it('a scroll that cuts a row of chips redraws every chip on it, across a wrap at a narrow width', async () => {
    const t = await lit(chipsDoc(LIGANDS, true), { width: 220, height: 60 });
    const L = t.layout(t.section).lines;
    // The chips wrap onto two rows or more.
    const rows = new Set(Array.from({ length: L.count }, (_, k) => (L.block[k] >= 1 ? L.top[k] : -1)));
    rows.delete(-1);
    expect(rows.size).toBeGreaterThan(1);
    const be = recording(COPY_BACK);
    t.frame(be, Dirty.all);
    let checked = 0;
    for (let i = 0; i < 40; i++) {
      t.vp?.panBy(0, 5, 'pointer');
      be.reset();
      t.frame(be, Dirty.tracks);
      const kept = surviving(be);
      // Each chip row the scroll cut was redrawn: all its chips' pills stay.
      const l = t.layout(t.section);
      for (const top of rows) {
        const on = Array.from({ length: L.count }, (_, k) => k).filter((k) => L.block[k] >= 1 && L.top[k] === top);
        // The pills' tops: 3 px into the row's letters band.
        let e = 0;
        while (!(l.entryRow[e] === on[0] && l.entryKind[e] === 0)) e++;
        const y = Math.round(l.scrolledTop + l.entryTop[e] - (t.vp?.live.current.scrollY ?? 0) + 3);
        // By x: a frame may draw a pill twice (its row redrawn whole after the strip that uncovered it).
        const drawn = new Set(kept.filter((r) => Math.abs(r[1] - y) <= 1).map((r) => r[0]));
        if (drawn.size === 0) continue;
        expect(drawn.size, `row at ${top}, scroll ${t.vp?.live.current.scrollY}`).toBe(on.length);
        checked++;
      }
    }
    expect(checked).toBeGreaterThan(0);
  });

  it("the hover and the selection on a chip are drawn at the chip's x, not its row's start", async () => {
    const t = await lit(chipsDoc(LIGANDS, true), { width: 500, height: 600 });
    const l = t.layout(t.section);
    const L = l.lines;
    const ax = t.ctx.layout.axis(t.section);
    for (let k = 0; k < L.count; k++) {
      if (L.block[k] < 1 || L.left[k] === 0) continue;
      const x0 = Math.round(l.dataLeft + L.left[k]);
      const w = Math.round(ax?.runs[L.block[k]]?.advance ?? 0);
      t.hover(x0 + 4, l.scrolledTop + L.top[k] + l.entryHeight[0] / 2 + 12);
      const hovered = passRects(t.overlay(), Pass.overlay);
      expect(hovered.length, `chip ${k}`).toBeGreaterThan(0);
      expect(hovered.every((r) => Math.abs(r[0] - x0) <= 1 && r[2] <= w + 1)).toBe(true);
      t.ctx.selection.set(t.ctx.hover.toLoci() as never, 'api');
      const selected = passRects(t.overlay(), Pass.overlay);
      expect(selected.every((r) => r[0] >= x0 - 1 && r[0] + r[2] <= x0 + w + 1)).toBe(true);
      t.ctx.selection.set(null as never, 'api');
    }
  });
});
