/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { hitView, pairsBinHover } from '@molstar/seqx-core/render';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { loadGenome } from '@molstar/seqx-core/resolve';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { advance, drawn, exported, openRegion } from './region-twin.js';

// Features/Region reads/Pairs, headless on a fake clock (the story's composition): pairs by region zoomed in (items as
// arcs, EC-208), zoomed out (bin pairs stroked by count, the hover naming the bin pair, EC-171, EC-208) and on a ring
// (bin-pair chords). Each in the canvas and in the SVG, light and dark.

let viewer: SeqViewer | undefined;
// the genome chunk landed first: its import is real I/O, never on the fake clock
beforeAll(async () => {
  await loadGenome();
});
afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

interface PairsLay {
  readonly n: number;
  readonly mode: string;
  readonly carets: { readonly length: number };
  readonly prep?: { readonly weight?: Float32Array; readonly level?: number };
  readonly curve: Float64Array;
  readonly width: Float32Array;
}

/** The pairs layers of the view's first section. */
const lays = (v: SeqViewer) => v.frame!.layout.sections[0].sectionLays.map((l) => l.lay as unknown as PairsLay);

describe('Features/Region reads/Pairs, headless', () => {
  it('zoomed in (EC-208, EC-210): the items as arcs between two assemblies’ chromosomes, each alignment block one, read from either end; waiting first, then drawn', async () => {
    const v = (viewer = await openRegion('pairs-zoomed-in'));
    // before the source answers: no arc, no issue
    await advance(v, 100);
    expect(lays(v)[0]?.n ?? 0).toBe(0);
    await advance(v, 7000);
    await advance(v, 20_000);
    const [lay] = lays(v);
    // under the gate: the items, no bins
    expect(lay.prep?.weight).toBeUndefined();
    expect(lay.n).toBeGreaterThan(30);
    expect(lay.n).toBeLessThan(4096);
    // each pair's ends on the two genomes' blocks (A's chrA2, B's chrB2), found once from either
    const prep = (
      lay as unknown as { prep: { endUnits: string[]; data: { aUnit: Int32Array; bUnit: Int32Array; n: number } } }
    ).prep;
    const ends = new Set<string>();
    for (let k = 0; k < prep.data.n; k++)
      ends.add(`${prep.endUnits[prep.data.aUnit[k]]} ${prep.endUnits[prep.data.bUnit[k]]}`);
    expect([...ends]).toEqual(['/data/a/(sequence :name chrA2) /data/b/(sequence :name chrB2)']);
    expect(lay.carets.length).toBe(0);
    // a block the alignment holds with its B end reversed is twisted: arcs of two kinds
    expect(drawn(v, 'path')).toBeGreaterThan(30);
    for (const dark of [false, true]) {
      const svg = await exported(v, dark);
      expect(svg.settled).toBe(true);
      expect(svg.issues.map((i) => `${i.code}: ${i.message}`)).toEqual([]);
      expect(svg.svg).toContain('<path');
    }
    expect((await exported(v, true)).svg).not.toBe((await exported(v, false)).svg);
  });

  it('zoomed out (EC-171, EC-208): bin pairs stroked by count, the BEDPE’s bin pairs and the held links’ one picture (the pairs level), the hover naming the bin pair', async () => {
    const v = (viewer = await openRegion('pairs-zoomed-out'));
    await advance(v, 8000);
    const [held, bedpe] = lays(v);
    // more than 4,096 links in view: both read as bin pairs, nothing to say (no arc left out)
    expect(v.layoutIssues.value.map((i) => i.code)).toEqual([]);
    // the pairs read by region (the BEDPE): bin pairs, their counts summing to the links
    expect(bedpe.prep?.weight).toBeDefined();
    expect(bedpe.n).toBeGreaterThan(100);
    expect(Array.from(bedpe.prep!.weight!).reduce((a, b) => a + b, 0)).toBe(12_000);
    expect(Math.max(...bedpe.width.subarray(0, bedpe.n))).toBeGreaterThan(
      Math.min(...bedpe.width.subarray(0, bedpe.n)),
    );
    // the held links past the gate: bin pairs at the pairs level the BEDPE is read at, one picture (the same bins, the
    // same counts)
    expect(held.prep?.weight).toBeDefined();
    expect(held.prep?.level).toBe(bedpe.prep?.level);
    const bins = (l: PairsLay) => {
      const d = (l.prep as unknown as { data: { n: number; aStart: Int32Array; bStart: Int32Array } }).data;
      const out = new Map<string, number>();
      for (let k = 0; k < d.n; k++) {
        const key = `${d.aStart[k]}|${d.bStart[k]}`;
        out.set(key, (out.get(key) ?? 0) + l.prep!.weight![k]);
      }
      return out;
    };
    expect(bins(held)).toEqual(bins(bedpe));
    // the strongest bin pair of the held track: hit at its arc's middle, its hover a bin pair, never a pair
    const widths = Array.from(bedpe.width.subarray(0, bedpe.n));
    const i = widths.indexOf(Math.max(...widths));
    const c = bedpe.curve.subarray(8 * i, 8 * i + 8);
    const x = (c[0] + 3 * c[2] + 3 * c[4] + c[6]) / 8;
    const y = (c[1] + 3 * c[3] + 3 * c[5] + c[7]) / 8;
    const h = hitView(v.frame!.layout, x, y);
    expect(String(h?.element?.key)).toMatch(/^bin:/);
    const card = pairsBinHover(bedpe as never, h!.element!.key, () => 'chr')!;
    expect(card.title).toMatch(/^chr:[\d,]+–[\d,]+ ↔ chr:[\d,]+–[\d,]+$/);
    expect(card.text).toMatch(/zoom in to see them$/);
    for (const dark of [false, true]) {
      const svg = await exported(v, dark);
      expect(svg.settled).toBe(true);
      expect(svg.issues.map((i) => `${i.code}: ${i.message}`)).toEqual([]);
    }
  });

  it('zoomed in to a few Mb of the 12,000 links: the items, the bins gone (the gate)', async () => {
    const v = (viewer = await openRegion('pairs-zoomed-out'));
    await advance(v, 8000);
    v.runCommand('address.focus', { address: '/data/chr/(range :from 0 :to 80000)' });
    await advance(v, 8000);
    const [, bedpe] = lays(v);
    expect(bedpe.prep?.weight).toBeUndefined();
    expect(bedpe.n).toBeGreaterThan(0);
    expect(bedpe.n).toBeLessThan(4096);
  });

  it('on a genome ring at the fit (region-reads Default 3): bin-pair chords between the chromosomes, stroked and shaded by count within the bound, the strongest kept', async () => {
    const v = (viewer = await openRegion('pairs-ring'));
    await advance(v, 8000);
    const [lay] = lays(v);
    const prep = (
      lay as unknown as {
        prep: { endUnits: string[]; weight?: Float32Array; data: { aUnit: Int32Array; bUnit: Int32Array; n: number } };
      }
    ).prep;
    // bin pairs (a count each), never the links themselves at the fit
    expect(prep.weight).toBeDefined();
    expect(lay.n).toBeGreaterThan(100);
    expect(lay.n).toBeLessThanOrEqual(4096);
    expect(lay.carets.length).toBe(0);
    // no chord's stroke past the bound (1 to 6 px), the strong ones thicker; alphas from a quarter to whole
    const widths = Array.from(lay.width.subarray(0, lay.n));
    expect(Math.max(...widths)).toBeLessThanOrEqual(6);
    expect(Math.min(...widths)).toBeGreaterThanOrEqual(1);
    expect(Math.max(...widths)).toBeGreaterThan(3);
    const color = (lay as unknown as { color: Uint32Array }).color;
    const alphas = Array.from(color.subarray(0, lay.n), (c) => c & 0xff);
    expect(Math.min(...alphas)).toBeLessThan(128);
    // chords join different chromosomes (and some stay within one), each end on one of the six
    let across = 0;
    for (let i = 0; i < lay.n; i++) {
      const k = (lay as unknown as { pair: Int32Array }).pair[i];
      if (prep.data.aUnit[k] !== prep.data.bUnit[k]) across++;
    }
    expect(across).toBeGreaterThan(50);
    expect(across).toBeLessThan(lay.n);
    expect(new Set(Array.from(prep.data.aUnit)).size).toBe(6);
    expect(drawn(v, 'path')).toBeGreaterThan(100);
    for (const dark of [false, true]) {
      const svg = await exported(v, dark);
      expect(svg.settled).toBe(true);
    }
    expect((await exported(v, true)).svg).not.toBe((await exported(v, false)).svg);
  });
});
