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

  it('zoomed out (EC-171, EC-208): bin pairs stroked by count, the BEDPE’s bin pairs and the held links binned by column, the hover naming the bin pair', async () => {
    const v = (viewer = await openRegion('pairs-zoomed-out'));
    await advance(v, 8000);
    const [held, bedpe] = lays(v);
    // more than 4,096 links in view: binned, said once
    expect(v.layoutIssues.value.map((i) => i.code)).toEqual(['seq.render.pairs.limit']);
    // the pairs read by region (the BEDPE): bin pairs, their counts summing to the links
    expect(bedpe.prep?.weight).toBeDefined();
    expect(bedpe.n).toBeGreaterThan(100);
    expect(Array.from(bedpe.prep!.weight!).reduce((a, b) => a + b, 0)).toBe(12_000);
    expect(Math.max(...bedpe.width.subarray(0, bedpe.n))).toBeGreaterThan(
      Math.min(...bedpe.width.subarray(0, bedpe.n)),
    );
    // the held links past the gate: arcs binned by their ends' device columns (the pairs level is the store's, a held
    // table is not read through it), a density strip of the rest
    expect(held.prep?.weight).toBeUndefined();
    expect(held.n).toBeGreaterThan(100);
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
      expect(svg.issues.map((i) => `${i.code}: ${i.message}`)).toEqual([
        'seq.render.pairs.limit: More than 4096 pairs in view: binned',
      ]);
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

  it('on a genome ring at the fit: chords between the chromosomes and within each, none dangling, in both themes', async () => {
    const v = (viewer = await openRegion('pairs-ring'));
    await advance(v, 8000);
    const [lay] = lays(v);
    // the toy assembly's 1,376 positions at about a position a pixel: the links themselves, past 4,096 binned by column
    expect(lay.n).toBeGreaterThan(100);
    expect(lay.carets.length).toBe(0);
    const prep = (
      lay as unknown as { prep: { endUnits: string[]; data: { aUnit: Int32Array; bUnit: Int32Array; n: number } } }
    ).prep;
    expect(prep.data.n).toBe(12_000);
    let across = 0;
    for (let k = 0; k < prep.data.n; k++) if (prep.data.aUnit[k] !== prep.data.bUnit[k]) across++;
    expect(across).toBeGreaterThan(100);
    expect(new Set(Array.from(prep.data.aUnit)).size).toBe(6);
    expect(drawn(v, 'path')).toBeGreaterThan(100);
    for (const dark of [false, true]) {
      const svg = await exported(v, dark);
      expect(svg.settled).toBe(true);
    }
    expect((await exported(v, true)).svg).not.toBe((await exported(v, false)).svg);
  });
});
