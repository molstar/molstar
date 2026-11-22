/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  DefaultOpenLimits,
  type PairsBins,
  type PairsList,
  type SourceEntry,
  binPairs,
  openSeq,
} from '@molstar/seqx-core/data';
import { DefaultRenderConfig } from '@molstar/seqx-core/render';
import { BedRuntime } from '@molstar/seqx-core/runtimes/bed';
import { describe, expect, it } from 'vitest';
import {
  LinkCount,
  LinkLevels,
  LinksChromosome,
  LinksChromosomes,
  RingChromosomes,
  linksBedpe,
  linksOnChromosome,
  linksSource,
} from '../regions/links.js';

const signal = new AbortController().signal;
const assets = {
  readText: async () => '',
  readBytes: async () => new Uint8Array(),
};

/** The BEDPE's pairs the BED runtime reads whole and serves with an end in a region (what the store's pyramid reads). */
async function bedpePairs(): Promise<PairsList> {
  const source: SourceEntry = { address: '/data/links', format: 'bedpe', value: linksBedpe() as never };
  const opened = await BedRuntime.open(source, { assets, signal, limits: DefaultOpenLimits } as never);
  const a = await (
    opened as unknown as {
      pairs(s: string, a: number, b: number, g: AbortSignal): Promise<{ ok: boolean; pairs: PairsList }>;
    }
  ).pairs(LinksChromosome.name, 0, LinksChromosome.length, signal);
  if (!a.ok) throw new Error('not read');
  return a.pairs;
}

const counts = (b: PairsBins) => {
  const m = new Map<string, number>();
  for (let i = 0; i < b.aBin.length; i++) {
    const lo = Math.min(b.aBin[i], b.bBin[i]);
    const hi = Math.max(b.aBin[i], b.bBin[i]);
    m.set(`${lo}:${hi}`, (m.get(`${lo}:${hi}`) ?? 0) + b.count[i]);
  }
  return m;
};

describe('the generated links (Features/Region reads)', () => {
  it('are 12,000, deterministic, every end inside its chromosome', () => {
    const a = linksOnChromosome();
    expect(a.aStart.length).toBe(LinkCount);
    expect(Array.from(linksOnChromosome().aStart.subarray(0, 50))).toEqual(Array.from(a.aStart.subarray(0, 50)));
    for (let i = 0; i < LinkCount; i++) {
      expect(a.aEnd[i]).toBeLessThanOrEqual(LinksChromosome.length);
      expect(a.bEnd[i]).toBeLessThanOrEqual(LinksChromosome.length);
    }
  });

  it('held and as a BEDPE are one picture: the same bin pairs with the same counts (EC-208)', async () => {
    // the file is served at its items' level (the store's summaries make its bins, as for held links)
    const items = await bedpePairs();
    expect(items.aStart.length).toBe(LinkCount);
    const heldCounts = counts(binPairs(linksOnChromosome(), 1024));
    expect(heldCounts.size).toBeGreaterThan(100);
    expect(counts(binPairs(items, 1024))).toEqual(heldCounts);
  });

  it('over a toy assembly are a source with bins at two levels, every link found from either end', async () => {
    const data = linksSource();
    expect(data.coverage.levels).toEqual([...LinkLevels]);
    const src = openSeq(data);
    let ends = 0;
    for (const s of LinksChromosomes) {
      const a = await src.read({ system: s.name, start: 0, end: s.length, level: 1, shapes: ['pairs'] }, signal);
      if (!a.ok) throw new Error(a.issue.message);
      const t = a.tables.pairs as PairsList;
      for (let i = 0; i < t.aStart.length; i++)
        ends += (t.aSystem!.get(i) === s.name ? 1 : 0) + (t.bSystem!.get(i) === s.name ? 1 : 0);
    }
    // each link counted once a chromosome end (twice when both ends are on one)
    expect(ends).toBe(2 * LinkCount);
    const bins = await src.read(
      { system: 'chr1', start: 0, end: LinksChromosomes[0].length, level: 16_384, shapes: ['pairs'] },
      signal,
    );
    expect(bins.ok && bins.tables.pairs?.form).toBe('bins');
  });

  it('on the genome ring they are structure, not a disc: a few strong blocks, local links, a thin scatter, in the arcs’ bound', async () => {
    const data = linksSource(RingChromosomes, LinkCount, 9, true);
    // deterministic
    const again = linksSource(RingChromosomes, LinkCount, 9, true);
    expect(again.pairs[0].answer).toEqual(data.pairs[0].answer);
    const src = openSeq(data);
    const level = LinkLevels[LinkLevels.length - 1];
    const bins = new Map<string, number>();
    for (const s of RingChromosomes) {
      const a = await src.read({ system: s.name, start: 0, end: s.length, level, shapes: ['pairs'] }, signal);
      if (!a.ok) throw new Error(a.issue.message);
      const t = a.tables.pairs as {
        aBin: Int32Array;
        bBin: Int32Array;
        count: Float32Array;
        aSystem: { get(i: number): string };
        bSystem: { get(i: number): string };
      };
      for (let i = 0; i < t.count.length; i++)
        bins.set(`${t.aSystem.get(i)}:${t.aBin[i]}|${t.bSystem.get(i)}:${t.bBin[i]}`, t.count[i]);
    }
    const counts = [...bins.values()].sort((x, y) => y - x);
    const total = counts.reduce((a, b) => a + b, 0);
    expect(total).toBe(LinkCount);
    // the strongest chords are bundles, an order above the thin ones (a disc of equal chords has no strongest)
    const median = counts[counts.length >> 1];
    expect(counts[9]).toBeGreaterThanOrEqual(6 * median);
    // every chord drawn: fewer bin pairs than the arcs' cap (no density strip), and links within a chromosome among them
    expect(bins.size).toBeLessThan(DefaultRenderConfig.pairs.maxArcs);
    expect([...bins.keys()].some((k) => k.split('|').every((e, _, a) => e.split(':')[0] === a[0].split(':')[0]))).toBe(
      true,
    );
  });
});
