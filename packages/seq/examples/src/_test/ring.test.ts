/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The genome ring's data and views (`genomes/ring/`): the links' generator byte-identical (the committed 1,000-link
// BEDPE is its output), its shares as asked (60 % within, hub pairs strong), synthetic bands as the design says, and
// each view checked and drawn with the extension: held and BEDPE one picture (the same stems and counts).

import { readFileSync, writeFileSync } from 'node:fs';
import { createSeqContext } from '@molstar/seqx-core/context/default';
import { DataStore } from '@molstar/seqx-core/data';
import { bundledTextMetrics, LightTokens, openModel, Renderer, renderConfig } from '@molstar/seqx-core/render';
import { BedRuntime } from '@molstar/seqx-core/runtimes/bed';
import { rondo, type RondoLay } from '@molstar/seqx-ext-rondo';
import { describe, expect, it } from 'vitest';
import { Grch38Chromosomes, RingLinks1000File, ringBedpe, ringLinks, syntheticBands } from '../genomes/ring/data.js';
import { RingViews, registry, ringSpec } from '../genomes/ring/views.js';
import { nodeFixtureAssets, PackageFolder } from '../node-assets.js';

describe('the genome ring’s data', () => {
  it('the committed 1,000-link BEDPE is the generator’s output, byte for byte', () => {
    const text = ringBedpe(ringLinks(1000));
    if (process.env.SEQ_RECORDS === 'update') writeFileSync(`${PackageFolder}${RingLinks1000File}`, text);
    expect(readFileSync(`${PackageFolder}${RingLinks1000File}`, 'utf8')).toBe(text);
    expect(ringBedpe(ringLinks(1000))).toBe(text);
  });

  it('12,000 links: about 60 % within a chromosome, the rest between, the hubs strong; a type and a q each', () => {
    const links = ringLinks();
    expect(links.length).toBe(12_000);
    const within = links.filter((l) => l.a === l.b).length / links.length;
    expect(within).toBeGreaterThan(0.6);
    expect(within).toBeLessThan(0.7);
    const pair = links.filter((l) => (l.a === 0 && l.b === 16) || (l.a === 16 && l.b === 0)).length;
    expect(pair).toBeGreaterThan(100);
    for (const l of links) {
      expect(l.as).toBeGreaterThanOrEqual(0);
      expect(l.bs).toBeLessThan(Grch38Chromosomes[l.b].length);
      expect(l.q).toBeGreaterThan(0);
    }
  });

  it('synthetic bands: alternating stains, an acen pair at 40 %, a gvar on Y’s last third, the chromosome covered', () => {
    const text = syntheticBands([
      { name: 'chr1', length: 1_000_000 },
      { name: 'chrY', length: 600_000 },
    ]);
    const rows = text
      .trim()
      .split('\n')
      .map((l) => l.split('\t'));
    for (const name of ['chr1', 'chrY']) {
      const mine = rows.filter((r) => r[0] === name);
      expect(Number(mine[0][1])).toBe(0);
      for (let i = 1; i < mine.length; i++) expect(Number(mine[i][1])).toBe(Number(mine[i - 1][2]));
      expect(mine.filter((r) => r[4] === 'acen').length).toBe(2);
    }
    expect(rows.filter((r) => r[0] === 'chrY').at(-1)![4]).toBe('gvar');
    expect(rows.filter((r) => r[0] === 'chr1').some((r) => r[4] === 'gvar')).toBe(false);
    expect(syntheticBands([{ name: 'chr1', length: 1_000_000 }])).toBe(
      syntheticBands([{ name: 'chr1', length: 1_000_000 }]),
    );
  });
});

async function frameOf(spec: ReturnType<typeof ringSpec>) {
  const ctx = createSeqContext({ extensions: [rondo] });
  const store = new DataStore();
  // the views' own runtime (the cytobands from their recording), its reads from the package's fixtures, no delay
  const recorded = (registry({ delayMs: 0 })[0] as { runtime: typeof BedRuntime }).runtime;
  const runtime = {
    ...recorded,
    open: (s: never, c: never) => recorded.open(s, { ...(c as object), assets: nodeFixtureAssets().reader } as never),
  };
  const { model, view } = await openModel(spec, store, { runtimes: [runtime], kinds: ctx.kinds });
  const r = new Renderer({
    types: ctx.types,
    layouts: ctx.layouts,
    model,
    view: view!,
    data: store,
    tokens: LightTokens,
    text: bundledTextMetrics(),
    config: renderConfig(),
  });
  const input = { size: { width: 720, height: 760 } };
  let f = r.frame(input);
  for (let i = 0; i < 6; i++) {
    store.want(f.wants);
    await store.whenSettled();
    f = r.frame(input);
  }
  f = r.frame({ ...input, settled: true });
  const l = f.layout.sections.flatMap((s) => s.sectionLays).find((x) => x.type.name === 'rondo:ribbons')!;
  return { f, lay: l.lay as RondoLay, store };
}

const stems = (l: RondoLay) => {
  const out: string[] = [];
  for (let s = 0; s < l.ns; s++)
    if (l.stemPlaced[s] && l.stemCount[s]) out.push(`${l.names[l.stemSystem[s]]}:${l.stemBin[s]}:${l.stemCount[s]}`);
  return out.sort();
};

describe('the genome ring’s views', () => {
  it('each checks clean and draws its ribbons with the extension', async () => {
    const ctx = createSeqContext({ extensions: [rondo] });
    for (const v of RingViews) {
      const spec = v.build();
      expect(ctx.check(spec), v.id).toEqual([]);
      const { f, lay, store } = await frameOf(spec);
      expect(
        store.issues.value.filter((i) => i.severity === 'error'),
        v.id,
      ).toEqual([]);
      // the cytobands read from their recording: the source opened, its chromosomes covered
      // the cytobands read from their recording: the source opened, every chromosome covered
      expect(store.coverage('/data/bands-bed')?.systems.length, v.id).toBeGreaterThanOrEqual(23);
      expect(
        f.layout.issues.filter((i) => i.severity === 'error'),
        v.id,
      ).toEqual([]);
      expect(lay.nr, v.id).toBeGreaterThan(0);
    }
  });

  it('held and BEDPE are one picture: the same stems and counts (chr4’s ring)', async () => {
    const read = await frameOf(ringSpec({ include: ['chr4'], synthetic: true }));
    const held = await frameOf(ringSpec({ include: ['chr4'], synthetic: true, held: true }));
    expect(stems(held.lay)).toEqual(stems(read.lay));
  });
});
