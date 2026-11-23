/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The pairs reader (_spec/README.md "Render", The pairs layer: the reader): one reader of bins at a level and items in
// a region for held and region pairs, the built-in pairs layer reading through it. Its draw lists are held to the
// recording made before the reader was split out of the layer (byte-equal: the same frames, the same bytes), over one
// system (bins, items), held pairs past the gate, pairs over an assembly (bins at the fit, items on one block) and a
// bridge of two genomes; and a reader of its own (an extension's) reads what the layer reads, with its own gate.

import { murmurHash3_32 } from '@molstar/core/data/util/hash-functions';
import { categoryColumn } from '@molstar/core/repr/tree';
import { binPairs, type PairsList, type SeqData, SeqRuntime, seqToValue } from '@molstar/seqx-core/data';
import { assembly, createSeqSpec, genome } from '@molstar/seqx-core/lang';
import { BedRuntime } from '@molstar/seqx-core/runtimes/bed';
import { describe, expect, it } from 'vitest';
import { BuiltInLayerTypes } from '../built-ins.js';
import { type FrameInput, type LayerType, recordDrawLists } from '../index.js';
import { MoleculeLettersType } from '../letters-molecule.js';
import { type PairsLay, type PairsPrep, PairsType } from '../pairs.js';
import { type PairsRead, PairsReader, type PairsReadOptions, type PairsSource } from '../pairs-reader.js';
import { sideBlocks, sidesBridgeSpec } from './bridge-sides-fixtures.js';
import { settle } from './helpers.js';
import { LinksLength, linksSpec } from './region-render-fixtures.js';

/** Three chromosomes of 1 Mb, `n` links between hotspots listed on both ends' tiles, bins at 1 kb and 16 kb. */
function assemblyLinks(layout: 'rings' | 'tracks' = 'rings', n = 3000) {
  const names = ['c1', 'c2', 'c3'];
  const length = 1_000_000;
  let x = 11;
  const r = () => (x = (x * 1103515245 + 12345) % 2 ** 31) / 2 ** 31;
  const ends = Array.from({ length: n }, () => {
    const a = Math.floor(r() * 3);
    const b = Math.floor(r() * 3);
    return { a, b, as: Math.floor(r() * 20) * 48_000 + Math.floor(r() * 900), bs: Math.floor(r() * 20) * 48_000 };
  });
  const listOf = (k: number): PairsList => {
    const mine = ends.map((e, i) => ({ ...e, i })).filter((e) => e.a === k || e.b === k);
    return {
      form: 'list',
      aSystem: categoryColumn(mine.map((e) => names[e.a])),
      bSystem: categoryColumn(mine.map((e) => names[e.b])),
      aStart: Int32Array.from(mine, (e) => e.as),
      aEnd: Int32Array.from(mine, (e) => e.as + 10),
      bStart: Int32Array.from(mine, (e) => e.bs),
      bEnd: Int32Array.from(mine, (e) => e.bs + 10),
      id: categoryColumn(mine.map((e) => `l${e.i}`)),
    };
  };
  const levels = [1, 1024, 16_384];
  const data: SeqData = {
    coverage: {
      systems: names.map((name) => ({ name, length })),
      shapes: names.map((system) => ({ system, shape: 'pairs' as const })),
      levels,
    },
    pairs: names.flatMap((system, k) =>
      levels.map((level) => ({
        request: { system, start: 0, end: length, level, shapes: ['pairs' as const] },
        answer: { ok: true as const, tables: { pairs: level === 1 ? listOf(k) : binPairs(listOf(k), level) } },
      })),
    ),
  };
  const b = createSeqSpec({ title: 'links over an assembly' });
  const g = b.data.assembly(
    'toy',
    assembly(
      'Toy',
      names.map((name) => ({ name, length })),
    ),
  );
  const src = b.data.source('src', { format: 'seq', value: seqToValue(data) as never });
  const links = b.data.pairs('links', { a: g, source: src } as never);
  b.views.view('v', (v) => {
    genome(v, 'genome', g, { layout }, (s) => {
      s.track('links').layer('links', { type: 'pairs', data: links });
    });
  });
  return b.build();
}

/** 6,000 held pairs between 60 hotspots of a 2 Mb chromosome: past the gate, binned by the pairs level. */
function heldSpec() {
  const n = 6000;
  let x = 7;
  const r = () => (x = (x * 1103515245 + 12345) % 2 ** 31) / 2 ** 31;
  const as = Int32Array.from({ length: n }, () => Math.floor(r() * 60) * 33_000 + Math.floor(r() * 500));
  const bs = Int32Array.from({ length: n }, () => Math.floor(r() * 60) * 33_000 + Math.floor(r() * 500));
  const b = createSeqSpec({ title: 'held' });
  const sizes = b.data.source('sizes', { format: 'chrom-sizes', value: 'chr\t2000000\n' as never });
  const chr = b.data.sequence('chr', { alphabet: 'dna', source: sizes, name: 'chr' });
  const held = b.data.pairs('held', {
    a: chr,
    rows: { a_start: as, a_end: as.map((v) => v + 10), b_start: bs, b_end: bs.map((v) => v + 10) },
  } as never);
  b.views.view('v', (v) => v.section('chr', { of: chr }).track('held').layer('held', { type: 'pairs', data: held }));
  return b.build();
}

/** A spec's frame after its reads land, its draw lists recorded and hashed (the recording's length beside). */
async function recorded(spec: ReturnType<typeof linksSpec>, input: FrameInput, types?: Record<string, LayerType>) {
  const s = await settle(spec, [SeqRuntime, BedRuntime], {}, types ? { types } : {});
  await s.frame(input);
  await s.store.whenSettled();
  await s.frame(input);
  await s.store.whenSettled();
  const f = await s.frame(input);
  const text = recordDrawLists(f.lists);
  return { s, f, sig: `${text.length}:${murmurHash3_32(text, 7)}` };
}

const Wide = { width: 720, height: 300 };
const Square = { width: 720, height: 720 };

/** The recording before the reader was split out (`regionPairs`, `assemblyPairs`, `heldBins` called by the layer). */
const Before = {
  binsOneSystem: '190884:3637016521',
  itemsOneSystem: '296570:4050371117',
  held: '72599:608758320',
  assemblyFit: '351770:3501582133',
  assemblyBlock: '17728:569659795',
  bridge: '14125:4255137094',
};

describe('the pairs reader', () => {
  it('the pairs layer reads through it: its draw lists byte-equal to the layer’s before the split', async () => {
    const got = {
      binsOneSystem: (await recorded(linksSpec(6000), { size: Wide, viewport: { start: 0, end: LinksLength } })).sig,
      itemsOneSystem: (await recorded(linksSpec(6000), { size: Wide, viewport: { start: 0, end: 400_000 } })).sig,
      held: (await recorded(heldSpec(), { size: Wide })).sig,
      assemblyFit: (await recorded(assemblyLinks('rings', 6000), { size: Square })).sig,
      assemblyBlock: (
        await recorded(assemblyLinks('tracks'), {
          size: Wide,
          viewports: { '/views/v/genome': { start: 0, end: 60_000 } },
        } as FrameInput)
      ).sig,
      bridge: (await recorded(sidesBridgeSpec(sideBlocks(50, 3)), { size: Wide })).sig,
    };
    expect(got).toEqual(Before);
  });

  it('a reader of its own reads what the layer reads, at its own gate (bins where the layer shows items)', async () => {
    const reads: (PairsRead | undefined)[] = [];
    /** The pairs layer with a reader beside it (an extension's type would hold its own): asked with `options`. */
    const beside = (options?: PairsReadOptions): Record<string, LayerType> => {
      const reader = new PairsReader();
      const Beside: LayerType<PairsPrep, PairsLay> = {
        ...PairsType,
        layout(prep, ctx, data, overlay, prev) {
          reads.push(reader.read(prep as PairsSource, ctx, data, options));
          return PairsType.layout(prep, ctx, data, overlay, prev);
        },
      };
      return { ...BuiltInLayerTypes, letters: MoleculeLettersType, pairs: Beside };
    };
    const zoomed = { size: Wide, viewport: { start: 0, end: 400_000 } };
    // the layer's own numbers: the same read, items
    const same = await recorded(linksSpec(6000), zoomed, beside());
    const lay = same.f.layout.sections[0].sectionLays[0].lay as PairsLay;
    const mine = reads[reads.length - 1]!;
    expect([lay.prep!.level, mine.level, mine.data?.n]).toEqual([undefined, undefined, lay.prep!.data!.n]);
    expect(mine.data!.keys).toEqual(lay.prep!.data!.keys);
    // a gate of 1,024 (the ribbons'): 3,884 pairs in view are past it, so bins, every pair counted once
    reads.length = 0;
    await recorded(linksSpec(6000), zoomed, beside({ gate: 1024 }));
    const gated = reads[reads.length - 1]!;
    expect(gated.level).toBe(1024);
    expect(Array.from(gated.weight!).reduce((a, b) => a + b, 0)).toBe(6000);
  });

  it('held pairs read alike: binned past its gate, themselves under it', async () => {
    const reads: (PairsRead | undefined)[] = [];
    const reader = new PairsReader();
    const types = (options: PairsReadOptions): Record<string, LayerType> => ({
      ...BuiltInLayerTypes,
      letters: MoleculeLettersType,
      pairs: {
        ...PairsType,
        layout(prep, ctx, data, overlay, prev) {
          reads.push(reader.read(prep as PairsSource, ctx, data, options));
          return PairsType.layout(prep, ctx, data, overlay, prev);
        },
      } as LayerType<PairsPrep, PairsLay>,
    });
    await recorded(heldSpec(), { size: Wide }, types({ gate: 10_000 }));
    expect([reads[reads.length - 1]!.level, reader.items]).toEqual([undefined, true]);
    await recorded(heldSpec(), { size: Wide }, types({ gate: 100 }));
    expect(reads[reads.length - 1]!.level).toBeGreaterThan(1);
    expect(reader.items).toBe(false);
  });
});
